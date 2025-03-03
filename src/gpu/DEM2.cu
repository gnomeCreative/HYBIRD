#include "DEM2.h"
#include <numeric>
#include <algorithm>

#include <cub/device/device_scan.cuh>
#include <cub/device/device_radix_sort.cuh>

#include "DEMParams.h"
#include "cub_temp_mem.h"
#include "lattice.h"
#include "LBParams.h"

__global__ void d_buildList_dual_bool(unsigned int* counter_a, unsigned int* buffer_a, const bool* test_a, const unsigned int count_a, 
                                      unsigned int* counter_b, unsigned int* buffer_b, const bool* test_b, const unsigned int count_b,
                                      const unsigned int threadCount) {    // Grid stride loop
    for (unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
        i < threadCount;
        i += blockDim.x * gridDim.x)
    {
        if (i < count_a && test_a[i]) {
            const unsigned int offset = atomicInc(counter_a, UINT_MAX);
            buffer_a[offset] = i;
        }
        if (i < count_b && test_b[i]) {
            const unsigned int offset = atomicInc(counter_b, UINT_MAX);
            buffer_b[offset] = i;
        }
    }
}
template<>
void DEM2::buildActiveLists<CUDA>() {
    // This is a simple implementation, there may be faster approaches
    // Alternate approach, stable pair-sort indices by type, then scan to identify boundaries

    // Ensure builder list is atleast min(19*hd_nodes.activeCount, count)
    auto& ctb = CubTempMem::GetBufferSingleton();
    ctb.resize((hd_elements.count + 1 + hd_particles.count + 1) * sizeof(unsigned int));
    unsigned int* builderI = reinterpret_cast<unsigned int*>(ctb.getPtr());
    // Init index 0 of each buffer to 0, this will be used as an atomic counter
    CUDA_CALL(cudaMemset(builderI, 0, sizeof(unsigned int)));
    CUDA_CALL(cudaMemset(builderI + hd_elements.count + 1, 0, sizeof(unsigned int)));
    // Launch kernel as grid stride loop
    int numSMs;
    cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, 0); // TODO Assumes device 0 in multi device system
    const unsigned int num_threads = std::max(hd_elements.count, hd_particles.count);
    d_buildList_dual_bool<<<32 * numSMs, 256>>>(builderI, &builderI[1], hd_elements.active, hd_elements.count,
                                                &builderI[hd_elements.count + 1], &builderI[hd_elements.count + 2], hd_particles.active, hd_particles.count,
                                      num_threads);
    CUDA_CHECK();
    // Copy back element result to main list
    unsigned int new_element_count = 0;
    CUDA_CALL(cudaMemcpy(&new_element_count, builderI, sizeof(unsigned int), cudaMemcpyDeviceToHost)); // illegal memory access?
    if (new_element_count > hd_elements.activeAlloc) {
        // Resize interface buffer (it doesn't currently ever scale back down)
        if (hd_elements.activeI) {
            CUDA_CALL(cudaFree(hd_elements.activeI));
        }
        CUDA_CALL(cudaMalloc(&hd_elements.activeI, new_element_count * sizeof(unsigned int)));
        hd_elements.activeAlloc = new_element_count;
    }
    hd_elements.activeCount = new_element_count;
    // Copy back element result to main list
    unsigned int new_particle_count = 0;
    CUDA_CALL(cudaMemcpy(&new_particle_count, builderI, sizeof(unsigned int), cudaMemcpyDeviceToHost)); // illegal memory access?
    if (new_particle_count > hd_particles.activeAlloc) {
        // Resize interface buffer (it doesn't currently ever scale back down)
        if (hd_particles.activeI) {
            CUDA_CALL(cudaFree(hd_particles.activeI));
        }
        CUDA_CALL(cudaMalloc(&hd_particles.activeI, new_particle_count * sizeof(unsigned int)));
        hd_particles.activeAlloc = new_particle_count;
    }
    hd_particles.activeCount = new_particle_count;
    // Sort list into it's new storage
    size_t temp_storage_bytes = 0;
    if (new_element_count > new_particle_count) {
        CUDA_CALL(cub::DeviceRadixSort::SortKeys(nullptr, temp_storage_bytes, &builderI[1], hd_elements.activeI, new_element_count));
    } else {
        CUDA_CALL(cub::DeviceRadixSort::SortKeys(nullptr, temp_storage_bytes, &builderI[hd_elements.count + 2], hd_particles.activeI, new_particle_count));
    }
    auto& ctm = CubTempMem::GetTempSingleton();
    ctm.resize(temp_storage_bytes);
    CUDA_CALL(cub::DeviceRadixSort::SortKeys(ctm.getPtr(), temp_storage_bytes, &builderI[1], hd_elements.activeI, new_element_count));
    CUDA_CALL(cub::DeviceRadixSort::SortKeys(ctm.getPtr(), temp_storage_bytes, &builderI[hd_elements.count + 2], hd_particles.activeI, new_particle_count));
    CUDA_CHECK()
    // Update device struct (new size and ptr, but whole struct because eh)
    CUDA_CALL(cudaMemcpy(d_elements, &hd_elements, sizeof(Element2), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_particles, &hd_particles, sizeof(Particle2), cudaMemcpyHostToDevice));
}
__global__ void d_atomicHistogram3D(
    Particle2 *d_particles, unsigned int *d_histogram) {
    const unsigned int i = (blockIdx.x * blockDim.x) + threadIdx.x;
    // Kill excess threads
    if (i >= d_particles->count) return;

    // Calculate this particle's bin within the neighbour grid
    const unsigned int hash = d_particles->x0[i].linearizePosition();
    if (hash >= (DEM_P.nCells[0] * DEM_P.nCells[1] * DEM_P.nCells[2]) + 1) {
        printf("Hash: %u >= %u at location (%f, %f, %f)\n", hash, (DEM_P.nCells[0] * DEM_P.nCells[1] * DEM_P.nCells[2]) + 1, d_particles->x0[i].x, d_particles->x0[i].y, d_particles->x0[i].z);
        
    }
    // Contribute to the histogram, and log our position within our neighbour grid bin
    const unsigned int bin_idx = atomicInc(&d_histogram[hash], 0xFFFFFFFF);
    d_particles->neighbour_index[i] = bin_idx;
}
__global__ void d_updateNeighbourIndex(
    Particle2* d_particles) {
    const unsigned int i = (blockIdx.x * blockDim.x) + threadIdx.x;
    // Kill excess threads
    if (i >= d_particles->count) return;

    // Calculate this particle's bin within the neighbour grid
    const unsigned int hash = d_particles->x0[i].linearizePosition();
    // Store our index within the neighbour grid so our data can be looked up
    const unsigned int bin_data_start = d_particles->PBM[hash];
    const unsigned int bin_sub_index = d_particles->neighbour_index[i];
    d_particles->neighbour_index[bin_data_start + bin_sub_index] = i;
}
template<>
void DEM2::evalNeighborTable<CUDA>() {
    // Build list of active elements
    buildActiveLists<CUDA>();
    // Allocate enough memory to store partition boundary matrix (PBM)
    // This will become an index to where each bin's data exists in memory.
    const unsigned int PBM_len = (DEM_P.nCells[0] * DEM_P.nCells[1] * DEM_P.nCells[2]) + 1;
    if (hd_particles.PBM_alloc < PBM_len) {
        if (hd_particles.PBM) {
            CUDA_CALL(cudaFree(hd_particles.PBM));
        }
        CUDA_CALL(cudaMalloc(&hd_particles.PBM, PBM_len * sizeof(unsigned int)));
        hd_particles.PBM_alloc = PBM_len;
        // Sync this change to device (just copy the whole struct)
        CUDA_CALL(cudaMemcpy(d_particles, &hd_particles, sizeof(Particle2), cudaMemcpyHostToDevice));
    }
    auto& ctb = CubTempMem::GetBufferSingleton();
    ctb.resize(PBM_len * sizeof(unsigned int));
    unsigned int *d_histogram = reinterpret_cast<unsigned int*>(ctb.getPtr());
    { // Build a histogram of particles using atomics
        CUDA_CALL(cudaMemset(d_histogram, 0x00000000, PBM_len * sizeof(unsigned int)));
        int blockSize = 0;  // The launch configurator returned block size
        int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
        int gridSize = 0;  // The actual grid size needed, based on input size
        cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_atomicHistogram3D, 0, hd_particles.count);
        // Round up to accommodate required threads
        gridSize = (hd_particles.count + blockSize - 1) / blockSize;
        d_atomicHistogram3D<<<gridSize, blockSize>>>(d_particles, d_histogram);
        CUDA_CHECK()
    }
    {  // Scan (sum) histogram, to finalise PBM
        size_t temp_storage_bytes = 0;
        CUDA_CALL(cub::DeviceScan::ExclusiveSum(nullptr, temp_storage_bytes, d_histogram, hd_particles.PBM, PBM_len));
        auto& ctm = CubTempMem::GetTempSingleton();
        ctm.resize(temp_storage_bytes);
        CUDA_CALL(cub::DeviceScan::ExclusiveSum(ctm.getPtr(), temp_storage_bytes, d_histogram, hd_particles.PBM, PBM_len));
    }
    { // Store particle indexes in the neighbour grid
        int blockSize = 0;  // The launch configurator returned block size
        int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
        int gridSize = 0;  // The actual grid size needed, based on input size
        cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_updateNeighbourIndex, 0, hd_particles.count);
        // Round up to accommodate required threads
        gridSize = (hd_particles.count + blockSize - 1) / blockSize;
        d_updateNeighbourIndex <<<gridSize, blockSize>>>(d_particles);
    }
}

/**
 * predictor()
 */
__global__ void d_predictor(Element2* d_elements) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    // Kill excess threads early
    if (i >= d_elements->activeCount) return;

    d_elements->predict(d_elements->activeI[i]);
}
template<>
void DEM2::predictor<CUDA>() {
    // Launch cuda kernel to update
    int blockSize = 0;  // The launch configurator returned block size
    int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
    int gridSize = 0;  // The actual grid size needed, based on input size
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_predictor, 0, hd_elements.activeCount);
    // Round up to accommodate required threads
    gridSize = (hd_elements.activeCount + blockSize - 1) / blockSize;
    d_predictor<<<gridSize, blockSize>>>(d_elements);
    CUDA_CHECK();
}

/**
 * updateParticlesPredicted()
 */
__global__ void d_updateParticlesPredicted(Particle2* d_particles, Element2* d_elements) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    // Kill excess threads early
    if (i >= d_particles->activeCount) return;

    d_particles->updatePredicted(d_particles->activeI[i], d_elements);
}
template<>
void DEM2::updateParticlesPredicted<CUDA>() {
    // Launch cuda kernel to update
    int blockSize = 0;  // The launch configurator returned block size
    int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
    int gridSize = 0;  // The actual grid size needed, based on input size
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_updateParticlesPredicted, 0, hd_particles.activeCount);
    // Round up to accommodate required threads
    gridSize = (hd_particles.activeCount + blockSize - 1) / blockSize;
    d_updateParticlesPredicted<<<gridSize, blockSize>>>(d_particles, d_elements);
    CUDA_CHECK();
}

/**
 * evaluateForces()
 */
__host__ __device__ __forceinline__ void common_resetElementForces(const unsigned int ae_i, Element2* elements) {
    elements->FParticle[ae_i].reset();
    elements->FWall[ae_i].reset();
    elements->MParticle[ae_i].reset();
    elements->MWall[ae_i].reset();
    elements->FSpringP[ae_i].reset();
    elements->FSpringW[ae_i].reset();
    elements->MRolling[ae_i].reset();
    // reset also apparent accellerations
    elements->ACoriolis[ae_i].reset();
    elements->ACentrifugal[ae_i].reset();
    // reset also maximum overlap (moved to IO function)
    elements->maxDtOverlap[ae_i] = std::max(elements->maxDtOverlap[ae_i], elements->maxOverlap[ae_i]);
    elements->maxOverlap[ae_i] = 0;
    elements->slippingCase[ae_i] = 0;
    // information about phase change and connectivity
    elements->solidIntensity[ae_i].reset();
    elements->coordination[ae_i] = 0;
}
__global__ void d_resetForces(Particle2* d_particles, Element2* d_elements, Wall2* d_walls, Object2 *d_objects) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (DEM_P.staticFrictionSolve && i < d_particles->activeCount) {
        // @todo springs
        //for (int t = 0; t < 4; ++t) {
        //    for (int s = 0; s < particles[p].springs[t].size(); ++s) {
        //        particles[p].springs[t][s].active = false;
        //    }
        //}
    }
    if (i < d_elements->activeCount) {
        common_resetElementForces(d_elements->activeI[i], d_elements);
    }
    if (i < d_walls->count) {
        d_walls->FParticle[i].reset();
    }
    if (i < d_objects->count) {
        d_objects->FParticle[i].reset();
    }
}
__device__ __forceinline__ double normalContact(const double overlap, const double vrelnnorm, const double effRad, const double effMass) {

    // total normal force
    double fn = 0.0;

    switch (DEM_P.sphereMat.contactModel) {
    case HERTZIAN:
    {
        // square root of effective radius
        const double sqrtEffRad = sqrt(effRad);
        // square root of overlap
        const double sqrtOverlap = sqrt(overlap);
        // normal stiffness = 2/3*(Reff^1/2)*Y/(1-nu²)
        const double kn = DEM_P.sphereMat.knConst * sqrtEffRad * sqrtOverlap;
        // damping
        const double gamman = 2.0 * DEM_P.sphereMat.dampCoeff * sqrt(kn * effMass); // 2 was missing
        const double dumpfn = -gamman * vrelnnorm;
        // elastic normal force (Hertzian contact)
        const double elasticfn = kn * overlap;
        // total normal force (limit to repulsive force)
        fn = std::max(elasticfn + dumpfn, 0.0);
        break;
    }
    case LINEAR:
    {
        // elastic normal force (linear contact)
        const double elasticfn = DEM_P.sphereMat.linearStiff * overlap;
        // damping
        const double gamman = 2.0 * DEM_P.sphereMat.dampCoeff * sqrt(DEM_P.sphereMat.linearStiff * effMass);
        const double dumpfn = -gamman * vrelnnorm;
        // total normal force (limit to repulsive force)
        fn = std::max(elasticfn + dumpfn, 0.0);
        break;
    }
    }

    return fn;

}
__device__ __forceinline__ tVect FRtangentialContact(const tVect& tangRelVelContact, const double fn, const double overlap, const double effRad,
        const double effMass, Elongation* elongation_new, const double friction, const double tangStiff, const double viscTang) {
    // tangential force
    tVect fs;
    // tangent stiffness
    double ks = 0.0;
    switch (DEM_P.sphereMat.contactModel) {
        case HERTZIAN:
        {
            // static const double power = 1.0 / 3.0;
            // square root of effective radius
            const double sqrtEffRad = sqrt(effRad);
            // square root of overlap
            const double sqrtOverlap = sqrt(overlap);
            // normal stiffness (from normal contac)
            const double kn = DEM_P.sphereMat.knConst * sqrtEffRad*sqrtOverlap;
            // tangent stiffness -> not physically sound, static friction is actually missing
            ks = 2.0 / 7.0 * kn; //sphereMat.ksConst * sqrtEffRad * std::pow(std::abs(fn), power);
            break;
        }
        case LINEAR:
        {
            ks = 2.0 / 7.0 * tangStiff;
            break;
        }
    }

    // maximum force due to static friction
    const double fsMax = friction*fn;
    const double fdMax = 0.9 * friction*fn;

    // viscous force
    const tVect viscousForce = 2.0 * viscTang * sqrt(effMass * ks) * tangRelVelContact;

    // model with static friction
    if (DEM_P.staticFrictionSolve) {
        /* @todo springs
        const tVect F_control = elongation_new->e * ks + viscousForce;

        const double F_control_norm = F_control.norm();

        tVect et = tVect(0.0, 0.0, 0.0);
        if (F_control_norm != 0) {
            et = F_control / F_control_norm;
        }

        if (elongation_new->slipping) {
            if (F_control_norm < fdMax) {
                // was slipping, goes back to static
                fs = F_control;
                elongation_new->e = elongation_new->e + tangRelVelContact* DEM_P.deltat;
                elongation_new->slipping = false;
                elongation_new->slippingCase = 1;
            } else {
                // keeps slipping
                fs = fdMax*et;
                const tVect f = fdMax * et - viscousForce;
                elongation_new->e = f / ks;
                elongation_new->slippingCase = 2;
            }

        } else {

            if (F_control_norm < fsMax) {
                // keeps being static
                fs = F_control;
                elongation_new->e = elongation_new->e + tangRelVelContact* DEM_P.deltat;
                elongation_new->slippingCase = 3;

            } else {
                // was static, begins slipping
                fs = fdMax*et;
                elongation_new->slipping = true;
                const tVect f = fdMax * et - viscousForce;
                elongation_new->e = f / ks;
                elongation_new->slippingCase = 4;
            }

        }
        */

    }// model without static friction
    else {
        const double viscousForceNorm = viscousForce.norm();
        const double viscousForceReal = std::min(viscousForceNorm, fsMax);
        fs = viscousForce / viscousForceNorm*viscousForceReal;
    }
    return fs;
}

__device__ __forceinline__ tVect rollingContact(const tVect& wI, const tVect& wJ, const double effRad, const double fn, const double rolling) {

    const tVect wrel = wI - wJ;
    const double wrel_norm = wrel.norm();

    tVect wIJ = tVect(0.0, 0.0, 0.0);
    if (wrel_norm != 0.0) {
        wIJ = wrel / wrel_norm;
    }

    const tVect fr = wIJ * rolling * fn * effRad;

    return fr;
}

__device__ __forceinline__ void d_particleParticleCollision(Particle2* d_particles, Element2* d_elements, const unsigned int p_i, const unsigned int p_j, const tVect& vectorDistance, Elongation* elongation_new) {
    /**
     * @todo This could potentially be improved by calculating the delta for many of the vectors (MParticle, FParticle, SolidIntensity)
     * @todo and then updating them globally once rather than with each force calc.
     * @todo This would increase register pressure for reduced atomic contention, which feels like a win
     */
    // pointers to elements
    const unsigned int e_i = d_particles->clusterIndex[p_i];
    const unsigned int e_j = d_particles->clusterIndex[p_j];

    // NORMAL FORCE ///////////////////////////////////////////////////////////////////

    // geometry /////////////////////////////
    const double radI = d_particles->r[p_i];
    const double radJ = d_particles->r[p_j];
    // distance norm
    const double distance = vectorDistance.norm();
    // overlap
    const double overlap = radI + radJ - distance;
    // relative velocity
    const tVect relVel = d_particles->x1[p_j] - d_particles->x1[p_i];
    // first local unit vector (normal)
    const tVect en = vectorDistance / distance;
    // relative normal velocity (modulus)
    const double normRelVel = relVel.dot(en);
    // relative normal velocity
    const tVect normalRelVel = en*normRelVel;
    // effective mass
    const double effMass = d_elements->m[e_i] * d_elements->m[e_j] / (d_elements->m[e_i] + d_elements->m[e_j]);
    // effective radius
    const double effRad = radI * radJ / (radI + radJ);
    /*cout<<"overlap "<<overlap<<endl;
    cout<<"effRad "<<effRad<<" effMass "<<effMass<<" normRelVel "<<normRelVel<<" en";
    en.show();
    vectorDistance.show();
    cout<<distance<<" "<<partI->particleIndex<<" "<<partJ->particleIndex;*/

    // force computation /////////////////////////////////
    const double normNormalForce = normalContact(overlap, normRelVel, effRad, effMass);
    assert(normNormalForce >= 0.0);

    //    // Overlap elastic potential energy
    //    energy.elastic+=0.5*fn*xi*xi;

    // force application ///////////////////////////////////
    // vectorial normal force
    const tVect normalForce = en*normNormalForce;

    // cluster geometry
    // vectorized radii
    const tVect vecRadI = radI*en;
    const tVect vecRadj = -radJ*en;
    // vectorized distance contactPoint-center of cluster
    const tVect centerDistI = vecRadI + d_particles->radiusVec[p_i];
    const tVect centerDistJ = vecRadj + d_particles->radiusVec[p_j];

    // @todo Particles are never marked as ghost, does this matter? (e.g. can these if's be removed as always true)
    if (!d_particles->isGhost[p_i]) {
        d_elements->FParticle[e_i]._atomicSub(normalForce);
        d_elements->solidIntensity[e_i]._atomicAdd(normalForce.abs());
        //  moment generated in non-spherical particles
        if (d_elements->size[e_i] > 1) {
            d_elements->MParticle[e_i]._atomicSub(centerDistI.cross(normalForce));
        }
    }
    if (!d_particles->isGhost[p_j]) {
        d_elements->FParticle[e_j]._atomicAdd(normalForce);
        d_elements->solidIntensity[e_j]._atomicAdd(normalForce.abs());
        //  moment generated in non-spherical particles
        if (d_elements->size[e_j] > 1) {
            d_elements->MParticle[e_j]._atomicAdd(centerDistJ.cross(normalForce));
        }
    }

    // TANGENTIAL FORCE ///////////////////////////////////////////////////////////////////

    // angular velocities (global reference frame)
    const tVect wI = d_elements->wpGlobal[e_i]; //2.0*quat2vec( elmtI->qp1.multiply( elmtI->qp0.adjoint() ) );
    const tVect wJ = d_elements->wpGlobal[e_j]; //2.0*quat2vec( elmtJ->qp1.multiply( elmtJ->qp0.adjoint() ) );
    // relative velocity at contact point
    const tVect relVelContact = relVel - wI.cross(vecRadI) + wJ.cross(vecRadj);
    // tangential component of relative velocity
    const tVect tangRelVelContact = relVelContact - normalRelVel;
    // norm of tangential velocity
    const double normTangRelVelContact = tangRelVelContact.norm();
    // checking if there is any tangential motion
    if (normTangRelVelContact != 0.0) {

        // update spring
        if (DEM_P.staticFrictionSolve) {
            /* @todo springs
            const double elong_old_norm = elongation_new->e.norm();

            double normElong = elongation_new->e.dot(en);
            const tVect normalElong = en*normElong;

            elongation_new->e = elongation_new->e - normalElong;

            if (elongation_new->e.norm() != 0) {
                const double scaling = elong_old_norm / elongation_new->e.norm();
                elongation_new->e = elongation_new->e*scaling;
            } else {
                const double scaling = 0.0;
                elongation_new->e = elongation_new->e*scaling;
            }
            */
        }

        //elongation.e.show();
        const tVect tangForce = FRtangentialContact(tangRelVelContact, normNormalForce, overlap, effRad, effMass, elongation_new, DEM_P.sphereMat.frictionCoefPart,
                DEM_P.sphereMat.linearStiff, DEM_P.sphereMat.viscTang);

        // torque updating
        // @todo Particles are never marked as ghost, does this matter? (e.g. can these if's be removed as always true)
        if (!d_particles->isGhost[p_i]) {
            d_elements->MParticle[e_i]._atomicAdd(centerDistI.cross(tangForce));
            d_elements->FSpringP[e_i]._atomicAdd(tangForce);
            d_elements->FParticle[e_i]._atomicAdd(tangForce);
            d_elements->solidIntensity[e_i]._atomicAdd(tangForce.abs());
            /* @todo springs
            if (DEM_P.staticFrictionSolve) {
                // @todo, why is this equals? couldn't there be multiple that take this path on same element?
                d_elements->slippingCase[e_i] = elongation_new->slippingCase;
            }
            */

        }
        if (!d_particles->isGhost[p_j]) {
            d_elements->MParticle[e_j]._atomicSub(centerDistJ.cross(tangForce));
            d_elements->FSpringP[e_j]._atomicSub(tangForce);
            d_elements->FParticle[e_j]._atomicSub(tangForce);
            d_elements->solidIntensity[e_j]._atomicAdd(tangForce.abs());
            /* @todo springs
            if (DEM_P.staticFrictionSolve) {
                // @todo, why is this equals? couldn't there be multiple that take this path on same element?
                d_elements->slippingCase[e_j] = elongation_new->slippingCase;
            }
            */
        }

    }


    //ROLLING
    const tVect rollingMoment = rollingContact(wI, wJ, effRad, normNormalForce,
        DEM_P.sphereMat.rollingCoefPart);

    // @todo Particles are never marked as ghost, does this matter? (e.g. can these if's be removed as always true)
    if (!d_particles->isGhost[p_i]) {
        d_elements->MRolling[e_i]._atomicSub(rollingMoment);
    }
    if (!d_particles->isGhost[p_j]) {
        d_elements->MRolling[e_j]._atomicAdd(rollingMoment);
    }

    // save overlap
    d_elements->atomicMaxOverlap(e_i, overlap);
    d_elements->atomicMaxOverlap(e_j, overlap);

    // updating connectivity
    atomicInc(&d_elements->coordination[e_i], std::numeric_limits<unsigned int>::max());
    atomicInc(&d_elements->coordination[e_j], std::numeric_limits<unsigned int>::max());
}
__global__ void d_particleParticleContacts(Particle2 *d_particles, Element2 *d_elements) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    // Kill excess threads early
    if (tid >= d_particles->activeCount)
        return;

    const unsigned int a_i = d_particles->activeI[tid];
    const unsigned int my_cluster = d_particles->clusterIndex[a_i];

    // Iterate a_i's Moore neighbourhood to find neighbours particles
    // @todo reduce bins checked/make forces apply recriprocally, what's best way to layout threads if we want to process 1 bin per thread?
    // a_i's 3D bin coordinate
    const tVect pos = d_particles->x0[a_i];
    const std::array<int, 3> bin_pos = {
    static_cast<int>(floor(pos.x / DEM_P.cellWidth[0]) + 1),
    static_cast<int>(floor(pos.y / DEM_P.cellWidth[1]) + 1),
    static_cast<int>(floor(pos.z / DEM_P.cellWidth[2]) + 1) };
    for (int dx = -1; dx <= 1; ++dx) {
        // Clamp/Wrap x coord
        int x = bin_pos[0] + dx;
        if (x < 0) {
            if (LB_P.boundary[0] == PERIODIC) {
                x = DEM_P.nCells[0] - 1;  // Note bugprone if less than 3 bins across
            } else {
                continue;
            }
        } else if (x >= DEM_P.nCells[0]) {
            if (LB_P.boundary[1] == PERIODIC) {
                x = 0;
            } else {
                continue;
            }
        }
        for (int dy = -1; dy <= 1; ++dy) {
            // Clamp/Wrap y coord
            int y = bin_pos[1] + dy;
            if (y < 0) {
                if (LB_P.boundary[2] == PERIODIC) {
                    y = DEM_P.nCells[1] - 1;  // Note bugprone if less than 3 bins across
                } else {
                    continue;
                }
            } else if (y >= DEM_P.nCells[1]) {
                if (LB_P.boundary[3] == PERIODIC) {
                    y = 0;
                } else {
                    continue;
                }
            }
            for (int dz = -1; dz <= 1; ++dz) {
                // Clamp/Wrap y coord
                int z = bin_pos[1] + dz;                
                if (z < 0) {
                    if (LB_P.boundary[4] == PERIODIC) {
                        z = DEM_P.nCells[2] - 1;  // Note bugprone if less than 3 bins across
                    } else {
                        continue;
                    }
                } else if (z >= DEM_P.nCells[2]) {
                    if (LB_P.boundary[5] == PERIODIC) {
                        z = 0;
                    } else {
                        continue;
                    }
                }
                // Iterate particles within the bin
                const unsigned int hash = x + DEM_P.nCells[0] * (y + DEM_P.nCells[1] * x);
                for (unsigned int i = d_particles->PBM[hash]; i < d_particles->PBM[hash+1]; ++i){
                    const unsigned int n_i = d_particles->neighbour_index[i];
                    // @temp Filter out particles that have lower index, we are operating naive reciprocity
                    if (n_i > a_i) {
                        // Filter out particles that are same cluster or not in contact
                        if (d_particles->clusterIndex[n_i] != my_cluster) {
                            // checking for overlap
                            const double ri = d_particles->r[a_i];
                            const double rj = d_particles->r[n_i];
                            const double sigij = ri + rj;
                            const double sigij2 = sigij * sigij;
                            // Convert the neighbour's location to it's closest ghost position in a wrapped environment
                            tVect n_x0 = d_particles->x0[n_i];
                            const tVect an_x0 = n_x0 - d_particles->x0[a_i];
                            n_x0 = 
                            {
                                abs(an_x0.x) > DEM_P.half_demSize.x ? n_x0.x - (an_x0.x / abs(an_x0.x) * DEM_P.half_demSize.x) : n_x0.x,
                                abs(an_x0.y) > DEM_P.half_demSize.y ? n_x0.y - (an_x0.y / abs(an_x0.y) * DEM_P.half_demSize.y) : n_x0.y,
                                abs(an_x0.z) > DEM_P.half_demSize.z ? n_x0.z - (an_x0.z / abs(an_x0.z) * DEM_P.half_demSize.z) : n_x0.z,
                            };
                            // recalculate distance between centers with virtual neighbour position
                            const tVect vectorDistance = d_particles->x0[a_i] - n_x0;
                            const double distance2 = an_x0.norm2();
                            // check for contact
                            if (distance2 < sigij2) {
                                // pointer to elongation, initially pointing to an empty spring
                                Elongation* elongation_here_new;
                                /* @TODO SPRINGS
                                if (DEM_P.staticFrictionSolve) {
                                    elongation_here_new = findSpring(0, a_i, n_i);
                                }*/
                                d_particleParticleCollision(d_particles, d_elements, a_i, n_i, vectorDistance, elongation_here_new);
                            }
                        }
                    }
                }
            }
        }
    }

}
template<>
void DEM2::particleParticleContacts<CUDA>() {
    // Launch cuda kernel to update
    int blockSize = 0;  // The launch configurator returned block size
    int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
    int gridSize = 0;  // The actual grid size needed, based on input size
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_particleParticleContacts, 0, hd_particles.activeCount);
    // Round up to accommodate required threads
    gridSize = (hd_particles.activeCount + blockSize - 1) / blockSize;
    d_particleParticleContacts<<<gridSize, blockSize>>>(d_particles, d_elements);
    CUDA_CHECK();
}

__device__ __forceinline__ void d_wallParticleCollision(Particle2 *d_particles, Wall2 *d_walls, Element2* d_elements, const unsigned int w_i, const unsigned int p_i, const double overlap, Elongation* elongation_new) {

    // pointers to element
    const unsigned int e_i = d_particles->clusterIndex[p_i];

    // NORMAL FORCE ///////////////////////////////////////////////////////////////////

    // geometry ///////////////
    // particle radius
    const double radJ = d_particles->r[p_i];
    // first local unit vector (normal)
    const tVect en = d_walls->n[w_i];
    // speed of the wall at contact point
    const tVect contactPointVelocity = d_walls->getSpeed(w_i, d_particles->x0[p_i]); // fix this, contact point not defined
    // relative velocity
    const tVect relVel = d_particles->x1[p_i] - contactPointVelocity;
    // relative normal velocity (modulus)
    const double normRelVel = relVel.dot(en);
    // relative normal velocity
    const tVect normalRelVel = en*normRelVel;

    // force computation /////////////////////////////////
    double normNormalForce = normalContact(overlap, normRelVel, radJ, d_elements->m[e_i]); // removed 2.0 *

    //    switch (problemName) {
    //        case TRIAXIAL:
    //        {
    //            if (wallI->index > 5) {
    //                normNormalForce /= 10;
    //            }
    //        }
    //    }


    // Overlap elastic potential energy
    //                    energy.elastic+=fn*xj/2.5;

    // force application ///////////////////////////////////
    // vectorial normal force
    const tVect normalForce = en*normNormalForce;

    // cluster geometry
    // vectorized radius
    const tVect vecRadJ = -radJ*en;
    // vectorized distance contactPoint-center of cluster
    tVect centerDistJ = vecRadJ;
    if (d_elements->size[e_i] > 1) {
        // vectorized distance contactPoint-center of cluster
        centerDistJ += (d_particles->x0[p_i] - d_elements->xp0[e_i]);
    }

    // force updating
    d_elements->FWall[e_i]._atomicAdd(normalForce);
    d_walls->FParticle[w_i]._atomicSub(normalForce);
    d_elements->solidIntensity[e_i]._atomicAdd(normalForce.abs());
    // torque updating
    if (d_elements->size[e_i] > 1) {
        d_elements->MWall[e_i]._atomicAdd(centerDistJ.cross(normalForce));
    }

    // TANGENTIAL FORCE ///////////////////////////////////////////////////////////////////

    // angular velocities (global reference frame)
    const tVect wJ = d_elements->wpGlobal[e_i]; //2.0*quat2vec( elmtJ->qp1.multiply( elmtJ->qp0.adjoint() ) );
    // relative velocity at contact point
    const tVect relVelContact = relVel + wJ.cross(vecRadJ);
    // tangential component of relative velocity
    const tVect tangRelVelContact = relVelContact - normalRelVel;
    // norm of tangential velocity
    const double normTangRelVelContact = tangRelVelContact.norm();


    // checking if there is any tangential motion
    if (normTangRelVelContact != 0.0) {

        // update spring
        /** @todo springs
        if (DEM_P.staticFrictionSolve) {
            const double elong_old_norm = elongation_new->e.norm();

            double normElong = elongation_new->e.dot(en);
            const tVect normalElong = en*normElong;

            elongation_new->e = elongation_new->e - normalElong;

            if (elongation_new->e.norm() != 0) {
                const double scaling = elong_old_norm / elongation_new->e.norm();
                elongation_new->e = elongation_new->e*scaling;
            } else {
                const double scaling = 0.0;
                elongation_new->e = elongation_new->e*scaling;
            }
        }
        */

        tVect tangForce = FRtangentialContact(tangRelVelContact, normNormalForce, overlap, radJ, d_elements->m[e_i], elongation_new, DEM_P.sphereMat.frictionCoefWall,
            DEM_P.sphereMat.linearStiff, DEM_P.sphereMat.viscTang);

        // torque updating
        d_elements->MWall[e_i]._atomicSub(centerDistJ.cross(tangForce));
        // force updating
        d_elements->FSpringW[e_i]._atomicSub(tangForce);
        d_elements->FWall[e_i]._atomicSub(tangForce);
        d_elements->solidIntensity[e_i]._atomicAdd(tangForce.abs());
        d_walls->FParticle[w_i]._atomicAdd(tangForce);

        /** @todo springs
        if (DEM_P.staticFrictionSolve) {
            // @todo, this is suspicious, race condition?
            elmtJ->slippingCase = elongation_new->slippingCase;
        }
        */

    }
    //ROLLING

    const tVect rollingMoment = rollingContact({0,0,0}, wJ, radJ, normNormalForce,
            2.0 * DEM_P.sphereMat.rollingCoefPart);

    // @todo Particles are never currently ghost, this if is always true, remove?
    if (!d_particles->isGhost[p_i]) {
        d_elements->MRolling[e_i]._atomicAdd(rollingMoment);
    }

    // save overlap
    d_elements->atomicMaxOverlap(e_i, overlap);

    // updating connectivity
    atomicInc(&d_elements->coordination[e_i], std::numeric_limits<unsigned int>::max());
}
__global__ void d_wallParticleContacts(Particle2* d_particles, Wall2* d_walls, Element2 *d_elements) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    // Kill excess threads early
    if (tid >= d_particles->activeCount)
        return;

    const unsigned int p_i = d_particles->activeI[tid];
    const tVect p_x0 = d_particles->x0[p_i];
    const double p_x0_Xp = p_x0.dot(Xp);
    const double p_x0_Yp = p_x0.dot(Yp);
    const double p_x0_Zp = p_x0.dot(Zp);

    // Check each particle against every wall
    // We don't (currently) have a nearWallTable
    for (unsigned int w_i = 0; w_i < d_walls->count; ++w_i) {
        

        // distance from wall (norm)
        const double distance = d_walls->dist(w_i, p_x0);

        if (d_walls->limited[w_i]) {
            if (p_x0_Xp < d_walls->xMin[w_i] ||
                p_x0_Xp > d_walls->xMax[w_i] ||
                p_x0_Yp < d_walls->yMin[w_i] ||
                p_x0_Yp > d_walls->yMax[w_i] ||
                p_x0_Zp < d_walls->zMin[w_i] ||
                p_x0_Zp > d_walls->zMax[w_i] ||
                // check if we are above the plane
                distance < 0.0
                ) {
                continue;
            }
        }
        
        // radius
        const double rj = d_particles->r[p_i];

        // distance before contact
        const double overlap = rj - distance;

        if (overlap > 0.0) {
            // pointer to elongation, initially pointing to an empty spring
            Elongation* elongation_here_new;
            /** @todo springs
            if (DEM_P.staticFrictionSolve) {
                const unsigned int indexI = d_walls->index[w_i];
                elongation_here_new = findSpring(1, indexI, p_i);
            }
            */
            d_wallParticleCollision(d_particles, d_walls, d_elements, w_i, p_i, overlap, elongation_here_new);
        }
    }
}
template<>
void DEM2::wallParticleContacts<CUDA>() {
    // Launch cuda kernel to update
    int blockSize = 0;  // The launch configurator returned block size
    int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
    int gridSize = 0;  // The actual grid size needed, based on input size
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_particleParticleContacts, 0, hd_particles.activeCount);
    // Round up to accommodate required threads
    gridSize = (hd_particles.activeCount + blockSize - 1) / blockSize;
    d_wallParticleContacts<<<gridSize, blockSize>>>(d_particles, d_walls, d_elements);
    CUDA_CHECK();
}
__device__ __forceinline__ void d_objectParticleCollision(Particle2 *d_particles, Object2 *d_objects, Element2* d_elements, const unsigned int o_i, const unsigned int p_i, const tVect& vectorDistance, Elongation* elongation_new) {
    // pointers to element
    const unsigned int e_i = d_particles->clusterIndex[p_i];

    
    // NORMAL FORCE ///////////////////////////////////////////////////////////////////

    // geometry ///////////////
    // particle radius
    const double radJ = d_particles->r[p_i];
    // distance from object (norm)
    const double distance = vectorDistance.norm();
    // distance before contact
    double overlap = radJ + d_objects->r[o_i] - distance;
    // first local unit vector (normal)
    const tVect en = (1.0 / distance) * vectorDistance;
    // speed of the wall at contact point
    const tVect contactPointVelocity = d_objects->x1[o_i]; // fix this, contact point not defined
    // relative velocity
    const tVect relVel = d_particles->x1[p_i] - contactPointVelocity;
    // relative normal velocity (modulus)
    const double normRelVel = relVel.dot(en);
    // relative normal velocity
    const tVect normalRelVel = en*normRelVel;

    // force computation /////////////////////////////////
    const double normNormalForce = normalContact(overlap, normRelVel, radJ, d_elements->m[e_i]); // was 2.0 * overlap

    // Overlap elastic potential energy
    //                    energy.elastic+=fn*xj/2.5;

    // force application ///////////////////////////////////
    // vectorial normal force
    const tVect normalForce = en*normNormalForce;

    // cluster geometry
    // vectorized radius
    const tVect vecRadJ = -radJ*en;
    // vectorized distance contactPoint-center of cluster
    const tVect centerDistJ = vecRadJ + (d_particles->x0[p_i] - d_elements->xp0[e_i]);

    // force updating
    d_elements->FWall[e_i]._atomicAdd(normalForce);
    d_objects->FParticle[o_i]._atomicSub(normalForce);
    d_elements->solidIntensity[e_i]._atomicAdd(normalForce.abs());

    // torque updating
    if (d_elements->size[e_i] > 1) {
        d_elements->MWall[e_i]._atomicAdd(centerDistJ.cross(normalForce));
    }

    // TANGENTIAL FORCE ///////////////////////////////////////////////////////////////////

    // angular velocities (global reference frame)
    const tVect wJ = d_elements->wpGlobal[e_i]; //2.0*quat2vec( elmtJ->qp1.multiply( elmtJ->qp0.adjoint() ) );
    // relative velocity at contact point
    const tVect relVelContact = relVel + wJ.cross(vecRadJ);
    // tangential component of relative velocity
    const tVect tangRelVelContact = relVelContact - normalRelVel;
    // norm of tangential velocity
    const double normTangRelVelContact = tangRelVelContact.norm();
    // checking if there is any tangential motion
    if (normTangRelVelContact != 0.0) {
        // update spring
        if (DEM_P.staticFrictionSolve) {

            const double elong_old_norm = elongation_new->e.norm();

            double normElong = elongation_new->e.dot(en);
            const tVect normalElong = en*normElong;

            elongation_new->e = elongation_new->e - normalElong;

            if (elongation_new->e.norm() != 0) {
                const double scaling = elong_old_norm / elongation_new->e.norm();
                elongation_new->e = elongation_new->e*scaling;
            } else {
                const double scaling = 0.0;
                elongation_new->e = elongation_new->e*scaling;
            }
        }

        tVect tangForce = FRtangentialContact(tangRelVelContact, normNormalForce, overlap, radJ, d_elements->m[e_i], elongation_new, DEM_P.sphereMat.frictionCoefObj,
                DEM_P.sphereMat.linearStiff, DEM_P.sphereMat.viscTang);

        // torque updating
        d_elements->MWall[e_i]._atomicSub(centerDistJ.cross(tangForce));
        // force updating
        d_elements->FSpringW[e_i]._atomicSub(tangForce);
        d_elements->FWall[e_i]._atomicSub(tangForce);
        d_elements->solidIntensity->_atomicAdd(tangForce.abs());
        d_objects->FParticle[o_i]._atomicAdd(tangForce);

    }
    //ROLLING

    tVect rollingMoment = rollingContact({ 0,0,0 }, wJ, radJ, normNormalForce,
            2.0 * DEM_P.sphereMat.rollingCoefPart);

    // @todo Remove redundant if(true)? there are no longer any ghosts
    if (!d_particles->isGhost[p_i]) {
        d_elements->MRolling[e_i]._atomicAdd(rollingMoment);
    }

    // save overlap
    d_elements->atomicMaxOverlap(e_i, overlap);

    // updating connectivity
    atomicInc(&d_elements->coordination[e_i], std::numeric_limits<unsigned int>::max());
}
__global__ void d_objectParticleContacts(Particle2* d_particles, Object2* d_objects, Element2* d_elements) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    // Kill excess threads early
    if (tid >= d_particles->activeCount)
        return;

    const unsigned int p_i = d_particles->activeI[tid];

    // Check each particle against every object
    // We don't (currently) have a nearObjectTable
    for (unsigned int o_i = 0; o_i < d_objects->count; ++o_i) {
        // radius
        const double rj = d_particles->r[p_i];
        // distance from object (vector)
        const tVect vectorDistance = d_particles->x0[p_i] - d_objects->x0[o_i];
        // distance from object (norm)
        const double distance = vectorDistance.norm();
        // distance before contact
        const double overlap = rj + d_objects->r[o_i] - distance;

        if (overlap > 0.0) {
            // pointer to elongation, initially pointing to an empty spring
            Elongation* elongation_here_new;
            /** @todo springs
            if (staticFrictionSolve) {
                const unsigned int indexI = d_objects->index[o_i];
                elongation_here_new = findSpring(3, indexI, p_i);
            }
            */
            d_objectParticleCollision(d_particles, d_objects, d_elements, o_i, p_i, vectorDistance, elongation_here_new);
        }
    }
}
template<>
void DEM2::objectParticleContacts<CUDA>() {
    // Launch cuda kernel to update
    int blockSize = 0;  // The launch configurator returned block size
    int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
    int gridSize = 0;  // The actual grid size needed, based on input size
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_objectParticleContacts, 0, hd_particles.activeCount);
    // Round up to accommodate required threads
    gridSize = (hd_particles.activeCount + blockSize - 1) / blockSize;
    d_objectParticleContacts<<<gridSize, blockSize>>>(d_particles, d_objects, d_elements);
    CUDA_CHECK();
}

__device__ __forceinline__ void d_cylinderParticleCollision(Particle2* d_particles, Cylinder2* d_cylinders, Element2* d_elements, const unsigned int c_i, const unsigned int p_i, const double overlap, Elongation* elongation_new) {
    // pointers to element
    const unsigned int e_i = d_particles->clusterIndex[p_i];

    // NORMAL FORCE ///////////////////////////////////////////////////////////////////

    // geometry ///////////////
    // particle radius
    const double radJ = d_particles->r[p_i];
    // vectorial distance
    const tVect vecDistance = d_cylinders->vecDist(c_i, d_particles->x0[p_i]);
    // contact point
    const tVect contactPoint = d_particles->x0[p_i] - vecDistance;
    // first local unit vector (normal)
    const tVect en = vecDistance / (radJ - overlap);
    // speed of the cylinder at contact point
    const tVect contactPointVelocity = d_cylinders->getSpeed(c_i, contactPoint);
    // relative velocity
    const tVect relVel = d_particles->x1[p_i] - contactPointVelocity;
    // relative normal velocity (modulus)
    const double normRelVel = relVel.dot(en);
    // relative normal velocity
    const tVect normalRelVel = en*normRelVel;

    // force computation /////////////////////////////////
    const double normNormalForce = normalContact(overlap, normRelVel, radJ, d_elements->m[e_i]); // was 2.0 * overlap

    // Overlap elastic potential energy
    //                    energy.elastic+=fn*xj/2.5;

    // force application ///////////////////////////////////
    // vectorial normal force
    const tVect normalForce = en*normNormalForce;

    // cluster geometry
    // vectorized radius
    const tVect vecRadJ = -radJ*en;
    // vectorized distance contactPoint-center of cluster
    const tVect centerDistJ = vecRadJ + (d_particles->x0[p_i] - d_elements->xp0[e_i]);

    // force updating
    d_elements->FWall[e_i]._atomicAdd(normalForce);
    d_elements->solidIntensity[e_i] += normalForce.abs();
    // wallI->FParticle=wallI->FParticle-fnv;
    // torque updating
    if (d_elements->size[e_i] > 1) {
        d_elements->MWall[e_i]._atomicAdd(centerDistJ.cross(normalForce));
    }

    // TANGENTIAL FORCE ///////////////////////////////////////////////////////////////////

    // angular velocities (global reference frame)
    const tVect wJ = d_elements->wpGlobal[e_i]; //2.0*quat2vec( elmtJ->qp1.multiply( elmtJ->qp0.adjoint() ) );
    // relative velocity at contact point
    const tVect relVelContact = relVel + wJ.cross(vecRadJ); // couldn't we just use elmtJ.w?
    // tangential component of relative velocity
    const tVect tangRelVelContact = relVelContact - normalRelVel;
    // norm of tangential velocity
    const double normTangRelVelContact = tangRelVelContact.norm();
    // checking if there is any tangential motion
    if (normTangRelVelContact != 0.0) {

        /** @todo springs
        // update spring
        if (DEM_P.staticFrictionSolve) {
            const double elong_old_norm = elongation_new->e.norm();

            double normElong = elongation_new->e.dot(en);
            const tVect normalElong = en*normElong;

            elongation_new->e = elongation_new->e - normalElong;

            if (elongation_new->e.norm() != 0) {
                const double scaling = elong_old_norm / elongation_new->e.norm();
                elongation_new->e = elongation_new->e*scaling;
            } else {
                const double scaling = 0.0;
                elongation_new->e = elongation_new->e*scaling;
            }
        }
        */

        tVect tangForce = FRtangentialContact(tangRelVelContact, normNormalForce, overlap, radJ, d_elements->m[e_i], elongation_new, DEM_P.sphereMat.frictionCoefWall,
            DEM_P.sphereMat.linearStiff, DEM_P.sphereMat.viscTang);

        // torque updating
        d_elements->MWall[e_i]._atomicSub(centerDistJ.cross(tangForce));
        // force updating
        d_elements->FWall[e_i]._atomicSub(tangForce);
        d_elements->solidIntensity[e_i] += tangForce.abs();
        //wallI->FParticle = wallI->FParticle+ftv;
    }
    //ROLLING

    const tVect rollingMoment = rollingContact({ 0,0,0 }, wJ, radJ, normNormalForce,
            2.0 * DEM_P.sphereMat.rollingCoefPart);


    // @todo There are not curerntly any ghosts, remove redundant if(true)?
    if (!d_particles->isGhost[p_i]) {
        d_elements->MRolling->_atomicAdd(rollingMoment);
    }

    // save overlap
    d_elements->atomicMaxOverlap(e_i, overlap);

    // updating connectivity
    atomicInc(&d_elements->coordination[e_i], std::numeric_limits<unsigned int>::max());
}
__global__ void d_cylinderParticleContacts(Particle2* d_particles, Cylinder2* d_cylinders, Element2* d_elements) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    // Kill excess threads early
    if (tid >= d_particles->activeCount)
        return;

    const unsigned int p_i = d_particles->activeI[tid];

    // Check each particle against every object
    // We don't (currently) have a nearCylinderTable
    for (unsigned int c_i = 0; c_i < d_cylinders->count; ++c_i) {
        // radius
        const double rj = d_particles->r[p_i];
        // distance from wall (norm)
        const double distance = d_cylinders->dist(c_i, d_particles->x0[p_i]);
        // distance before contact
        const double overlap = rj - distance;

        /* @todo unused?
        // distance to point 1 of axis
        const tVect p1dist = d_particles->x0[p_i] - d_cylinders->p1[c_i];
        // same but projected on the axis
        const tVect p1distax = (p1dist.dot(d_cylinders->naxes[c_i])) * d_cylinders->naxes[c_i];
        // distance of point from cylinder axis
        const tVect p1distcylinder = p1distax - p1dist;
        */

        // check for contact
        if (overlap > 0.0) {
            // pointer to elongation, initially pointing to an empty spring
            Elongation* elongation_here_new;
            /** @todo springs
            if (DEM_P.staticFrictionSolve) {
                const unsigned int indexI = d_cylinders->index[c_i];
                elongation_here_new = findSpring(2, indexI, partJ);
            }
            */
            d_cylinderParticleCollision(d_particles, d_cylinders, d_elements, c_i, p_i, overlap, elongation_here_new);
        }
    }
}
template<>
void DEM2::cylinderParticleContacts<CUDA>() {
    // Launch cuda kernel to update
    int blockSize = 0;  // The launch configurator returned block size
    int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
    int gridSize = 0;  // The actual grid size needed, based on input size
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_cylinderParticleContacts, 0, hd_particles.activeCount);
    // Round up to accommodate required threads
    gridSize = (hd_particles.activeCount + blockSize - 1) / blockSize;
    d_cylinderParticleContacts<<<gridSize, blockSize>>>(d_particles, d_cylinders, d_elements);
    CUDA_CHECK();
}
__global__ void d_computeApparentForces(Element2* d_elements) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int e_i = blockIdx.x * blockDim.x + threadIdx.x;
    // Kill excess threads early
    if (e_i >= d_elements->count)
        return;
    //local coordinates and velocity

    // Calculate the centrifugal acceleration
    if (DEM_P.solveCentrifugal) {
        // centrifugal acceleration
        d_elements->ACentrifugal[e_i] = computeCentrifugal(d_elements->xp0[e_i], DEM_P.demRotCenter, DEM_P.demRot);
    } else {
        d_elements->ACentrifugal[e_i].reset();
    }

    // Calculate Coriolis acceleration
    if (DEM_P.solveCoriolis) {
        // Get Coriolis acceleration: 
        d_elements->ACoriolis[e_i] = computeCoriolis(d_elements->xp1[e_i], DEM_P.demRot);
    } else {
        d_elements->ACoriolis[e_i].reset();
    }
}
template<>
void DEM2::computeApparentForces<CUDA>() {
    if (DEM_P.solveCentrifugal || DEM_P.solveCoriolis) {
        // Launch cuda kernel to update
        int blockSize = 0;  // The launch configurator returned block size
        int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
        int gridSize = 0;  // The actual grid size needed, based on input size
        cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_computeApparentForces, 0, hd_elements.count);
        // Round up to accommodate required threads
        gridSize = (hd_elements.count + blockSize - 1) / blockSize;
        d_computeApparentForces << <gridSize, blockSize >> > (d_elements);
        CUDA_CHECK();
    } else {
        CUDA_CALL(cudaMemset(hd_elements.ACentrifugal, 0, sizeof(tVect) * hd_elements.count));
        CUDA_CALL(cudaMemset(hd_elements.ACoriolis, 0, sizeof(tVect) * hd_elements.count));
    }
}
__global__ void d_saveObjectForces(Object2 *d_objects, tVect *d_totalForce) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int o_i = blockIdx.x * blockDim.x + threadIdx.x;
    // Kill excess threads early
    if (o_i >= d_objects->count)
        return;

    const tVect direction = Xp;

    // save maximum local force
    d_objects->updateMax(o_i, direction, DEM_P.demTime);

    // compute total force on objects
    d_totalForce->_atomicAdd(d_objects->FHydro[o_i] + d_objects->FParticle[o_i]);
}
template<>
void DEM2::saveObjectForces<CUDA>() {
    //@todo If there aren't many objects, this might be faster on the host
    // Get some memory for the reduction
    auto& ctm = CubTempMem::GetBufferSingleton();
    ctm.resize(sizeof(tVect));
    CUDA_CALL(cudaMemset(ctm.getPtr(), 0, sizeof(tVect)));
    tVect* d_totalForce = static_cast<tVect*>(ctm.getPtr());
    // Launch cuda kernel to update
    int blockSize = 0;  // The launch configurator returned block size
    int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
    int gridSize = 0;  // The actual grid size needed, based on input size
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_saveObjectForces, 0, hd_elements.count);
    // Round up to accommodate required threads
    gridSize = (hd_elements.count + blockSize - 1) / blockSize;
    d_saveObjectForces<<<gridSize, blockSize>>>(d_objects, d_totalForce);
    CUDA_CHECK();
    // Copy back the result of the reduction
    tVect totalForce;
    CUDA_CALL(cudaMemcpy(&totalForce, d_totalForce, sizeof(tVect), cudaMemcpyDeviceToHost));

    const tVect direction = Xp;
    if (DEM_P.objMaxTotalForce.dot(direction) < totalForce.dot(direction)) {
        DEM_P.objMaxTotalForce = totalForce;
        // Copy the updated parameter to device (Covered by evolveBoundaries() which syncs param before objMaxTotalForce is used)
        // CUDA_CALL(cudaMemcpyToSymbol(d_DEM_P, &h_DEM_P, sizeof(DEMParams)));
        // Save forces
        CUDA_CALL(cudaMemcpy(hd_objects.savedFParticle, hd_objects.FParticle, sizeof(tVect) * hd_objects.count, cudaMemcpyDeviceToDevice));
    }
}
__global__ void d_newtonEquationsSolution(Element2* d_elements) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    // Kill excess threads early
    if (tid >= d_elements->activeCount)
        return;

    const unsigned int e_i = d_elements->activeI[tid];

    // numerical viscosity for stability
    // see "Viscous torque on a sphere under arbitrary rotation" by Lei,  Yang, and Wu, Applied Physics Letters 89, 181908 (2006)
    const tVect FVisc = -6.0 * M_PI * DEM_P.numVisc * d_elements->radius[e_i] * d_elements->xp1[e_i];
    const tVect MVisc = -8.0 * M_PI * DEM_P.numVisc * d_elements->radius[e_i] * d_elements->radius[e_i] * d_elements->radius[e_i] * d_elements->wpGlobal[e_i];

    // translational motion
    // acceleration (sum of forces / mass + sum of accelerations)
    d_elements->x2[e_i] = (FVisc + d_elements->FHydro[e_i] + d_elements->FParticle[e_i] + d_elements->FWall[e_i]) / d_elements->m[e_i] + DEM_P.demF + d_elements->ACoriolis[e_i] + d_elements->ACentrifugal[e_i];

    // rotational motion
    // adjoint of orientation quaternion
    //const tQuat q0adj=d_elements->qp0.adjoint();
    // rotational velocity (body-fixed reference frame)
    //const tVect wBf=2.0*quat2vec( q0adj.multiply( d_elements->qp1 ) );
    // moment in global reference frame
    const tVect moment = MVisc + d_elements->MHydro[e_i] + d_elements->MParticle[e_i] + d_elements->MWall[e_i] + d_elements->MRolling[e_i];

    // moment in body-fixed reference frame
    //if (d_elements->size
    const tVect momentBf = project(moment, d_elements->qp0[e_i].adjoint());
    // rotational acceleration (body-fixed reference frame) (Newton equation for principal system)
    const tVect waBf = newtonAcc(momentBf, d_elements->I[e_i], d_elements->wpLocal[e_i]);
    // rotational acceleration (vector)
    d_elements->w1[e_i] = project(waBf, d_elements->qp0[e_i]);
    // rotational acceleration (quaternion)
    if (d_elements->size[e_i] > 1) {
        const tQuat waQuat = quatAcc(waBf, d_elements->qp1[e_i]);
        d_elements->q2[e_i] = 0.5 * d_elements->qp0[e_i].multiply(waQuat);
    }
}
template<>
void DEM2::newtonEquationsSolution<CUDA>() {
    // Launch cuda kernel to update
    int blockSize = 0;  // The launch configurator returned block size
    int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
    int gridSize = 0;  // The actual grid size needed, based on input size
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_newtonEquationsSolution, 0, hd_elements.activeCount);
    // Round up to accommodate required threads
    gridSize = (hd_elements.activeCount + blockSize - 1) / blockSize;
    d_newtonEquationsSolution<<<gridSize, blockSize>>>(d_elements);
    CUDA_CHECK();
}

template<>
void DEM2::evaluateForces<CUDA>() {
    { // Reset forces (this implementation feels grim)
        // Launch cuda kernel to update
        int blockSize = 0;  // The launch configurator returned block size
        int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
        int gridSize = 0;  // The actual grid size needed, based on input size
        const unsigned int num_threads = std::max(std::max(hd_particles.activeCount, hd_elements.activeCount),
                                                  std::max(hd_walls.count, hd_objects.count));
        cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_resetForces, 0, num_threads);
        // Round up to accommodate required threads
        gridSize = (num_threads + blockSize - 1) / blockSize;
        d_resetForces<<<gridSize, blockSize>>>(d_particles, d_elements, d_walls, d_objects);
        CUDA_CHECK();
    }
    //@todo should these contacts kernels be fused? Or ran in separate streams?
    // forces due to particle overlap and lubrication
    particleParticleContacts<CUDA>();

    // forces due to contact with plane walls and lubrication
    wallParticleContacts<CUDA>();

    // forces due to contact with stationary spheres (objects)
    objectParticleContacts<CUDA>();

    // forces due to contact with cylinders
    cylinderParticleContacts<CUDA>();


    /** @todo springs
    totSprings = 0;
    if (staticFrictionSolve) {
        // erase springs that are still inactive
        for (int a = 0; a < activeParticles.size(); ++a) {
            unsigned int p = activeParticles[a];
            for (int t = 0; t < 4; ++t) {
                for (int s = 0; s < particles[p].springs[t].size(); ++s) {
                    if (particles[p].springs[t][s].active == false) {
                        particles[p].springs[t].erase(particles[p].springs[t].begin() + s);
                        --s;
                    } else {
                        ++totSprings;
                    }
                }
            }
        }
    }
    */

    // compute apparent accelerations
    if (DEM_P.solveCentrifugal || DEM_P.solveCoriolis) {
        computeApparentForces<CUDA>();
    }

    // save info on maximum local force on the particles
    // @todo This could be in a stream parallel with apparent forces/newtonEquations
    saveObjectForces<CUDA>();
    
    //  Newton equations solution
    newtonEquationsSolution<CUDA>();
}

__global__ void d_corrector(Element2* d_elements, const std::array<double, 6> coeff1ord, const std::array<double, 6> coeff2ord) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    // Kill excess threads early
    if (i >= d_elements->activeCount)
        return;
    const unsigned int ae_i = d_elements->activeI[i];
    d_elements->correct(ae_i, coeff1ord, coeff2ord);
}
template<>
void DEM2::corrector<CUDA>() {
    //static double gear[6] = {3.0/16.0, 251.0/360.0, 1.0, 11.0/18.0, 1.0/6.0, 1.0/60.0};
//static const double c[5]={deltat, deltat*deltat/2.0, deltat*deltat*deltat/6.0, deltat*deltat*deltat*deltat/24.0, deltat*deltat*deltat*deltat*deltat/120.0};
//static const double coeff[6]={gear[0]*c[1], gear[1]*c[1]/c[0], gear[2]*c[1]/c[1], gear[3]*c[1]/c[2], gear[4]*c[1]/c[3], gear[5]*c[1]/c[4]};


// see: Computer Simulation of Liquids By M. P. Allen, D. J. Tildesle 
    constexpr std::array<double, 6> gear1ord = { 95.0 / 288.0, 1.0, 25.0 / 24.0, 35.0 / 72.0, 5.0 / 48.0, 1.0 / 120.0 }; // coefficients for a first-order equation
    constexpr std::array<double, 6> gear2ord = { 3.0 / 16.0, 251.0 / 360.0, 1.0, 11.0 / 18.0, 1.0 / 6.0, 1.0 / 60.0 }; // coefficients for a second-order equation

    static const std::array<double, 6> c = { DEM_P.deltat, DEM_P.deltat * DEM_P.deltat / 2.0,DEM_P.deltat * DEM_P.deltat * DEM_P.deltat / 6.0, DEM_P.deltat * DEM_P.deltat * DEM_P.deltat * DEM_P.deltat / 24.0, DEM_P.deltat * DEM_P.deltat * DEM_P.deltat * DEM_P.deltat * DEM_P.deltat / 120.0 };
    // Careful the c's must be changed for eq of 1st order
    static const std::array<double, 6> coeff1ord = { gear1ord[0] * c[0], gear1ord[1] * c[0] / c[0], gear1ord[2] * c[0] / c[1], gear1ord[3] * c[0] / c[2], gear1ord[4] * c[0] / c[3], gear1ord[5] * c[0] / c[4] };
    // coeff2ord has different c's sequence 
    static const std::array<double, 6> coeff2ord = { gear2ord[0] * c[1], gear2ord[1] * c[1] / c[0], gear2ord[2] * c[1] / c[1], gear2ord[3] * c[1] / c[2], gear2ord[4] * c[1] / c[3], gear2ord[5] * c[1] / c[4] };


    //    doubleList coeff;
    //    coeff.resize(6);
    //    coeff[0]=gear[0]*c[1];
    //    coeff[1]=gear[1]*c[1]/c[0];
    //    coeff[2]=gear[2]*c[1]/c[1];
    //    coeff[3]=gear[3]*c[1]/c[2];
    //    coeff[4]=gear[4]*c[1]/c[3];
    //    coeff[5]=gear[5]*c[1]/c[4];
    
    // Launch cuda kernel to update
    int blockSize = 0;  // The launch configurator returned block size
    int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
    int gridSize = 0;  // The actual grid size needed, based on input size
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_newtonEquationsSolution, 0, hd_elements.activeCount);
    // Round up to accommodate required threads
    gridSize = (hd_elements.activeCount + blockSize - 1) / blockSize;
    d_corrector<<<gridSize, blockSize>>>(d_elements, coeff1ord, coeff2ord);
    CUDA_CHECK();
}

__global__ void d_updateParticlesCorrected(Particle2 *d_particles, Element2 *d_elements) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    // Kill excess threads early
    if (i >= d_particles->activeCount)
        return;

    const unsigned int ap_i = d_particles->activeI[i];

    //getting belonging element index
    const unsigned int clusterIndex = d_particles->clusterIndex[ap_i];
    d_particles->updateCorrected(ap_i, d_elements, clusterIndex);
}
template<>
void DEM2::updateParticlesCorrected<CUDA>() {
    // this subroutine assumes that there is no change in the structure of the particles
    // any creation, deletion or reorganization of particles requires a different subroutine

    // Launch cuda kernel to update
    int blockSize = 0;  // The launch configurator returned block size
    int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
    int gridSize = 0;  // The actual grid size needed, based on input size
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_newtonEquationsSolution, 0, hd_particles.activeCount);
    // Round up to accommodate required threads
    gridSize = (hd_particles.activeCount + blockSize - 1) / blockSize;
    d_updateParticlesCorrected<<<gridSize, blockSize>>>(d_particles, d_elements);
    CUDA_CHECK();
    /** @todo Not currently any ghosts
    if (ghosts.size() != 0) {
        // updating ghost particles
        for (int g = 0; g < ghosts.size(); ++g) {
            //getting original particle index
            const unsigned int originParticleIndex = ghosts[g].ghostIndex;
            const unsigned int ghostParticleIndex = stdPartNumber + g;
            // reconstructing ghost particle
            particles[ghostParticleIndex].ghostUpdate(particles[originParticleIndex], ghosts[g].pbcVector);
        }
    }
    */
}

void DEM2::discreteElementInit(const std::array<types, 6>& externalBoundary, const std::array<double, 3>& externalSize, const std::array<tVect, 6>& externalBoundaryLocation,
                               const tVect& externalAccel, const tVect& externalRotation, const tVect& externalRotationCenter, const bool externalSolveCoriolis, const bool externalSolveCentrifugal, const double externalTimeStep) {
    // initializing DEM parameters from external parameters (LES or LBM))

    // domain size is the same of LBM. This is expressed in lattice coordinates
    DEM_P.demSize[0] = externalSize[0];
    DEM_P.demSize[1] = externalSize[1];
    DEM_P.demSize[2] = externalSize[2];
    // @todo do lattice coords need to be converted to dem coords?
    h_DEM_P.half_demSize.x = externalSize[0] / 2;
    h_DEM_P.half_demSize.y = externalSize[1] / 2;
    h_DEM_P.half_demSize.z = externalSize[2] / 2;

    // switchers for apparent accelerations
    DEM_P.solveCoriolis = externalSolveCoriolis;
    DEM_P.solveCentrifugal = externalSolveCentrifugal;

    // acceleration field
    DEM_P.demF = externalAccel;

    // roation of the reference frame
    DEM_P.demRotCenter = externalRotationCenter;
    DEM_P.demRot = externalRotation;

    // initializing particles
    const double partDensity = DEM_P.sphereMat.density;
    
    unsigned int particleCount = 0;
    for (unsigned int i = 0; i < h_elements.count; ++i) {
        // initialize element
        h_elements.initialize(i, partDensity);
        // Count required particles
        particleCount += h_elements.size[i];

    }
    // Allocate memory for particles
    h_particles.memoryAlloc<CPU>(particleCount);
    // Allocate memory for element components
    h_elements.allocComponentsData();

    for (unsigned int i = 0; i < h_elements.count; ++i) {
        h_elements.generateParticles(i, h_particles.count, h_particles);
    }

    // the number of standard particles (=not ghosts) is now fixed, and WILL NOT BE CHANGED
    stdPartNumber = h_particles.count;
    actvPartNumber = stdPartNumber;

    // initializing wall for DEM
    h_walls.initialize(externalBoundary, externalBoundaryLocation);
    // initializing cylinders for DEM
    h_cylinders.initialize();
    // initializing periodic boundary conditions for DEM
    // @todo initializePbcs(externalBoundary, externalBoundaryLocation);
    // initializing destroy planes for DEM
    //initializeDestroy(externalBoundary, externalBoundaryLocation);

    // used to get tot number of objects for CG (coarse graining)
    stdObjects = h_objects.count;

    // Copy all initialise DEM data to device
    syncElementsToDevice();
    syncWallsToDevice();
    syncObjectsToDevice();
    syncCylindersToDevice();
    syncParticlesToDevice();

    // initialize neighbor list parameters (also generate particles and set ghosts)
    if (h_elements.count) {
        // Init params for neighbours, and make them available on device
        initNeighborParameters();
        syncParamsToDevice();
        // @todo periodicObjects();
        evalNeighborTable<IMPL>();
    }


    const double totMass = std::reduce(h_elements.m, h_elements.m + h_elements.count, 0.0);
    const double minPartRadius = *std::min_element(h_elements.radius, h_elements.radius + h_elements.count);
    const double maxPartRadius = *std::max_element(h_elements.radius, h_elements.radius + h_elements.count);
    const double meanPartRadius = h_elements.count ? std::reduce(h_elements.radius, h_elements.radius + h_elements.count) / h_elements.count : 0;

    const double minObjRadius = *std::min_element(h_objects.r, h_objects.r + h_objects.count);
    const double maxObjRadius = *std::max_element(h_objects.r, h_objects.r + h_objects.count);
    const double meanObjRadius = h_objects.count ? std::reduce(h_objects.r, h_objects.r + h_objects.count) / h_objects.count : 0;

    // DEM time step
    // if multistep is 0, it should be calculated by the program here
    determineTimeStep(externalTimeStep);

    cout << "DEM parameters\n";
    cout << "domain size: xdim =" << DEM_P.demSize[0] << "; ydim= " << DEM_P.demSize[1] << "; zdim= " << DEM_P.demSize[2] << ";" << endl;
    cout << "Tot elements: " << h_elements.count << ";\t";
    cout << "Tot objects: " << h_objects.count << ";\t";
    cout << "Object radius: Mean=" << meanObjRadius << " Min=" << minObjRadius << " Max=" << maxObjRadius << ";" << endl;
    cout << "Tot standard particles: " << stdPartNumber << endl;
    cout << "Particle radius: Mean=" << meanPartRadius << " Min=" << minPartRadius << " Max=" << maxPartRadius << ";" << endl;
    cout << "Total particle mass=" << totMass << endl;
    if (DEM_P.multiStep > 0) {
        cout << "Deltat =" << DEM_P.deltat << ", therefore " << DEM_P.multiStep << " multiple steps, as a " << DEM_P.criticalRatio << ":1 ratio to estimated collision time" << endl;
    } else {
        cout << "Deltat =" << DEM_P.deltat << ", by imposing " << DEM_P.multiStep << " substeps" << endl;
    }

    switch (DEM_P.sphereMat.contactModel) {
        case LINEAR:
        {
            cout << "Contact model: linear dashpot" << endl;
            cout << "Normal stiffness = " << DEM_P.sphereMat.linearStiff << endl;
            break;
        }
        case HERTZIAN:
        {
            cout << "Contact model: damped Hertzian contact" << endl;
            break;
        }
    }
    cout << "Damping ratio = " << DEM_P.sphereMat.dampCoeff << ", equivalent to a coefficient of restitution of c=" << DEM_P.sphereMat.restitution << endl;
    cout << "Tangential viscosity = " << DEM_P.sphereMat.viscTang << endl;
    cout << "Particle-particle friction = " << DEM_P.sphereMat.frictionCoefPart
            << ", wall-particle friction = " << DEM_P.sphereMat.frictionCoefWall
            << ", object-particle friction = " << DEM_P.sphereMat.frictionCoefObj << endl;
    cout << "Rolling coefficient = " << DEM_P.sphereMat.rollingCoefPart << endl;
    cout << "Numerical viscosity =" << DEM_P.numVisc << endl;
}

void DEM2::initNeighborParameters() {
    // initializes all parameters useful for neighbor list algorithm

    cout << "Initialize neighbor list\n";

    const double radiusMultiplier = 2.5;

    // maximum radius
    const double maxRad = std::max(*std::max_element(h_elements.radius, h_elements.radius + h_elements.count),
                            *std::max_element(h_objects.r, h_objects.r + h_objects.count));
    // minimum radius
    const double minRad = *std::min_element(h_elements.radius, h_elements.radius + h_elements.count);
    
    cout << "Max radius=" << maxRad << endl;
    cout << "Min radius=" << minRad << endl;
    // if there are no particles, just to avoid nonsense numbers, we use an only cell
    // therefore the cell width is the same as the simulation box (which comes from the LB settings)
    if (h_elements.count == 0) {
        DEM_P.cellWidth[0] = DEM_P.demSize[0];
        DEM_P.cellWidth[1] = DEM_P.demSize[1];
        DEM_P.cellWidth[2] = DEM_P.demSize[2];
    }// otherwise the cell size is determined by the size of the particles
    else {
        DEM_P.cellWidth[0] = std::min(maxRad * radiusMultiplier, DEM_P.demSize[0]);
        DEM_P.cellWidth[1] = std::min(maxRad * radiusMultiplier, DEM_P.demSize[1]);
        DEM_P.cellWidth[2] = std::min(maxRad * radiusMultiplier, DEM_P.demSize[2]);
    }

    for (int k = 0; k < 3; ++k) {
        // ncells = numbers of cells for the linked cell algorithm
        DEM_P.nCells[k] = (int)floor(DEM_P.demSize[k] / DEM_P.cellWidth[k]);
        // width of the cells (actual)
        DEM_P.cellWidth[k] = DEM_P.demSize[k] / (double)DEM_P.nCells[k];
        // increase by two to give a container for ghost cells, in case of periodicity
        DEM_P.nCells[k] += 2;
    }

    // may need a revision
    DEM_P.nebrRange = std::max(maxRad * radiusMultiplier, 0.5 * std::min(DEM_P.cellWidth[0], std::min(DEM_P.cellWidth[1], DEM_P.cellWidth[2])));
    DEM_P.maxDisp = 100.0 * DEM_P.nebrRange;
    cout << "Neighbor list parameters\n";
    cout << "DEM size " << DEM_P.demSize[0] << " " << DEM_P.demSize[1] << " " << DEM_P.demSize[2] << "\n";
    cout << "Number of Cells " << DEM_P.nCells[0] << " " << DEM_P.nCells[1] << " " << DEM_P.nCells[2] << "\n";
    cout << "Cell width " << DEM_P.cellWidth[0] << " " << DEM_P.cellWidth[1] << " " << DEM_P.cellWidth[2] << "\n";
    cout << "Range " << DEM_P.nebrRange << "\n";
}
void DEM2::determineTimeStep(const double externalTimeStep) {
    // if there are particles, compute a deltat
    if (h_elements.count) {
        // if multistep is 0, it should be calculated by the program here
        if (DEM_P.multiStep == 0) {
            // find critical deltaT
            const double crit = DEM_P.criticalRatio * criticalTimeStep();
            // if the critical time is bigger than the external time step, then just use the LBM time step
            if (crit >= externalTimeStep) {
                DEM_P.multiStep = 1;
                DEM_P.deltat = externalTimeStep;
            }// if it is lower, calculate the number of substeps
            else {
                const double ratio = externalTimeStep / crit;
                ASSERT(ratio >= 1);
                DEM_P.multiStep = std::floor(ratio) + 1;
                DEM_P.deltat = externalTimeStep / (double)DEM_P.multiStep;
            }
        }// multistep can also be imposed by the user
        else {
            DEM_P.deltat = externalTimeStep / (double)DEM_P.multiStep;
        }
    }// otherwise just take the fluid time step
    else {
        DEM_P.deltat = externalTimeStep;
        DEM_P.multiStep = 1;
    }
    h_DEM_P.init_prototypeC1C2(); // Init dependant constants
}

double DEM2::criticalTimeStep() const {
    // determines the critical time step based on the stiffness and mass of the elements
    // we use YADE documentation, see https://yade-dem.org/doc/formulation.html

    // maximum damping coefficient to avoid having too large deltat
    constexpr double maxDampCoef = 0.9;

    const double minRad = *std::min_element(h_elements.radius, h_elements.radius + h_elements.count);
    const double minMass = *std::min_element(h_elements.m, h_elements.m + h_elements.count);

    // double const k=8.0/15.0*sphereMat.youngMod/(1-sphereMat.poisson*sphereMat.poisson)*sqrt(minRad);

    double deltaTCrit = 0.0;
    switch (DEM_P.sphereMat.contactModel) {
    case HERTZIAN:
    {
        // maximum length scale
        const double maxDist = std::max(DEM_P.demSize[0], std::max(DEM_P.demSize[1], DEM_P.demSize[2]));
        // modulus of acceleration
        const double maxAccel = DEM_P.demF.norm();
        // estimate maximum velocity (assuming governed by acceleration field)
        const double maxVel = std::sqrt(2.0 * maxAccel * maxDist);
        // see Landau & Lifshitz, or better Antypov & Elliott
        deltaTCrit = 2.214 * (2.0 * minRad) * std::pow(DEM_P.sphereMat.density / DEM_P.sphereMat.youngMod, 2.0 / 5.0) * std::pow(maxVel, 1.0 / 5.0);
        cout << "Computed duration of collisions: E=" << DEM_P.sphereMat.youngMod << " r=" << minRad << ", then t_c=" << deltaTCrit << endl;
        // deltaTCrit=minRad*sqrt(sphereMat.density/sphereMat.youngMod);
        break;
    }
    case LINEAR:
    {
        // compute dampCoef to use here
        double dampCoefDeltaT = 0.0;
        if (DEM_P.sphereMat.dampCoeff > maxDampCoef) {
            cout << "reducing dampCoeff (alpha_n) from " << DEM_P.sphereMat.dampCoeff << " to " << maxDampCoef << " for the purpose of deltaT computation" << endl;
            dampCoefDeltaT = maxDampCoef;
        }
        else {
            dampCoefDeltaT = DEM_P.sphereMat.dampCoeff;
        }
        deltaTCrit = M_PI / sqrt(DEM_P.sphereMat.linearStiff / minMass * (1.0 - dampCoefDeltaT * dampCoefDeltaT));
        cout << "Computed duration of collisions: k=" << DEM_P.sphereMat.linearStiff << " alpha_n=" << dampCoefDeltaT << " m=" << minMass << ", then t_c=" << deltaTCrit << endl;
        break;
    }
    }

    return deltaTCrit;
}
extern ProblemName problemName;
void DEM2::evolveBoundaries() {
#ifdef USE_CUDA
    // @todo Does anything else need syncing from device
    CUDA_CALL(cudaMemcpy(h_walls.FParticle, hd_walls.FParticle, sizeof(tVect) * h_walls.count, cudaMemcpyDeviceToHost));
#endif

    switch (problemName) {

    case TRIAXIAL:
    {
        // double frac = 1.0 / DEM_P.demInitialRepeat / ((double)DEM_P.multiStep);
        const double reduceFactor = 10.0 / (DEM_P.deltat * (double)DEM_P.multiStep);

        // wall in X
        {
            const double sideY = h_walls.p[7].dot(Yp) - h_walls.p[2].dot(Yp);
            const double sideZ = h_walls.p[8].dot(Zp) - h_walls.p[4].dot(Zp);
            const double area = sideY * sideZ;

            const double expectedForce = DEM_P.triIsopressure * area;
            const double measuredForce = h_walls.FParticle[6].dot(-1.0 * h_walls.n[6]);
            DEM_P.pressureX = measuredForce / area;
            const double deltaForce = expectedForce - measuredForce;

            const double displacement = deltaForce / DEM_P.sphereMat.linearStiff / reduceFactor;

            h_walls.p[6] = h_walls.p[6] + displacement * h_walls.n[6];
            //cout << "X: deltaForce " << deltaForce << " displacement " << displacement << " position "<< h_walls.p[6].dot(Xp)<< endl;

        }
        {
            // wall in Y
            const double sideZ = h_walls.p[8].dot(Zp) - h_walls.p[4].dot(Zp);
            const double sideX = h_walls.p[6].dot(Xp) - h_walls.p[0].dot(Xp);
            const double area = sideX * sideZ;

            const double expectedForce = DEM_P.triIsopressure * area;
            const double measuredForce = h_walls.FParticle[7].dot(-1.0 * h_walls.n[7]);
            DEM_P.pressureY = measuredForce / area;
            const double deltaForce = expectedForce - measuredForce;

            const double displacement = deltaForce / DEM_P.sphereMat.linearStiff / reduceFactor;

            h_walls.p[7] = h_walls.p[7] + displacement * h_walls.n[7];
            //cout << "Y: deltaForce " << deltaForce << " displacement " << displacement << endl;
        }
        {
            // wall in Z
            const double sideX = h_walls.p[6].dot(Xp) - h_walls.p[0].dot(Xp);
            const double sideY = h_walls.p[7].dot(Yp) - h_walls.p[2].dot(Yp);
            const double area = sideX * sideY;

            const double expectedForce = DEM_P.triIsopressure * area;
            const double measuredForce = h_walls.FParticle[8].dot(-1.0 * h_walls.n[8]);
            DEM_P.pressureZ = measuredForce / area;


            if (DEM_P.triDefSpeed == 0.0) {
                const double deltaForce = expectedForce - measuredForce;
                const double displacement = deltaForce / DEM_P.sphereMat.linearStiff / reduceFactor;
                h_walls.p[8] = h_walls.p[8] + displacement * h_walls.n[8];
            } else {
                h_walls.p[8] = h_walls.p[8] + DEM_P.triDefSpeed * DEM_P.deltat * (double)DEM_P.multiStep * h_walls.n[8];
            }
            //cout << "Z: deltaForce " << deltaForce << " displacement " << displacement << endl;
        }

        break;
    }
    }

    // update the position of the walls (rigid body)
    for (unsigned int n = 0; n < h_walls.count; ++n) {
        if (h_walls.translating[n]) {
            h_walls.p[n] = h_walls.p[n] + DEM_P.deltat * h_walls.trans[n];
        }
    }
    // update the position of the cylinders (rigid body)
    for (unsigned int c = 0; c < h_cylinders.count; c++) {
        if (h_cylinders.translating[c] == true) {
            // update p1 and p2
            h_cylinders.p1[c] = h_cylinders.p1[c] + DEM_P.deltat * h_cylinders.trans[c];
            h_cylinders.p2[c] = h_cylinders.p2[c] + DEM_P.deltat * h_cylinders.trans[c];
        }
    }
    // update the position of the object (rigid body)
    for (unsigned int o = 0; o < h_objects.count; o++) {
        // Check if the velocity vector is != 0
        const double objTrans = h_objects.x1[o].x + h_objects.x1[o].y + h_objects.x1[o].z;
        // update object position
        if (h_objects.translating[o] == true && objTrans != 0.0) {
            h_objects.x0[o] = h_objects.x0[o] + DEM_P.deltat * h_objects.x1[o];
        }
    }

#ifdef USE_CUDA
    // Sync changes back to device
    CUDA_CALL(cudaMemcpy(hd_walls.p, h_walls.p, h_walls.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_cylinders.p1, h_cylinders.p1, h_cylinders.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_cylinders.p2, h_cylinders.p2, h_cylinders.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_objects.x0, h_objects.x0, h_objects.count * sizeof(tVect), cudaMemcpyHostToDevice));
    syncParamsToDevice(); // DEM.pressureX-Z (DEM.objMaxTotalForce)
#endif
}

void DEM2::discreteElementStep() {
    // set trigger for new neighbor list
    static const double neighListTrigger = 0.25 * DEM_P.nebrRange;

    for (int demIter = 0; demIter < DEM_P.multiStep; ++demIter) {
        DEM_P.demTimeStep++;
        DEM_P.demTime += DEM_P.deltat;

        syncParamsToDevice(); // Update time on device
        // neighbor management

        evalNeighborTable<IMPL>();
        /* At this stage by rebuilding neighbour table every step
        //        cout<<"1.0 ("<<demIter<<") "<<"a="<<activeElmts[214]<<endl;
        evalMaxDisp();
        //        cout<<"1.1 ("<<demIter<<") "<<"a="<<activeElmts[214]<<endl;
        if (maxDisp > neighListTrigger) {
            maxDisp = 0.0;
            //cout<<"new neighbor list"<<endl;
            evalNeighborTable<IMPL>();
        }
        //cout<<"1.2"<<endl;
        //        cout<<"1.2 ("<<demIter<<") "<<"a="<<activeElmts[214]<<endl;
        */
        // predictor step
        predictor<IMPL>();
        //        cout<<"1.3 ("<<demIter<<") "<<"a="<<activeElmts[214]<<endl;
        // particles generation
        updateParticlesPredicted<IMPL>();

        // force evaluation
        evaluateForces<IMPL>();

        // corrector step
        corrector<IMPL>();

        // particles re-generation
        updateParticlesCorrected<IMPL>();
    }
}


Particle2& DEM2::getParticles() {
    syncParticlesFromDevice();
    return h_particles;
}
Element2& DEM2::getElements() {
    syncElementsFromDevice();
    return h_elements;
}
DEMParams& DEM2::getParams() {
    return DEM_P;
}
void DEM2::syncParamsToDevice() {
#ifdef USE_CUDA
    CUDA_CALL(cudaMemcpyToSymbol(d_DEM_P, &h_DEM_P, sizeof(DEMParams)));
#endif
}
void DEM2::syncParticlesToDevice() {
#ifdef USE_CUDA
    // note, does not currently sync neighbour index, this is created on the device
    if (!d_particles) {
        CUDA_CALL(cudaMalloc(&d_particles, sizeof(Particle2)));
    }
    // Copy particle data from h_particles to d_particles/hd_particles
    bool updateDeviceStruct = false;
    if (hd_particles.alloc < h_particles.count) {
        // Grow device buffers
        if (hd_particles.particleIndex) {
            CUDA_CALL(cudaFree(hd_particles.particleIndex));
            CUDA_CALL(cudaFree(hd_particles.clusterIndex));
            CUDA_CALL(cudaFree(hd_particles.protoIndex));
            CUDA_CALL(cudaFree(hd_particles.active));
            CUDA_CALL(cudaFree(hd_particles.isGhost));
            CUDA_CALL(cudaFree(hd_particles.r));
            CUDA_CALL(cudaFree(hd_particles.x0));
            CUDA_CALL(cudaFree(hd_particles.x1));
            CUDA_CALL(cudaFree(hd_particles.radiusVec));
            //todo springs
            CUDA_CALL(cudaFree(hd_particles.neighbour_index));
        }
        hd_particles.alloc = h_particles.count;
        CUDA_CALL(cudaMalloc(&hd_particles.particleIndex, h_particles.count * sizeof(unsigned int)));
        CUDA_CALL(cudaMalloc(&hd_particles.clusterIndex, h_particles.count * sizeof(unsigned int)));
        CUDA_CALL(cudaMalloc(&hd_particles.protoIndex, h_particles.count * sizeof(unsigned int)));
        CUDA_CALL(cudaMalloc(&hd_particles.active, h_particles.count * sizeof(bool)));
        CUDA_CALL(cudaMalloc(&hd_particles.isGhost, h_particles.count * sizeof(bool)));
        CUDA_CALL(cudaMalloc(&hd_particles.r, h_particles.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_particles.x0, h_particles.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_particles.x1, h_particles.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_particles.radiusVec, h_particles.count * sizeof(tVect)));
        //todo springs
        CUDA_CALL(cudaMalloc(&hd_particles.neighbour_index, h_particles.count * sizeof(unsigned int)));
        updateDeviceStruct = true;
    }
    hd_particles.count = h_particles.count;
    if (updateDeviceStruct) {
        // Copy updated device pointers to device
        CUDA_CALL(cudaMemcpy(d_particles, &hd_particles, sizeof(Particle2), cudaMemcpyHostToDevice));
    } else {
        // Copy updated particle count to device
        CUDA_CALL(cudaMemcpy(&d_particles->count, &hd_particles.count, sizeof(unsigned int), cudaMemcpyHostToDevice));
    }
    // Copy data to device buffers
    CUDA_CALL(cudaMemcpy(hd_particles.particleIndex, h_particles.particleIndex, h_particles.count * sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_particles.clusterIndex, h_particles.clusterIndex, h_particles.count * sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_particles.protoIndex, h_particles.protoIndex, h_particles.count * sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_particles.active, h_particles.active, h_particles.count * sizeof(bool), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_particles.r, h_particles.r, h_particles.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_particles.isGhost, h_particles.isGhost, h_particles.count * sizeof(bool), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_particles.x0, h_particles.x0, h_particles.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_particles.x1, h_particles.x1, h_particles.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_particles.radiusVec, h_particles.radiusVec, h_particles.count * sizeof(tVect), cudaMemcpyHostToDevice));
    //todo springs
#endif
}

void DEM2::syncWallsToDevice() {
#ifdef USE_CUDA
    if (!d_walls) {
        CUDA_CALL(cudaMalloc(&d_walls, sizeof(Wall2)));
    }
    // Copy latest wall data from HOST DEM to the device
    if (hd_walls.count < h_walls.count) {
        // Grow device buffers
        if (hd_walls.n) {
            CUDA_CALL(cudaFree(hd_walls.index));
            CUDA_CALL(cudaFree(hd_walls.n));
            CUDA_CALL(cudaFree(hd_walls.p));
            CUDA_CALL(cudaFree(hd_walls.FHydro));
            CUDA_CALL(cudaFree(hd_walls.FParticle));
            CUDA_CALL(cudaFree(hd_walls.rotCenter));
            CUDA_CALL(cudaFree(hd_walls.omega));
            CUDA_CALL(cudaFree(hd_walls.vel));
            CUDA_CALL(cudaFree(hd_walls.moving));
            CUDA_CALL(cudaFree(hd_walls.slip));
            CUDA_CALL(cudaFree(hd_walls.translating));
            CUDA_CALL(cudaFree(hd_walls.trans));
            CUDA_CALL(cudaFree(hd_walls.limited));
            CUDA_CALL(cudaFree(hd_walls.xMin));
            CUDA_CALL(cudaFree(hd_walls.xMax));
            CUDA_CALL(cudaFree(hd_walls.yMin));
            CUDA_CALL(cudaFree(hd_walls.yMax));
            CUDA_CALL(cudaFree(hd_walls.zMin));
            CUDA_CALL(cudaFree(hd_walls.zMax));
        }
        CUDA_CALL(cudaMalloc(&hd_walls.index, h_walls.count * sizeof(unsigned int)));
        CUDA_CALL(cudaMalloc(&hd_walls.n, h_walls.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_walls.p, h_walls.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_walls.FHydro, h_walls.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_walls.FParticle, h_walls.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_walls.rotCenter, h_walls.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_walls.omega, h_walls.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_walls.vel, h_walls.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_walls.moving, h_walls.count * sizeof(bool)));
        CUDA_CALL(cudaMalloc(&hd_walls.slip, h_walls.count * sizeof(bool)));
        CUDA_CALL(cudaMalloc(&hd_walls.translating, h_walls.count * sizeof(bool)));
        CUDA_CALL(cudaMalloc(&hd_walls.trans, h_walls.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_walls.limited, h_walls.count * sizeof(bool)));
        CUDA_CALL(cudaMalloc(&hd_walls.xMin, h_walls.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_walls.xMax, h_walls.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_walls.yMin, h_walls.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_walls.yMax, h_walls.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_walls.zMin, h_walls.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_walls.zMax, h_walls.count * sizeof(double)));
        hd_walls.count = h_walls.count;
        // Copy updated device pointers to device
        CUDA_CALL(cudaMemcpy(d_walls, &hd_walls, sizeof(Wall2), cudaMemcpyHostToDevice));
    } else if(hd_walls.count != h_walls.count) {
        // Buffer has shrunk, so just update size
        hd_walls.count = h_walls.count;
        // Copy updated particle count to device
        CUDA_CALL(cudaMemcpy(&d_walls->count, &h_walls.count, sizeof(unsigned int), cudaMemcpyHostToDevice));
    }
    // Copy data to device buffers
    CUDA_CALL(cudaMemcpy(hd_walls.index, h_walls.index, h_walls.count * sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_walls.n, h_walls.n, h_walls.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_walls.p, h_walls.p, h_walls.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_walls.FHydro, h_walls.FHydro, h_walls.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_walls.FParticle, h_walls.FParticle, h_walls.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_walls.rotCenter, h_walls.rotCenter, h_walls.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_walls.omega, h_walls.omega, h_walls.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_walls.vel, h_walls.vel, h_walls.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_walls.moving, h_walls.moving, h_walls.count * sizeof(bool), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_walls.slip, h_walls.slip, h_walls.count * sizeof(bool), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_walls.translating, h_walls.translating, h_walls.count * sizeof(bool), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_walls.trans, h_walls.trans, h_walls.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_walls.limited, h_walls.limited, h_walls.count * sizeof(bool), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_walls.xMin, h_walls.xMin, h_walls.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_walls.xMax, h_walls.xMax, h_walls.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_walls.yMin, h_walls.yMin, h_walls.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_walls.yMax, h_walls.yMax, h_walls.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_walls.zMin, h_walls.zMin, h_walls.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_walls.zMax, h_walls.zMax, h_walls.count * sizeof(double), cudaMemcpyHostToDevice));
#endif
}

void DEM2::syncCylindersToDevice() {
#ifdef USE_CUDA
    if (!d_cylinders) {
        CUDA_CALL(cudaMalloc(&d_cylinders, sizeof(Cylinder2)));
    }
    // Copy latest wall data from HOST DEM to the device
    if (hd_cylinders.count < h_cylinders.count) {
        // Grow device buffers
        if (hd_cylinders.p1) {
            CUDA_CALL(cudaFree(hd_cylinders.p1));
            CUDA_CALL(cudaFree(hd_cylinders.p2));
            CUDA_CALL(cudaFree(hd_cylinders.R));
            CUDA_CALL(cudaFree(hd_cylinders.naxes));
            CUDA_CALL(cudaFree(hd_cylinders.omega));
            CUDA_CALL(cudaFree(hd_cylinders.moving));
            CUDA_CALL(cudaFree(hd_cylinders.slip));
            CUDA_CALL(cudaFree(hd_cylinders.type));
            CUDA_CALL(cudaFree(hd_cylinders.limited));
            CUDA_CALL(cudaFree(hd_walls.xMin));
            CUDA_CALL(cudaFree(hd_walls.xMax));
            CUDA_CALL(cudaFree(hd_walls.yMin));
            CUDA_CALL(cudaFree(hd_walls.yMax));
            CUDA_CALL(cudaFree(hd_walls.zMin));
            CUDA_CALL(cudaFree(hd_walls.zMax));
        }
        CUDA_CALL(cudaMalloc(&hd_cylinders.p1, h_cylinders.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_cylinders.p2, h_cylinders.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_cylinders.R, h_cylinders.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_cylinders.naxes, h_cylinders.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_cylinders.omega, h_cylinders.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_cylinders.moving, h_cylinders.count * sizeof(bool)));
        CUDA_CALL(cudaMalloc(&hd_cylinders.slip, h_cylinders.count * sizeof(bool)));
        CUDA_CALL(cudaMalloc(&hd_cylinders.type, h_cylinders.count * sizeof(cylinderType)));
        CUDA_CALL(cudaMalloc(&hd_cylinders.limited, h_cylinders.count * sizeof(bool)));
        CUDA_CALL(cudaMalloc(&hd_cylinders.xMin, h_cylinders.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_cylinders.xMax, h_cylinders.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_cylinders.yMin, h_cylinders.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_cylinders.yMax, h_cylinders.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_cylinders.zMin, h_cylinders.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_cylinders.zMax, h_cylinders.count * sizeof(double)));
        hd_cylinders.count = h_cylinders.count;
        // Copy updated device pointers to device
        CUDA_CALL(cudaMemcpy(d_cylinders, &hd_cylinders, sizeof(Cylinder2), cudaMemcpyHostToDevice));
    } else if(hd_cylinders.count != h_cylinders.count) {
        // Buffer has shrunk, so just update size
        hd_cylinders.count = h_cylinders.count;
        // Copy updated particle count to device
        CUDA_CALL(cudaMemcpy(&d_cylinders->count, &h_cylinders.count, sizeof(unsigned int), cudaMemcpyHostToDevice));
    }
    // Copy data to device buffers
    CUDA_CALL(cudaMemcpy(hd_cylinders.p1, h_cylinders.p1, h_cylinders.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_cylinders.p2, h_cylinders.p2, h_cylinders.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_cylinders.R, h_cylinders.R, h_cylinders.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_cylinders.naxes, h_cylinders.naxes, h_cylinders.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_cylinders.omega, h_cylinders.omega, h_cylinders.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_cylinders.moving, h_cylinders.moving, h_cylinders.count * sizeof(bool), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_cylinders.slip, h_cylinders.slip, h_cylinders.count * sizeof(bool), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_cylinders.type, h_cylinders.type, h_cylinders.count * sizeof(cylinderType), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_cylinders.limited, h_cylinders.limited, h_cylinders.count * sizeof(bool), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_cylinders.xMin, h_cylinders.xMin, h_cylinders.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_cylinders.xMax, h_cylinders.xMax, h_cylinders.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_cylinders.yMin, h_cylinders.yMin, h_cylinders.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_cylinders.yMax, h_cylinders.yMax, h_cylinders.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_cylinders.zMin, h_cylinders.zMin, h_cylinders.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_cylinders.zMax, h_cylinders.zMax, h_cylinders.count * sizeof(double), cudaMemcpyHostToDevice));
#endif
}

void DEM2::syncParticlesFromDevice() {
#ifdef USE_CUDA
    // If using CUDA, data is on device by default, so sync back.
    if (hd_particles.count > h_particles.count) {
        h_particles.memoryAlloc<CPU>(hd_particles.count);
    }
    h_particles.count = hd_particles.count;
    // Copy main buffers back to host
    CUDA_CALL(cudaMemcpy(h_particles.particleIndex, hd_particles.particleIndex, hd_particles.count * sizeof(unsigned int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_particles.clusterIndex, hd_particles.clusterIndex, hd_particles.count * sizeof(unsigned int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_particles.protoIndex, hd_particles.protoIndex, hd_particles.count * sizeof(unsigned int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_particles.active, hd_particles.active, hd_particles.count * sizeof(bool), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_particles.isGhost, hd_particles.isGhost, hd_particles.count * sizeof(bool), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_particles.neighbour_index, hd_particles.neighbour_index, hd_particles.count * sizeof(unsigned int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_particles.r, hd_particles.r, hd_particles.count * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_particles.x0, hd_particles.x0, hd_particles.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_particles.x1, hd_particles.x1, hd_particles.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_particles.radiusVec, hd_particles.radiusVec, hd_particles.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    // Copy tertiary buffers back to host
    h_particles.activeCount = hd_particles.activeCount;
    if (hd_particles.activeCount > h_particles.activeAlloc) {
        if (h_particles.activeI) {
            free(h_particles.activeI);
        }
        h_particles.activeI = (unsigned int*)malloc(hd_particles.activeCount * sizeof(unsigned int));
        h_particles.activeAlloc = hd_particles.activeCount;
    }
    CUDA_CALL(cudaMemcpy(h_particles.activeI, hd_particles.activeI, hd_particles.activeCount * sizeof(unsigned int), cudaMemcpyDeviceToHost));

#endif
}
void DEM2::syncElementsToDevice() {
#ifdef USE_CUDA
    if (!d_elements) {
        CUDA_CALL(cudaMalloc(&d_elements, sizeof(Element2)));
    }
    // Copy particle data from h_elements to d_elements/hd_elements
    bool updateDeviceStruct = false;
    if (hd_elements.alloc < h_elements.count) {
        // Grow device buffers
        if (hd_elements.wSolver) {
            CUDA_CALL(cudaFree(h_elements.wSolver));
            CUDA_CALL(cudaFree(h_elements.index));
            CUDA_CALL(cudaFree(h_elements.active));
            CUDA_CALL(cudaFree(h_elements.size));
            CUDA_CALL(cudaFree(h_elements.radius));
            CUDA_CALL(cudaFree(h_elements.m));
            CUDA_CALL(cudaFree(h_elements.I));
            CUDA_CALL(cudaFree(h_elements.x0));
            CUDA_CALL(cudaFree(h_elements.x1));
            CUDA_CALL(cudaFree(h_elements.x2));
            CUDA_CALL(cudaFree(h_elements.x3));
            CUDA_CALL(cudaFree(h_elements.x4));
            CUDA_CALL(cudaFree(h_elements.x5));
            CUDA_CALL(cudaFree(h_elements.xp0));
            CUDA_CALL(cudaFree(h_elements.xp1));
            CUDA_CALL(cudaFree(h_elements.xp2));
            CUDA_CALL(cudaFree(h_elements.xp3));
            CUDA_CALL(cudaFree(h_elements.xp4));
            CUDA_CALL(cudaFree(h_elements.xp5));
            CUDA_CALL(cudaFree(h_elements.x0history));
            CUDA_CALL(cudaFree(h_elements.q0));
            CUDA_CALL(cudaFree(h_elements.q1));
            CUDA_CALL(cudaFree(h_elements.q2));
            CUDA_CALL(cudaFree(h_elements.q3));
            CUDA_CALL(cudaFree(h_elements.q4));
            CUDA_CALL(cudaFree(h_elements.q5));
            CUDA_CALL(cudaFree(h_elements.qp0));
            CUDA_CALL(cudaFree(h_elements.qp1));
            CUDA_CALL(cudaFree(h_elements.qp2));
            CUDA_CALL(cudaFree(h_elements.qp3));
            CUDA_CALL(cudaFree(h_elements.qp4));
            CUDA_CALL(cudaFree(h_elements.qp5));
            CUDA_CALL(cudaFree(h_elements.wGlobal));
            CUDA_CALL(cudaFree(h_elements.wLocal));
            CUDA_CALL(cudaFree(h_elements.wpLocal));
            CUDA_CALL(cudaFree(h_elements.wpGlobal));
            CUDA_CALL(cudaFree(h_elements.w0));
            CUDA_CALL(cudaFree(h_elements.w1));
            CUDA_CALL(cudaFree(h_elements.w2));
            CUDA_CALL(cudaFree(h_elements.w3));
            CUDA_CALL(cudaFree(h_elements.w4));
            CUDA_CALL(cudaFree(h_elements.w5));
            CUDA_CALL(cudaFree(h_elements.wp0));
            CUDA_CALL(cudaFree(h_elements.wp1));
            CUDA_CALL(cudaFree(h_elements.wp2));
            CUDA_CALL(cudaFree(h_elements.wp3));
            CUDA_CALL(cudaFree(h_elements.wp4));
            CUDA_CALL(cudaFree(h_elements.wp5));
            CUDA_CALL(cudaFree(h_elements.FHydro));
            CUDA_CALL(cudaFree(h_elements.FParticle));
            CUDA_CALL(cudaFree(h_elements.FWall));
            CUDA_CALL(cudaFree(h_elements.FGrav));
            CUDA_CALL(cudaFree(h_elements.FSpringP));
            CUDA_CALL(cudaFree(h_elements.FSpringW));
            CUDA_CALL(cudaFree(h_elements.MHydro));
            CUDA_CALL(cudaFree(h_elements.MParticle));
            CUDA_CALL(cudaFree(h_elements.MWall));
            CUDA_CALL(cudaFree(h_elements.MRolling));
            CUDA_CALL(cudaFree(h_elements.solidIntensity));
            CUDA_CALL(cudaFree(h_elements.coordination));
            CUDA_CALL(cudaFree(h_elements.ACoriolis));
            CUDA_CALL(cudaFree(h_elements.ACentrifugal));
            CUDA_CALL(cudaFree(h_elements.fluidVolume));
            CUDA_CALL(cudaFree(h_elements.maxOverlap));
            CUDA_CALL(cudaFree(h_elements.maxDtOverlap));
            CUDA_CALL(cudaFree(h_elements.slippingCase));
            CUDA_CALL(cudaFree(h_elements.componentsIndex));
            CUDA_CALL(cudaFree(h_elements.componentsData));
        }
        hd_elements.alloc = h_elements.count;
        CUDA_CALL(cudaMalloc(&hd_elements.wSolver, h_elements.count * sizeof(bool)));
        CUDA_CALL(cudaMalloc(&hd_elements.index, h_elements.count * sizeof(unsigned int)));
        CUDA_CALL(cudaMalloc(&hd_elements.active, h_elements.count * sizeof(bool)));
        CUDA_CALL(cudaMalloc(&hd_elements.size, h_elements.count * sizeof(unsigned int)));
        CUDA_CALL(cudaMalloc(&hd_elements.radius, h_elements.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_elements.m, h_elements.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_elements.I, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.x0, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.x1, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.x2, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.x3, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.x4, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.x5, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.xp0, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.xp1, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.xp2, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.xp3, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.xp4, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.xp5, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.x0history, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.q0, h_elements.count * sizeof(tQuat)));
        CUDA_CALL(cudaMalloc(&hd_elements.q1, h_elements.count * sizeof(tQuat)));
        CUDA_CALL(cudaMalloc(&hd_elements.q2, h_elements.count * sizeof(tQuat)));
        CUDA_CALL(cudaMalloc(&hd_elements.q3, h_elements.count * sizeof(tQuat)));
        CUDA_CALL(cudaMalloc(&hd_elements.q4, h_elements.count * sizeof(tQuat)));
        CUDA_CALL(cudaMalloc(&hd_elements.q5, h_elements.count * sizeof(tQuat)));
        CUDA_CALL(cudaMalloc(&hd_elements.qp0, h_elements.count * sizeof(tQuat)));
        CUDA_CALL(cudaMalloc(&hd_elements.qp1, h_elements.count * sizeof(tQuat)));
        CUDA_CALL(cudaMalloc(&hd_elements.qp2, h_elements.count * sizeof(tQuat)));
        CUDA_CALL(cudaMalloc(&hd_elements.qp3, h_elements.count * sizeof(tQuat)));
        CUDA_CALL(cudaMalloc(&hd_elements.qp4, h_elements.count * sizeof(tQuat)));
        CUDA_CALL(cudaMalloc(&hd_elements.qp5, h_elements.count * sizeof(tQuat)));
        CUDA_CALL(cudaMalloc(&hd_elements.wGlobal, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.wLocal, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.wpLocal, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.wpGlobal, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.w0, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.w1, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.w2, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.w3, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.w4, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.w5, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.wp0, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.wp1, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.wp2, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.wp3, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.wp4, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.wp5, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.FHydro, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.FParticle, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.FWall, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.FGrav, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.FSpringP, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.FSpringW, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.MHydro, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.MParticle, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.MWall, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.MRolling, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.solidIntensity, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.coordination, h_elements.count * sizeof(unsigned int)));
        CUDA_CALL(cudaMalloc(&hd_elements.ACoriolis, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.ACentrifugal, h_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.fluidVolume, h_elements.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_elements.maxOverlap, h_elements.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_elements.maxDtOverlap, h_elements.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_elements.slippingCase, h_elements.count * sizeof(int)));
        CUDA_CALL(cudaMalloc(&hd_elements.componentsIndex, (h_elements.count + 1) * sizeof(unsigned int)));
        CUDA_CALL(cudaMalloc(&hd_elements.componentsData, h_elements.componentsIndex[h_elements.count] * sizeof(unsigned int)));
        //todo springs
        updateDeviceStruct = true;
    }
    hd_elements.count = h_elements.count;
    if (updateDeviceStruct) {
        // Copy updated device pointers to device
        CUDA_CALL(cudaMemcpy(d_elements, &hd_elements, sizeof(Element2), cudaMemcpyHostToDevice));
    } else {
        // Copy updated particle count to device
        CUDA_CALL(cudaMemcpy(&d_elements->count, &hd_elements.count, sizeof(unsigned int), cudaMemcpyHostToDevice));
    }
    // Copy data to device buffers
    CUDA_CALL(cudaMemcpy(hd_elements.wSolver, h_elements.wSolver, hd_elements.count * sizeof(bool), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.index, h_elements.index, hd_elements.count * sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.active, h_elements.active, hd_elements.count * sizeof(bool), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.size, h_elements.size, hd_elements.count * sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.radius, h_elements.radius, hd_elements.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.m, h_elements.m, hd_elements.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.I, h_elements.I, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.x0, h_elements.x0, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.x1, h_elements.x1, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.x2, h_elements.x2, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.x3, h_elements.x3, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.x4, h_elements.x4, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.x5, h_elements.x5, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.xp0, h_elements.xp0, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.xp1, h_elements.xp1, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.xp2, h_elements.xp2, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.xp3, h_elements.xp3, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.xp4, h_elements.xp4, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.xp5, h_elements.xp5, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.x0history, h_elements.x0history, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.q0, h_elements.q0, hd_elements.count * sizeof(tQuat), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.q1, h_elements.q1, hd_elements.count * sizeof(tQuat), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.q2, h_elements.q2, hd_elements.count * sizeof(tQuat), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.q3, h_elements.q3, hd_elements.count * sizeof(tQuat), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.q4, h_elements.q4, hd_elements.count * sizeof(tQuat), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.q5, h_elements.q5, hd_elements.count * sizeof(tQuat), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.qp0, h_elements.qp0, hd_elements.count * sizeof(tQuat), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.qp1, h_elements.qp1, hd_elements.count * sizeof(tQuat), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.qp2, h_elements.qp2, hd_elements.count * sizeof(tQuat), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.qp3, h_elements.qp3, hd_elements.count * sizeof(tQuat), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.qp4, h_elements.qp4, hd_elements.count * sizeof(tQuat), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.qp5, h_elements.qp5, hd_elements.count * sizeof(tQuat), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.wGlobal, h_elements.wGlobal, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.wLocal, h_elements.wLocal, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.wpLocal, h_elements.wpLocal, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.wpGlobal, h_elements.wpGlobal, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.w0, h_elements.w0, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.w1, h_elements.w1, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.w2, h_elements.w2, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.w3, h_elements.w3, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.w4, h_elements.w4, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.w5, h_elements.w5, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.wp0, h_elements.wp0, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.wp1, h_elements.wp1, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.wp2, h_elements.wp2, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.wp3, h_elements.wp3, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.wp4, h_elements.wp4, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.wp5, h_elements.wp5, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.FHydro, h_elements.FHydro, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.FParticle, h_elements.FParticle, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.FWall, h_elements.FWall, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.FGrav, h_elements.FGrav, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.FSpringP, h_elements.FSpringP, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.FSpringW, h_elements.FSpringW, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.MHydro, h_elements.MHydro, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.MParticle, h_elements.MParticle, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.MWall, h_elements.MWall, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.MRolling, h_elements.MRolling, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.solidIntensity, h_elements.solidIntensity, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.coordination, h_elements.coordination, hd_elements.count * sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.ACoriolis, h_elements.ACoriolis, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.ACentrifugal, h_elements.ACentrifugal, hd_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.fluidVolume, h_elements.fluidVolume, hd_elements.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.maxOverlap, h_elements.maxOverlap, hd_elements.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.maxDtOverlap, h_elements.maxDtOverlap, hd_elements.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.slippingCase, h_elements.slippingCase, hd_elements.count * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.componentsIndex, h_elements.componentsIndex, (hd_elements.count + 1) * sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.componentsData, h_elements.componentsData, h_elements.componentsIndex[h_elements.count] * sizeof(unsigned int), cudaMemcpyHostToDevice));
    // Copy tertiary buffers to device
    h_elements.activeCount = hd_elements.activeCount;
    if (h_elements.activeCount > hd_elements.activeAlloc) {
        if (hd_elements.activeI) {
            CUDA_CALL(cudaFree(hd_elements.activeI));
        }
        CUDA_CALL(cudaMalloc(&hd_elements.activeI, hd_elements.activeCount * sizeof(unsigned int)));
        hd_elements.activeAlloc = hd_elements.activeCount;
    }
    CUDA_CALL(cudaMemcpy(hd_elements.activeI, h_elements.activeI, hd_elements.activeCount * sizeof(unsigned int), cudaMemcpyHostToDevice));
#endif
}
void DEM2::syncElementsFromDevice() {
#ifdef USE_CUDA
    // If using CUDA, data is on device by default, so sync back.
    if (hd_elements.count > h_elements.count) {
        h_elements.memoryAlloc<CPU>(hd_elements.count);
    }
    h_elements.count = hd_elements.count;
    // Copy main buffers back to host
    CUDA_CALL(cudaMemcpy(h_elements.wSolver, hd_elements.wSolver, hd_elements.count * sizeof(bool), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.index, hd_elements.index, hd_elements.count * sizeof(unsigned int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.active, hd_elements.active, hd_elements.count * sizeof(bool), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.componentsIndex, hd_elements.componentsIndex, hd_elements.count * sizeof(unsigned int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.size, hd_elements.size, hd_elements.count * sizeof(unsigned int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.radius, hd_elements.radius, hd_elements.count * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.m, hd_elements.m, hd_elements.count * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.I, hd_elements.I, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.x0, hd_elements.x0, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.x1, hd_elements.x1, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.x2, hd_elements.x2, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.x3, hd_elements.x3, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.x4, hd_elements.x4, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.x5, hd_elements.x5, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.xp0, hd_elements.xp0, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.xp1, hd_elements.xp1, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.xp2, hd_elements.xp2, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.xp3, hd_elements.xp3, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.xp4, hd_elements.xp4, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.xp5, hd_elements.xp5, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.x0history, hd_elements.x0history, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.q0, hd_elements.q0, hd_elements.count * sizeof(tQuat), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.q1, hd_elements.q1, hd_elements.count * sizeof(tQuat), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.q2, hd_elements.q2, hd_elements.count * sizeof(tQuat), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.q3, hd_elements.q3, hd_elements.count * sizeof(tQuat), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.q4, hd_elements.q4, hd_elements.count * sizeof(tQuat), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.q5, hd_elements.q5, hd_elements.count * sizeof(tQuat), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.qp0, hd_elements.qp0, hd_elements.count * sizeof(tQuat), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.qp1, hd_elements.qp1, hd_elements.count * sizeof(tQuat), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.qp2, hd_elements.qp2, hd_elements.count * sizeof(tQuat), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.qp3, hd_elements.qp3, hd_elements.count * sizeof(tQuat), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.qp4, hd_elements.qp4, hd_elements.count * sizeof(tQuat), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.qp5, hd_elements.qp5, hd_elements.count * sizeof(tQuat), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.wGlobal, hd_elements.wGlobal, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.wLocal, hd_elements.wLocal, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.wpLocal, hd_elements.wpLocal, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.wpGlobal, hd_elements.wpGlobal, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.w0, hd_elements.w0, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.w1, hd_elements.w1, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.w2, hd_elements.w2, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.w3, hd_elements.w3, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.w4, hd_elements.w4, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.w5, hd_elements.w5, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.wp0, hd_elements.wp0, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.wp1, hd_elements.wp1, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.wp2, hd_elements.wp2, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.wp3, hd_elements.wp3, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.wp4, hd_elements.wp4, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.wp5, hd_elements.wp5, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.FHydro, hd_elements.FHydro, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.FParticle, hd_elements.FParticle, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.FWall, hd_elements.FWall, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.FGrav, hd_elements.FGrav, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.FSpringP, hd_elements.FSpringP, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.FSpringW, hd_elements.FSpringW, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.MHydro, hd_elements.MHydro, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.MParticle, hd_elements.MParticle, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.MWall, hd_elements.MWall, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.MRolling, hd_elements.MRolling, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.solidIntensity, hd_elements.solidIntensity, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.coordination, hd_elements.coordination, hd_elements.count * sizeof(unsigned int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.ACoriolis, hd_elements.ACoriolis, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.ACentrifugal, hd_elements.ACentrifugal, hd_elements.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.fluidVolume, hd_elements.fluidVolume, hd_elements.count * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.maxOverlap, hd_elements.maxOverlap, hd_elements.count * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.maxDtOverlap, hd_elements.maxDtOverlap, hd_elements.count * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_elements.slippingCase, hd_elements.slippingCase, hd_elements.count * sizeof(int), cudaMemcpyDeviceToHost));
    // Copy tertiary buffers back to host
    h_elements.activeCount = hd_elements.activeCount;
    if (hd_elements.activeCount > h_elements.activeAlloc) {
        if (h_elements.activeI) {
            CUDA_CALL(cudaFree(h_elements.activeI));
        }
        CUDA_CALL(cudaMalloc(&h_elements.activeI, hd_elements.activeCount * sizeof(unsigned int)));
        h_elements.activeAlloc = hd_elements.activeCount;
    }
    CUDA_CALL(cudaMemcpy(h_elements.activeI, hd_elements.activeI, hd_elements.activeCount * sizeof(unsigned int), cudaMemcpyDeviceToHost));
#endif
}

void DEM2::syncObjectsToDevice() {
#ifdef USE_CUDA
    if (!d_objects) {
        CUDA_CALL(cudaMalloc(&d_objects, sizeof(Object2)));
    }
    // Copy latest wall data from HOST DEM to the device
    if (hd_objects.count < h_objects.count) {
        // Grow device buffers
        if (hd_objects.r) {
            CUDA_CALL(cudaFree(hd_objects.ElID));
            CUDA_CALL(cudaFree(hd_objects.r));
            CUDA_CALL(cudaFree(hd_objects.x0));
            CUDA_CALL(cudaFree(hd_objects.x1));
            CUDA_CALL(cudaFree(hd_objects.FParticle));
            CUDA_CALL(cudaFree(hd_objects.FHydro));
            CUDA_CALL(cudaFree(hd_objects.maxFParticle));
            CUDA_CALL(cudaFree(hd_objects.timeMaxFParticle));
            CUDA_CALL(cudaFree(hd_objects.savedFParticle));
            CUDA_CALL(cudaFree(hd_objects.translating));
            CUDA_CALL(cudaFree(hd_objects.trans));
        }
        CUDA_CALL(cudaMalloc(&hd_objects.ElID, h_objects.count * sizeof(unsigned int)));
        CUDA_CALL(cudaMalloc(&hd_objects.r, h_objects.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_objects.x0, h_objects.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_objects.x1, h_objects.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_objects.FParticle, h_objects.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_objects.FHydro, h_objects.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_objects.maxFParticle, h_objects.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_objects.timeMaxFParticle, h_objects.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_objects.savedFParticle, h_objects.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_objects.translating, h_objects.count * sizeof(bool)));
        CUDA_CALL(cudaMalloc(&hd_objects.trans, h_objects.count * sizeof(tVect)));
        hd_objects.count = h_objects.count;
        // Copy updated device pointers to device
        CUDA_CALL(cudaMemcpy(d_objects, &hd_objects, sizeof(Object2), cudaMemcpyHostToDevice));
    } else if(hd_objects.count != h_objects.count) {
        // Buffer has shrunk, so just update size
        hd_objects.count = h_objects.count;
        // Copy updated particle count to device
        CUDA_CALL(cudaMemcpy(&d_objects->count, &h_objects.count, sizeof(unsigned int), cudaMemcpyHostToDevice));
    }
    // Copy data to device buffers
    CUDA_CALL(cudaMemcpy(hd_objects.ElID, h_objects.ElID, h_objects.count * sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_objects.r, h_objects.r, h_objects.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_objects.x0, h_objects.x0, h_objects.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_objects.x1, h_objects.x1, h_objects.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_objects.FParticle, h_objects.FParticle, h_objects.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_objects.FHydro, h_objects.FHydro, h_objects.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_objects.maxFParticle, h_objects.maxFParticle, h_objects.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_objects.timeMaxFParticle, h_objects.timeMaxFParticle, h_objects.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_objects.savedFParticle, h_objects.savedFParticle, h_objects.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_objects.translating, h_objects.translating, h_objects.count * sizeof(bool), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_objects.trans, h_objects.trans, h_objects.count * sizeof(tVect), cudaMemcpyHostToDevice));
#endif
}

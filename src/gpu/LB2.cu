#include "LB2.h"

#include <cstdlib>

#include "cuda_helper.h"

#include "DEM.h"

/**
 * Storage for static members must be defined
 */
std::unique_ptr<CubTempMem> CubTempMem::_singletonT;
std::unique_ptr<CubTempMem> CubTempMem::_singletonB;

/**
 * (Temporary) DEM data synchronisation
 * Reformat DEM data to structure of arrays (for CPU), and copy it to device (for CUDA)
 */
template<>
bool LB2::syncElements<CPU>(const elmtList &elements) {
    bool componentsHasGrown = false;
    if (h_elements.count < elements.size()) {
        // Grow host buffers
         if (h_elements.x1) {
             free(h_elements.x1);
             free(h_elements.wGlobal);
             free(h_elements.FHydro);
             free(h_elements.MHydro);
             free(h_elements.fluidVolume);
         }
         h_elements.x1 = (tVect*)malloc(elements.size() * sizeof(tVect));
         h_elements.wGlobal = (tVect*)malloc(elements.size() * sizeof(tVect));
         h_elements.FHydro = (tVect*)malloc(elements.size() * sizeof(tVect));
         h_elements.MHydro = (tVect*)malloc(elements.size() * sizeof(tVect));
         h_elements.fluidVolume = (double*)malloc(elements.size() * sizeof(double));
    }
    // Update size
    h_elements.count = static_cast<unsigned int>(elements.size());
    // Repackage host particle data from array of structures, to structure of arrays
     for (unsigned int i = 0; i < h_elements.count; ++i) {
         h_elements.x1[i] = elements[i].x1;
         h_elements.wGlobal[i] = elements[i].wGlobal;
         h_elements.FHydro[i] = elements[i].FHydro;
         // h_elements.MHydro[i] = elements[i].MHydro; // This is zero'd before use in latticeBoltzmannStep()
         // h_elements.fluidVolume[i] = elements[i].fluidVolume; // This is zero'd before use in latticeBoltzmannStep()
     }
    // Construct the components storage
    {
        // Allocate memory for componentsData
        unsigned int totalComponents = 0;
        for (const auto& e : elements)
            totalComponents += static_cast<unsigned int>(e.components.size());
        if (!h_elements.componentsIndex || totalComponents >= h_elements.componentsIndex[elements.size()]) {
            if (h_elements.componentsData)
                free(h_elements.componentsData);
            h_elements.componentsData = (unsigned int*)malloc(totalComponents * sizeof(unsigned int));
            componentsHasGrown = true;
        }
        // Allocate componentsIndex if first pass
        if (!h_elements.componentsIndex)
            h_elements.componentsIndex = (unsigned int*)malloc((elements.size() + 1) * sizeof(unsigned int));
        // Fill componentsIndex and componentsData
        totalComponents = 0;
        for (int i = 0; i < elements.size(); ++i) {
            h_elements.componentsIndex[i] = totalComponents;
            if (!elements[i].components.empty()) {
                memcpy(h_elements.componentsData + totalComponents, elements[i].components.data(), elements[i].components.size() * sizeof(unsigned int));
                totalComponents += static_cast<unsigned int>(elements[i].components.size());
            }
        }
        h_elements.componentsIndex[elements.size()] = totalComponents;
    }
    return componentsHasGrown;
}
template<>
void LB2::syncParticles<CPU>(const particleList &particles) {
    if (h_particles.count < particles.size()) {
        // Grow host buffers
        if (h_particles.clusterIndex) {
            free(h_particles.clusterIndex);
            free(h_particles.r);
            free(h_particles.x0);
            free(h_particles.radiusVec);
        }
        h_particles.clusterIndex = (unsigned int*)malloc(particles.size() * sizeof(unsigned int));
        h_particles.r = (double*)malloc(particles.size() * sizeof(double));
        h_particles.x0 = (tVect*)malloc(particles.size() * sizeof(tVect));
        h_particles.radiusVec = (tVect*)malloc(particles.size() * sizeof(tVect));
    }
    // Update size
    h_particles.count = static_cast<unsigned int>(particles.size());
    // Repackage host particle data from array of structures, to structure of arrays
    for (int i = 0; i < h_particles.count; ++i) {
        h_particles.clusterIndex[i] = particles[i].clusterIndex;
        h_particles.r[i] = particles[i].r;
        h_particles.x0[i] = particles[i].x0;
        h_particles.radiusVec[i] = particles[i].radiusVec;
    }
}
template<>
void LB2::syncCylinders<CPU>(const cylinderList &cylinders) {
    if (h_cylinders.count < cylinders.size()) {
        // Grow host buffers
        if (h_cylinders.p1) {
            free(h_cylinders.p1);
            free(h_cylinders.p2);
            free(h_cylinders.R);
            free(h_cylinders.naxes);
            free(h_cylinders.omega);
            free(h_cylinders.moving);
        }
        h_cylinders.p1 = (tVect*)malloc(cylinders.size() * sizeof(tVect));
        h_cylinders.p2 = (tVect*)malloc(cylinders.size() * sizeof(tVect));
        h_cylinders.R = (double*)malloc(cylinders.size() * sizeof(double));
        h_cylinders.naxes = (tVect*)malloc(cylinders.size() * sizeof(tVect));
        h_cylinders.omega = (tVect*)malloc(cylinders.size() * sizeof(tVect));
        h_cylinders.moving = (bool*)malloc(cylinders.size() * sizeof(bool));
    }
    // Update size
    h_cylinders.count = static_cast<unsigned int>(cylinders.size());
    // Repackage host particle data from array of structures, to structure of arrays
    for (unsigned int i = 0; i < h_cylinders.count; ++i) {
        h_cylinders.p1[i] = cylinders[i].p1;
        h_cylinders.p2[i] = cylinders[i].p2;
        h_cylinders.R[i] = cylinders[i].R;
        h_cylinders.naxes[i] = cylinders[i].naxes;
        h_cylinders.omega[i] = cylinders[i].omega;
        h_cylinders.moving[i] = cylinders[i].moving;
    }
}
template<>
void LB2::syncWalls<CPU>(const wallList &walls) {
    if (h_walls.count < walls.size()) {
        // Grow host buffers
        if (h_walls.n) {
            free(h_walls.n);
            free(h_walls.p);
            free(h_walls.rotCenter);
            free(h_walls.omega);
            free(h_walls.vel);
            free(h_walls.FHydro);
        }
        h_walls.n = (tVect*)malloc(walls.size() * sizeof(tVect));
        h_walls.p = (tVect*)malloc(walls.size() * sizeof(tVect));
        h_walls.rotCenter = (tVect*)malloc(walls.size() * sizeof(tVect));
        h_walls.omega = (tVect*)malloc(walls.size() * sizeof(tVect));
        h_walls.vel = (tVect*)malloc(walls.size() * sizeof(tVect));
        h_walls.FHydro = (tVect*)malloc(walls.size() * sizeof(tVect));
    }
    // Update size
    h_walls.count = static_cast<unsigned int>(walls.size());
    // Repackage host particle data from array of structures, to structure of arrays
    for (unsigned int i = 0; i < h_walls.count; ++i) {
        h_walls.n[i] = walls[i].n;
        h_walls.p[i] = walls[i].p;
        h_walls.rotCenter[i] = walls[i].rotCenter;
        h_walls.omega[i] = walls[i].omega;
        h_walls.vel[i] = walls[i].vel;
        // h_walls.FHydro[i] = walls[i].FHydro; // Zero'd before use in streaming()
    }
}
template<>
void LB2::syncObjects<CPU>(const objectList &objects) {
    if (h_objects.count < objects.size()) {
        // Grow host buffers
        if (h_objects.r) {
            free(h_objects.r);
            free(h_objects.x0);
            free(h_objects.x1);
            free(h_objects.FHydro);
        }
        h_objects.r = (double*)malloc(objects.size() * sizeof(double));
        h_objects.x0 = (tVect*)malloc(objects.size() * sizeof(tVect));
        h_objects.x1 = (tVect*)malloc(objects.size() * sizeof(tVect));
        h_objects.FHydro = (tVect*)malloc(objects.size() * sizeof(tVect));
    }
    // Update size
    h_objects.count = static_cast<unsigned int>(objects.size());
    // Repackage host particle data from array of structures, to structure of arrays
    for (unsigned int i = 0; i < h_objects.count; ++i) {
        h_objects.r[i] = objects[i].r;
        h_objects.x0[i] = objects[i].x0;
        h_objects.x1[i] = objects[i].x1;
        // h_objects.FHydro[i] = objects[i].FHydro; // Zero'd before use in streaming()
    }
}
#ifdef USE_CUDA
template<>
bool LB2::syncElements<CUDA>(const elmtList &elements) {
    if (!d_elements) {
        CUDA_CALL(cudaMalloc(&d_elements, sizeof(Element2)));
    }
    // @todo copy hd_elements to d_elements
    // Copy latest particle data from HOST DEM to the device
    // @todo Can these copies be done ahead of time async?
    // @todo These copies will be redundant when DEM is moved to CUDA
    bool componentsHasGrown = this->syncElements<CPU>(elements);
    bool updateDeviceStruct = false;
    if (hd_elements.count < elements.size()) {
        if (hd_elements.x1) {
            CUDA_CALL(cudaFree(hd_elements.x1));
            CUDA_CALL(cudaFree(hd_elements.wGlobal));
            CUDA_CALL(cudaFree(hd_elements.FHydro));
            CUDA_CALL(cudaFree(hd_elements.MHydro));
            CUDA_CALL(cudaFree(hd_elements.fluidVolume));
        }
        // Initially allocate device buffers except components
        CUDA_CALL(cudaMalloc(&hd_elements.x1, hd_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.wGlobal, hd_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.FHydro, hd_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.MHydro, hd_elements.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_elements.fluidVolume, hd_elements.count * sizeof(double)));
        updateDeviceStruct = true;
    }
    if (componentsHasGrown || !hd_elements.componentsIndex) {
        // Allocate components
        if (hd_elements.componentsIndex)
            CUDA_CALL(cudaFree(hd_elements.componentsIndex));
        if (hd_elements.componentsData)
            CUDA_CALL(cudaFree(hd_elements.componentsData));
        // Allocate components
        CUDA_CALL(cudaMalloc(&hd_elements.componentsIndex, (h_elements.count + 1) * sizeof(unsigned int)));
        CUDA_CALL(cudaMalloc(&hd_elements.componentsData, h_elements.componentsIndex[h_elements.count] * sizeof(unsigned int)));
        updateDeviceStruct = true;
    }
    // Update size
    hd_elements.count = elements.size();
    if (updateDeviceStruct) {
        // Copy updated device pointers to device (@todo When/where is d_elements allocated??)
        CUDA_CALL(cudaMemcpy(d_elements, &hd_elements, sizeof(Element2), cudaMemcpyHostToDevice));
    } else {
        // Copy updated device pointers to device (@todo When/where is d_elements allocated??)
        CUDA_CALL(cudaMemcpy(&d_elements->count, &hd_elements.count, sizeof(unsigned int), cudaMemcpyHostToDevice));
    }
    // Copy data to device buffers
    CUDA_CALL(cudaMemcpy(hd_elements.x1, &h_elements.x1, h_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.wGlobal, &h_elements.wGlobal, h_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.FHydro, &h_elements.FHydro, h_elements.count * sizeof(tVect), cudaMemcpyHostToDevice));
    // CUDA_CALL(cudaMemcpy(hd_elements.MHydro, &h_elements.MHydro, h_elements.count * sizeof(tVect), cudaMemcpyHostToDevice)); // This is zero'd before use in latticeBoltzmannStep()
    // CUDA_CALL(cudaMemcpy(hd_elements.fluidVolume, &h_elements.fluidVolume, h_elements.count * sizeof(double), cudaMemcpyHostToDevice)); // This is zero'd before use in latticeBoltzmannStep()
    CUDA_CALL(cudaMemcpy(hd_elements.componentsIndex, h_elements.componentsIndex, (h_elements.count + 1) * sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.componentsData, h_elements.componentsData, h_elements.componentsIndex[h_elements.count] * sizeof(unsigned int), cudaMemcpyHostToDevice));
    return componentsHasGrown;
}
template<>
void LB2::syncParticles<CUDA>(const particleList &particles) {
    if (!d_particles) {
        CUDA_CALL(cudaMalloc(&d_particles, sizeof(Particle2)));
    }
    // Copy latest particle data from HOST DEM to the device
    // @todo Can these copies be done ahead of time async?
    // @todo These copies will be redundant when DEM is moved to CUDA
    this->syncParticles<CPU>(particles);
    if (hd_particles.count < particles.size()) {
        // Grow device buffers
        if (hd_particles.clusterIndex) {
            CUDA_CALL(cudaFree(hd_particles.clusterIndex));
            CUDA_CALL(cudaFree(hd_particles.r));
            CUDA_CALL(cudaFree(hd_particles.x0));
            CUDA_CALL(cudaFree(hd_particles.radiusVec));
        }
        CUDA_CALL(cudaMalloc(&hd_particles.clusterIndex, h_particles.count * sizeof(unsigned int)));
        CUDA_CALL(cudaMalloc(&hd_particles.r, h_particles.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_particles.x0, h_particles.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_particles.radiusVec, h_particles.count * sizeof(tVect)));
        hd_particles.count = h_particles.count;
        // Copy updated device pointers to device (@todo When/where is d_particles allocated??)
        CUDA_CALL(cudaMemcpy(d_particles, &h_particles, sizeof(Particle2), cudaMemcpyHostToDevice));
    } else if(hd_particles.count != particles.size()) {
        // Buffer has shrunk, so just update size
        hd_particles.count = static_cast<unsigned int>(particles.size());
        // Copy updated particle count to device
        CUDA_CALL(cudaMemcpy(&d_particles->count, &h_particles.count, sizeof(unsigned int), cudaMemcpyHostToDevice));
    }
    // Copy data to device buffers
    CUDA_CALL(cudaMemcpy(hd_particles.clusterIndex, h_particles.clusterIndex, h_particles.count * sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_particles.r, h_particles.r, h_particles.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_particles.x0, h_particles.x0, h_particles.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_particles.radiusVec, h_particles.radiusVec, h_particles.count * sizeof(tVect), cudaMemcpyHostToDevice));
}
template<>
void LB2::syncCylinders<CUDA>(const cylinderList &cylinders) {
    if (!d_cylinders) {
        CUDA_CALL(cudaMalloc(&d_cylinders, sizeof(Cylinder2)));
    }
    // Copy latest particle data from HOST DEM to the device
    // @todo Can these copies be done ahead of time async?
    // @todo These copies will be redundant when DEM is moved to CUDA
    this->syncCylinders<CPU>(cylinders);
    if (hd_cylinders.count < cylinders.size()) {
        // Grow device buffers
        if (hd_cylinders.p1) {
            CUDA_CALL(cudaFree(hd_cylinders.p1));
            CUDA_CALL(cudaFree(hd_cylinders.p2));
            CUDA_CALL(cudaFree(hd_cylinders.R));
            CUDA_CALL(cudaFree(hd_cylinders.naxes));
            CUDA_CALL(cudaFree(hd_cylinders.omega));
            CUDA_CALL(cudaFree(hd_cylinders.moving));
        }
        CUDA_CALL(cudaMalloc(&hd_cylinders.p1, h_cylinders.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_cylinders.p2, h_cylinders.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_cylinders.R, h_cylinders.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_cylinders.naxes, h_cylinders.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_cylinders.omega, h_cylinders.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_cylinders.moving, h_cylinders.count * sizeof(bool)));
        hd_cylinders.count = h_cylinders.count;
        // Copy updated device pointers to device
        CUDA_CALL(cudaMemcpy(d_cylinders, &hd_cylinders, sizeof(Cylinder2), cudaMemcpyHostToDevice));
    } else if(hd_cylinders.count != cylinders.size()) {
        // Buffer has shrunk, so just update size
        hd_cylinders.count = static_cast<unsigned int>(cylinders.size());
        // Copy updated particle count to device (@todo When/where is d_elements allocated??)
        CUDA_CALL(cudaMemcpy(&d_cylinders->count, &h_walls.count, sizeof(unsigned int), cudaMemcpyHostToDevice));
    }
    // Copy data to device buffers
    CUDA_CALL(cudaMemcpy(hd_cylinders.p1, h_cylinders.p1, h_cylinders.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_cylinders.p2, h_cylinders.p2, h_cylinders.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_cylinders.R, h_cylinders.R, h_cylinders.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_cylinders.naxes, h_cylinders.naxes, h_cylinders.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_cylinders.omega, h_cylinders.omega, h_cylinders.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_cylinders.moving, h_cylinders.moving, h_cylinders.count * sizeof(bool), cudaMemcpyHostToDevice));
}
template<>
void LB2::syncWalls<CUDA>(const wallList &walls) {
    if (!d_walls) {
        CUDA_CALL(cudaMalloc(&d_walls, sizeof(Wall2)));
    }
    // Copy latest particle data from HOST DEM to the device
    // @todo Can these copies be done ahead of time async?
    // @todo These copies will be redundant when DEM is moved to CUDA
    this->syncWalls<CPU>(walls);
    if (hd_walls.count < walls.size()) {
        // Grow device buffers
        if (hd_walls.n) {
            CUDA_CALL(cudaFree(hd_walls.n));
            CUDA_CALL(cudaFree(hd_walls.p));
            CUDA_CALL(cudaFree(hd_walls.rotCenter));
            CUDA_CALL(cudaFree(hd_walls.omega));
            CUDA_CALL(cudaFree(hd_walls.vel));
            CUDA_CALL(cudaFree(hd_walls.FHydro));
        }
        CUDA_CALL(cudaMalloc(&hd_walls.n, h_walls.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_walls.p, h_walls.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_walls.rotCenter, h_walls.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_walls.omega, h_walls.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_walls.vel, h_walls.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_walls.FHydro, h_walls.count * sizeof(tVect)));
        hd_walls.count = h_walls.count;
        // Copy updated device pointers to device
        CUDA_CALL(cudaMemcpy(d_walls, &hd_walls, sizeof(Wall2), cudaMemcpyHostToDevice));
    } else if(hd_walls.count != walls.size()) {
        // Buffer has shrunk, so just update size
        hd_walls.count = static_cast<unsigned int>(walls.size());
        // Copy updated particle count to device (@todo When/where is d_walls allocated??)
        CUDA_CALL(cudaMemcpy(&d_walls->count, &h_walls.count, sizeof(unsigned int), cudaMemcpyHostToDevice));
    }
    // Copy data to device buffers
    CUDA_CALL(cudaMemcpy(hd_walls.n, h_walls.n, h_walls.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_walls.p, h_walls.p, h_walls.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_walls.rotCenter, h_walls.rotCenter, h_walls.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_walls.omega, h_walls.omega, h_walls.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_walls.vel, h_walls.vel, h_walls.count * sizeof(tVect), cudaMemcpyHostToDevice));
    // CUDA_CALL(cudaMemcpy(hd_walls.FHydro, h_walls.FHydro, h_walls.count * sizeof(tVect), cudaMemcpyHostToDevice)); // Zero'd before use in streaming()
}
template<>
void LB2::syncObjects<CUDA>(const objectList &walls) {
    if (!d_objects) {
        CUDA_CALL(cudaMalloc(&d_objects, sizeof(Object2)));
    }
    // Copy latest particle data from HOST DEM to the device
    // @todo Can these copies be done ahead of time async?
    // @todo These copies will be redundant when DEM is moved to CUDA
    this->syncObjects<CPU>(walls);
    if (hd_objects.count < walls.size()) {
        // Grow device buffers
        if (hd_objects.r) {
            CUDA_CALL(cudaFree(hd_objects.r));
            CUDA_CALL(cudaFree(hd_objects.x0));
            CUDA_CALL(cudaFree(hd_objects.x1));
            CUDA_CALL(cudaFree(hd_objects.FHydro));
        }
        CUDA_CALL(cudaMalloc(&hd_objects.r, h_objects.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_objects.x0, h_objects.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_objects.x1, h_objects.count * sizeof(tVect)));
        CUDA_CALL(cudaMalloc(&hd_objects.FHydro, h_objects.count * sizeof(tVect)));
        hd_objects.count = h_objects.count;
        // Copy updated device pointers to device
        CUDA_CALL(cudaMemcpy(d_objects, &hd_objects, sizeof(Object2), cudaMemcpyHostToDevice));
    } else if(hd_objects.count != walls.size()) {
        // Buffer has shrunk, so just update size
        hd_objects.count = static_cast<unsigned int>(walls.size());
        // Copy updated particle count to device (@todo When/where is d_objects allocated??)
        CUDA_CALL(cudaMemcpy(&d_objects->count, &h_objects.count, sizeof(unsigned int), cudaMemcpyHostToDevice));
    }
    // Copy data to device buffers
    CUDA_CALL(cudaMemcpy(hd_objects.r, h_objects.r, h_objects.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_objects.x0, h_objects.x0, h_objects.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_objects.x1, h_objects.x1, h_objects.count * sizeof(tVect), cudaMemcpyHostToDevice));
    // CUDA_CALL(cudaMemcpy(hd_objects.FHydro, h_objects.FHydro, h_objects.count * sizeof(tVect), cudaMemcpyHostToDevice)); // Zero'd before use in streaming()
}
#endif


/**
 * initializeParticleBoundaries()
 */
__host__ __device__ __forceinline__ inline double common_initializeParticleBoundaries(const unsigned int i, Node2* nodes, Particle2* particles) {
    // Fetch the index of the (active) node being processed
    const unsigned int an_i = nodes->activeI[i];
    const tVect node_position = nodes->getPosition(an_i);
    for (unsigned int p_i = 0; p_i < particles->count; ++p_i) {
        const tVect convertedPosition = particles->x0[p_i] / PARAMS.unit.Length;
        // @todo pre-compute PARAMS.hydrodynamicRadius / PARAMS.unit.Length ?
        const double convertedRadius = particles->r[p_i] * PARAMS.hydrodynamicRadius / PARAMS.unit.Length;
        if (node_position.insideSphere(convertedPosition, convertedRadius)) { //-0.5?
            nodes->setInsideParticle(an_i, true);
            nodes->solidIndex[an_i] = p_i;
            return nodes->mass[an_i];  // @todo in original code it doesn't break after setting
        }
    }
    return 0.0;
}
template<>
double LB2::initializeParticleBoundaries<CPU>() {
    // Reset all nodes to outside
    memset(hd_nodes.p, 0, h_nodes.count * sizeof(bool));

    // @todo can we parallelise at a higher level?
    double totalParticleMass = 0;
#pragma omp parallel for reduction(+:totalParticleMass) 
    for (unsigned int i = 0; i < d_nodes->activeCount; ++i) {
        // Pass the active node index to the common implementation
        totalParticleMass += common_initializeParticleBoundaries(i, d_nodes, d_particles);
    }
    return totalParticleMass;
}
#ifdef USE_CUDA
__global__ void d_initializeParticleBoundaries(Node2* d_nodes, Particle2* d_particles, double *node_in_particle_mass) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    // Kill excess threads early
    if (i >= d_nodes->activeCount) return;
    // Pass the active node index to the common implementation
    const double t = common_initializeParticleBoundaries(i, d_nodes, d_particles);

    if (t != 0.0) {
        atomicAdd(node_in_particle_mass, t);
    }
}
template<>
double LB2::initializeParticleBoundaries<CUDA>() {
    // Reset all nodes to outside
    CUDA_CALL(cudaMemset(hd_nodes.p, 0, h_nodes.count * sizeof(bool)));
    // Initialise reduction variable
    auto &t = CubTempMem::GetTempSingleton();
    t.resize(sizeof(double));
    double *d_return = static_cast<double*>(t.getPtr());
    double h_return = 0;
    CUDA_CALL(cudaMemcpy(d_return, &h_return, sizeof(double), cudaMemcpyHostToDevice));

    // Launch cuda kernel to update
    // @todo Try unrolling this, so 1 thread per node+particle combination (2D launch?)
    int blockSize = 0;  // The launch configurator returned block size
    int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
    int gridSize = 0;  // The actual grid size needed, based on input size
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_initializeParticleBoundaries, 0, h_nodes.activeCount);
    // Round up to accommodate required threads
    gridSize = (h_nodes.activeCount + blockSize - 1) / blockSize;
    d_initializeParticleBoundaries << <gridSize, blockSize >> > (d_nodes, d_particles, d_return);
    CUDA_CHECK();

    // Copy back return value
    CUDA_CALL(cudaMemcpy(&h_return, d_return, sizeof(double), cudaMemcpyDeviceToHost));
    return h_return;
}
#endif

/**
 * findNewActive()
 */
__host__ __device__ __forceinline__ inline void common_findNewActive(const unsigned int i, Node2* nodes, Particle2* particles, Element2* elements) {
    // Fetch the index of the (active) node being processed
    const unsigned int an_i = nodes->activeI[i];
    if (nodes->p[an_i]) {
        const tVect nodePosition = nodes->getPosition(an_i);
        // solid index to identify cluster
        const unsigned int particleIndex = nodes->solidIndex[an_i];
        const unsigned int clusterIndex = particles->clusterIndex[particleIndex];
        // in this case check if it has been uncovered (must be out of all particles of the cluster) - we start with a true hypothesis
        // cycling through component particles
        const unsigned int first_component = elements->componentsIndex[clusterIndex];
        const unsigned int last_component = elements->componentsIndex[clusterIndex + 1];
        for (unsigned int j = first_component; j < last_component; ++j) {
            // getting indexes from particle composing the cluster
            const unsigned int componentIndex = elements->componentsData[j];
            // checking if it has been uncovered in component j of the cluster
            // radius need to be increased by half a lattice unit
            // this is because solid boundaries are located halfway between solid and fluid nodes
            const tVect convertedPosition = particles->x0[componentIndex] / PARAMS.unit.Length;
            // @todo pre-compute PARAMS.hydrodynamicRadius / PARAMS.unit.Length ?
            const double convertedRadius = particles->r[componentIndex] * PARAMS.hydrodynamicRadius / PARAMS.unit.Length;
            if (nodePosition.insideSphere(convertedPosition, convertedRadius)) { //-0.5?
                // if the node is still inside the element, the hypothesis of new active is not true anymore
                // and we can get out of the cycle
                return;
            }
        }
        // turning up the cell as we didn't exit early
        nodes->setInsideParticle(an_i, false);
    }
}
template<>
void LB2::findNewActive<CPU>() {
    // @todo can we parallelise at a higher level?
#pragma omp parallel for
    for (unsigned int i = 0; i < d_nodes->activeCount; ++i) {
        // Pass the active node index to the common implementation
        common_findNewActive(i, d_nodes, d_particles, d_elements);
    }
}
#ifdef USE_CUDA
__global__ void d_findNewActive(Node2* d_nodes, Particle2* d_particles, Element2* d_elements) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    // Kill excess threads early
    if (i >= d_nodes->activeCount) return;
    // Pass the active node index to the common implementation
    common_findNewActive(i, d_nodes, d_particles, d_elements);
}
template<>
void LB2::findNewActive<CUDA>() {
    // Launch cuda kernel to update
    int blockSize = 0;  // The launch configurator returned block size
    int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
    int gridSize = 0;  // The actual grid size needed, based on input size
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_initializeParticleBoundaries, 0, h_nodes.activeCount);
    // Round up to accommodate required threads
    gridSize = (h_nodes.activeCount + blockSize - 1) / blockSize;
    d_findNewActive << <gridSize, blockSize >> > (d_nodes, d_particles, d_elements);
    CUDA_CHECK();
}
#endif

/**
 * findNewSolid()
 */
__host__ __device__ __forceinline__ inline void common_findNewSolid(const unsigned int i, Node2* nodes, Particle2* particles, Element2* elements) {
    const unsigned int an_i = nodes->activeI[i];
    if (nodes->isInsideParticle(an_i)) {  // If node is inside particle
        // solid index to identify cluster
        const unsigned int particleIndex = nodes->solidIndex[an_i];
        const unsigned int clusterIndex = particles->clusterIndex[particleIndex];
        // cycle through first neighbors
        const unsigned int nodeCount = nodes->count;
        for (int k = 1; k < lbmMainDirec; ++k) {
            const unsigned int l_i = nodes->d[nodeCount * k + an_i];
            if (l_i != std::numeric_limits<unsigned int>::max()) {
                // checking if solid particle is close to an active one -> we have an active node to check
                if (!nodes->isInsideParticle(l_i) && nodes->isActive(l_i)) {
                    const tVect linkPosition = nodes->getPosition(l_i);
                    // check if neighbors has been covered (by any of the particles of the cluster) - we start with a false hypothesis
                    // cycling through all components of the cluster
                    const unsigned int first_component = elements->componentsIndex[clusterIndex];
                    const unsigned int last_component = elements->componentsIndex[clusterIndex + 1];
                    for (unsigned int j = first_component; j < last_component; ++j) {
                        // getting component particle index
                        const unsigned int componentIndex = elements->componentsData[j];
                        // check if it getting inside
                        // radius need to be increased by half a lattice unit
                        // this is because solid boundaries are located halfway between solid and fluid nodes
                        // @todo pre-compute PARAMS.hydrodynamicRadius / PARAMS.unit.Length ?
                        if (linkPosition.insideSphere(particles->x0[componentIndex] / PARAMS.unit.Length, particles->r[componentIndex] * PARAMS.hydrodynamicRadius / PARAMS.unit.Length)) { //-0.5?
                            // if so, then the false hypothesis does not hold true anymore
                            nodes->solidIndex[l_i] = componentIndex;
                            // By setting particle to inside, it won't be checked again, newSolidNodes hence becomes redundant
                            nodes->setInsideParticle(l_i, true);  // @todo Is this a race condition? Multiple nodes may share a link node?
                            // and we exit the cycle
                            break;
                        }
                    }
                }
            }
        }
    }
}
template<>
void LB2::findNewSolid<CPU>() {
    // @todo can we parallelise at a higher level?
#pragma omp parallel for
    for (unsigned int i = 0; i < d_nodes->activeCount; ++i) {
        // Pass the active node index to the common implementation
        common_findNewSolid(i, d_nodes, d_particles, d_elements);
    }
}

#ifdef USE_CUDA
__global__ void d_findNewSolid(Node2* d_nodes, Particle2* d_particles, Element2* d_elements) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    // Kill excess threads early
    if (i >= d_nodes->activeCount) return;
    // Pass the active node index to the common implementation
    common_findNewSolid(i, d_nodes, d_particles, d_elements);
}
template<>
void LB2::findNewSolid<CUDA>() {
    // Launch cuda kernel to update
    int blockSize = 0;  // The launch configurator returned block size
    int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
    int gridSize = 0;  // The actual grid size needed, based on input size
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_findNewSolid, 0, h_nodes.activeCount);
    // Round up to accommodate required threads
    gridSize = (h_nodes.activeCount + blockSize - 1) / blockSize;
    d_findNewSolid << <gridSize, blockSize >> > (d_nodes, d_particles, d_elements);
    CUDA_CHECK();
}
#endif

/**
 * checkNewInterfaceParticles()
 */
__host__ __device__ __forceinline__ inline void common_checkNewInterfaceParticles(const unsigned int e_i, Node2* nodes, Particle2* particles, Element2* elements) {
    // INITIAL PARTICLE POSITION ////////////////////////
    if (elements->FHydro[e_i].norm2() == 0.0) {
        const unsigned int first_component = elements->componentsIndex[e_i];
        const unsigned int last_component = elements->componentsIndex[e_i + 1];
        for (unsigned int n = first_component; n < last_component; ++n) {
            const unsigned int componentIndex = elements->componentsData[n];
            const tVect convertedPosition = particles->x0[componentIndex] / PARAMS.unit.Length;
            // @todo pre-compute PARAMS.hydrodynamicRadius / PARAMS.unit.Length ?
            const double convertedRadius = particles->r[componentIndex] * PARAMS.hydrodynamicRadius / PARAMS.unit.Length;
            for (unsigned int i_i = 0; i_i < nodes->interfaceCount; ++i_i) {
                const unsigned int nodeHere = nodes->interfaceI[i_i];
                if (!nodes->isInsideParticle(nodeHere)) {
                    // checking if node is inside a particle
                    const tVect nodePosition = nodes->getPosition(nodeHere);
                    if (nodePosition.insideSphere(convertedPosition, convertedRadius)) { //-0.5?
                        nodes->setInsideParticle(nodeHere, true);
                        nodes->solidIndex[nodeHere] = componentIndex;
                    }
                }
            }
        }
    }
}
template<>
void LB2::checkNewInterfaceParticles<CPU>() {
#pragma omp parallel for
    for (unsigned int e_i = 0; e_i < d_elements->count; ++e_i) {
        common_checkNewInterfaceParticles(e_i, d_nodes, d_particles, d_elements);
    }
}
#ifdef USE_CUDA
__global__ void d_checkNewInterfaceParticles(Node2* d_nodes, Particle2* d_particles, Element2* d_elements) {
    // Get unique CUDA thread index, which corresponds to element 
    const unsigned int e_i = blockIdx.x * blockDim.x + threadIdx.x;
    // Kill excess threads early
    if (e_i >= d_elements->count) return;
    // Pass the active node index to the common implementation
    common_checkNewInterfaceParticles(e_i, d_nodes, d_particles, d_elements);
}
template<>
void LB2::checkNewInterfaceParticles<CUDA>() {
    // Launch cuda kernel to update
    int blockSize = 0;  // The launch configurator returned block size
    int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
    int gridSize = 0;  // The actual grid size needed, based on input size
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_checkNewInterfaceParticles, 0, h_elements.count);
    // Round up to accommodate required threads
    gridSize = (h_elements.count + blockSize - 1) / blockSize;
    // @todo Are there more elements or particles? This may want to be inverted, and we can go straight to particles rather than components?
    d_checkNewInterfaceParticles << <gridSize, blockSize >> > (d_nodes, d_particles, d_elements);
    CUDA_CHECK();
}
#endif

/**
 * reconstruct()
 * computeHydroForces()
 * collision()
 */
__host__ __device__ __forceinline__ inline void common_computeHydroForces(const unsigned int an_i, Node2* nodes, Particle2* particles, Element2* elements) {
    // resetting hydrodynamic forces on nodes
    nodes->hydroForce[an_i].reset();
    if (nodes->isInsideParticle(an_i)) {
        // getting the index of the particle to compute force in the right object
        const unsigned int index = nodes->coord[an_i];
        const unsigned int particleIndex = nodes->solidIndex[an_i];
        const unsigned int clusterIndex = particles->clusterIndex[particleIndex];
        // calculating velocity of the solid boundary at the node (due to rotation of particles)
        // vectorized radius (real units)
        const tVect radius = nodes->getPosition(index) - particles->x0[particleIndex] / PARAMS.unit.Length + particles->radiusVec[particleIndex] / PARAMS.unit.Length;
        // update velocity of the particle node (u=v_center+omega x radius) (real units)
        const tVect localVel = elements->x1[clusterIndex] / PARAMS.unit.Speed + (elements->wGlobal[clusterIndex].cross(radius)) / PARAMS.unit.AngVel;

        // calculate differential velocity
        const tVect diffVel = nodes->age[an_i] * nodes->age[an_i] * nodes->liquidFraction(an_i) * (nodes->u[an_i] - localVel);

        // force on fluid
        nodes->hydroForce[an_i] += -1.0 * diffVel;

        // force on particle
#ifdef __CUDA_ARCH__
        // CUDA atomics
        atomicAdd(&elements->fluidVolume[clusterIndex], nodes->mass[an_i]);
        atomicAdd(&elements->FHydro[clusterIndex].x, 1.0 * diffVel.x);
        atomicAdd(&elements->FHydro[clusterIndex].y, 1.0 * diffVel.y);
        atomicAdd(&elements->FHydro[clusterIndex].z, 1.0 * diffVel.z);
        const tVect t = 1.0 * radius.cross(diffVel);
        atomicAdd(&elements->MHydro[clusterIndex].x, t.x);
        atomicAdd(&elements->MHydro[clusterIndex].y, t.y);
        atomicAdd(&elements->MHydro[clusterIndex].z, t.z);
#else
        // CPU atomics
#pragma omp atomic update
        elements->fluidVolume[clusterIndex] += nodes->mass[an_i];
#pragma omp atomic update
        elements->FHydro[clusterIndex] += 1.0 * diffVel;
#pragma omp atomic update
        elements->MHydro[clusterIndex] += 1.0 * radius.cross(diffVel);
#endif
    }
}
template<>
void LB2::reconstructHydroCollide<CPU>() {
    // @todo the inside of this loop could be merged with d_reconstructHydroCollide()
#pragma omp parallel for
    for (unsigned int i = 0; i < d_nodes->activeCount; ++i) {
        // Convert index to active node index
        const unsigned int an_i = d_nodes->activeI[i];

        // reconstruction of macroscopic variables from microscopic distribution
        // this step is necessary to proceed to the collision step
        d_nodes->reconstruct(an_i);

        // compute interaction forces
        if (d_elements->count) {
            common_computeHydroForces(an_i, d_nodes, d_particles, d_elements);
        }

        //collision operator
        d_nodes->collision(an_i);
    }
}
#ifdef USE_CUDA
__global__ void d_reconstructHydroCollide(Node2* d_nodes, Particle2* d_particles, Element2* d_elements) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    // Kill excess threads early
    if (i >= d_nodes->activeCount) return;

    // Convert index to active node index
    const unsigned int an_i = d_nodes->activeI[i];

    // reconstruction of macroscopic variables from microscopic distribution
    // this step is necessary to proceed to the collision step
    d_nodes->reconstruct(an_i);

    // compute interaction forces
    if (d_elements->count) {
        common_computeHydroForces(an_i, d_nodes, d_particles, d_elements);
    }

    //collision operator
    d_nodes->collision(an_i);
}
template<>
void LB2::reconstructHydroCollide<CUDA>() {
    // Launch cuda kernel to update
    int blockSize = 0;  // The launch configurator returned block size
    int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
    int gridSize = 0;  // The actual grid size needed, based on input size
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_reconstructHydroCollide, 0, h_nodes.activeCount);
    // Round up to accommodate required threads
    gridSize = (h_nodes.activeCount + blockSize - 1) / blockSize;
    d_reconstructHydroCollide << <gridSize, blockSize >> > (d_nodes, d_particles, d_elements);
    CUDA_CHECK();
}

__host__ __device__ void common_streaming(const unsigned int i, Node2* nodes, Wall2* walls) {
    // Convert index to active node index
    const unsigned int an_i = nodes->activeI[i];

    // coefficient for free-surface
    constexpr double C2x2 = 9.0;
    constexpr double C3x2 = 3.0;
    // coefficient for slip conditions
    const double S1 = PARAMS.slipCoefficient;
    const double S2 = (1.0 - PARAMS.slipCoefficient);
    // creating list for collision function @todo can this be precomputed, rather than once per node?
    std::array<double, lbmDirec> staticPres;
    for (int j = 0; j < lbmDirec; j++) {
        staticPres[j] = PARAMS.fluidMaterial.initDensity * coeff[j];
    }

    // coefficient for bounce-back
    constexpr double BBCoeff = 2.0 * 3.0;

    const unsigned int A_OFFSET = an_i * lbmDirec;
    // cycling through neighbours
    for (unsigned int j = 1; j < lbmDirec; ++j) {
        // getting neighbour index
        const unsigned int ln_i = nodes->d[nodes->count * j + an_i];
        // if neighbour is normal fluid cell what follows is true

        if (ln_i == std::numeric_limits<unsigned int>::max()) { // is gas
            // additional variables for equilibrium f computation
            const double usq = nodes->u[an_i].norm2();
            const double vuj = nodes->u[an_i].dot(v[j]);
            // streaming with constant pressure interface
            nodes->f[A_OFFSET + opp[j]] = -nodes->fs[A_OFFSET + j] + coeff[j] * PARAMS.fluidMaterial.initDensity * (2.0 + C2x2 * (vuj * vuj) - C3x2 * usq);
        } else {
            const unsigned int L_OFFSET = ln_i * lbmDirec;
            // @todo this could be greatly improved by stacking matching cases to reduce divergence
            switch (nodes->type[ln_i]) {
            case INTERFACE:
#ifdef DEBUG
            {
                // TEST USING AGE //////////////////////////////////////
                const double usq = nodes->u[an_i].norm2();
                const double vuj = nodes->u[an_i].dot(v[j]);
                nodes->f[A_OFFSET + opp[j]] = nodes->age[ln_i] * nodes->fs[L_OFFSET + opp[j]] +
                    (1.0 - nodes->age[ln_i]) * (-nodes->fs[A_OFFSET + j] + coeff[j] * PARAMS.fluidMaterial.initDensity * (2.0 + C2x2 * (vuj * vuj) - C3x2 * usq));
                break;
            }
#endif // INTERFACE falls through to LIQUID if DEBUG not defined
            case LIQUID:
            {
                nodes->f[A_OFFSET + opp[j]] = nodes->fs[L_OFFSET + opp[j]];
                break;
            }
            // for walls there is simple bounce-back
            case STAT_WALL:
            {
#ifndef DEBUG 
                if (nodes->type[an_i] == INTERFACE) {
                    // additional variables for equilibrium f computation
                    const double usq = nodes->u[an_i].norm2();
                    const double vuj = nodes->u[an_i].dot(v[j]);
                    //streaming with constant pressure interface
                    nodes->f[A_OFFSET + opp[j]] = -nodes->fs[A_OFFSET + j] + coeff[j] * PARAMS.fluidMaterial.initDensity * (2.0 + C2x2 * (vuj * vuj) - C3x2 * usq);
                    break;
                }
#endif      
                // getting the index of the wall to compute force in the right object
                const unsigned int solidIndex = nodes->solidIndex[ln_i];

                // static pressure is subtracted in order to correctly compute buoyancy for floating objects
                const tVect BBforce = nodes->bounceBackForce(an_i, j, staticPres, 0.0);
                // updating force and torque on the object (lattice units). This point is critical since many nodes update the force on the same object (lattice units)
#ifdef __CUDA_ARCH__
                    // CUDA atomics
                atomicAdd(&walls->FHydro[solidIndex].x, BBforce.x);
                atomicAdd(&walls->FHydro[solidIndex].y, BBforce.y);
                atomicAdd(&walls->FHydro[solidIndex].z, BBforce.z);
#else
                    // CPU atomics
#pragma omp atomic update
                walls->FHydro[solidIndex] += BBforce;
#endif
                // Fall through to TOPO
            }
            // for curved walls there is the rule of Mei-Luo-Shyy
            case TOPO:
            {
                nodes->f[A_OFFSET + opp[j]] = nodes->fs[A_OFFSET + j];
                break;
            }
            case OUTLET:
            {
                nodes->f[A_OFFSET + opp[j]] = std::min(nodes->fs[A_OFFSET + opp[j]], nodes->fs[A_OFFSET + j]);
                break;
            }
            // for moving walls there is simple bounce-back with velocity correction
            case DYN_WALL:
            {
                // getting the index of the wall to compute force in the right object
                const unsigned int solidIndex = nodes->solidIndex[ln_i];
                // velocity of the wall
                const tVect vel = nodes->u[ln_i];
                // variation in Bounce-Back due to moving object
                const double BBi = BBCoeff * nodes->n[an_i] * coeff[j] * vel.dot(v[j]); // mass!!!!!

                // static pressure is subtracted in order to correctly compute buoyancy for floating objects
                const tVect BBforce = nodes->bounceBackForce(an_i, j, staticPres, BBi);
                // updating force and torque on the object (lattice units). This point is critical since many nodes update the force on the same object (lattice units)
#ifdef __CUDA_ARCH__
                    // CUDA atomics
                atomicAdd(&walls->FHydro[solidIndex].x, BBforce.x);
                atomicAdd(&walls->FHydro[solidIndex].y, BBforce.y);
                atomicAdd(&walls->FHydro[solidIndex].z, BBforce.z);
#else
                    // CPU atomics
#pragma omp atomic update
                walls->FHydro[solidIndex] += BBforce;
#endif
                nodes->f[A_OFFSET + opp[j]] = nodes->fs[A_OFFSET + j] - BBi;
                // adding the extra mass to the surplus //@todo extraMass required for parity
                // extraMass = BBi * nodes->mass[an_i];  // redistributeMass() currently not used, so this isn't implemented properly
                break;
            }// for walls there is simple bounce-back
            case OBJ:
            {
                // getting the index of the wall to compute force in the right object
                const unsigned int solidIndex = nodes->solidIndex[ln_i];
                // static pressure is subtracted in order to correctly compute buoyancy for floating objects
                const tVect BBforce = nodes->bounceBackForce(an_i, j, staticPres, 0.0);
                // updating force and torque on the object (lattice units). This point is critical since many nodes update the force on the same object (lattice units)
#ifdef __CUDA_ARCH__
                    // CUDA atomics
                atomicAdd(&walls->FHydro[solidIndex].x, BBforce.x);
                atomicAdd(&walls->FHydro[solidIndex].y, BBforce.y);
                atomicAdd(&walls->FHydro[solidIndex].z, BBforce.z);
#else
                    // CPU atomics
#pragma omp atomic update
                walls->FHydro[solidIndex] += BBforce;
#endif
                nodes->f[A_OFFSET + opp[j]] = nodes->fs[A_OFFSET + j];
                break;
            }
            case SLIP_STAT_WALL:
            {
                if (j > 6) {
                    const unsigned int nodeCheck1 = nodes->d[slip1Check[j] * nodes->count + an_i];
                    const unsigned int nodeCheck2 = nodes->d[slip2Check[j] * nodes->count + an_i];
                    // check for the environment
                    const bool active1 = nodeCheck1 != std::numeric_limits<unsigned int>::max() && nodes->isActive(nodeCheck1);
                    const bool active2 = nodeCheck2 != std::numeric_limits<unsigned int>::max() && nodes->isActive(nodeCheck2);
                    // given the environment, perform the right operation
                    if (active1 && !active2) {
                        // first
                        nodes->f[A_OFFSET + opp[j]] = S1 * nodes->fs[nodeCheck1 * lbmDirec + slip1[j]] + S2 * nodes->fs[A_OFFSET + j];
                    }
                    else if (!active1 && active2) {
                        // second
                        nodes->f[A_OFFSET + opp[j]] = S1 * nodes->fs[nodeCheck2 * lbmDirec + slip2[j]] + S2 * nodes->fs[A_OFFSET + j];
                    }
                    else {
                        // standard BB
                        nodes->f[A_OFFSET + opp[j]] = nodes->fs[A_OFFSET + j];
                    }
                }
                else {
                    // standard BB
                    nodes->f[A_OFFSET + opp[j]] = nodes->fs[A_OFFSET + j];
                }
                break;
            }
            case SLIP_DYN_WALL:
            {
                // velocity of the wall
                const tVect vel = nodes->u[ln_i];
                // variation in Bounce-Back due to moving object
                const double BBi = BBCoeff * nodes->n[an_i] * coeff[j] * vel.dot(v[j]);
                if (j > 6) {
                    const unsigned int nodeCheck1 = nodes->d[slip1Check[j] * nodes->count + an_i];
                    const unsigned int nodeCheck2 = nodes->d[slip2Check[j] * nodes->count + an_i];
                    // check for the environment
                    const bool active1 = nodeCheck1 != std::numeric_limits<unsigned int>::max() && nodes->isActive(nodeCheck1);
                    const bool active2 = nodeCheck2 != std::numeric_limits<unsigned int>::max() && nodes->isActive(nodeCheck2);
                    // given the environment, perform the right operation
                    if (active1 && !active2) {
                        // first
                        nodes->f[A_OFFSET + opp[j]] = S1 * nodes->fs[nodeCheck1 * lbmDirec + slip1[j]] + S2 * (nodes->fs[A_OFFSET + j] - BBi);
                        // adding the extra mass to the surplus //@todo extraMass required for parity
                        // extraMass += S2 * nodes->mass[an_i] * BBi;  // redistributeMass() currently not used, so this isn't implemented properly
                    }
                    else if (!active1 && active2) {
                        // second
                        nodes->f[A_OFFSET + opp[j]] = S1 * nodes->fs[nodeCheck2 * lbmDirec + slip2[j]] + S2 * (nodes->fs[A_OFFSET + j] - BBi);
                        // adding the extra mass to the surplus //@todo extraMass required for parity
                        // extraMass += S2 * nodes->mass[an_i] * BBi;  // redistributeMass() currently not used, so this isn't implemented properly
                    }
                    else {
                        // standard BB
                        nodes->f[A_OFFSET + opp[j]] = nodes->fs[A_OFFSET + j] - BBi;
                        // adding the extra mass to the surplus //@todo extraMass required for parity
                        // extraMass += nodes->mass[an_i] * BBi;  // redistributeMass() currently not used, so this isn't implemented properly
                    }
                }
                else {
                    // standard BB
                    nodes->f[A_OFFSET + opp[j]] = nodes->fs[A_OFFSET + j] - BBi;
                    // adding the extra mass to the surplus //@todo extraMass required for parity
                    // extraMass += nodes->mass[an_i] * BBi;  // redistributeMass() currently not used, so this isn't implemented properly
                }
                break;
            }
            case UNUSED:
            case GAS:
            case PERIODIC:
            case CYL:
            default:
            {
                {
                    // @todo This may print out of order if multiple threads break in parallel
                    tVect pos = nodes->getPosition(an_i);
                    printf("%u(%f, %f, %f) %s TYPE ERROR:\n", an_i, pos.x, pos.y, pos.z, typeString(nodes->type[an_i]));
                    for (unsigned int k = 1; k < lbmDirec; ++k) {
                        printf("before error: j=%u link=%u\n", k, nodes->coord[nodes->d[k * nodes->count + an_i]]);
                    }
                    pos = nodes->getPosition(ln_i);
                    printf("(%f, %f, %f) %s TYPE ERROR\n", pos.x, pos.y, pos.z, typeString(nodes->type[ln_i]));
                    // @todo aborting from CUDA is harder, especially if the printf() is to be saved
#ifndef __CUDA_ARCH__
                    std::abort();
#endif
                    return;
                }
                break;

            }
            }
        }
    }
}
template<>
void LB2::streaming<CPU>() {
    // STREAMING STEP
    // Init forces to zero
    hd_walls.initForces<CPU>();
    hd_objects.initForces<CPU>();
    // Init streaming support vector
    hd_nodes.store<CPU>();

#pragma omp parallel for // @note extraMass reduction is not currently implemented
    for (unsigned int i = 0; i < d_nodes->activeCount; ++i) {
        common_streaming(i, d_nodes, d_walls);
    }

    // redistributing extra mass due to bounce back to interface cells
    // redistributeMass(extraMass);  // extraMass hasn't been implemented properly
}
#ifdef USE_CUDA
__global__ void d_streaming(Node2* d_nodes, Wall2* d_walls) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    // Kill excess threads early
    if (i >= d_nodes->activeCount) return;

    common_streaming(i, d_nodes, d_walls);
}
template<>
void LB2::streaming<CUDA>() {
    // STREAMING STEP
    // Init forces to zero
    hd_walls.initForces<CUDA>();
    hd_objects.initForces<CUDA>();
    // Init streaming support vector
    hd_nodes.store<CUDA>();

    // Launch cuda kernel to update
    int blockSize = 0;  // The launch configurator returned block size
    int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
    int gridSize = 0;  // The actual grid size needed, based on input size
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_initializeParticleBoundaries, 0, h_nodes.activeCount);
    // Round up to accommodate required threads
    gridSize = (h_nodes.activeCount + blockSize - 1) / blockSize;
    d_streaming << <gridSize, blockSize >> > (d_nodes, d_walls);
    CUDA_CHECK();

#ifdef _DEBUG
    CUDA_CALL(cudaMemcpy(h_nodes.f, hd_nodes.f, sizeof(double) * h_nodes.count * lbmDirec, cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_nodes.activeI, hd_nodes.activeI, sizeof(unsigned int) * h_nodes.activeCount, cudaMemcpyDeviceToHost));
    for (unsigned int in = 0; in < h_nodes.activeCount; ++in) {
        const unsigned int a_i = h_nodes.activeI[in];
        for (unsigned int j = 1; j < lbmDirec; ++j) {
            if (h_nodes.f[a_i * lbmDirec + j] == 0) {
                cout << "Error!" << endl;
            }
        }
    }
#endif

    // redistributing extra mass due to bounce back to interface cells
    // redistributeMass(extraMass);  // extraMass hasn't been implemented properly
}
#endif
template<>
void LB2::shiftToPhysical<CPU>() {
    for (unsigned int i = 0; i < d_elements->count; ++i) {
        d_elements->FHydro[i] *= PARAMS.unit.Force;
        d_elements->MHydro[i] *= PARAMS.unit.Torque;
        d_elements->fluidVolume[i] *= PARAMS.unit.Volume;
    }
    for (unsigned int i = 0; i < d_walls->count; ++i) {
        d_walls->FHydro[i] *= PARAMS.unit.Force;
    }
    for (unsigned int i = 0; i < d_objects->count; ++i) {
        d_objects->FHydro[i] *= PARAMS.unit.Force;
    }
}
#ifdef USE_CUDA
__global__ void d_shiftToPhysical(Element2* d_elements, Wall2* d_walls, Object2* d_objects) {
    // Get unique CUDA thread index
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < d_elements->count) {
        d_elements->FHydro[i] *= PARAMS.unit.Force;
        d_elements->MHydro[i] *= PARAMS.unit.Torque;
        d_elements->fluidVolume[i] *= PARAMS.unit.Volume;
    }
    if (i < d_walls->count) {
        d_walls->FHydro[i] *= PARAMS.unit.Force;
    }
    if (i < d_objects->count) {
        d_objects->FHydro[i] *= PARAMS.unit.Force;
    }
}
template<>
void LB2::shiftToPhysical<CUDA>() {
    // Launch enough threads to accomodate everything
    const unsigned int maxCount = std::max(std::max(h_elements.count, h_walls.count), h_objects.count);
    // Launch cuda kernel to update
    int blockSize = 0;  // The launch configurator returned block size
    int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
    int gridSize = 0;  // The actual grid size needed, based on input size
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_initializeParticleBoundaries, 0, maxCount);
    // Round up to accommodate required threads
    gridSize = (maxCount + blockSize - 1) / blockSize;
    d_shiftToPhysical << <gridSize, blockSize >> > (d_elements, d_walls, d_objects);
    CUDA_CHECK();
}
#endif


void LB2::latticeBoltzmannCouplingStep(bool& newNeighbourList) {
    // identifies which nodes need to have an update due to particle movement
    // the complexity arises from trying to keep the scaling as close to linear as possible
    // maybe the best idea is to do this in two step:
    // 1) the first is to check for new active nodes and initialise them
    // 2) the second is checking for new solid nodes.
    // this automatically check also for the hideous case of particle to particle double transition

    /**
     * @todo The parallelisation of each of these methods should be reviewed
     *       Most are 2D loops, the range of each being unclear
     *       Likewise, can OpenMP parallel block be moved outside of each member?
     */

     // first we check if a new neighbour table has been defined. In that case, the indexing needs to be reinitialised
    if (newNeighbourList) {
        cout << endl << "New neighbour list" << endl;
        this->initializeParticleBoundaries<IMPL>();
        newNeighbourList = false;
    }
    else {
        // SOLID TO ACTIVE CHECK
        // @note Calling this directly after initializeParticleBoundaries() is redundant, hence else
        this->findNewActive<IMPL>();
    }

    // ACTIVE TO SOLID CHECK
    this->findNewSolid<IMPL>();

    if (PARAMS.freeSurface) {
        this->checkNewInterfaceParticles<IMPL>();
    }
}
void LB2::latticeBoltzmannStep() {
    // Reconstruct active list
    hd_nodes.cleanLists<IMPL>();

    // Initializing the elements forces (lattice units)
    hd_elements.initElements<IMPL>();

    // Initialise lattice boltzmann force vector
    if (!h_PARAMS.forceField) {
        h_PARAMS.lbF.reset();
        syncParams();
    }

    // reconstruct(), computeHydroForces(), collision()
    // Reconstruct macroscopic variables from microscopic distribution
    // Compute interaction forces with DEM elmts
    // Collision step
    this->reconstructHydroCollide<IMPL>();

    // Streaming operator
    this->streaming<IMPL>();

    // Shift element/wall/object forces and torques to physical units
    this->shiftToPhysical<IMPL>();
}
#endif

Node2& LB2::getNodes() {
#ifdef USE_CUDA
    // If using CUDA, data is on device by default, so sync back.
    if (hd_nodes.count > h_nodes.count) {
        // Resize main buffers
        if (h_nodes.f) free(h_nodes.f);
        h_nodes.f = static_cast<double*>(malloc(hd_nodes.count * lbmDirec * sizeof(double)));
        if (h_nodes.fs) free(h_nodes.fs);
        h_nodes.fs = static_cast<double*>(malloc(hd_nodes.count * lbmDirec * sizeof(double)));
        if (h_nodes.n) free(h_nodes.n);
        h_nodes.n = static_cast<double*>(malloc(hd_nodes.count * sizeof(double)));
        if (h_nodes.u) free(h_nodes.u);
        h_nodes.u = static_cast<tVect*>(malloc(hd_nodes.count * sizeof(tVect)));
        if (h_nodes.hydroForce) free(h_nodes.hydroForce);
        h_nodes.hydroForce = static_cast<tVect*>(malloc(hd_nodes.count * sizeof(tVect)));
        if (h_nodes.centrifugalForce) free(h_nodes.centrifugalForce);
        h_nodes.centrifugalForce = static_cast<tVect*>(malloc(hd_nodes.count * sizeof(tVect)));
        if (h_nodes.mass) free(h_nodes.mass);
        h_nodes.mass = static_cast<double*>(malloc(hd_nodes.count * sizeof(double)));
        if (h_nodes.visc) free(h_nodes.visc);
        h_nodes.visc = static_cast<double*>(malloc(hd_nodes.count * sizeof(double)));
        if (h_nodes.basal) free(h_nodes.basal);
        h_nodes.basal = static_cast<bool*>(malloc(hd_nodes.count * sizeof(bool)));
        if (h_nodes.friction) free(h_nodes.friction);
        h_nodes.friction = static_cast<double*>(malloc(hd_nodes.count * sizeof(double)));
        if (h_nodes.age) free(h_nodes.age);
        h_nodes.age = static_cast<float*>(malloc(hd_nodes.count * sizeof(float)));
        if (h_nodes.solidIndex) free(h_nodes.solidIndex);
        h_nodes.solidIndex = static_cast<unsigned int*>(malloc(hd_nodes.count * sizeof(unsigned int)));
        if (h_nodes.d) free(h_nodes.d);
        h_nodes.d = static_cast<unsigned int*>(malloc(hd_nodes.count * lbmDirec * sizeof(unsigned int)));
        if (h_nodes.curved) free(h_nodes.curved);
        h_nodes.curved = static_cast<unsigned int*>(malloc(hd_nodes.count * sizeof(unsigned int)));
        if (h_nodes.type) free(h_nodes.type);
        h_nodes.type = static_cast<types*>(malloc(hd_nodes.count * sizeof(types)));
        if (h_nodes.p) free(h_nodes.p);
        h_nodes.p = static_cast<bool*>(malloc(hd_nodes.count * sizeof(bool)));
    }
    h_nodes.count = hd_nodes.count;
    // Resize misc buffers
    if (hd_nodes.activeCount > h_nodes.activeAlloc) {
        if (h_nodes.activeI) free(h_nodes.activeI);
        h_nodes.activeI = static_cast<unsigned int*>(malloc(hd_nodes.activeCount * sizeof(unsigned int)));
    }
    h_nodes.activeCount = hd_nodes.activeCount;
    if (hd_nodes.interfaceCount > h_nodes.interfaceCount) {
        if (h_nodes.interfaceI) free(h_nodes.interfaceI);
        h_nodes.interfaceI = static_cast<unsigned int*>(malloc(hd_nodes.interfaceCount * sizeof(unsigned int)));
    }
    h_nodes.interfaceCount = hd_nodes.interfaceCount;
    if (hd_nodes.fluidCount > h_nodes.fluidCount) {
        if (h_nodes.fluidI) free(h_nodes.fluidI);
        h_nodes.fluidI = static_cast<unsigned int*>(malloc(hd_nodes.fluidCount * sizeof(unsigned int)));
    }
    h_nodes.fluidCount = hd_nodes.fluidCount;
    if (hd_nodes.wallCount > h_nodes.wallCount) {
        if (h_nodes.wallI) free(h_nodes.wallI);
        h_nodes.wallI = static_cast<unsigned int*>(malloc(hd_nodes.wallCount * sizeof(unsigned int)));
    }
    h_nodes.wallCount = hd_nodes.wallCount;
    if (hd_nodes.curveCount > h_nodes.curveCount) {
        if (h_nodes.curves) free(h_nodes.curves);
        h_nodes.curves = static_cast<curve*>(malloc(hd_nodes.curveCount * sizeof(curve)));
    }
    h_nodes.curveCount = hd_nodes.curveCount;
    // Copy main buffers back to host
    CUDA_CALL(cudaMemcpy(h_nodes.coord, hd_nodes.coord, h_nodes.count * sizeof(unsigned int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_nodes.f, hd_nodes.f, h_nodes.count * lbmDirec * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_nodes.fs, hd_nodes.fs, h_nodes.count * lbmDirec * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_nodes.n, hd_nodes.n, h_nodes.count * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_nodes.u, hd_nodes.u, h_nodes.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_nodes.hydroForce, hd_nodes.hydroForce, h_nodes.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_nodes.centrifugalForce, hd_nodes.centrifugalForce, h_nodes.count * sizeof(tVect), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_nodes.mass, hd_nodes.mass, h_nodes.count * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_nodes.visc, hd_nodes.visc, h_nodes.count * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_nodes.basal, hd_nodes.basal, h_nodes.count * sizeof(bool), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_nodes.friction, hd_nodes.friction, h_nodes.count * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_nodes.age, hd_nodes.age, h_nodes.count * sizeof(float), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_nodes.solidIndex, hd_nodes.solidIndex, h_nodes.count * sizeof(unsigned int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_nodes.d, hd_nodes.d, h_nodes.count * lbmDirec * sizeof(unsigned int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_nodes.curved, hd_nodes.curved, h_nodes.count * sizeof(unsigned int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_nodes.type, hd_nodes.type, h_nodes.count * sizeof(types), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_nodes.p, hd_nodes.p, h_nodes.count * sizeof(bool), cudaMemcpyDeviceToHost));
    // Copy misc buffers back to host
    CUDA_CALL(cudaMemcpy(h_nodes.activeI, hd_nodes.activeI, h_nodes.activeCount * sizeof(unsigned int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_nodes.interfaceI, hd_nodes.interfaceI, h_nodes.interfaceCount * sizeof(unsigned int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_nodes.fluidI, hd_nodes.fluidI, h_nodes.fluidCount * sizeof(unsigned int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_nodes.wallI, hd_nodes.wallI, h_nodes.wallCount * sizeof(unsigned int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_nodes.curves, hd_nodes.curves, h_nodes.curveCount * sizeof(curve), cudaMemcpyDeviceToHost));
#endif
    return h_nodes;
}
void LB2::initDeviceNodes() {
#ifdef USE_CUDA
    // Allocate the main storage
    if (d_nodes) {
        fprintf(stderr, "LB2::initDeviceNodes() should only be called once.");
        throw std::exception();
    }
    CUDA_CALL(cudaMalloc(&d_nodes, sizeof(Node2)));
    // Build HD struct
    hd_nodes.activeCount = h_nodes.activeCount;
    hd_nodes.activeAlloc = h_nodes.activeCount;
    CUDA_CALL(cudaMalloc(&hd_nodes.activeI, hd_nodes.activeCount * sizeof(unsigned int)));
    CUDA_CALL(cudaMemcpy(hd_nodes.activeI, h_nodes.activeI, hd_nodes.activeCount * sizeof(unsigned int), cudaMemcpyHostToDevice));
    hd_nodes.interfaceCount = h_nodes.interfaceCount;
    CUDA_CALL(cudaMalloc(&hd_nodes.interfaceI, hd_nodes.interfaceCount * sizeof(unsigned int)));
    CUDA_CALL(cudaMemcpy(hd_nodes.interfaceI, h_nodes.interfaceI, hd_nodes.interfaceCount * sizeof(unsigned int), cudaMemcpyHostToDevice));
    hd_nodes.fluidCount = h_nodes.fluidCount;
    CUDA_CALL(cudaMalloc(&hd_nodes.fluidI, hd_nodes.fluidCount * sizeof(unsigned int)));
    CUDA_CALL(cudaMemcpy(hd_nodes.fluidI, h_nodes.fluidI, hd_nodes.fluidCount * sizeof(unsigned int), cudaMemcpyHostToDevice));
    hd_nodes.wallCount = h_nodes.wallCount;
    CUDA_CALL(cudaMalloc(&hd_nodes.wallI, hd_nodes.wallCount * sizeof(unsigned int)));
    CUDA_CALL(cudaMemcpy(hd_nodes.wallI, h_nodes.wallI, hd_nodes.wallCount * sizeof(unsigned int), cudaMemcpyHostToDevice));
    hd_nodes.curveCount = h_nodes.curveCount;
    CUDA_CALL(cudaMalloc(&hd_nodes.curves, hd_nodes.curveCount * sizeof(curve)));
    CUDA_CALL(cudaMemcpy(hd_nodes.curves, h_nodes.curves, hd_nodes.curveCount * sizeof(curve), cudaMemcpyHostToDevice));
    hd_nodes.count = h_nodes.count;
    CUDA_CALL(cudaMalloc(&hd_nodes.coord, hd_nodes.count * sizeof(unsigned int)));
    CUDA_CALL(cudaMemcpy(hd_nodes.coord, h_nodes.coord, hd_nodes.count * sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMalloc(&hd_nodes.f, hd_nodes.count * lbmDirec * sizeof(double)));
    CUDA_CALL(cudaMemcpy(hd_nodes.f, h_nodes.f, hd_nodes.count * lbmDirec * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMalloc(&hd_nodes.fs, hd_nodes.count * lbmDirec * sizeof(double)));
    CUDA_CALL(cudaMemcpy(hd_nodes.fs, h_nodes.fs, hd_nodes.count * lbmDirec * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMalloc(&hd_nodes.n, hd_nodes.count * sizeof(double)));
    CUDA_CALL(cudaMemcpy(hd_nodes.n, h_nodes.n, hd_nodes.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMalloc(&hd_nodes.u, hd_nodes.count * sizeof(tVect)));
    CUDA_CALL(cudaMemcpy(hd_nodes.u, h_nodes.u, hd_nodes.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMalloc(&hd_nodes.hydroForce, hd_nodes.count * sizeof(tVect)));
    CUDA_CALL(cudaMemcpy(hd_nodes.hydroForce, h_nodes.hydroForce, hd_nodes.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMalloc(&hd_nodes.centrifugalForce, hd_nodes.count * sizeof(tVect)));
    CUDA_CALL(cudaMemcpy(hd_nodes.centrifugalForce, h_nodes.centrifugalForce, hd_nodes.count * sizeof(tVect), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMalloc(&hd_nodes.mass, hd_nodes.count * sizeof(double)));
    CUDA_CALL(cudaMemcpy(hd_nodes.mass, h_nodes.mass, hd_nodes.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMalloc(&hd_nodes.visc, hd_nodes.count * sizeof(double)));
    CUDA_CALL(cudaMemcpy(hd_nodes.visc, h_nodes.visc, hd_nodes.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMalloc(&hd_nodes.basal, hd_nodes.count * sizeof(bool)));
    CUDA_CALL(cudaMemcpy(hd_nodes.basal, h_nodes.basal, hd_nodes.count * sizeof(bool), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMalloc(&hd_nodes.friction, hd_nodes.count * sizeof(double)));
    CUDA_CALL(cudaMemcpy(hd_nodes.friction, h_nodes.friction, hd_nodes.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMalloc(&hd_nodes.age, hd_nodes.count * sizeof(float)));
    CUDA_CALL(cudaMemcpy(hd_nodes.age, h_nodes.age, hd_nodes.count * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMalloc(&hd_nodes.solidIndex, hd_nodes.count * sizeof(unsigned int)));
    CUDA_CALL(cudaMemcpy(hd_nodes.solidIndex, h_nodes.solidIndex, hd_nodes.count * sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMalloc(&hd_nodes.d, hd_nodes.count * lbmDirec * sizeof(unsigned int)));
    CUDA_CALL(cudaMemcpy(hd_nodes.d, h_nodes.d, hd_nodes.count * lbmDirec * sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMalloc(&hd_nodes.curved, hd_nodes.count * sizeof(unsigned int)));
    CUDA_CALL(cudaMemcpy(hd_nodes.curved, h_nodes.curved, hd_nodes.count * sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMalloc(&hd_nodes.type, hd_nodes.count * sizeof(unsigned int)));
    CUDA_CALL(cudaMemcpy(hd_nodes.type, h_nodes.type, hd_nodes.count * sizeof(types), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMalloc(&hd_nodes.p, hd_nodes.count * sizeof(bool)));
    CUDA_CALL(cudaMemcpy(hd_nodes.p, h_nodes.p, hd_nodes.count * sizeof(bool), cudaMemcpyHostToDevice));
    // Copy struct containing device pointers and counts
    CUDA_CALL(cudaMemcpy(d_nodes, &hd_nodes, sizeof(Node2), cudaMemcpyHostToDevice));
#endif
    
}

void LB2::init(cylinderList& cylinders, wallList& walls, particleList& particles, objectList& objects, bool externalSolveCoriolis, bool externalSolveCentrifugal) {
    // Convert from AoS format to SoA and copy to device
    syncCylinders<IMPL>(cylinders);
    syncWalls<IMPL>(walls);
    syncParticles<IMPL>(particles);
    syncObjects<IMPL>(objects);

    //  Lattice Boltzmann initialization steps

    // switchers for apparent accelerations
    h_PARAMS.solveCoriolis = externalSolveCoriolis;
    h_PARAMS.solveCentrifugal = externalSolveCentrifugal;

    // first comes the initialization of the data structures
    cout << "Initializing nodes containers and types" << endl;
    // total number of nodes
    h_PARAMS.totPossibleNodes = h_PARAMS.lbSize[0] * h_PARAMS.lbSize[1] * h_PARAMS.lbSize[2];

    // Count the number of nodes to be created, so memory can be pre-allocated
    std::map<unsigned int, NewNode> newNodes;
    // application of lattice boundaries
    countLatticeBoundaries(newNodes);
    // then the initial node type must be identified for every node (if not specified, it is already Fluid)
    countTypes(newNodes, walls, cylinders, objects);

    // Build a temporary buffer of curves
    std::vector<curve> curves;

    ifstream fluidFileID;
    if (h_PARAMS.lbRestart) {
        // open fluid restart file
        fluidFileID.open(init_params.lbRestartFile.c_str(), ios::in);
        ASSERT(fluidFileID.is_open());
        // check if the restart file size is ok
        unsigned int restartX, restartY, restartZ, restartNodes;
        fluidFileID >> restartX;
        fluidFileID >> restartY;
        fluidFileID >> restartZ;
        fluidFileID >> restartNodes;
        ASSERT(restartX == h_PARAMS.lbSize[0]);
        ASSERT(restartY == h_PARAMS.lbSize[1]);
        ASSERT(restartZ == h_PARAMS.lbSize[2]);
        // read active nodes from file and generate
        // @todo
        throw std::exception("lbRestart is not yet supported.");
        // restartInterface(fluidFileID, restartNodes);
    } else {
        // initialize interface
        countInterface(newNodes);
        // Create and initialize active nodes
        generateInitialNodes(newNodes, curves);
    }

    // initialize variables for wall nodes
    initializeWalls();

    // initializing curved properties
    initializeCurved(curves);

    // Allocate the curves storage
    h_nodes.curveCount = static_cast<unsigned int>(curves.size());
    h_nodes.curves = static_cast<curve*>(malloc(h_nodes.curveCount * sizeof(curve)));
    memcpy(h_nodes.curves, curves.data(), h_nodes.curveCount * sizeof(curve));

    // initialize h_nodes's fluid, interface and active lists
    initializeLists();

    // Setup hd_nodes, copy it to d_nodes
    initDeviceNodes();

    // application of particle initial position
    const double inside_mass = initializeParticleBoundaries<IMPL>();

    // in case mass needs to be kept constant, compute it here
    h_PARAMS.totalMass = 0.0;
    if (h_PARAMS.imposeFluidVolume) {
        // volume and mass is the same in lattice units
        h_PARAMS.totalMass = h_PARAMS.imposedFluidVolume / h_PARAMS.unit.Volume;
    } else {
        switch (problemName) {
        case DRUM:
        {
            h_PARAMS.totalMass = h_PARAMS.fluidMass / h_PARAMS.unit.Mass;
            break;
        }
        case STAVA:
        {
            h_PARAMS.totalMass = 200000.0 / h_PARAMS.unit.Volume;
            break;
        }
        default:
        {
            // @todo This needs to be calculated on device if using CUDA, h_nodes is out of date follow initializeParticleBoundaries
            for (int i = 0; i < h_nodes.activeCount; ++i) {
                const unsigned int a_i = h_nodes.activeI[i];
                if (!h_nodes.isInsideParticle(a_i)) {
                    h_PARAMS.totalMass += h_nodes.mass[a_i];
                }
            }
            break;
        }
        }
    }
    if (h_PARAMS.increaseVolume) {
        h_PARAMS.deltaVolume /= h_PARAMS.unit.Volume;
        h_PARAMS.deltaTime /= h_PARAMS.unit.Time;
    }

    syncParams();

    cout << "Done with initialization" << endl;
}

void LB2::countLatticeBoundaries(std::map<unsigned int, NewNode> &newNodes) {
    // Based on initialiseLatticeBoundaries()
    // XY
    for (unsigned int x = 0; x < h_PARAMS.lbSize[0]; ++x) {
        for (unsigned int y = 0; y < h_PARAMS.lbSize[1]; ++y) {
            // bottom
            newNodes.emplace(h_PARAMS.getIndex(x, y, 0), NewNode{ h_PARAMS.boundary[4] });
            // top
            newNodes.emplace(h_PARAMS.getIndex(x, y, h_PARAMS.lbSize[2] - 1), NewNode{ h_PARAMS.boundary[5] });
        }
    }

    // YZ
    for (unsigned int y = 0; y < h_PARAMS.lbSize[1]; ++y) {
        for (unsigned int z = 0; z < h_PARAMS.lbSize[2]; ++z) {
            // bottom
            newNodes.emplace(h_PARAMS.getIndex(0, y, z), NewNode{ h_PARAMS.boundary[0] });
            // top
            newNodes.emplace(h_PARAMS.getIndex(h_PARAMS.lbSize[0] - 1, y, z), NewNode{ h_PARAMS.boundary[1] });
        }
    }

    // ZX
    for (unsigned int z = 0; z < h_PARAMS.lbSize[2]; ++z) {
        for (unsigned int x = 0; x < h_PARAMS.lbSize[0]; ++x) {
            // bottom
            newNodes.emplace(h_PARAMS.getIndex(x, 0, z), NewNode{ h_PARAMS.boundary[2] });
            // top
            newNodes.emplace(h_PARAMS.getIndex(x, h_PARAMS.lbSize[1] - 1, z), NewNode{ h_PARAMS.boundary[3] });
        }
    }
}
void LB2::countTypes(std::map<unsigned int, NewNode> &newNodes, const wallList& walls, const cylinderList& cylinders, const objectList& objects) {
    countWallBoundaries(newNodes, walls);
    // application of solid cylinders
    countCylinderBoundaries(newNodes, cylinders);
    // application of objects
    countObjectBoundaries(newNodes, objects);
    // initializing topography if one is present
    countTopography(newNodes);
}
void LB2::countWallBoundaries(std::map<unsigned int, NewNode> &newNodes, const wallList& walls) {
    // Based on initializeWallBoundaries()
    // const double wallThickness = 2.0 * h_PARAMS.unit.Length;
    // SOLID WALLS ////////////////////////
    for (unsigned int iw = 0; iw < walls.size(); ++iw) {
        const tVect convertedWallp = walls[iw].p / h_PARAMS.unit.Length;
        const tVect normHere = walls[iw].n;
        const unsigned int indexHere = walls[iw].index;
        const bool slipHere = walls[iw].slip;
        const bool movingHere = walls[iw].moving;
        for (unsigned int it = 0; it < h_PARAMS.totPossibleNodes; ++it) {
            // check if the node is solid
            // all walls have max thickness 2 nodes
            const tVect pos = h_PARAMS.getPosition(it);
            const double wallDistance = pos.distance2Plane(convertedWallp, normHere);
            if (wallDistance > -2.0 && wallDistance < 0.0) {
                //check for borders in limted walls
                if (walls[iw].limited) {
                    const double xHere = pos.x * h_PARAMS.unit.Length;
                    const double yHere = pos.y * h_PARAMS.unit.Length;
                    const double zHere = pos.z * h_PARAMS.unit.Length;
                    // check if beyond limits
                    if (xHere < walls[iw].xMin || xHere > walls[iw].xMax ||
                        yHere < walls[iw].yMin || yHere > walls[iw].yMax ||
                        zHere < walls[iw].zMin || zHere > walls[iw].zMax) {
                        continue;
                    }
                }
                // setting type: 5-6=slip, 7-8=no-slip
                if (slipHere) {
                    // setting type for slip: 5=static, 6=moving
                    if (movingHere) {
                        const auto nn = newNodes.emplace(it, NewNode{ SLIP_DYN_WALL, iw });
                        nn.first->second.solidIndex = iw; // Update solidIndex, even if node creation already requested
                    } else {
                        const auto nn = newNodes.emplace(it, NewNode{ SLIP_STAT_WALL, iw });
                        nn.first->second.solidIndex = iw; // Update solidIndex, even if node creation already requested
                    }
                } else {
                    // setting type for no-slip: 7=static, 8=moving
                    if (movingHere) {
                        const auto nn = newNodes.emplace(it, NewNode{ DYN_WALL, iw });
                        nn.first->second.solidIndex = iw; // Update solidIndex, even if node creation already requested
                    } else {
                        const auto nn = newNodes.emplace(it, NewNode{ STAT_WALL, iw });
                        nn.first->second.solidIndex = iw; // Update solidIndex, even if node creation already requested
                    }
                }
            }
        }
    }

}
void LB2::countObjectBoundaries(std::map<unsigned int, NewNode> &newNodes, const objectList& objects) {
    // Based on initializeObjectBoundaries()
    // SOLID WALLS ////////////////////////
    for (int io = 0; io < objects.size(); ++io) {
        const tVect convertedPosition = objects[io].x0 / h_PARAMS.unit.Length;
        const double convertedRadius = objects[io].r / h_PARAMS.unit.Length;
        const unsigned int indexHere = objects[io].index;
        for (unsigned int it = 0; it < h_PARAMS.totPossibleNodes; ++it) {
            const tVect nodePosition = h_PARAMS.getPosition(it);
            if (nodePosition.insideSphere(convertedPosition, convertedRadius)) {
                newNodes.emplace(it, NewNode{ OBJ, indexHere });
            }
        }
    }
}
void LB2::countCylinderBoundaries(std::map<unsigned int, NewNode> &newNodes, const cylinderList& cylinders) {
    // Based on initializeCylinderBoundaries()
    // SOLID CYLINDERS ////////////////////////
    for (int ic = 0; ic < cylinders.size(); ++ic) {

        const tVect convertedCylinderp1 = cylinders[ic].p1 / h_PARAMS.unit.Length;
        const tVect naxesHere = cylinders[ic].naxes;
        const double convertedRadius = cylinders[ic].R / h_PARAMS.unit.Length;
        const unsigned int indexHere = cylinders[ic].index;
        const bool slipHere = cylinders[ic].slip;
        const bool movingHere = cylinders[ic].moving;
        for (unsigned int it = 0; it < h_PARAMS.totPossibleNodes; ++it) {
            // creating solid cells
            const bool isOutside = h_PARAMS.getPosition(it).insideCylinder(convertedCylinderp1, naxesHere, convertedRadius, convertedRadius + 3.0);
            const bool isInside = h_PARAMS.getPosition(it).insideCylinder(convertedCylinderp1, naxesHere, max(convertedRadius - 3.0, 0.0), convertedRadius);
            if ((cylinders[ic].type == FULL && isInside) ||
                (cylinders[ic].type == EMPTY && isOutside)) {
                //check for borders in limted walls
                if (cylinders[ic].limited) {
                    const tVect here = h_PARAMS.getPosition(it) * h_PARAMS.unit.Length;
                    // check if beyond limits
                    if (here.x < cylinders[ic].xMin || here.x > cylinders[ic].xMax ||
                        here.y < cylinders[ic].yMin || here.y > cylinders[ic].yMax ||
                        here.z < cylinders[ic].zMin || here.z > cylinders[ic].zMax) {
                        continue;
                    }
                }
                // setting type: 5-6=slip, 7-8=no-slip
                if (slipHere) {
                    // setting type for slip: 5=static, 6=moving
                    if (movingHere) {
                        newNodes.emplace(it, NewNode{ SLIP_DYN_WALL, indexHere });
                    } else {
                        newNodes.emplace(it, NewNode{ SLIP_STAT_WALL, indexHere });
                    }
                } else {
                    // setting type for no-slip: 7=static, 8=moving
                    if (movingHere) {
                        newNodes.emplace(it, NewNode{ DYN_WALL, indexHere });
                    } else {
                        newNodes.emplace(it, NewNode{ STAT_WALL, indexHere });
                    }
                }
            }
        }
    }
}
void LB2::countTopography(std::map<unsigned int, NewNode> &newNodes) {
    // Based on initializeTopography()
    
    const double surfaceThickness = 1.75 * h_PARAMS.unit.Length;

    // TOPOGRAPHY ////////////////////////
    if (h_PARAMS.lbTopography) {
        lbTop.readFromFile(init_params.lbTopographyFile, h_PARAMS.translateTopographyX, h_PARAMS.translateTopographyY, h_PARAMS.translateTopographyZ);
        lbTop.show();
        // check if topography grid contains the fluid domain
        ASSERT(lbTop.coordX[0] < h_PARAMS.unit.Length);
        ASSERT(lbTop.coordY[0] < h_PARAMS.unit.Length);

        cout << "lbTop.coordX[lbTop.sizeX - 1]=" << lbTop.coordX[lbTop.sizeX - 1] << endl;
        cout << "lbSize[0]) * unit.Length=" << h_PARAMS.lbSize[0] * h_PARAMS.unit.Length << endl;
        ASSERT(lbTop.coordX[lbTop.sizeX - 1] > h_PARAMS.lbSize[0] * h_PARAMS.unit.Length);
        cout << "lbTop.coordY[lbTop.sizeY - 1]=" << lbTop.coordY[lbTop.sizeY - 1] << endl;
        cout << "lbSize[1]) * unit.Length=" << h_PARAMS.lbSize[1] * h_PARAMS.unit.Length << endl;
        ASSERT(lbTop.coordY[lbTop.sizeY - 1] > h_PARAMS.lbSize[1] * h_PARAMS.unit.Length);


        for (unsigned int ix = 1; ix < h_PARAMS.lbSize[0] - 1; ++ix) {
            for (unsigned int iy = 1; iy < h_PARAMS.lbSize[1] - 1; ++iy) {
                for (unsigned int iz = 1; iz < h_PARAMS.lbSize[2] - 1; ++iz) {
                    const tVect nodePosition = tVect(ix, iy, iz) * h_PARAMS.unit.Length;
                    const double distanceFromTopography = lbTop.distance(nodePosition);
                    
                    if (distanceFromTopography < 0.0 && distanceFromTopography>-1.0 * surfaceThickness) {
                        const unsigned int it = ix + iy * h_PARAMS.lbSize[0] + iz * h_PARAMS.lbSize[0] * h_PARAMS.lbSize[1];
                        newNodes.emplace(it, NewNode{ TOPO });
                    }
                }
            }
        }
    }
}
void LB2::countInterface(std::map<unsigned int, NewNode> &newNodes) {
    // Based on initializeInterface()
    // @Currently only default case is supported
    // creates an interface electing interface cells from active cells
    if (h_PARAMS.lbTopographySurface) {
        // Formerly setTopographySurface()
        for (unsigned int it = 0; it < h_PARAMS.totPossibleNodes; ++it) {
            if (newNodes.find(it) == newNodes.end()) {
                // control is done in real coordinates
                const tVect nodePosition = h_PARAMS.getPosition(it) * h_PARAMS.unit.Length;
                const double surfaceIsoparameterHere = lbTop.surfaceIsoparameter(nodePosition);
                if (surfaceIsoparameterHere > 0.0 && surfaceIsoparameterHere <= 1.0) {// setting solidIndex
                    newNodes.emplace(it, NewNode{ LIQUID });
                }
            }
        }
    } else {
        switch (problemName) {
        case NONE:
        case SHEARCELL:
        case AVALANCHE:
        case DRUM:
        case NET:
        case BARRIER:
        case ZHOU:
        case OPENBARRIER:
        case HONGKONG:
        case STVINCENT:
        case STAVA:
        case NIGRO:
        case CAROLINE:
        case DAMBREAK:
        case GRAY_DAMBREAK:
        case GRAY_DAMBREAK_2D:
        case INCLINEFLOW:
        case HOURGLASS:
        case IERVOLINO:
        case IERVOLINO_2D:
        case IERVOLINO_CYLINDERTEST:
        case HEAP:
        case TRIAXIAL:
        case JOP:
        case WILL:
        case WILL_SETTLING:
        case MANGENEY:
        case GRAY:
        case ESERCITAZIONE:
        case FILIPPO_SILOS:
        case HK_SMALL:
        case HK_LARGE:
        case KELVIN:
        case SHEARCELL2023:
        case INTRUDER:
        case OBJMOVING:
        default:
            {
                cout << "X=(" << double(h_PARAMS.freeSurfaceBorders[0]) * h_PARAMS.unit.Length << ", " << double(h_PARAMS.freeSurfaceBorders[1]) * h_PARAMS.unit.Length << ")" << endl;
                cout << "Y=(" << double(h_PARAMS.freeSurfaceBorders[2]) * h_PARAMS.unit.Length << ", " << double(h_PARAMS.freeSurfaceBorders[3]) * h_PARAMS.unit.Length << ")" << endl;
                cout << "Z=(" << double(h_PARAMS.freeSurfaceBorders[4]) * h_PARAMS.unit.Length << ", " << double(h_PARAMS.freeSurfaceBorders[5]) * h_PARAMS.unit.Length << ")" << endl;
                for (unsigned int it = 0; it < h_PARAMS.totPossibleNodes; ++it) {
                    if (newNodes.find(it) == newNodes.end()) {
                        // creating fluid cells
                        const tVect pos = h_PARAMS.getPosition(it);
                        if ((pos.x > h_PARAMS.freeSurfaceBorders[0]) &&
                            (pos.x < h_PARAMS.freeSurfaceBorders[1]) &&
                            (pos.y > h_PARAMS.freeSurfaceBorders[2]) &&
                            (pos.y < h_PARAMS.freeSurfaceBorders[3]) &&
                            (pos.z > h_PARAMS.freeSurfaceBorders[4]) &&
                            (pos.z < h_PARAMS.freeSurfaceBorders[5])) {
                            newNodes.emplace(it, NewNode{ LIQUID });
                        }
                    }
                }
                break;
            }
        }
    }
}
void LB2::generateInitialNodes(const std::map<unsigned int, NewNode> &newNodes, std::vector<curve> &curves) {
    // Allocate enough memory for these nodes
    assert(h_nodes.count == 0);  // No nodes should exist at the time this is called
    h_nodes.count = newNodes.size();
    // Allocate host buffers
    h_nodes.coord = static_cast<unsigned int*>(malloc(h_nodes.count * sizeof(unsigned int)));
    h_nodes.f = static_cast<double*>(malloc(h_nodes.count * lbmDirec * sizeof(double)));
    h_nodes.fs = static_cast<double*>(malloc(h_nodes.count * lbmDirec * sizeof(double)));
    h_nodes.n = static_cast<double*>(malloc(h_nodes.count * sizeof(double)));
    h_nodes.u = static_cast<tVect*>(malloc(h_nodes.count * sizeof(tVect)));
    h_nodes.hydroForce = static_cast<tVect*>(malloc(h_nodes.count * sizeof(tVect)));
    h_nodes.centrifugalForce = static_cast<tVect*>(malloc(h_nodes.count * sizeof(tVect)));
    h_nodes.mass = static_cast<double*>(malloc(h_nodes.count * sizeof(double)));
    h_nodes.visc = static_cast<double*>(malloc(h_nodes.count * sizeof(double)));
    h_nodes.basal = static_cast<bool*>(malloc(h_nodes.count * sizeof(bool)));
    h_nodes.friction = static_cast<double*>(malloc(h_nodes.count * sizeof(double)));
    h_nodes.age = static_cast<float*>(malloc(h_nodes.count * sizeof(float)));
    h_nodes.solidIndex = static_cast<unsigned int*>(malloc(h_nodes.count * sizeof(unsigned int)));
    h_nodes.d = static_cast<unsigned int*>(malloc(h_nodes.count * lbmDirec * sizeof(unsigned int)));
    h_nodes.curved = static_cast<unsigned int*>(malloc(h_nodes.count * sizeof(unsigned int)));
    h_nodes.type = static_cast<types*>(malloc(h_nodes.count * sizeof(types)));
    h_nodes.p = static_cast<bool*>(malloc(h_nodes.count * sizeof(bool)));
    // Zero initialisation
    memset(h_nodes.f, 0, h_nodes.count * lbmDirec * sizeof(double));
    memset(h_nodes.fs, 0, h_nodes.count * lbmDirec * sizeof(double));
    memset(h_nodes.n, 0, h_nodes.count * sizeof(double));
    memset(h_nodes.u, 0, h_nodes.count * sizeof(tVect));
    memset(h_nodes.hydroForce, 0, h_nodes.count * sizeof(tVect));
    memset(h_nodes.mass, 0, h_nodes.count * sizeof(double));
    // h_nodes.visc is instead init to 1 below
    memset(h_nodes.basal, 0, h_nodes.count * sizeof(bool));
    memset(h_nodes.friction, 0, h_nodes.count * sizeof(double));
    memset(h_nodes.age, 0, h_nodes.count * sizeof(float));
    memset(h_nodes.curved, std::numeric_limits<unsigned int>::max(), h_nodes.count * sizeof(unsigned int));
    memset(h_nodes.p, 0, h_nodes.count * sizeof(bool));
    // Perform the generateNode() loop for each item in newNodes
    std::map<unsigned int, unsigned int> idIndexMap;
    {
        unsigned int i = 0;
        for (const auto& [id, nn] : newNodes) {
            idIndexMap.emplace(id, i);
            h_nodes.visc[i] = 1;
            h_nodes.coord[i] = id;
            h_nodes.type[i] = nn.type;
            h_nodes.solidIndex[i] = nn.solidIndex;
            // set centrifugal acceleration
            h_nodes.centrifugalForce[i] = computeCentrifugal(h_nodes.getPosition(i), PARAMS.rotationCenter, PARAMS.rotationSpeed);
            ++i;
        }
    }
    // Perform a second pass for handling neighbours
    for (unsigned int i = 0; i < h_nodes.count; ++i) {
        // findNeighbors()
        std::array<unsigned int, lbmDirec> neighbourCoord;
        for (int j = 1; j < lbmDirec; ++j) {
            neighbourCoord[j] = h_nodes.coord[i] + PARAMS.ne[j];
        }
        // Boundary conditions
        if (h_nodes.isWall(i)) {
            const std::array<int, 3> pos = h_nodes.getGridPosition(i);// getPosition() adds + 0.5x
            if (pos[0] == 0) {
                for (unsigned int j = 1; j < lbmDirec; ++j) {
                    if (v[j].dot(Xp) < 0.0) {
                        neighbourCoord[j] = h_nodes.coord[i];
                    }
                }
            } else if (pos[0] == static_cast<int>(PARAMS.lbSize[0] - 1)) {
                for (unsigned int j = 1; j < lbmDirec; ++j) {
                    if (v[j].dot(Xp) > 0.0) {
                        neighbourCoord[j] = h_nodes.coord[i];
                    }
                }
            }
            if (pos[1] == 0) {
                for (unsigned int j = 1; j < lbmDirec; ++j) {
                    if (v[j].dot(Yp) < 0.0) {
                        neighbourCoord[j] = h_nodes.coord[i];
                    }
                }
            } else if (pos[1] == static_cast<int>(PARAMS.lbSize[1] - 1)) {
                for (unsigned int j = 1; j < lbmDirec; ++j) {
                    if (v[j].dot(Yp) > 0.0) {
                        neighbourCoord[j] = h_nodes.coord[i];
                    }
                }
            }
            if (pos[2] == 0) {
                for (unsigned int j = 1; j < lbmDirec; ++j) {
                    if (v[j].dot(Zp) < 0.0) {
                        neighbourCoord[j] = h_nodes.coord[i];
                    }
                }
            } else if (pos[2] == static_cast<int>(PARAMS.lbSize[2] - 1)) {
                for (unsigned int j = 1; j < lbmDirec; ++j) {
                    if (v[j].dot(Zp) > 0.0) {
                        neighbourCoord[j] = h_nodes.coord[i];
                    }
                }
            }
        } else if (h_nodes.isActive(i)) {
            /*
            // PERIODICITY ////////////////////////////////////////////////
            // assigning periodicity conditions (this needs to be done after applying boundary conditions)
            // runs through free cells and identifies neighboring cells. If neighbor cell is
            // a special cell (periodic) then the proper neighboring condition is applied
            // calculates the effect of periodicity
            */
            // 
            // neighboring and periodicity vector for boundary update
            std::array<unsigned int, lbmDirec> pbc = {};

            auto f = idIndexMap.find(neighbourCoord[1]);
            if (f != idIndexMap.end()) {
                if (h_nodes.type[f->second] == PERIODIC) {
                    for (int j = 1; j < lbmDirec; ++j) {
                        if (v[j].dot(Xp) > 0.0) {
                            pbc[j] -= PARAMS.domain[0];
                        }
                    }
                }
            }
            f = idIndexMap.find(neighbourCoord[2]);
            if (f != idIndexMap.end()) {
                if (h_nodes.type[f->second] == PERIODIC) {
                    for (int j = 1; j < lbmDirec; ++j) {
                        if (v[j].dot(Xp) < 0.0) {
                            pbc[j] += PARAMS.domain[0];
                        }
                    }
                }
            }
            f = idIndexMap.find(neighbourCoord[3]);
            if (f != idIndexMap.end()) {
                if (h_nodes.type[f->second] == PERIODIC) {
                    for (int j = 1; j < lbmDirec; ++j) {
                        if (v[j].dot(Yp) > 0.0) {
                            pbc[j] -= PARAMS.domain[1];
                        }
                    }
                }
            }
            f = idIndexMap.find(neighbourCoord[4]);
            if (f != idIndexMap.end()) {
                if (h_nodes.type[f->second] == PERIODIC) {
                    for (int j = 1; j < lbmDirec; ++j) {
                        if (v[j].dot(Yp) < 0.0) {
                            pbc[j] += PARAMS.domain[1];
                        }
                    }
                }
            }
            f = idIndexMap.find(neighbourCoord[5]);
            if (f != idIndexMap.end()) {
                if (h_nodes.type[f->second] == PERIODIC) {
                    for (int j = 1; j < lbmDirec; ++j) {
                        if (v[j].dot(Zp) > 0.0) {
                            pbc[j] -= PARAMS.domain[2];
                        }
                    }
                }
            }
            f = idIndexMap.find(neighbourCoord[6]);
            if (f != idIndexMap.end()) {
                if (h_nodes.type[f->second] == PERIODIC) {
                    for (int j = 1; j < lbmDirec; ++j) {
                        if (v[j].dot(Zp) < 0.0) {
                            pbc[j] += PARAMS.domain[2];
                        }
                    }
                }
            }

            // apply periodicity
            for (int j = 1; j < lbmDirec; ++j) {
                neighbourCoord[j] += pbc[j];
            }
        }
        // assign neighbour nodes
        for (int j = 1; j < lbmDirec; ++j) {
            auto f = idIndexMap.find(neighbourCoord[j]);
            // check if node at that location exists
            if (f != idIndexMap.end()) {
                const unsigned int l_i = f->second;
                // assign neighbor for local node
                h_nodes.d[j * h_nodes.count + i] = l_i;
                // if neighbor node is also active, link it to local node
                if (h_nodes.isActive(i)) {
                    h_nodes.d[opp[j] * h_nodes.count + l_i] = i;
                    // if the neighbor is a curved wall, set parameters accordingly
                    if (h_nodes.type[l_i] == TOPO) {
                        if (h_nodes.curved[i] == std::numeric_limits<unsigned int>::max()) {
                            h_nodes.curved[i] = static_cast<unsigned int>(curves.size());
                            curves.emplace_back();
                        }
                        // set curved
                        const tVect nodePosHere = PARAMS.unit.Length * h_nodes.getPosition(i);
                        // xf - xw
                        const double topographyDistance = 1.0 * lbTop.directionalDistance(nodePosHere, vDirec[j]) / PARAMS.unit.Length;
                        // wall normal
                        curves.back().wallNormal = lbTop.surfaceNormal(nodePosHere);
                        //cout << topographyDistance << endl;
                        const double deltaHere = topographyDistance / vNorm[j];
                        curves.back().delta[j] = std::min(0.99, std::max(0.01, deltaHere));
                        curves.back().computeCoefficients();
                    }
                    if (h_nodes.isWall(l_i)) {
                        h_nodes.basal[i] = true;
                    }
                }
            } else {
                // Neighbour is gas
                h_nodes.d[j * h_nodes.count + i] = std::numeric_limits<unsigned int>::max();
            }
        }
    }

    // Formerly initializeVariables()    
    cout << "Initializing variables" << endl;
    // note that interface is not defined here. All fluid, interface and gas cells are uninitialized at the moment
    // calculate maximum height of the fluid

    // find "taller" and "deepest" points
    double minProjection = std::numeric_limits<double>::max();
    double maxProjection = -std::numeric_limits<double>::max();
        
    if (!PARAMS.solveCentrifugal) {
        for (unsigned int i = 0; i < h_nodes.count; ++i) {
            if (h_nodes.isActive(i)) {
                const tVect position = h_nodes.getPosition(i);
                const double projection = position.dot(PARAMS.lbF);
                minProjection = std::min(minProjection, projection);
                maxProjection = std::max(maxProjection, projection);
            }
        }
        cout << "minProjection = " << minProjection << endl;
    } else {
        for (unsigned int i = 0; i < h_nodes.count; ++i) {
            if (h_nodes.isActive(i)) {
                const tVect position = h_nodes.getPosition(i);
                const double projection = position.dot(h_nodes.centrifugalForce[i]);
                minProjection = std::min(minProjection, projection);
                maxProjection = std::max(maxProjection, projection);
            }
        }
        cout << "minProjection = " << minProjection << endl;
    }

    // checking for boundary between gas and fluid and assigning interface properties
    // at this point fluid cells contain actual fluid cells and potential interface cells, so we create the node anyway
    double massFluid = 0.0;
    double massInterface = 0.0;
    for (unsigned int i = 0; i < h_nodes.count; ++i) {
        if (h_nodes.type[i] == LIQUID) {
            // check if it is interface
            for (int j = 1; j < lbmDirec; ++j) {
                unsigned int linkNode = h_nodes.d[j * h_nodes.count + i];
                if (linkNode == std::numeric_limits<unsigned int>::max()) {
                    h_nodes.type[i] = INTERFACE;
                    break;
                }
            }
        }
        // now assign macroscopic quantities accordingly
        // FLUID NODES ////
        if (h_nodes.type[i] == LIQUID) {
            massFluid += 1.0;
            // setting macroscopic variables
            // density is calculated using hydrostatic profile
            const tVect position = h_nodes.getPosition(i);
            if (!PARAMS.solveCentrifugal) {
                const double projection = position.dot(PARAMS.lbF);
                h_nodes.initialize(i, PARAMS.fluidMaterial.initDensity + 3.0 * PARAMS.fluidMaterial.initDensity * (projection-minProjection), PARAMS.initVelocity, PARAMS.fluidMaterial.initDensity, PARAMS.fluidMaterial.initDynVisc, PARAMS.lbF, 1.0, Zero);
            } else {
                const double projection = position.dot(h_nodes.centrifugalForce[i]);
                h_nodes.initialize(i, PARAMS.fluidMaterial.initDensity + 3.0 * PARAMS.fluidMaterial.initDensity * (projection-minProjection), PARAMS.initVelocity, PARAMS.fluidMaterial.initDensity, PARAMS.fluidMaterial.initDynVisc, PARAMS.lbF, 1.0, PARAMS.rotationSpeed);
            }
        }// INTERFACE NODES ////
        else if (h_nodes.type[i] == INTERFACE) {
            massInterface += 0.5;
            // setting macroscopic variables
            h_nodes.initialize(i, PARAMS.fluidMaterial.initDensity, PARAMS.initVelocity, 0.5 * PARAMS.fluidMaterial.initDensity, PARAMS.fluidMaterial.initDynVisc, PARAMS.lbF, 1.0, PARAMS.rotationSpeed);
        }

    }
    cout << "Approximate volume = " << massFluid * PARAMS.unit.Volume << " (fluid body), " << massInterface * PARAMS.unit.Volume << " (interface), " << (massFluid + massInterface) * PARAMS.unit.Volume << " (tot), " << endl;
}
void LB2::initializeWalls() {
    cout << "Initializing wall nodes" << endl;
    const double zero = 0.0;

    std::vector<unsigned int> wallNodes;

    // initializing wall nodes
    // note that, in the hypothesis that these walls are not evolving, only nodes at the interface need creation
    for (unsigned int i = 0; i < h_nodes.count; ++i) {
        if (h_nodes.isWall(i)) {
            // initialize node
            // STATIC WALL NODES ////
            if (h_nodes.type[i] == STAT_WALL ||
                h_nodes.type[i] == SLIP_STAT_WALL ||
                h_nodes.type[i] ==  OBJ ||
                h_nodes.type[i] == TOPO) {
                // reset velocity and mass (useful for plotting)
                // density=0.0; velocity=(0.0,0.0,0.0), mass=0.0; viscosity=0.0; force=(0.0,0.0,0.0)
                h_nodes.initialize(i, PARAMS.fluidMaterial.initDensity, Zero, zero, zero, Zero, 1.0, Zero);
            }// DYNAMIC WALL NODES ////
            else if (h_nodes.type[i] == DYN_WALL || 
                     h_nodes.type[i] == SLIP_DYN_WALL || 
                     h_nodes.type[i] == CYL) {
                // need to define velocity. It could be part of a cylinder or wall, we check both
                tVect solidVelocity;
                const tVect nodePosition = h_nodes.getPosition(i);
                unsigned int solidIndex = h_nodes.solidIndex[i];
                // wall
                if (solidIndex < h_walls.count && nodePosition.insidePlane(h_walls.p[solidIndex] / PARAMS.unit.Length, h_walls.n[solidIndex])) {
                    solidVelocity = h_walls.getSpeed(solidIndex, nodePosition * PARAMS.unit.Length) / PARAMS.unit.Speed;
                }// cylinder
                else if (solidIndex < h_cylinders.count && !nodePosition.insideCylinder(h_cylinders.p1[solidIndex] / PARAMS.unit.Length, h_cylinders.naxes[solidIndex], 0.0, h_cylinders.R[solidIndex] / PARAMS.unit.Length)) {
                    solidVelocity = h_cylinders.getSpeed(solidIndex, nodePosition * PARAMS.unit.Length) / PARAMS.unit.Speed;
                }// objects
                else if (solidIndex < h_objects.count && nodePosition.insideSphere(h_objects.x0[solidIndex] / PARAMS.unit.Length, h_objects.r[solidIndex] / PARAMS.unit.Length)) {
                    solidVelocity = h_objects.x1[solidIndex] / PARAMS.unit.Speed;
                }
                // reset velocity and mass (useful for plotting)
                // density=0.0; velocity=solidVelocity, mass=0.0; viscosity=0.0; force=(0.0,0.0,0.0)
                h_nodes.initialize(i, PARAMS.fluidMaterial.initDensity, solidVelocity, zero, zero, Zero, 1.0, PARAMS.rotationSpeed);
            }
            // add node to list
            wallNodes.push_back(i);
        }
    }
    // Allocate the wall nodes storage
    h_nodes.wallCount = static_cast<unsigned int>(wallNodes.size());
    h_nodes.wallI = static_cast<unsigned int*>(malloc(h_nodes.wallCount * sizeof(unsigned int)));
    memcpy(h_nodes.wallI, wallNodes.data(), h_nodes.wallCount * sizeof(unsigned int));    
}
void LB2::initializeCurved(std::vector<curve> &curves) {
    cout << "Initializing curved boundaries" << endl;
    for (unsigned int i = 0; i < h_nodes.wallCount; ++i) {
        const unsigned int w_i = h_nodes.wallI[i];
        if (h_nodes.type[w_i] == CYL) {
            assert(h_nodes.curved[w_i] == std::numeric_limits<unsigned int>::max());
            h_nodes.curved[w_i] = static_cast<unsigned int>(curves.size());
            curves.emplace_back();
            const tVect nodePos = PARAMS.unit.Length * h_nodes.getPosition(w_i);
            for (int j = 1; j < lbmDirec; ++j) {
                curves.back().delta[j] = 1.0 - h_cylinders.segmentIntercept(0, nodePos, PARAMS.unit.Length * v[j]);
                curves.back().computeCoefficients();
            }
        }
    }
}
void LB2::initializeLists() {
    cout << "Resetting lists ...";

    // note that interface is not defined here. All fluid, interface and gas cells are 0 at the moment
    std::vector<unsigned int> fluidNodes;
    std::vector<unsigned int> interfaceNodes;

    // creating list and initialize macroscopic variables for all nodes except walls
    for (unsigned int i = 0; i < h_nodes.count; ++i) {
        if (h_nodes.type[i] == LIQUID) {
            fluidNodes.push_back(i);
        } else if (h_nodes.type[i] == INTERFACE) {
            interfaceNodes.push_back(i);
        }

    }

    // Array to Buffer
    assert(!h_nodes.fluidI);
    h_nodes.fluidCount = static_cast<unsigned int>(fluidNodes.size());
    h_nodes.fluidI = static_cast<unsigned int*>(malloc(h_nodes.fluidCount * sizeof(unsigned int)));
    memcpy(h_nodes.fluidI, fluidNodes.data(), h_nodes.fluidCount * sizeof(unsigned int));

    assert(!h_nodes.interfaceI);
    h_nodes.interfaceCount = static_cast<unsigned int>(interfaceNodes.size());
    h_nodes.interfaceI = static_cast<unsigned int*>(malloc(h_nodes.interfaceCount * sizeof(unsigned int)));
    memcpy(h_nodes.interfaceI, interfaceNodes.data(), h_nodes.interfaceCount * sizeof(unsigned int));

    // Build a sorted active nodes list
    fluidNodes.insert(fluidNodes.end(), interfaceNodes.begin(), interfaceNodes.end());
    std::sort(fluidNodes.begin(), fluidNodes.end());

    // Array to buffer
    assert(!h_nodes.activeI);
    h_nodes.activeCount = static_cast<unsigned int>(fluidNodes.size());
    h_nodes.activeI = static_cast<unsigned int*>(malloc(h_nodes.activeCount * sizeof(unsigned int)));
    memcpy(h_nodes.activeI, fluidNodes.data(), h_nodes.activeCount * sizeof(unsigned int));
    
    cout << " done" << endl;
}
void LB2::step(const DEM &dem, bool io_demSolver) {
    this->syncDEM(dem.elmts, dem.particles, dem.walls, dem.objects);

    if (io_demSolver) {
        this->latticeBoltzmannCouplingStep(dem.newNeighborList);
    }

    if (dem.demTime >= dem.demInitialRepeat) {
        this->latticeBoltzmannStep();

        // Lattice Boltzmann core steps @todo after latticeBoltzmannStep() has been tested
        // if (this->freeSurface) {
        //     this->latticeBoltzmannFreeSurfaceStep();
        // }
    }
}

void LB2::syncDEM(const elmtList &elmts, const particleList &particles, const wallList &walls, const objectList &objects) {
    // Sync DEM data to structure of arrays format (and device memory)
    syncElements<IMPL>(elmts);
    syncParticles<IMPL>(particles);
    syncWalls<IMPL>(walls);
    syncObjects<IMPL>(objects);
}
void LB2::setParams(const LBParams& params, const LBInitParams& initParams, bool skip_sync) {
    // CPU
    h_PARAMS = params;
    init_params = initParams;
    // CUDA
    if (!skip_sync)
        syncParams();
}
void LB2::syncParams() {
#ifdef USE_CUDA
    CUDA_CALL(cudaMemcpyToSymbol(d_PARAMS, &h_PARAMS, sizeof(LBParams)));
#endif
}

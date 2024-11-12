#include "LB2.h"

#include <cstdlib>

#include "cuda_helper.h"

#include "DEM.h"

void LB2::step(const DEM &dem, bool io_demSolver) {
    this->syncDEM(dem.elmts, dem.particles, dem.walls, dem.objects);

    if (io_demSolver) {
        this->latticeBoltzmannCouplingStep(dem.newNeighborList);
    }

    if (dem.demTime >= dem.demInitialRepeat) {
        // @todo elmts/particles sync only occurs if io_demSolver = true
        this->latticeBoltzmannStep();

        // Lattice Boltzmann core steps
        if (this->freeSurface) {
            this->latticeBoltzmannFreeSurfaceStep();
        }
    }
}

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
    h_elements.count = elements.size();
    // Repackage host particle data from array of structures, to structure of arrays
     for (int i = 0; i < h_elements.count; ++i) {
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
            totalComponents += e.components.size();
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
                totalComponents += elements[i].components.size();
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
    h_particles.count = particles.size();
    // Repackage host particle data from array of structures, to structure of arrays
    for (int i = 0; i < h_particles.count; ++i) {
        h_particles.clusterIndex[i] = particles[i].clusterIndex;
        h_particles.r[i] = particles[i].r;
        h_particles.x0[i] = particles[i].x0;
        h_particles.radiusVec[i] = particles[i].radiusVec;
    }
}
template<>
void LB2::syncWalls<CPU>(const wallList &walls) {
    if (h_walls.count < walls.size()) {
        // Grow host buffers
        if (h_walls.FHydro) {
            free(h_walls.FHydro);
        }
        h_walls.FHydro = (tVect*)malloc(walls.size() * sizeof(tVect));
    }
    // Update size
    h_walls.count = walls.size();
    // Repackage host particle data from array of structures, to structure of arrays
    // for (int i = 0; i < h_walls.count; ++i) {
        // h_walls.FHydro[i] = walls[i].FHydro; // Zero'd before use in streaming()
    // }
}
template<>
void LB2::syncObjects<CPU>(const objectList &objects) {
    if (h_objects.count < objects.size()) {
        // Grow host buffers
        if (h_objects.FHydro) {
            free(h_objects.FHydro);
        }
        h_objects.FHydro = (tVect*)malloc(objects.size() * sizeof(tVect));
    }
    // Update size
    h_objects.count = objects.size();
    // Repackage host particle data from array of structures, to structure of arrays
    // for (int i = 0; i < h_objects.count; ++i) {
        // h_objects.FHydro[i] = objects[i].FHydro; // Zero'd before use in streaming()
    // }
}
#ifdef USE_CUDA
template<>
bool LB2::syncElements<CUDA>(const elmtList &elements) {
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
        CUDA_CALL(cudaMemcpy(d_elements, &h_elements, sizeof(Element2), cudaMemcpyHostToDevice));
    } else {
        // Copy updated device pointers to device (@todo When/where is d_elements allocated??)
        CUDA_CALL(cudaMemcpy(&d_elements->count, &h_elements.count, sizeof(unsigned int), cudaMemcpyHostToDevice));
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
        hd_particles.count = particles.size();
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
void LB2::syncWalls<CUDA>(const wallList &walls) {
    // Copy latest particle data from HOST DEM to the device
    // @todo Can these copies be done ahead of time async?
    // @todo These copies will be redundant when DEM is moved to CUDA
    this->syncWalls<CPU>(walls);
    if (hd_walls.count < walls.size()) {
        // Grow device buffers
        if (hd_walls.FHydro) {
            CUDA_CALL(cudaFree(hd_walls.FHydro));
        }
        CUDA_CALL(cudaMalloc(&hd_walls.FHydro, h_walls.count * sizeof(tVect)));
        hd_walls.count = h_walls.count;
        // Copy updated device pointers to device
        CUDA_CALL(cudaMemcpy(d_walls, &h_walls, sizeof(Wall2), cudaMemcpyHostToDevice));
    } else if(hd_walls.count != walls.size()) {
        // Buffer has shrunk, so just update size
        hd_walls.count = walls.size();
        // Copy updated particle count to device (@todo When/where is d_elements allocated??)
        CUDA_CALL(cudaMemcpy(&d_walls->count, &h_walls.count, sizeof(unsigned int), cudaMemcpyHostToDevice));
    }
    // Copy data to device buffers
    // CUDA_CALL(cudaMemcpy(hd_walls.FHydro, h_walls.FHydro, h_walls.count * sizeof(tVect), cudaMemcpyHostToDevice)); // Zero'd before use in streaming()
}
template<>
void LB2::syncObjects<CUDA>(const objectList &walls) {
    // Copy latest particle data from HOST DEM to the device
    // @todo Can these copies be done ahead of time async?
    // @todo These copies will be redundant when DEM is moved to CUDA
    this->syncObjects<CPU>(walls);
    if (hd_objects.count < walls.size()) {
        // Grow device buffers
        if (hd_objects.FHydro) {
            CUDA_CALL(cudaFree(hd_objects.FHydro));
        }
        CUDA_CALL(cudaMalloc(&hd_objects.FHydro, h_objects.count * sizeof(tVect)));
        hd_objects.count = h_objects.count;
        // Copy updated device pointers to device
        CUDA_CALL(cudaMemcpy(d_objects, &h_objects, sizeof(Wall2), cudaMemcpyHostToDevice));
    } else if(hd_objects.count != walls.size()) {
        // Buffer has shrunk, so just update size
        hd_objects.count = walls.size();
        // Copy updated particle count to device (@todo When/where is d_elements allocated??)
        CUDA_CALL(cudaMemcpy(&d_objects->count, &h_objects.count, sizeof(unsigned int), cudaMemcpyHostToDevice));
    }
    // Copy data to device buffers
    // CUDA_CALL(cudaMemcpy(hd_objects.FHydro, h_objects.FHydro, h_objects.count * sizeof(tVect), cudaMemcpyHostToDevice)); // Zero'd before use in streaming()
}
#endif

void LB2::syncDEM(const elmtList &elmts, const particleList &particles, const wallList &walls, const objectList &objects) {
    // Sync DEM data to structure of arrays format (and device memory)
    syncElements<IMPL>(elmts);
    syncParticles<IMPL>(particles);
    syncWalls<IMPL>(walls);
    syncObjects<IMPL>(objects);
}

/**
 * initializeParticleBoundaries()
 */
__host__ __device__ __forceinline__ inline void common_initializeParticleBoundaries(const unsigned int an_i, Node2 *nodes, Particle2 *particles) {
    // Fetch the index of the (active) node being processed
    const unsigned int n_i = nodes->activeI[an_i];
    const tVect node_position = nodes->getPosition(n_i);
    for (int p_i = 0; p_i < particles->count; ++p_i) {
        const tVect convertedPosition = particles->x0[p_i] / PARAMS.unit.Length;
        // @todo pre-compute PARAMS.hydrodynamicRadius / PARAMS.unit.Length ?
        const double convertedRadius = particles->r[p_i] * PARAMS.hydrodynamicRadius / PARAMS.unit.Length;
        if (node_position.insideSphere(convertedPosition, convertedRadius)) { //-0.5?
            nodes->setInsideParticle(n_i, true);
            nodes->solidIndex[n_i] = p_i;
            break;
        }
    }
}
template<>
void LB2::initializeParticleBoundaries<CPU>() {
    // Reset all nodes to outside (e.g. to std::numeric_limits<unsigned int>::max())
    memset(hd_nodes.solidIndex, 0xffffffff, h_nodes.count * sizeof(unsigned int));

    // @todo can we parallelise at a higher level?
    #pragma omp parallel for
    for (size_t an_i = 0; an_i < d_nodes->activeCount; ++an_i) {
        // Pass the active node index to the common implementation
        common_initializeParticleBoundaries(an_i, d_nodes, d_particles);
    }
}
#ifdef USE_CUDA
__global__ void d_initializeParticleBoundaries(Node2 *d_nodes, Particle2 *d_particles) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int an_i = blockIdx.x * blockDim.x + threadIdx.x;
    // Kill excess threads early
    if (an_i >= d_nodes->activeCount) return;
    // Pass the active node index to the common implementation
    common_initializeParticleBoundaries(an_i, d_nodes, d_particles);
}
template<>
void LB2::initializeParticleBoundaries<CUDA>() {
    // Reset all nodes to outside (e.g. to std::numeric_limits<unsigned int>::max())
    CUDA_CALL(cudaMemset(hd_nodes.solidIndex, 0xffffffff, h_nodes.count * sizeof(unsigned int)));
    
    // Launch cuda kernel to update
    // @todo Try unrolling this, so 1 thread per node+particle combination (2D launch?)
    int blockSize = 0;  // The launch configurator returned block size
    int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
    int gridSize = 0;  // The actual grid size needed, based on input size
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_initializeParticleBoundaries, 0, h_nodes.activeCount);
    // Round up to accommodate required threads
    gridSize = (h_nodes.activeCount + blockSize - 1) / blockSize;
    d_initializeParticleBoundaries<<<gridSize, blockSize>>>(d_nodes, d_particles);
    CUDA_CHECK();
}
#endif

/**
 * findNewActive()
 */
__host__ __device__ __forceinline__ inline void common_findNewActive(const unsigned int an_i, Node2 *nodes, Particle2 *particles, Element2 *elements) {
    // Fetch the index of the (active) node being processed
    const unsigned int n_i = nodes->activeI[an_i];
    if (nodes->p[n_i]) {
        const tVect nodePosition = nodes->getPosition(n_i);
        // solid index to identify cluster
        const unsigned int particleIndex = nodes->solidIndex[n_i];
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
        nodes->setInsideParticle(n_i, false);
    }
}
template<>
void LB2::findNewActive<CPU>() {
    // @todo can we parallelise at a higher level?
#pragma omp parallel for
    for (unsigned int an_i = 0; an_i < d_nodes->activeCount; ++an_i) {
        // Pass the active node index to the common implementation
        common_findNewActive(an_i, d_nodes, d_particles, d_elements);
    }
}
#ifdef USE_CUDA
__global__ void d_findNewActive(Node2* d_nodes, Particle2* d_particles, Element2* d_elements) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int an_i = blockIdx.x * blockDim.x + threadIdx.x;
    // Kill excess threads early
    if (an_i >= d_nodes->activeCount) return;
    // Pass the active node index to the common implementation
    common_findNewActive(an_i, d_nodes, d_particles, d_elements);
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
    d_findNewActive<<<gridSize, blockSize>>>(d_nodes, d_particles, d_elements);
    CUDA_CHECK();
}
#endif

/**
 * findNewSolid()
 */
__host__ __device__ __forceinline__ inline void common_findNewSolid(const unsigned int an_i, Node2 *nodes, Particle2 *particles, Element2 *elements) {
    const unsigned int a_i = nodes->activeI[an_i];
    if (nodes->isInsideParticle(a_i)) {  // If node is inside particle
        // solid index to identify cluster
        const unsigned int particleIndex = nodes->solidIndex[a_i];
        const unsigned int clusterIndex = particles->clusterIndex[particleIndex];
        // cycle through first neighbors
        const unsigned int nodeCount = nodes->count;
        for (int k = 1; k < lbmMainDirec; ++k) {
            const unsigned int l_i = nodes->d[nodeCount * k + a_i];
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
                            nodes->setInsideParticle(l_i, true);
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
    for (unsigned int an_i = 0; an_i < d_nodes->activeCount; ++an_i) {
        // Pass the active node index to the common implementation
        common_findNewSolid(an_i, d_nodes, d_particles, d_elements);
    }
}

#ifdef USE_CUDA
__global__ void d_findNewSolid(Node2 *d_nodes, Particle2 *d_particles, Element2 *d_elements) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int an_i = blockIdx.x * blockDim.x + threadIdx.x;
    // Kill excess threads early
    if (an_i >= d_nodes->activeCount) return;
    // Pass the active node index to the common implementation
    common_findNewSolid(an_i, d_nodes, d_particles, d_elements);
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
    d_findNewSolid<<<gridSize, blockSize>>>(d_nodes, d_particles, d_elements);
    CUDA_CHECK();
}
#endif

/**
 * checkNewInterfaceParticles()
 */
__host__ __device__ __forceinline__ inline void common_checkNewInterfaceParticles(const unsigned int e_i, Node2 *nodes, Particle2 *particles, Element2 *elements) {
    // INITIAL PARTICLE POSITION ////////////////////////
    if (elements->FHydro[e_i].norm2() == 0.0) {
        const unsigned int first_component = elements->componentsIndex[e_i];
        const unsigned int last_component = elements->componentsIndex[e_i + 1];
        for (unsigned int n = first_component; n < last_component; ++n) {
            const unsigned int componentIndex = elements->componentsData[n];
            const tVect convertedPosition = particles->x0[componentIndex] / PARAMS.unit.Length;
            // @todo pre-compute PARAMS.hydrodynamicRadius / PARAMS.unit.Length ?
            const double convertedRadius = particles->r[componentIndex] * PARAMS.hydrodynamicRadius / PARAMS.unit.Length;
            for (int i_i = 0; i_i < nodes->interfaceCount; ++i_i) {
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
__global__ void d_checkNewInterfaceParticles(Node2 *d_nodes, Particle2 *d_particles, Element2 *d_elements) {
    // Get unique CUDA thread index, which corresponds to element 
    const unsigned int e_i = blockIdx.x * blockDim.x + threadIdx.x;
    // Kill excess threads early
    if (e_i >= d_elements->count) return;
    // Pass the active node index to the common implementation
    common_findNewSolid(e_i, d_nodes, d_particles, d_elements);
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
    d_checkNewInterfaceParticles<<<gridSize, blockSize>>>(d_nodes, d_particles, d_elements);
    CUDA_CHECK();
}
#endif

/**
 * reconstruct()
 * computeHydroForces()
 * collision()
 */
__host__ __device__ __forceinline__ inline void common_computeHydroForces(const unsigned int an_i, Node2 *nodes, Particle2 *particles, Element2 *elements) {
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
        const tVect diffVel = nodes->age[an_i] * nodes->age[an_i] * nodes->liquidFraction(an_i)*(nodes->u[an_i] - localVel);

        // force on fluid
        nodes->hydroForce[an_i] += -1.0 * diffVel;

        // force on particle
#ifdef __CUDA_ARCH__
        // CUDA atomics
        atomicAdd(&elements->fluidVolume[clusterIndex], nodes->mass[an_i]);
        atomicAdd(&elements->FHydro[clusterIndex], 1.0 * diffVel);
        atomicAdd(&elements->MHydro[clusterIndex], 1.0 * radius.cross(diffVel));        
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
#pragma omp parallel for
    for (unsigned int e_i = 0; e_i < d_elements->count; ++e_i) {
        common_checkNewInterfaceParticles(e_i, d_nodes, d_particles, d_elements);
    }
}
#ifdef USE_CUDA
__global__ void d_reconstructHydroCollide(Node2 *d_nodes, Particle2 *d_particles, Element2 *d_elements) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int an_i = blockIdx.x * blockDim.x + threadIdx.x;

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
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_initializeParticleBoundaries, 0, h_nodes.activeCount);
    // Round up to accommodate required threads
    gridSize = (h_nodes.activeCount + blockSize - 1) / blockSize;
    d_reconstructHydroCollide<<<gridSize, blockSize>>>(d_nodes, d_particles, d_elements);
    CUDA_CHECK();
}

__host__ __device__ void common_streaming(const unsigned int an_i, Node2 *nodes, Wall2 *walls) {

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
    for (int j = 1; j < lbmDirec; ++j) {
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
                    atomicAdd(&walls->FHydro[solidIndex], BBforce);
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
                    nodes->f[A_OFFSET + opp[j]] = -nodes->fs[A_OFFSET + j];
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
                    atomicAdd(&walls->FHydro[solidIndex], BBforce);
#else
                    // CPU atomics
                    #pragma omp atomic update
                    walls->FHydro[solidIndex] += BBforce;
#endif
                    nodes->f[A_OFFSET + opp[j]] = nodes->fs[A_OFFSET + j] - BBi;
                    // adding the extra mass to the surplus
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
                    atomicAdd(&walls->FHydro[solidIndex], BBforce);
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
                        } else if (!active1 && active2) {
                            // second
                            nodes->f[A_OFFSET + opp[j]] = S1 * nodes->fs[nodeCheck2 * lbmDirec + slip2[j]] + S2 * nodes->fs[A_OFFSET + j];
                        } else {
                            // standard BB
                            nodes->f[A_OFFSET + opp[j]] = nodes->fs[A_OFFSET + j];
                        }
                    } else {
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
                            // adding the extra mass to the surplus
                            // extraMass += S2 * nodes->mass[an_i] * BBi;  // redistributeMass() currently not used, so this isn't implemented properly
                        } else if (!active1 && active2) {
                            // second
                            nodes->f[A_OFFSET + opp[j]] = S1 * nodes->fs[nodeCheck2 * lbmDirec + slip2[j]] + S2 * (nodes->fs[A_OFFSET + j] - BBi);
                            // adding the extra mass to the surplus
                            // extraMass += S2 * nodes->mass[an_i] * BBi;  // redistributeMass() currently not used, so this isn't implemented properly
                        } else {
                            // standard BB
                            nodes->f[A_OFFSET + opp[j]] = nodes->fs[A_OFFSET + j] - BBi;
                            // adding the extra mass to the surplus
                            // extraMass += nodes->mass[an_i] * BBi;  // redistributeMass() currently not used, so this isn't implemented properly
                        }
                    } else {
                        // standard BB
                        nodes->f[A_OFFSET + opp[j]] = nodes->fs[A_OFFSET + j] - BBi;
                        // adding the extra mass to the surplus
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
                    printf("(%f, %f, %f) %s TYPE ERROR:\n", pos.x, pos.y, pos.z, typeString(nodes->type[an_i]));
                    for (int j = 1; j < lbmDirec; ++j) {
                        printf("before error: j=%d link=%u\n", j, nodes->coord[nodes->d[j * nodes->count + an_i]]);
                    }
                    pos = nodes->getPosition(ln_i);
                    printf("(%f, %f, %f) %s TYPE ERROR\n", pos.x, pos.y, pos.z, typeString(nodes->type[ln_i]));
                    // @todo aborting from CUDA is harder, especially if the printf() is to be saved
#ifdef __CUDA_ARCH__
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
    for (unsigned int an_i = 0; an_i < d_nodes->activeCount; ++an_i) {
        common_streaming(an_i, d_nodes, d_walls);
    }

    // redistributing extra mass due to bounce back to interface cells
    // redistributeMass(extraMass);  // extraMass hasn't been implemented properly
}
#ifdef USE_CUDA
__global__ void d_streaming(Node2 *d_nodes, Wall2 *d_walls) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int an_i = blockIdx.x * blockDim.x + threadIdx.x;
    
    common_streaming(an_i, d_nodes, d_walls);
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
    d_streaming<<<gridSize, blockSize>>>(d_nodes, d_walls);
    CUDA_CHECK();

#ifdef _DEBUG
    CUDA_CALL(cudaMemcpy(h_nodes.f, hd_nodes.f, sizeof(double) * h_nodes.activeCount * lbmDirec, cudaMemcpyDeviceToHost));
    for (unsigned int in = 0; in < h_nodes.activeCount; ++in) {        
        for (unsigned int j = 1; j < lbmDirec; ++j) {                
            if (h_nodes.f[in * lbmDirec + j] == 0) {
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
    for (int i = 0; i < d_elements->count; ++i) {
        d_elements->FHydro[i] *= PARAMS.unit.Force;
        d_elements->MHydro[i] *= PARAMS.unit.Torque;
        d_elements->fluidVolume[i] *= PARAMS.unit.Volume;
    }
    for (int i = 0; i < d_walls->count; ++i) {
        d_walls->FHydro[i] *= PARAMS.unit.Force;
    }
    for (int i = 0; i < d_objects->count; ++i) {
        d_objects->FHydro[i] *= PARAMS.unit.Force;
    }
}
#ifdef USE_CUDA
__global__ void d_shiftToPhysical(Element2 *d_elements, Wall2 *d_walls, Object2* d_objects) {
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
    d_shiftToPhysical<<<gridSize, blockSize>>>(d_elements, d_walls, d_objects);
    CUDA_CHECK();
}
#endif


void LB2::latticeBoltzmannCouplingStep(bool &newNeighbourList) {
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
    } else {
        // SOLID TO ACTIVE CHECK
        // @note Calling this directly after initializeParticleBoundaries() is redundant, hence else
        this->findNewActive<IMPL>();
    }

    // ACTIVE TO SOLID CHECK
    this->findNewSolid<IMPL>();

    if (freeSurface) {
        this->checkNewInterfaceParticles<IMPL>();
    }
}
void LB2::latticeBoltzmannStep() {
    // Reconstruct active list
    hd_nodes.cleanLists<IMPL>();

    // Initializing the elements forces (lattice units)
    hd_elements.initElements<IMPL>();

    // Initialise lattice boltzmann force vector
    // @note currently ignored, doesn't seem fully integrated with model
    //if (!h_params.forceField) {
    //    lbF.reset(); // @todo, does this exist on host or device
    //}

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

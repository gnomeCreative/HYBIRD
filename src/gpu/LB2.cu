#include "LB2.h"

#include <cstdlib>

#include "cuda_helper.h"

#include "DEM.h"

void LB2::step(const DEM &dem, bool io_demSolver) {
    if (io_demSolver) {
        this->latticeBoltzmannCouplingStep(dem.newNeighborList, dem.elmts, dem.particles);
    }

    if (dem.demTime >= dem.demInitialRepeat) {
        this->latticeBolzmannStep(dem.elmts, dem.particles, dem.walls, dem.objects);

        // Lattice Boltzmann core steps
        if (this->freeSurface) {
            this->latticeBoltzmannFreeSurfaceStep();
        }
    }
}

template<>
void LB2::syncParticles<CPU>(const particleList& particles) {
    if (h_particles.count < particles.size()) {
        // Grow host buffers
        if (h_particles.clusterIndex) {
            free(h_particles.clusterIndex);
            free(h_particles.r);
            free(h_particles.x0);
        }
        h_particles.clusterIndex = (unsigned int*)malloc(particles.size() * sizeof(unsigned int));
        h_particles.r = (double*)malloc(particles.size() * sizeof(double));
        h_particles.x0 = (tVect*)malloc(particles.size() * sizeof(tVect));
    }
    // Update size
    h_particles.count = particles.size();
    // Repackage host particle data from array of structures, to structure of arrays
    for (int i = 0; i < h_particles.count; ++i) {
        h_particles.clusterIndex[i] = particles[i].clusterIndex;
        h_particles.r[i] = particles[i].r;
        h_particles.x0[i] = particles[i].x0;
    }
}
template<>
bool LB2::syncElements<CPU>(const elmtList& elements) {
    bool componentsHasGrown = false;
    if (h_elements.count < elements.size()) {
        // Grow host buffers
         if (h_elements.FHydro) {
             free(h_elements.FHydro);
         }
         h_elements.FHydro = (tVect*)malloc(elements.size() * sizeof(tVect));
    }
    // Update size
    h_elements.count = elements.size();
    // Repackage host particle data from array of structures, to structure of arrays
     for (int i = 0; i < h_elements.count; ++i) {
         h_elements.FHydro[i] = elements[i].FHydro;
     }
    // Construct the components storage
    {
        // Allocate memory for componentsData
        size_t totalComponents = 0;
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
#ifdef USE_CUDA
template<>
void LB2::syncParticles<CUDA>(const particleList& particles) {
    // Copy latest particle data from HOST DEM to the device
    // @todo Can these copies be done ahead of time async?
    // @todo These copies will be redundant when DEM is moved to CUDA
    this->syncParticles<CPU>(particles);
    if (hd_particles.count < particles.size()) {
        // Grow device buffers
        if (hd_particles.r) {
            CUDA_CALL(cudaFree(hd_particles.clusterIndex));
            CUDA_CALL(cudaFree(hd_particles.r));
            CUDA_CALL(cudaFree(hd_particles.x0));
        }
        CUDA_CALL(cudaMalloc(&hd_particles.clusterIndex, h_particles.count * sizeof(unsigned int)));
        CUDA_CALL(cudaMalloc(&hd_particles.r, h_particles.count * sizeof(double)));
        CUDA_CALL(cudaMalloc(&hd_particles.x0, h_particles.count * sizeof(tVect)));
        hd_particles.count = h_particles.count;
        // Copy updated device pointers to device (@todo When/where is d_particles allocated??)
        CUDA_CALL(cudaMemcpy(d_particles, &h_particles, sizeof(Particle2), cudaMemcpyHostToDevice));
    } else if(hd_particles.count != particles.size()) {
        // Buffer has shrunk, so just update size
        hd_particles.count = particles.size();
        // Copy updated particle count to device
        CUDA_CALL(cudaMemcpy(&d_particles->count, &h_particles.count, sizeof(size_t), cudaMemcpyHostToDevice));
    }
    // Copy data to device buffers
    CUDA_CALL(cudaMemcpy(hd_particles.clusterIndex, h_particles.clusterIndex, h_particles.count * sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_particles.r, h_particles.r, h_particles.count * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_particles.x0, h_particles.x0, h_particles.count * sizeof(tVect), cudaMemcpyHostToDevice));
}
template<>
bool LB2::syncElements<CUDA>(const elmtList& elements) {
    // @todo copy h_elements to d_elements
    // Copy latest particle data from HOST DEM to the device
    // @todo Can these copies be done ahead of time async?
    // @todo These copies will be redundant when DEM is moved to CUDA
    bool componentsHasGrown = this->syncElements<CPU>(elements);
    bool updateDeviceStruct = false;
    if (hd_elements.count < elements.size()) {
        if (hd_elements.FHydro) {
            CUDA_CALL(cudaFree(hd_elements.FHydro));
        }
        // Initially allocate device buffers except components
        CUDA_CALL(cudaMalloc(&hd_elements.FHydro, hd_elements.count * sizeof(tVect)));
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
    if (updateDeviceStruct) {
        // Buffer has shrunk, so just update size
        hd_elements.count = elements.size();
        // Copy updated device pointers to device (@todo When/where is d_elements allocated??)
        CUDA_CALL(cudaMemcpy(d_elements, &h_elements, sizeof(Element2), cudaMemcpyHostToDevice));
    }
    // Copy data to device buffers
    CUDA_CALL(cudaMemcpy(hd_elements.FHydro, &h_elements.FHydro, h_elements.count * sizeof(tVect), cudaMemcpyHostToDevice)); // @todo does this value ever change?
    CUDA_CALL(cudaMemcpy(hd_elements.componentsIndex, h_elements.componentsIndex, (h_elements.count + 1) * sizeof(unsigned int), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(hd_elements.componentsData, h_elements.componentsData, h_elements.componentsIndex[h_elements.count] * sizeof(unsigned int), cudaMemcpyHostToDevice));
    return componentsHasGrown;
}
#endif

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
__global__ void d_checkNewInterfaceParticles(Node2* d_nodes, Particle2* d_particles, Element2* d_elements) {
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

void LB2::latticeBoltzmannCouplingStep(bool &newNeighbourList, const elmtList& elmts, const particleList& particles) {
    // Sync DEM data to structure of arrays format (and device memory)
    syncParticles<IMPL>(particles);
    syncElements<IMPL>(elmts);
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

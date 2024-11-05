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
#ifndef USE_CUDA
template<>
void LB2::initializeParticleBoundaries<CPU>() {
    // Reset all nodes to outside (std::numeric_limits<unsigned short>::max())
    // (The type short is 2 bytes long, but memset requires 4 bytes, so we pass 0xff..)
    memset(d_nodes->solidIndex, 0xffffffff, h_nodes.count * sizeof(unsigned short));
    
    for (int p_i = 0; p_i < d_particles.count; ++i) {
        const tVect convertedPosition = d_particles->x0[p_i] / PARAMS.unit.Length;
        const double convertedRadius = PARAMS.hydrodynamicRadius * d_particles->r[p_i] / PARAMS.unit.Length;
#pragma omp parallel for
        for (size_t i = 0; i < d_nodes.activeCount; ++i) {
            // checking if node is inside a particle
            const unsigned int a_i = d_nodes->activeI[i];
            const tVect node_position = nodes->getPosition(a_i);
            if (node_position.insideSphere(convertedPosition, convertedRadius)) { //-0.5?
                d_nodes->setInsideParticle(a_i, true);
                d_nodes->solidIndex[a_i] = [p_i];
            }
        }
    }
}
#else
__global__ void d_initializeParticleBoundaries(Node2 *d_nodes, Particle2 *d_particles) {
    const unsigned int n_i = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int a_i = d_nodes->activeI[n_i];
    const tVect node_position = d_nodes->getPosition(a_i);
    for (int p_i = 0; p_i < d_particles->count; ++p_i) {
        const tVect convertedPosition = d_particles->x0[p_i] / PARAMS.unit.Length;
        // @todo pre-compute PARAMS.hydrodynamicRadius / PARAMS.unit.Length ?
        const double convertedRadius = d_particles->r[p_i] * PARAMS.hydrodynamicRadius / PARAMS.unit.Length;
        if (node_position.insideSphere(convertedPosition, convertedRadius)) { //-0.5?
            d_nodes->setInsideParticle(a_i, true);
            d_nodes->solidIndex[a_i] = p_i;
            break;
        }
    }
}
template<>
void LB2::initializeParticleBoundaries<CUDA>() {
    // Reset all nodes to outside (std::numeric_limits<unsigned short>::max())
    // (The type short is 2 bytes long, but memset requires 4 bytes, so we pass 0xff..)
    CUDA_CALL(cudaMemset(hd_nodes.solidIndex, 0xffffffff, h_nodes.count * sizeof(unsigned short)));
    
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
#ifndef USE_CUDA
template<>
void LB2::findNewActive<CPU>() {
    // SOLID TO ACTIVE CHECK
    // cycling through particle nodes
    for (size_t i = 0; i < h_nodes.activeCount; ++i) {
        const unsigned int a_i = d_nodes->activeI[i];
        if (d_nodes->isInsideParticle(a_i)) {  // If node is inside particle
            const tVect nodePosition = d_nodes->getPosition(a_i);
            // solid index to identify cluster
            const unsigned int particleIndex = d_nodes->solidIndex[a_i];
            const unsigned int clusterIndex = d_particles->clusterIndex[particleIndex];
            // in this case check if it has been uncovered (must be out of all particles of the cluster) - we start with a true hypothesis
            // cycling through component particles
            const unsigned int first_component = d_elements->componentsIndex[clusterIndex];
            const unsigned int last_component = d_elements->componentsIndex[clusterIndex + 1];
            for (unsigned int j = first_component; j < last_component; ++j) {
                // getting indexes from particle composing the cluster
                const unsigned int componentIndex = d_elements->componentsData[j];
                // checking if it has been uncovered in component j of the cluster
                // radius need to be increased by half a lattice unit
                // this is because solid boundaries are located halfway between solid and fluid nodes
                // @todo pre-compute PARAMS.hydrodynamicRadius / PARAMS.unit.Length?
                if (nodePosition.insideSphere(particles[componentIndex].x0 / PARAMS.unit.Length, particles[componentIndex].r * PARAMS.hydrodynamicRadius / PARAMS.unit.Length)) { //-0.5?
                    // if the node is still inside the element, the hypothesis of new active is not true anymore
                    d_nodes->isInsideParticle(a_i) = false;
                    // and we can get out of the cycle
                    break;
                }
            }
        }
    }
}
#else
__global__ void d_findNewActive(unsigned int threadCount, Node2* d_nodes, Particle2* d_particles, Element2* d_elements) {
    const unsigned int n_i = blockIdx.x * blockDim.x + threadIdx.x;
    if (n_i >= threadCount) return;
    const unsigned int a_i = d_nodes->activeI[n_i];
    if (d_nodes->p[a_i]) {
        const tVect nodePosition = d_nodes->getPosition(a_i);
        // solid index to identify cluster
        const unsigned int particleIndex = d_nodes->solidIndex[a_i];
        const unsigned int clusterIndex = d_particles->clusterIndex[particleIndex];
        // in this case check if it has been uncovered (must be out of all particles of the cluster) - we start with a true hypothesis
        // cycling through component particles
        const unsigned int first_component = d_elements->componentsIndex[clusterIndex];
        const unsigned int last_component = d_elements->componentsIndex[clusterIndex + 1];
        for (unsigned int j = first_component; j < last_component; ++j) {
            // getting indexes from particle composing the cluster
            const unsigned int componentIndex = d_elements->componentsData[j];
            // checking if it has been uncovered in component j of the cluster
            // radius need to be increased by half a lattice unit
            // this is because solid boundaries are located halfway between solid and fluid nodes
            const tVect convertedPosition = d_particles->x0[componentIndex] / PARAMS.unit.Length;
            // @todo pre-compute PARAMS.hydrodynamicRadius / PARAMS.unit.Length ?
            const double convertedRadius = d_particles->r[componentIndex] * PARAMS.hydrodynamicRadius / PARAMS.unit.Length;
            if (nodePosition.insideSphere(convertedPosition, convertedRadius)){ //-0.5?
                // if the node is still inside the element, the hypothesis of new active is not true anymore
                // and we can get out of the cycle
                return;
            }
        }
        // turning up the cell as we didn't exit early
        d_nodes->setInsideParticle(a_i, false);
    }
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
    d_findNewActive<<<gridSize, blockSize>>>(h_nodes.activeCount, d_nodes, d_particles, d_elements);
    CUDA_CHECK();
}
#endif
#ifndef USE_CUDA
template<>
void LB2::findNewSolid<CPU>() {
    // ACTIVE TO SOLID CHECK
    // we check only first order neighbors of particle nodes. This is done in order to avoid cycling through all active cells
    for (size_t i = 0; i < h_nodes.activeCount; ++i) {
        const unsigned int a_i = d_nodes->activeI[i];
        if (d_nodes->isInsideParticle(a_i)) {  // If node is inside particle
            // solid index to identify cluster
            const unsigned int particleIndex = d_nodes->solidIndex[a_i];
            const unsigned int clusterIndex = d_particles->clusterIndex[particleIndex];
            // cycle through first neighbors
            const unsigned int nodeCount = d_nodes->count;
            for (int k = 1; k < lbmMainDirec; ++k) {
                const unsigned int l_i = d_nodes->d[nodeCount * k + a_i];
                if (l_i != std::numeric_limits<unsigned int>::max()) {
                    // checking if solid particle is close to an active one -> we have an active node to check
                    if (!d_nodes->isInsideParticle(l_i) && d_nodes->isActive(l_i)) {
                        const tVect linkPosition = d_nodes->getPosition(l_i);
                        // check if neighbors has been covered (by any of the particles of the cluster) - we start with a false hypothesis
                        // cycling through all components of the cluster
                        const unsigned int first_component = d_elements->componentsIndex[clusterIndex];
                        const unsigned int last_component = d_elements->componentsIndex[clusterIndex + 1];
                        for (unsigned int j = first_component; j < last_component; ++j) {
                            // getting component particle index
                            unsigned short componentIndex = d_elements->componentsData[j];
                            // check if it getting inside
                            // radius need to be increased by half a lattice unit
                            // this is because solid boundaries are located halfway between solid and fluid nodes
                            // @todo pre-compute PARAMS.hydrodynamicRadius / PARAMS.unit.Length ?
                            if (linkPosition.insideSphere(d_particles->x0[componentIndex] / PARAMS.unit.Length, d_particles->r[componentIndex] * PARAMS.hydrodynamicRadius / PARAMS.unit.Length)) { //-0.5?
                                // if so, then the false hypothesis does not hold true anymore
                                d_nodes->solidIndex[l_i] = componentIndex;
                                // By setting particle to inside, it won't be checked again, newSolidNodes hence becomes redundant
                                d_nodes->setInsideParticle(l_i, true);
                                // and we exit the cycle
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
}
#else
__global__ void d_findNewSolid(unsigned int threadCount, Node2* d_nodes, Particle2* d_particles, Element2* d_elements) {
    const unsigned int n_i = blockIdx.x * blockDim.x + threadIdx.x;
    if (n_i >= threadCount) return;
    const unsigned int a_i = d_nodes->activeI[n_i];
    if (d_nodes->isInsideParticle(a_i)) {  // If node is inside particle
        // solid index to identify cluster
        const unsigned int particleIndex = d_nodes->solidIndex[a_i];
        const unsigned int clusterIndex = d_particles->clusterIndex[particleIndex];
        // cycle through first neighbors
        const unsigned int nodeCount = d_nodes->count;
        for (int k = 1; k < lbmMainDirec; ++k) {
            const unsigned int l_i = d_nodes->d[nodeCount * k + a_i];
            if (l_i != std::numeric_limits<unsigned int>::max()) {
                // checking if solid particle is close to an active one -> we have an active node to check
                if (!d_nodes->isInsideParticle(l_i) && d_nodes->isActive(l_i)) {
                    const tVect linkPosition = d_nodes->getPosition(l_i);
                    // check if neighbors has been covered (by any of the particles of the cluster) - we start with a false hypothesis
                    // cycling through all components of the cluster
                    const unsigned int first_component = d_elements->componentsIndex[clusterIndex];
                    const unsigned int last_component = d_elements->componentsIndex[clusterIndex + 1];
                    for (unsigned int j = first_component; j < last_component; ++j) {
                        // getting component particle index
                        unsigned short componentIndex = d_elements->componentsData[j];
                        // check if it getting inside
                        // radius need to be increased by half a lattice unit
                        // this is because solid boundaries are located halfway between solid and fluid nodes
                        // @todo pre-compute PARAMS.hydrodynamicRadius / PARAMS.unit.Length ?
                        if (linkPosition.insideSphere(d_particles->x0[componentIndex] / PARAMS.unit.Length, d_particles->r[componentIndex] * PARAMS.hydrodynamicRadius / PARAMS.unit.Length)) { //-0.5?
                            // if so, then the false hypothesis does not hold true anymore
                            d_nodes->solidIndex[l_i] = componentIndex;
                            // By setting particle to inside, it won't be checked again, newSolidNodes hence becomes redundant
                            d_nodes->setInsideParticle(l_i, true);
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
void LB2::findNewSolid<CUDA>() {
    // Launch cuda kernel to update
    int blockSize = 0;  // The launch configurator returned block size
    int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
    int gridSize = 0;  // The actual grid size needed, based on input size
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_findNewSolid, 0, h_nodes.activeCount);
    // Round up to accommodate required threads
    gridSize = (h_nodes.activeCount + blockSize - 1) / blockSize;
    d_findNewSolid<<<gridSize, blockSize>>>(h_nodes.activeCount, d_nodes, d_particles, d_elements);
    CUDA_CHECK();
}
#endif
#ifndef USE_CUDA
template<>
void LB2::checkNewInterfaceParticles<CPU>() {
    // INITIAL PARTICLE POSITION ////////////////////////
    for (int m = 0; m < d_elements->count; ++m) {
        if (d_elements->FHydro[m].norm2() == 0.0) {
            const unsigned int first_component = d_elements->componentsIndex[m];
            const unsigned int last_component = d_elements->componentsIndex[m + 1];
            for (unsigned int n = first_component; n < last_component; ++n) {
                const unsigned short componentIndex = d_elements->componentsData[n];
                const tVect convertedPosition = d_particles->x0[componentIndex] / PARAMS.unit.Length;
                // @todo pre-compute PARAMS.hydrodynamicRadius / PARAMS.unit.Length ?
                const double convertedRadius = d_particles->r[componentIndex] * PARAMS.hydrodynamicRadius / PARAMS.unit.Length;
#pragma omp parallel for
                for (int it = 0; it < d_nodes->interfaceCount; ++it) {
                    const unsigned int nodeHere = d_nodes->interfaceI[it];
                    if (!d_nodes->isInsideParticle(nodeHere)) {
                        // checking if node is inside a particle
                        const tVect nodePosition = d_nodes->getPosition(nodeHere);
                        if (nodePosition.insideSphere(convertedPosition, convertedRadius)) { //-0.5?
                            d_nodes->setInsideParticle(nodeHere, true);
                            d_nodes->solidIndex[nodeHere] = componentIndex;
                        }
                    }
                }
            }
        }
    }
}
#else
__global__ void d_checkNewInterfaceParticles(unsigned int threadCount, Node2* d_nodes, Particle2* d_particles, Element2* d_elements) {
    const unsigned int m = blockIdx.x * blockDim.x + threadIdx.x;
    if (m >= threadCount) return;

    // INITIAL PARTICLE POSITION ////////////////////////
    if (d_elements->FHydro[m].norm2() == 0.0) {
        const unsigned int first_component = d_elements->componentsIndex[m];
        const unsigned int last_component = d_elements->componentsIndex[m + 1];
        for (unsigned int n = first_component; n < last_component; ++n) {
            const unsigned short componentIndex = d_elements->componentsData[n];
            const tVect convertedPosition = d_particles->x0[componentIndex] / PARAMS.unit.Length;
            // @todo pre-compute PARAMS.hydrodynamicRadius / PARAMS.unit.Length ?
            const double convertedRadius = d_particles->r[componentIndex] * PARAMS.hydrodynamicRadius / PARAMS.unit.Length;
#pragma omp parallel for
            for (int it = 0; it < d_nodes->interfaceCount; ++it) {
                const unsigned int nodeHere = d_nodes->interfaceI[it];
                if (!d_nodes->isInsideParticle(nodeHere)) {
                    // checking if node is inside a particle
                    const tVect nodePosition = d_nodes->getPosition(nodeHere);
                    if (nodePosition.insideSphere(convertedPosition, convertedRadius)) { //-0.5?
                        d_nodes->setInsideParticle(nodeHere, true);
                        d_nodes->solidIndex[nodeHere] = componentIndex;
                    }
                }
            }
        }
    }
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
    d_checkNewInterfaceParticles<<<gridSize, blockSize>>>(h_elements.count, d_nodes, d_particles, d_elements);
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

    // @todo Surely d_nodes->activeI needs to be rebuilt here, after we've updated which particles are inside spheres?

    // ACTIVE TO SOLID CHECK
    this->findNewSolid<IMPL>();

    if (freeSurface) {
        this->checkNewInterfaceParticles<IMPL>();
    }
}

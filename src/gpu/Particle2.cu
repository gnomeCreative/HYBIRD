#include "Particle2.h"

#include "Element2.h"


template<>
inline void Particle2::memoryAlloc<CPU>(unsigned int num) {
    if (!num) return;
    alloc = num;

    if (particleIndex) {
        CUDA_CALL(cudaFree(particleIndex));
        CUDA_CALL(cudaFree(clusterIndex));
        CUDA_CALL(cudaFree(protoIndex));
        CUDA_CALL(cudaFree(active));
        CUDA_CALL(cudaFree(isGhost));
        CUDA_CALL(cudaFree(neighbour_index));
        CUDA_CALL(cudaFree(r));
        CUDA_CALL(cudaFree(x0));
        CUDA_CALL(cudaFree(x1));
        CUDA_CALL(cudaFree(radiusVec));
        // @todo springs
    }
    
    particleIndex = (unsigned int*)malloc(alloc * sizeof(unsigned int));
    clusterIndex = (unsigned int*)malloc(alloc * sizeof(unsigned int));
    protoIndex = (unsigned int*)malloc(alloc * sizeof(unsigned int));
    active = (bool*)malloc(alloc * sizeof(bool));
    isGhost = (bool*)malloc(alloc * sizeof(bool));
    neighbour_index = (unsigned int*)malloc(alloc * sizeof(unsigned int));
    r = (double*)malloc(alloc * sizeof(double));
    x0 = (tVect*)malloc(alloc * sizeof(tVect));
    x1 = (tVect*)malloc(alloc * sizeof(tVect));
    radiusVec = (tVect*)malloc(alloc * sizeof(tVect));
    // @todo springs
    // Init
    memset(particleIndex, 0, alloc * sizeof(unsigned int));
    memset(clusterIndex, 0, alloc * sizeof(unsigned int));
    memset(protoIndex, 0, alloc * sizeof(unsigned int));
    // active?
    std::fill(neighbour_index, neighbour_index + alloc, std::numeric_limits<unsigned int>::max());
    std::fill(isGhost, isGhost + alloc, false);
    memset(r, 0, alloc * sizeof(double));
    memset(x0, 0, alloc * sizeof(tVect));
    memset(x1, 0, alloc * sizeof(tVect));
    memset(radiusVec, 0, alloc * sizeof(tVect));
    // @todo springs
}
template<>
inline void Particle2::memoryAlloc<CUDA>(unsigned int num) {
    if (!num) return;
    alloc = num;

    if (particleIndex) {
        CUDA_CALL(cudaFree(particleIndex));
        CUDA_CALL(cudaFree(clusterIndex));
        CUDA_CALL(cudaFree(protoIndex));
        CUDA_CALL(cudaFree(active));
        CUDA_CALL(cudaFree(isGhost));
        CUDA_CALL(cudaFree(r));
        CUDA_CALL(cudaFree(x0));
        CUDA_CALL(cudaFree(x1));
        CUDA_CALL(cudaFree(radiusVec));
        // @todo springs
    }

    CUDA_CALL(cudaMalloc(&particleIndex, alloc * sizeof(unsigned int)));
    CUDA_CALL(cudaMalloc(&clusterIndex, alloc * sizeof(unsigned int)));
    CUDA_CALL(cudaMalloc(&protoIndex, alloc * sizeof(unsigned int)));
    CUDA_CALL(cudaMalloc(&active, alloc * sizeof(bool)));
    CUDA_CALL(cudaMalloc(&isGhost, alloc * sizeof(bool)));
    CUDA_CALL(cudaMalloc(&r, alloc * sizeof(double)));
    CUDA_CALL(cudaMalloc(&x0, alloc * sizeof(tVect)));
    CUDA_CALL(cudaMalloc(&x1, alloc * sizeof(tVect)));
    CUDA_CALL(cudaMalloc(&radiusVec, alloc * sizeof(tVect)));
    // @todo springs
}

void Particle2::updateCorrected(const unsigned int p_i, const Element2 &elements, const unsigned int e_i) {
    // updating position and velocity for simple case
    x0[p_i] = elements.x0[e_i];
    radiusVec[p_i] = Zero;
    x1[p_i] = elements.x1[e_i];

    if (elements.size[e_i] > 1) {
        x0[p_i] = x0[p_i] + r[p_i] * project(DEM_P.prototypes[elements.size[e_i]][protoIndex[p_i]], elements.q0[e_i]);
        // updating radius (distance of particle center of mass to element center of mass)
        radiusVec[p_i] = x0[p_i] - elements.x0[e_i];
        x1[p_i] = x1[p_i] + elements.wGlobal[e_i].cross(radiusVec[p_i]);
    }

}
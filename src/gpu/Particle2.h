#ifndef PARTICLE2_H
#define PARTICLE2_H
#include "DEMParams.h"
#include "Element2.h"
#include "LBParams.h"

/**
 * \brief Elements within the DEM are decomposed into spherical particles.
 *
 * This leads to simpler maths for elements that would contain hard edges.
 * @todo/note original member vars coord and d are fixed at runtime
 */
struct Particle2 {
    // The total number of active particles
    unsigned int activeCount = 0;
    // The allocated size of activeI
    // We don't shrink the buffer if the number of active particles decreases
    unsigned int activeAlloc = 0;
    // Index of nodes marked as active
    unsigned int* activeI = nullptr;

    // The neighbour grid index
    unsigned int PBM_alloc = 0;
    unsigned int* PBM = nullptr;
    // Index of particle found within neighbour grid
    // This buffer's length matches the number of particles
    unsigned int *neighbour_index = nullptr;

    // The total number of particles
    // @note This is variable at runtime
    unsigned int count = 0;
    // Size of allocated buffer, as it never shrinks
    unsigned int alloc = 0;
    
    // belonging element index
    unsigned  *particleIndex = nullptr;
    // belonging element index
    unsigned int *clusterIndex = nullptr;
    // belonging element index
    unsigned int *protoIndex = nullptr;
    // is it active? (false=destroyed)
    bool *active = nullptr;
    // is it a ghost?
    bool *isGhost = nullptr;
    // belonging cell for neighbor list
    //unsigned int *tableCell = nullptr;
    // particle radius
    double *r = nullptr;
    // position of the particle
    tVect *x0 = nullptr;
    // velocity of the particle
    tVect *x1 = nullptr;
    // vector connecting center of element to center of particle
    tVect *radiusVec = nullptr;
    //@todo how do springs look? How will they exist on gpu
    // springs[4]??
    /**
     * Does not update count, or allocate componentsIndex/componentsData
     * Initialises based on original constructor
     */
    template<int impl>
    void memoryAlloc(unsigned int num);

    /**
     * \brief 
     * @param i Index of the particle to operate on
     * @param elements elements structure
     * @param e_i Index of mother element
     */
    __host__ __device__ void updateCorrected(unsigned int i, const Element2 *elements);


    __host__ __device__ void updatePredicted(unsigned int i, const Element2 *elements);

};


__host__ __device__ __forceinline__ void Particle2::updateCorrected(const unsigned int i, const Element2 *elements) {
    // Mother element index
    const unsigned int e_i = clusterIndex[i];
    // updating position and velocity for simple case
    x0[i] = elements->x0[e_i];
    radiusVec[i] = {0,0,0};
    x1[i] = elements->x1[e_i];

    if (elements->size[e_i] > 1) {
        x0[i] += r[i] * project(prototypes[elements->size[e_i]][protoIndex[i]], elements->q0[e_i]);
        // Virtual distance, incase element is wrapped before particle
        tVect e_x0 = elements->x0[e_i];
        tVect ae_x0 = e_x0 - x0[i];
        e_x0 =
        {
            abs(ae_x0.x) > DEM_P.half_demSize.x ? e_x0.x - (ae_x0.x / abs(ae_x0.x) * DEM_P.demSize.x) : e_x0.x,
            abs(ae_x0.y) > DEM_P.half_demSize.y ? e_x0.y - (ae_x0.y / abs(ae_x0.y) * DEM_P.demSize.y) : e_x0.y,
            abs(ae_x0.z) > DEM_P.half_demSize.z ? e_x0.z - (ae_x0.z / abs(ae_x0.z) * DEM_P.demSize.z) : e_x0.z,
        };

        // updating radius (distance of particle center of mass to element center of mass)
        radiusVec[i] = x0[i] - e_x0;
        x1[i] += elements->wGlobal[e_i].cross(radiusVec[i]);
    }

    // Wrap particle location if periodic boundaries
    if (LB_P.boundary[0] == PERIODIC && x0[i].x < 0) {
        x0[i].x += DEM_P.demSize.x;
    } else if (LB_P.boundary[1] == PERIODIC && x0[i].x >= DEM_P.demSize.x) {
        x0[i].x -= DEM_P.demSize.x;
    }
    if (LB_P.boundary[2] == PERIODIC && x0[i].y < 0) {
        x0[i].y += DEM_P.demSize.y;
    } else if (LB_P.boundary[3] == PERIODIC && x0[i].y >= DEM_P.demSize.y) {
        x0[i].y -= DEM_P.demSize.y;
    }
    if (LB_P.boundary[4] == PERIODIC && x0[i].z < 0) {
        x0[i].z += DEM_P.demSize.z;
    } else if (LB_P.boundary[5] == PERIODIC && x0[i].z >= DEM_P.demSize.z) {
        x0[i].z -= DEM_P.demSize.z;
    }
}
__host__ __device__ __forceinline__ void Particle2::updatePredicted(const unsigned int i, const Element2 *elements) {
    // Mother element index
    const unsigned int me_i = clusterIndex[i];
    // updating position and velocity for simple case
    x0[i] = elements->xp0[me_i];
    radiusVec[i] = {0, 0, 0};
    x1[i] = elements->xp1[me_i];

    if (elements->size[me_i] > 1) {
        x0[i] = x0[i] + r[i] * project(prototypes[elements->size[me_i]][protoIndex[i]], elements->qp0[me_i]);
        // Virtual distance, incase element is wrapped before particle
        tVect e_x0 = elements->x0[me_i];
        tVect ae_x0 = e_x0 - x0[i];
        e_x0 =
        {
            abs(ae_x0.x) > DEM_P.half_demSize.x ? e_x0.x - (ae_x0.x / abs(ae_x0.x) * DEM_P.demSize.x) : e_x0.x,
            abs(ae_x0.y) > DEM_P.half_demSize.y ? e_x0.y - (ae_x0.y / abs(ae_x0.y) * DEM_P.demSize.y) : e_x0.y,
            abs(ae_x0.z) > DEM_P.half_demSize.z ? e_x0.z - (ae_x0.z / abs(ae_x0.z) * DEM_P.demSize.z) : e_x0.z,
        };

        // updating radius (distance of particle center of mass to element center of mass)
        radiusVec[i] = x0[i] - e_x0;
        x1[i] = x1[i] + elements->wpGlobal[me_i].cross(radiusVec[i]);
    }
}


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
#endif  // PARTICLE2_H

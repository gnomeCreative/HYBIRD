#ifndef PARTICLE2_H
#define PARTICLE2_H
#include "DEMParams.h"
#include "Element2.h"

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
     * @param p_i Index of the particle to operate on
     * @param elements elements structure
     * @param e_i Index of mother element
     */
    void updateCorrected(unsigned int p_i, const Element2& elements, unsigned int e_i);


    __host__ __device__ __forceinline__ void updatePredicted(unsigned int i, const Element2* elements);

};
__host__ __device__ __forceinline__ void Particle2::updatePredicted(const unsigned int i, const Element2 *elements) {
    // Mother element index
    const unsigned int me_i = clusterIndex[i];
    // updating position and velocity for simple case
    x0[i] = elements->xp0[me_i];
    radiusVec[i] = Zero;
    x1[i] = elements->xp1[me_i];

    if (elements->size[me_i] > 1) {
        x0[i] = x0[i] + r[i] * project(DEM_P.prototypes[elements->size[me_i]][protoIndex[i]], elements->qp0[me_i]);
        // updating radius (distance of particle center of mass to element center of mass)
        radiusVec[i] = x0[i] - elements->xp0[me_i];
        x1[i] = x1[i] + elements->wpGlobal[me_i].cross(radiusVec[i]);
    }
}
#endif  // PARTICLE2_H

#ifndef PARTICLE2_H
#define PARTICLE2_H
#include "DEMParams.h"

/**
 * \brief Elements within the DEM are decomposed into spherical particles.
 *
 * This leads to simpler maths for elements that would contain hard edges.
 * @todo/note original member vars coord and d are fixed at runtime
 */
struct Particle2 {
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
    unsigned int *tableCell = nullptr;
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
     * @param dem_p DEM parameters
     */
    void updateCorrected(const unsigned int p_i, const Element2& elements, const unsigned int e_i, const DEMParams& dem_p);

};

template<>
inline void Particle2::memoryAlloc<CPU>(unsigned int num) {
    alloc = num;
    if (!num) return;
    
    particleIndex = (unsigned int*)malloc(alloc * sizeof(unsigned int));
    clusterIndex = (unsigned int*)malloc(alloc * sizeof(unsigned int));
    protoIndex = (unsigned int*)malloc(alloc * sizeof(unsigned int));
    active = (bool*)malloc(alloc * sizeof(bool));
    isGhost = (bool*)malloc(alloc * sizeof(bool));
    tableCell = (unsigned int*)malloc(alloc * sizeof(unsigned int));
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
    memset(tableCell, 0, alloc * sizeof(unsigned int));
    std::fill(isGhost, isGhost + alloc, false);
    memset(r, 0, alloc * sizeof(double));
    memset(x0, 0, alloc * sizeof(tVect));
    memset(x1, 0, alloc * sizeof(tVect));
    memset(radiusVec, 0, alloc * sizeof(tVect));
    // @todo springs
}
#endif  // PARTICLE2_H

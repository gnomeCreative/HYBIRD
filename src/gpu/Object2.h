#ifndef OBJECT2_H
#define OBJECT2_H
#include "myvector.h"

struct Object2 {
    unsigned int alloc;
    unsigned int count = 0;
    
    // object index
    unsigned int *index = nullptr;
    // original object index in case of copies
    unsigned int *originalIndex = nullptr;
    // element index
    unsigned int *ElID = nullptr;
    // particle radius
    double *r = nullptr;
    // position of the object
    tVect *x0 = nullptr;
    // velocity of the object
    tVect *x1 = nullptr;
    // force on the object
    tVect *FParticle = nullptr, *FHydro = nullptr;
    // force on the object, maximum over simulation time, and time of occurrence
    tVect *maxFParticle = nullptr;
    double *timeMaxFParticle = nullptr;
    // force on the object, saved when necessary
    tVect *savedFParticle = nullptr;
    // is it translating?
    bool *translating = nullptr;
    // translation vector
    tVect *trans = nullptr;

    /**
     * @todo zero init all allocated buffers
     */
    template<int impl>
    void allocObjects(unsigned int num);
    template<int impl>
    void initForces();
};

template<int impl>
void Object2::allocObjects(unsigned int num) {
    alloc = num;

    index = (unsigned int*)malloc(alloc * sizeof(unsigned int));
    originalIndex = (unsigned int*)malloc(alloc * sizeof(unsigned int));
    ElID = (unsigned int*)malloc(alloc * sizeof(unsigned int));
    r = (double*)malloc(alloc * sizeof(double));
    x0 = (tVect*)malloc(alloc * sizeof(tVect));
    x1 = (tVect*)malloc(alloc * sizeof(tVect));
    FParticle = (tVect*)malloc(alloc * sizeof(tVect));
    FHydro = (tVect*)malloc(alloc * sizeof(tVect));
    maxFParticle = (tVect*)malloc(alloc * sizeof(tVect));
    timeMaxFParticle = (double*)malloc(alloc * sizeof(double));
    savedFParticle = (tVect*)malloc(alloc * sizeof(tVect));
    translating = (bool*)malloc(alloc * sizeof(bool));
    trans = (tVect*)malloc(alloc * sizeof(tVect));
}
template <>
inline void Object2::initForces<CPU>() {
    // initializing object forces
    memset(FHydro, 0, sizeof(tVect) * count);
}
#ifdef USE_CUDA
template <>
inline void Object2::initForces<CUDA>() {
    // initializing object forces
    CUDA_CALL(cudaMemset(FHydro, 0, sizeof(tVect) * count));
}
#endif

#endif  // OBJECT2_H

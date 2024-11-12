#ifndef OBJECT2_H
#define OBJECT2_H
#include "myvector.h"

struct Object2 {
    unsigned int count;

    tVect *FHydro = nullptr;

    template<int impl>
    void initForces();
};
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

#ifndef WALL2_H
#define WALL2_H
#include "myvector.h"

struct Wall2 {
    unsigned int count;

    tVect *FHydro = nullptr;

    template<int impl>
    void initForces();
};
template <>
inline void Wall2::initForces<CPU>() {
    // initializing wall forces
    memset(FHydro, 0, sizeof(tVect) * count);
}
#ifdef USE_CUDA
template <>
inline void Wall2::initForces<CUDA>() {
    // initializing wall forces
    CUDA_CALL(cudaMemset(FHydro, 0, sizeof(tVect) * count));
}
#endif

#endif  // WALL2_H

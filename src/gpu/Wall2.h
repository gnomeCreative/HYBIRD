#ifndef WALL2_H
#define WALL2_H
#include "myvector.h"

struct Wall2 {
    unsigned int count;

    // normal vector
    tVect *n = nullptr;
    // basic point
    tVect *p = nullptr;
    // hydraulic force on wall
    tVect *FHydro = nullptr;
    // center of rotation
    tVect *rotCenter = nullptr;
    // rotational speed
    tVect *omega = nullptr;
    // translational speed
    tVect *vel = nullptr;
    // is it a moving wall? (true = moving; false = fixed);
    bool *moving = nullptr;

    template<int impl>
    void initForces();

    // computes the speed of a point of the wall
    __host__ __device__ tVect getSpeed(const unsigned int i, const tVect pt) const {
        if (moving[i]) {
            // distance rotation center
            const tVect distFromCenter = pt - rotCenter[i];
            // distance from axes
            const tVect distFromAxes = distFromCenter - (distFromCenter.dot(n[i])) * n[i];
            // tangential velocity
            return vel[i] + omega[i].cross(distFromAxes);
        } else {
            return tVect(0.0, 0.0, 0.0);
        }
    }
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

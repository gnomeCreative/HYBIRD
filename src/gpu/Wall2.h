#ifndef WALL2_H
#define WALL2_H
#include "myvector.h"

struct DEMParams;

struct Wall2 {
    unsigned int alloc = 0;
    unsigned int count = 0;

    // index
    unsigned int *index;
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
    // is it a slipping wall?
    bool *slip = nullptr;
    // is it translating?
    bool* translating = nullptr;
    // translation vector
    tVect *trans = nullptr;
    // limits for walls
    bool* limited = nullptr;
    double *xMin = nullptr, *xMax = nullptr;
    double *yMin = nullptr, *yMax = nullptr;
    double *zMin = nullptr, *zMax = nullptr;

    /**
     * Does not update count
     * Initialises based on original constructor
     */
    template<int impl>
    void memoryAlloc(unsigned int num);
    void initialize(const DEMParams& dem_p, const std::array<types, 6>& externalBoundary, const std::array<tVect, 6>& boundaryLocation);
    // show wall characteristics
    void wallShow(unsigned int i) const;

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
template<>
inline void Wall2::memoryAlloc<CPU>(unsigned int num) {
    alloc = num;

    index = (unsigned int*)malloc(alloc * sizeof(unsigned int));
    n = (tVect*)malloc(alloc * sizeof(tVect));
    p = (tVect*)malloc(alloc * sizeof(tVect));
    FHydro = (tVect*)malloc(alloc * sizeof(tVect));
    rotCenter = (tVect*)malloc(alloc * sizeof(tVect));
    omega = (tVect*)malloc(alloc * sizeof(tVect));
    vel = (tVect*)malloc(alloc * sizeof(tVect));
    moving = (bool*)malloc(alloc * sizeof(bool));
    slip = (bool*)malloc(alloc * sizeof(bool));
    translating = (bool*)malloc(alloc * sizeof(bool));
    trans = (tVect*)malloc(alloc * sizeof(tVect));
    limited = (bool*)malloc(alloc * sizeof(bool));
    xMin = (double*)malloc(alloc * sizeof(double));
    xMax = (double*)malloc(alloc * sizeof(double));
    yMin = (double*)malloc(alloc * sizeof(double));
    yMax = (double*)malloc(alloc * sizeof(double));
    zMin = (double*)malloc(alloc * sizeof(double));
    zMax = (double*)malloc(alloc * sizeof(double));
    // Init
    memset(index, 0, alloc * sizeof(unsigned int));
    std::fill(n, n + alloc, tVect(1.0, 0.0, 0.0));
    memset(p, 0, alloc * sizeof(tVect));
    memset(FHydro, 0, alloc * sizeof(tVect));
    //memset(FParticle, 0, alloc * sizeof(tVect));
    //memset(maxFHydro, 0, alloc * sizeof(tVect));
    //memset(maxFParticle, 0, alloc * sizeof(tVect));
    memset(rotCenter, 0, alloc * sizeof(tVect));
    memset(omega, 0, alloc * sizeof(tVect));
    memset(vel, 0, alloc * sizeof(tVect));
    std::fill(moving, moving + alloc, false);
    std::fill(slip, slip + alloc, false);
    std::fill(translating, translating + alloc, false);
    memset(trans, 0, alloc * sizeof(tVect));
    std::fill(limited, limited + alloc, false);
    std::fill(xMin, xMin + alloc, numeric_limits<double>::max());
    std::fill(xMax, xMax + alloc, numeric_limits<double>::max());
    std::fill(yMin, yMin + alloc, numeric_limits<double>::max());
    std::fill(yMax, yMax + alloc, numeric_limits<double>::max());
    std::fill(zMin, zMin + alloc, numeric_limits<double>::max());
    std::fill(zMax, zMax + alloc, numeric_limits<double>::max());    
}
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

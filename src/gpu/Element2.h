#ifndef ELEMENT2_H
#define ELEMENT2_H
#include "cuda_helper.h"
#include "myvector.h"

struct Particle2;
struct DEMParams;
/**
 * \brief Representation of objects, such as walls, within the DEM.
 * @todo/note original member vars wsolver, index, active are fixed at runtime
 */
struct Element2 {
    // The total number of active elements
    unsigned int activeCount = 0;
    // The allocated size of activeI
    // We don't shrink the buffer if the number of active elements decreases
    unsigned int activeAlloc = 0;
    // Index of nodes marked as active
    unsigned int* activeI = nullptr;

    // The total number of elements
    // @note This is currently fixed at runtime, but may become variable in future
    unsigned int count = 0;
    // Size of allocated buffer, as it never shrinks
    unsigned int alloc = 0;

    bool *wSolver = nullptr;
    // element index
    unsigned int *index = nullptr;
    // is it active? (false=destroyed)
    bool *active = nullptr;
    // constitutive particles indexes
    intList *components = nullptr;
    // number of constitutive particles
    unsigned int *size = nullptr;
    // radius of constitutive particles (supposed, for the moment, of being constant)
    double *radius = nullptr;
    // mass stored in the element
    double *m = nullptr;
    // inertia tensor on principal axis (diagonal)
    tVect *I = nullptr;
    // Position and derivatives
    // position of the center of mass of the element
    tVect *x0 = nullptr;
    // velocity of the center of the element
    tVect *x1 = nullptr;
    // acceleration of the center of the element
    tVect *x2 = nullptr;
    // other velocity derivative for Gear scheme
    tVect *x3 = nullptr, *x4 = nullptr, *x5 = nullptr;
    // predicted quantities for Gear scheme
    tVect *xp0 = nullptr, *xp1 = nullptr, *xp2 = nullptr, *xp3 = nullptr, *xp4 = nullptr, *xp5 = nullptr;
    // total displacement
    tVect *x0history = nullptr;
    // Orientation and derivatives
    // orientation of the center of mass of the element
    tQuat *q0 = nullptr;
    // quaternion rates (related to angular velocity, acceleration...)
    tQuat *q1 = nullptr, *q2 = nullptr, *q3 = nullptr, *q4 = nullptr, *q5 = nullptr;
    // predicted quantities for Gear scheme
    tQuat *qp0 = nullptr, *qp1 = nullptr, *qp2 = nullptr, *qp3 = nullptr, *qp4 = nullptr, *qp5 = nullptr;
    // angular velocity in global and local reference frame
    tVect *wGlobal = nullptr, *wLocal = nullptr;
    // predicted angular velocities
    tVect *wpGlobal = nullptr, *wpLocal = nullptr;
    // angular velocity rates (global)
    tVect *w0 = nullptr, *w1 = nullptr, *w2 = nullptr, *w3 = nullptr, *w4 = nullptr, *w5 = nullptr;
    // predicted angular velocity rates (global)
    tVect *wp0 = nullptr, *wp1 = nullptr, *wp2 = nullptr, *wp3 = nullptr, *wp4 = nullptr, *wp5 = nullptr;
    // forces and moments
    tVect *FHydro = nullptr, *FParticle = nullptr, *FWall = nullptr, *FGrav = nullptr, *FSpringP = nullptr, *FSpringW = nullptr;
    tVect *MHydro = nullptr, *MParticle = nullptr, *MWall = nullptr, *MRolling = nullptr;
    // force intensities
    tVect *solidIntensity = nullptr;
    // connectivity (coordination number)
    unsigned int *coordination = nullptr;
    // apparent accelerations
    tVect *ACoriolis = nullptr, *ACentrifugal = nullptr;
    // fluid mass entrained by the particle
    double *fluidVolume = nullptr;
    // container for maximum overlap, used for output
    double *maxOverlap = nullptr;
    // container for maximum overlap between time steps, used for output
    double *maxDtOverlap = nullptr;
    // whether the element is slipping on a wall
    int *slippingCase = nullptr;
    /**
     * elmt::components is a std::vector, whereby each elmt may have a different length
     *
     * As such to represent it in CUDA we store all of these vectors compactly in order
     * inside Element2::componentsData
     *
     * In order to access a particular element's components that element's componentsIndex is read
     * This value is the first index within componentsData containing component data for the element.
     * The following elements componentsIndex points to the item after this elements last component
     * e.g.
     * auto e_i = 1  // element index
     * for (auto i = componentsIndex[e_i]; i < componentsIndex[e_i+1]; ++i)
     *     auto component = componentsData[i];
     */
    unsigned int *componentsIndex = nullptr; // length == count+1
    unsigned int *componentsData = nullptr;  // length == componentsIndex[count]

    /**
     * Does not update count, or allocate componentsIndex/componentsData
     * @todo zero init all allocated buffers
     */
    template<int impl>
    void memoryAlloc(unsigned int num);
    /**
     * Alloc the componentsData buffer
     * ComponentsIndex must be filled before this can occur
     */
    void allocComponentsData();
    /**
     * Initialise elements (on the host)
     * @note Large parts of this init could be bulk handled with memset during alloc
     */
    void initialize(unsigned int index, double partDensity);

    /**
     * Initialise particles for element at index, beginning with the particle index globalIndex
     * @param index Index of the element to generate particles for
     * @param globalIndex Index within particles to init the first particle
     * @param particles The particles structure to init particles into
     * @param dem_p The DEM parameters
     */
    void generateParticles(unsigned int index, unsigned int& globalIndex, Particle2& particles);

    template<int impl>
    void initElements();

    __host__ __device__ void predict(unsigned int index);
    __host__ __device__ void correct(unsigned int index, const std::array<double, 6>& coeff1ord, const std::array<double, 6>& coeff2ord);

    __host__ __device__ void atomicMaxOverlap(unsigned int index, double overlap);

};


__host__ __device__ __forceinline__ void Element2::predict(const unsigned int i) {
    tVect t = xp0[i];
    xp0[i] = x0[i] + x1[i] * DEM_P.c2[0] + x2[i] * DEM_P.c2[1] + x3[i] * DEM_P.c2[2] + x4[i] * DEM_P.c2[3] + x5[i] * DEM_P.c2[4];

    if (xp0[i].x < 0 || xp0[i].y < 0 || xp0[i].z < 0 || isnan(xp0[i].x) || isnan(xp0[i].y) || isnan(xp0[i].z)) {
        printf("Bad element predict xp0 (%f, %f, %f)->(%f, %f, %f) from x0(%f, %f, %f), x1(%f, %f, %f), x2(%f, %f, %f), x3(%f, %f, %f), x4 (%f, %f, %f), x5(%f, %f, %f)\n",
            t.x, t.y, t.z, xp0[i].x, xp0[i].y, xp0[i].z, x0[i].x, x0[i].y, x0[i].z, x1[i].x, x1[i].y, x1[i].z, x2[i].x, x2[i].y, x2[i].z, x3[i].x, x3[i].y, x3[i].z, x4[i].x, x4[i].y, x4[i].z, x5[i].x, x5[i].y, x5[i].z);
    }
    xp1[i] = x1[i] + x2[i] * DEM_P.c2[0] + x3[i] * DEM_P.c2[1] + x4[i] * DEM_P.c2[2] + x5[i] * DEM_P.c2[3];
    xp2[i] = x2[i] + x3[i] * DEM_P.c2[0] + x4[i] * DEM_P.c2[1] + x5[i] * DEM_P.c2[2];
    xp3[i] = x3[i] + x4[i] * DEM_P.c2[0] + x5[i] * DEM_P.c2[1];
    xp4[i] = x4[i] + x5[i] * DEM_P.c2[0];
    xp5[i] = x5[i];

    qp0[i] = q0[i] + q1[i] * DEM_P.c2[0] + q2[i] * DEM_P.c2[1] + q3[i] * DEM_P.c2[2] + q4[i] * DEM_P.c2[3] + q5[i] * DEM_P.c2[4];
    qp1[i] = q1[i] + q2[i] * DEM_P.c2[0] + q3[i] * DEM_P.c2[1] + q4[i] * DEM_P.c2[2] + q5[i] * DEM_P.c2[3];
    qp2[i] = q2[i] + q3[i] * DEM_P.c2[0] + q4[i] * DEM_P.c2[1] + q5[i] * DEM_P.c2[2];
    qp3[i] = q3[i] + q4[i] * DEM_P.c2[0] + q5[i] * DEM_P.c2[1];
    qp4[i] = q4[i] + q5[i] * DEM_P.c2[0];
    qp5[i] = q5[i];

    qp0[i].normalize();

    wp0[i] = w0[i] + w1[i] * DEM_P.c1[0] + w2[i] * DEM_P.c1[1] + w3[i] * DEM_P.c1[2] + w4[i] * DEM_P.c1[3] + w5[i] * DEM_P.c1[4];
    wp1[i] = w1[i] + w2[i] * DEM_P.c1[0] + w3[i] * DEM_P.c1[1] + w4[i] * DEM_P.c1[2] + w5[i] * DEM_P.c1[3];
    wp2[i] = w2[i] + w3[i] * DEM_P.c1[0] + w4[i] * DEM_P.c1[1] + w5[i] * DEM_P.c1[2];
    wp3[i] = w3[i] + w4[i] * DEM_P.c1[0] + w5[i] * DEM_P.c1[1];
    wp4[i] = w4[i] + w5[i] * DEM_P.c1[0];
    wp5[i] = w5[i];

    if (wSolver[i]) {
        wpGlobal[i] = wp0[i];
        wpLocal[i] = project(wpGlobal[i], qp0[i].adjoint());
    } else {
        //q1[i].forceStability(q0[i]);
        const tQuat qp0adj = qp0[i].adjoint();
        wpGlobal[i] = 2.0 * quat2vec(qp1[i].multiply(qp0adj));
        wpLocal[i] = 2.0 * quat2vec(qp0adj.multiply(qp1[i]));
    }
}

__host__ __device__ __forceinline__ void Element2::correct(const unsigned int i, const std::array<double, 6> &coeff1ord, const std::array<double, 6> &coeff2ord) {

    const tVect x2Corr = x2[i] - xp2[i];
    tVect t = x0[i];
    x0[i] = xp0[i] + x2Corr * coeff2ord[0];
    if (x0[i].x < 0 || x0[i].y < 0 || x0[i].z < 0 || isnan(x0[i].x) || isnan(x0[i].y) || isnan(x0[i].z)) {
        printf("Bad element correct x0 (%f, %f, %f)->(%f, %f, %f) from x2Corr: (%f, %f, %f). coeff2ord[0]: %f\n", t.x, t.y, t.z, x0[i].x, x0[i].y, x0[i].z, x2Corr.x, x2Corr.y, x2Corr.z, coeff2ord[0]);
    }
    x1[i] = xp1[i] + x2Corr * coeff2ord[1];
    // x2 calculated directly at the end of force routine
    x3[i] = xp3[i] + x2Corr * coeff2ord[3];
    x4[i] = xp4[i] + x2Corr * coeff2ord[4];
    x5[i] = xp5[i] + x2Corr * coeff2ord[5];

    xp0[i] = x0[i];
    xp1[i] = x1[i];
    xp2[i] = x2[i];
    xp3[i] = x3[i];
    xp4[i] = x4[i];
    xp5[i] = x5[i];

    const tVect w1Corr = w1[i] - wp1[i];

    w0[i] = wp0[i] + w1Corr * coeff1ord[0];
    // w1 calculated directly at the end of force routine
    w2[i] = wp2[i] + w1Corr * coeff1ord[2];
    w3[i] = wp3[i] + w1Corr * coeff1ord[3];
    w4[i] = wp4[i] + w1Corr * coeff1ord[4];
    w5[i] = wp5[i] + w1Corr * coeff1ord[5];

    wp0[i] = w0[i];
    wp1[i] = w1[i];
    wp2[i] = w2[i];
    wp3[i] = w3[i];
    wp4[i] = w4[i];
    wp5[i] = w5[i];

    const tQuat q2Corr = q2[i] - qp2[i];

    q0[i] = qp0[i] + q2Corr * coeff2ord[0];
    q1[i] = qp1[i] + q2Corr * coeff2ord[1];
    // q2 calculated directly at the end of force routine
    q3[i] = qp3[i] + q2Corr * coeff2ord[3];
    q4[i] = qp4[i] + q2Corr * coeff2ord[4];
    q5[i] = qp5[i] + q2Corr * coeff2ord[5];
    //normalization of q0
    q0[i].normalize();

    qp0[i] = q0[i];
    qp1[i] = q1[i];
    qp2[i] = q2[i];
    qp3[i] = q3[i];
    qp4[i] = q4[i];
    qp5[i] = q5[i];


    if (wSolver[i]) {
        wGlobal[i] = w0[i];
        wLocal[i] = project(wGlobal[i], q0[i].adjoint());
    } else {
        //q1[i].forceStability(q0[i]);
        const tQuat q0adj = q0[i].adjoint();
        wGlobal[i] = 2.0 * quat2vec(q1[i].multiply(q0adj));
        wLocal[i] = 2.0 * quat2vec(q0adj.multiply(q1[i]));
    }
}
__host__ __device__ inline void Element2::atomicMaxOverlap(unsigned int index, const double overlap) {
#ifdef __CUDA_ARCH__
    // CUDA does not support native atomicMin/atomicMax for float/double
    // So simulate it with a more expensive CAS loop
    unsigned long long int* addr_as_ull = (unsigned long long int*)&maxOverlap[index];
    unsigned long long int old = *addr_as_ull, assumed;

    do {
        assumed = old;
        if (__longlong_as_double(assumed) >= overlap)
            break;  // Exit if the current value is already greater or equal
        old = atomicCAS(addr_as_ull, assumed, __double_as_longlong(overlap));
    } while (assumed != old);  // Repeat until successful update
#endif
}

template<>
inline void Element2::memoryAlloc<CPU>(unsigned int num) {
    alloc = num;

    if(wSolver) {
        free(wSolver);
        free(index);
        free(size);
        free(m);
        free(I);
        free(x0);
        free(x1);
        free(x2);
        free(x3);
        free(x4);
        free(x5);
        free(xp5);
        free(xp0);
        free(xp1);
        free(xp2);
        free(xp3);
        free(xp4);
        free(xp5);
        free(x0history);
        free(q0);
        free(q1);
        free(q2);
        free(q3);
        free(q4);
        free(q5);
        free(qp5);
        free(qp0);
        free(qp1);
        free(qp2);
        free(qp3);
        free(qp4);
        free(qp5);
        free(wGlobal);
        free(wLocal);
        free(wpLocal);
        free(wpGlobal);
        free(w0);
        free(w1);
        free(w2);
        free(w3);
        free(w4);
        free(w5);
        free(wp5);
        free(wp0);
        free(wp1);
        free(wp2);
        free(wp3);
        free(wp4);
        free(wp5);
        free(FHydro);
        free(FParticle);
        free(FWall);
        free(FGrav);
        free(FSpringP);
        free(FSpringW);
        free(MHydro);
        free(MParticle);
        free(MWall);
        free(MRolling);
        free(solidIntensity);
        free(coordination);
        free(ACoriolis);
        free(ACentrifugal);
        free(fluidVolume);
        free(maxOverlap);
        free(maxDtOverlap);
        free(slippingCase);
        free(componentsIndex);
    }

    wSolver = (bool*)malloc(alloc * sizeof(bool));
    index = (unsigned int*)malloc(alloc * sizeof(unsigned int));
    active = (bool*)malloc(alloc * sizeof(bool));
    componentsIndex = (unsigned int*)malloc((alloc + 1) * sizeof(unsigned int));    
    // componentsData = (unsigned int*)malloc(? * sizeof(unsigned int)); // This is handled later once index is full
    size = (unsigned int*)malloc(alloc * sizeof(unsigned int));
    radius = (double*)malloc(alloc * sizeof(double));
    m = (double*)malloc(alloc * sizeof(double));
    I = (tVect*)malloc(alloc * sizeof(tVect));
    x0 = (tVect*)malloc(alloc * sizeof(tVect));
    x1 = (tVect*)malloc(alloc * sizeof(tVect));
    x2 = (tVect*)malloc(alloc * sizeof(tVect));
    x3 = (tVect*)malloc(alloc * sizeof(tVect));
    x4 = (tVect*)malloc(alloc * sizeof(tVect));
    x5 = (tVect*)malloc(alloc * sizeof(tVect));
    xp0 = (tVect*)malloc(alloc * sizeof(tVect));
    xp1 = (tVect*)malloc(alloc * sizeof(tVect));
    xp2 = (tVect*)malloc(alloc * sizeof(tVect));
    xp3 = (tVect*)malloc(alloc * sizeof(tVect));
    xp4 = (tVect*)malloc(alloc * sizeof(tVect));
    xp5 = (tVect*)malloc(alloc * sizeof(tVect));
    x0history = (tVect*)malloc(alloc * sizeof(tVect));
    q0 = (tQuat*)malloc(alloc * sizeof(tQuat));
    q1 = (tQuat*)malloc(alloc * sizeof(tQuat));
    q2 = (tQuat*)malloc(alloc * sizeof(tQuat));
    q3 = (tQuat*)malloc(alloc * sizeof(tQuat));
    q4 = (tQuat*)malloc(alloc * sizeof(tQuat));
    q5 = (tQuat*)malloc(alloc * sizeof(tQuat));
    qp0 = (tQuat*)malloc(alloc * sizeof(tQuat));
    qp1 = (tQuat*)malloc(alloc * sizeof(tQuat));
    qp2 = (tQuat*)malloc(alloc * sizeof(tQuat));
    qp3 = (tQuat*)malloc(alloc * sizeof(tQuat));
    qp4 = (tQuat*)malloc(alloc * sizeof(tQuat));
    qp5 = (tQuat*)malloc(alloc * sizeof(tQuat));
    wGlobal = (tVect*)malloc(alloc * sizeof(tVect));
    wLocal = (tVect*)malloc(alloc * sizeof(tVect));
    wpLocal = (tVect*)malloc(alloc * sizeof(tVect));
    wpGlobal = (tVect*)malloc(alloc * sizeof(tVect));
    w0 = (tVect*)malloc(alloc * sizeof(tVect));
    w1 = (tVect*)malloc(alloc * sizeof(tVect));
    w2 = (tVect*)malloc(alloc * sizeof(tVect));
    w3 = (tVect*)malloc(alloc * sizeof(tVect));
    w4 = (tVect*)malloc(alloc * sizeof(tVect));
    w5 = (tVect*)malloc(alloc * sizeof(tVect));
    wp0 = (tVect*)malloc(alloc * sizeof(tVect));
    wp1 = (tVect*)malloc(alloc * sizeof(tVect));
    wp2 = (tVect*)malloc(alloc * sizeof(tVect));
    wp3 = (tVect*)malloc(alloc * sizeof(tVect));
    wp4 = (tVect*)malloc(alloc * sizeof(tVect));
    wp5 = (tVect*)malloc(alloc * sizeof(tVect));
    FHydro = (tVect*)malloc(alloc * sizeof(tVect));
    FParticle = (tVect*)malloc(alloc * sizeof(tVect));
    FWall = (tVect*)malloc(alloc * sizeof(tVect));
    FGrav = (tVect*)malloc(alloc * sizeof(tVect));
    FSpringP = (tVect*)malloc(alloc * sizeof(tVect));
    FSpringW = (tVect*)malloc(alloc * sizeof(tVect));
    MHydro = (tVect*)malloc(alloc * sizeof(tVect));
    MParticle = (tVect*)malloc(alloc * sizeof(tVect));
    MWall = (tVect*)malloc(alloc * sizeof(tVect));
    MRolling = (tVect*)malloc(alloc * sizeof(tVect));
    solidIntensity = (tVect*)malloc(alloc * sizeof(tVect));
    coordination = (unsigned int*)malloc(alloc * sizeof(unsigned int));
    ACoriolis = (tVect*)malloc(alloc * sizeof(tVect));
    ACentrifugal = (tVect*)malloc(alloc * sizeof(tVect));
    fluidVolume = (double*)malloc(alloc * sizeof(double));
    maxOverlap = (double*)malloc(alloc * sizeof(double));
    maxDtOverlap = (double*)malloc(alloc * sizeof(double));
    slippingCase = (int*)malloc(alloc * sizeof(int));
}
template <>
inline void Element2::initElements<CPU>() {
    //initializing this time step hydrodynamic force
    memset(FHydro, 0, sizeof(tVect) * count);
    memset(MHydro, 0, sizeof(tVect) * count);
    // initializing the fluid mass for buoyancy
    memset(fluidVolume, 0, sizeof(double) * count);
}
#ifdef USE_CUDA
template <>
inline void Element2::initElements<CUDA>() {
    //initializing this time step hydrodynamic force
    CUDA_CALL(cudaMemset(FHydro, 0, sizeof(tVect) * count));
    CUDA_CALL(cudaMemset(MHydro, 0, sizeof(tVect) * count));
    // initializing the fluid mass for buoyancy
    CUDA_CALL(cudaMemset(fluidVolume, 0, sizeof(double) * count));
}
#endif
#endif  // ELEMENT2_H

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
    void initialize(unsigned int index, double partDensity, const DEMParams& dem_p);

    /**
     * Initialise particles for element at index, beginning with the particle index globalIndex
     * @param index Index of the element to generate particles for
     * @param globalIndex Index within particles to init the first particle
     * @param particles The particles structure to init particles into
     * @param dem_p The DEM parameters
     */
    void generateParticles(unsigned int index, unsigned int& globalIndex, Particle2& particles, const DEMParams& dem_p);

    template<int impl>
    void initElements();
};

template<>
inline void Element2::memoryAlloc<CPU>(unsigned int num) {
    alloc = num;

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

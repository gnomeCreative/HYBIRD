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
};


__host__ __device__ __forceinline__ void Element2::predict(const unsigned int i) {
    xp0[i] = x0[i] + x1[i] * DEM_P.c2[0] + x2[i] * DEM_P.c2[1] + x3[i] * DEM_P.c2[2] + x4[i] * DEM_P.c2[3] + x5[i] * DEM_P.c2[4];
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

#endif  // ELEMENT2_H

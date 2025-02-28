#ifndef DEMPARAMS_H
#define DEMPARAMS_H

#include "elmt.h"

struct Element2;
struct Object2;
class GetPot;

/**
 * Discrete Element model parameters
 * These should not change after initialisation
 * @note This struct has been moved into a separate file to avoid a circular dependency
 *       whilst maximising how much code can be inlined
 * @note Members marked t have not been checked whether redundant
 */
struct DEMParams {
    // neighbor list (here COMMENTS are needed)
    // neighbor parameters
    // Used in the transformation of wrapped positions to ghost positions
    tVect half_demSize = {};
    /// ?
    double maxDisp = {};
    /// ?
    double nebrRange = {};
    //double nebrShell;
    // With of neighbour grid cells in each dimension
    tVect cellWidth = {};
    // Number of neighbour grid cells in each dimension
    unsigned int nCells[3] = {};
    // numerical viscosity to avoid energy problems
    double numVisc = {};
    /**
     * Initialise struct from run args/input file
     */
    void discreteElementGet(GetPot& configFile, GetPot& commandLine, Element2& elmts, Object2& objects);
    // domain size: DEM grid is orthogonal
    tVect demSize = { 1.0, 1.0, 1.0 };
    // switchers for rotating local system
    bool solveCoriolis = {};
    bool solveCentrifugal = {};
    // switcher for static friction
    bool staticFrictionSolve = {};
    // number of dem time steps before activating LB
    double demInitialRepeat = 0;
    // time step
    unsigned int demTimeStep = 0;
    // time step duration
    double deltat = 1.0;
    // C1 and C2 used by DEM::predictor()
    std::array<double, 5> c1 = {};
    std::array<double, 5> c2 = {};
    void init_prototypeC1C2();
    // total time duration
    double demTime = 0.0;
    // force field (if not constant need to be calculated inside cycle) // should not be here!!!!!
    tVect demF = { 0,0,0 };
    // rotation of the local reference system
    tVect demRot = { 0,0,0 };
    // center of rotation of the local reference system
    tVect demRotCenter = { 0,0,0 };
    // material for elements // should not be here!!!
    material sphereMat = {};
    // relative difference between time step in DEM and time step in LBM
    unsigned int multiStep = 1;
    // ratio between time step and estimated duration of contacts
    double criticalRatio = 0.1;

    // total force on obstacles
    tVect objMaxTotalForce = {0, 0, 0};

    // PROBLEM-SPECIFIC
    // stuff for drum
    double drumSpeed = 0.0;
    // stuff for hong kong
    double hongkongSlitSize = {};
    bool hongkongSmoothWall = {};
    // stuff for hourglass
    double hourglassOutletSize = {};
    double hourglassOutletHeight = {};
    // stuff for continuum heap
    double heapBaseLevel = {};
    // stuff for triaxial tests
    double triIsopressure = {};
    double triBeginX = {}, triBeginY = {}, triBeginZ = {};
    double triDefSpeed = {};
    double pressureX = {}, pressureY = {}, pressureZ = {};
    // stuff for Usman
    bool depositArea = false;
};

// sqrt() is not constexpr until C++26
// Result of sqrt(2.0), sqrt(3.0), sqrt(6.0) (also used by lattice.h)
constexpr double sqrt_2 = 1.414213562373095;
constexpr double sqrt_3 = 1.732050807568877;
constexpr double sqrt_6 = 2.449489742783178;

__host__device__ constexpr tVect prototypes[5][4] = { {
        // prototypes[0]
    }, {
        // prototypes[1]
        tVect(0.0,0.0,0.0)
    }, {
        // prototypes[2]
        tVect(0.5, 0.0, 0.0),
        tVect(-0.5, 0.0, 0.0)
    }, {
        // prototypes[3]
        tVect(0.0, 1.0, 0.0),
        tVect(-sqrt_3 / 2, -1 / 2, 0.0), // @bug, see issue #19
        tVect(sqrt_3 / 2, -1 / 2, 0.0)   // @bug, see issue #19
    }, {
        // prototypes[4]
        tVect(0.0, 0.0, 1.0),
        tVect(0.0, 2.0 * sqrt_2 / 3.0, -1.0 / 3.0),
        tVect(2.0 * sqrt_6 / 6.0, -2.0 * sqrt_2 / 6.0, -1.0 / 3.0),
        tVect(-2.0 * sqrt_6 / 6.0, -2.0 * sqrt_2 / 6.0, -1.0 / 3.0)
    } };
/**
 *
 * Global scope Discrete Element Model parameters structure
 * It is defined (e.g. sans 'extern') inside DEMParams.cu
 * CUDA constant's can only be defined in global (or file) scope
 * Therefore, the host replacement is also defined this way
 *
 * In order for them to share a name within shared host/device code
 * we selectively define the macro DEM_P, which maps to h_DEM_P
 * or d_DEM_P subject to whether __CUDA_ARCH__ has been defined.
 * __CUDA_ARCH__ is automatically defined by the CUDA compiler whilst
 * it compiles device code, to specify the target architecture,
 * it is not defined whilst the CUDA compiler builds host code.
 * However d_params is behind USE_CUDA, as it must exist in both
 * host and device code if compiled by the CUDA compiler.
 *
 * If you wish to specifically access the host or device versions
 * h_DEM_P and d_DEM_P can be used respectively.
 *
 * It may be possible to replace DEM_P with a reference,
 * however this would require careful CUDA testing.
 */
 // Declare the global storage variables
#ifdef USE_CUDA
extern __constant__ DEMParams d_DEM_P;
#endif
extern DEMParams h_DEM_P;
// Define the DEM_P macro which maps to the corresponding storage variable
#ifdef __CUDA_ARCH__
#define DEM_P d_DEM_P
#else
#define DEM_P h_DEM_P
#endif

// Implemented here, rather than myvector.inl, to avoid circular dependency
// (whilst keeping it inlined)
inline __host__ __device__ int tVect::linearizePosition() const {
    const int xc = static_cast<int>(floor(x / DEM_P.cellWidth[0]) + 1);
    const int yc = static_cast<int>(floor(y / DEM_P.cellWidth[1]) + 1);
    const int zc = static_cast<int>(floor(z / DEM_P.cellWidth[2]) + 1);
    const int index = static_cast<int>(xc + DEM_P.nCells[0] * (yc + DEM_P.nCells[1] * zc));
    return index;
}

#endif  // DEMPARAMS_H
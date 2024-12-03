#ifndef LBPARAMS_H
#define LBPARAMS_H

#include <array>

#include "lattice.h"
#include "MeasureUnits.h"
#include "node.h"
#include "utils.h"

class GetPot;

/**
 * Lattice Boltzmann model parameters
 * These should not change during after initialisation
 * @note This struct has been moved into a seperate file to avoid a circular dependency
 *       whilst maximising how much code can be inlined
 * @note Members marked t have not been checked whether redundant
 */
struct LBParams {
    /**
     * Initialise struct from run args/input file
     */
    void latticeBoltzmannGet(GetPot &configFile, GetPot &commandLine);
    /**
     * Initialise dynamic lattice params that scale with the model configuration
     * @see lattice.h for the static lattice params
     */
    void latticeDefinition();
    /**
     * Output lattice characteristics to stdout
     */
    void LBShow() const;
    // model for the fluid
    FluidMaterial fluidMaterial = {};
    // slip coefficient
    double slipCoefficient = 0.0;
    // hydrodynamic radius (see Kumnar et al., Mechanics of granular column collapse in fluid at varying slope angles)
    double hydrodynamicRadius = 1;
    // initial velocity for the fluid
    tVect initVelocity = {};  // t
    // force field (if not constant need to be calculated inside cycle)
    tVect lbF = {};  // @todo this needs to be implemented properly
    // switchers for rotating local system
    bool solveCoriolis = false;
    bool solveCentrifugal = false;
    // rotation speed of the local coordinate system
    tVect rotationSpeed = { 0,0,0 };
    // cebnter of rotation of the local coordinate system
    tVect rotationCenter = {};  // t
    // total number of nodes
    unsigned int totPossibleNodes = 0; // t
    // total mass (initial)
    double totalMass = 0;
    // standard neighbors shifting
    std::array<int, lbmDirec> ne = {};
    std::array<unsigned int, 3> shift = { 1,1,1 };
    std::array<int, 3> domain = { 1,1,1 };
    // boundary conditions
    std::array<types, 6> boundary = {};  // t
    /////////////////
    // config file //
    /////////////////
    // switcher for restart
    bool lbRestart = false;  // t
    // restart file
    //std::string lbRestartFile = "";  // t //@todo
    // switcher for imposed volume, and imposed volume
    bool imposeFluidVolume = false;  // t
    double imposedFluidVolume = 0.0;  // t
    // switcher for increasing volume, and extra volume and time to add it
    bool increaseVolume = false;  // t
    double deltaVolume = 0.0;  // t
    double deltaTime = 0.0;  // t
    // switcher for topography
    bool lbTopography = false;  // t
    // switcher for topographic initial level
    bool lbTopographySurface = false;  // t
    // shifts for topography file
    double translateTopographyX = 0.0;  // t
    double translateTopographyY = 0.0;  // t
    double translateTopographyZ = 0.0;  // t
    // topography file
    //string lbTopographyFile = "";  // t //@todo
    // topography container
    //topography lbTop = {};  // t
    // switchers for force field, non-Newtonian and everything
    bool freeSurface = false;
    bool forceField = false;
    // two-relaxation time (TRT) solver
    bool TRTsolver = false;
    double magicNumber = 0.25;
    // number of LBM step before activating the free surface
    unsigned int lbmInitialRepeat = 0;  // t
    // absolute time
    unsigned int time = 0;
    // lbm size in cell units
    std::array<unsigned int, 3> lbSize = { 1,1,1 };
    // lbm size in physical units (with boundaries)
    std::array<double, 3> lbPhysicalSize = { -1,-1,-1 };  // t
    // lbm size in physical units (without boundaries)
    std::array<double, 3> lbInnerPhysicalSize = { 1,1,1 };  // t
    // lbm box boundary locations
    std::array<tVect, 6> lbBoundaryLocation = {};  // t
    // standard borders for free surface
    std::array<int, 6> freeSurfaceBorders = { 0,1,0,1,0,1 };  // t
    //////////////////////
    // conversion units //
    //////////////////////
    // fluid and granular matter are solved in different measure units
    // for a reference, check Feng, Han, Owen, 2007
    MeasureUnits unit;
    ////////////////////////////
    // problem-specific stuff //
    ////////////////////////////
    // DRUM
    double fluidMass = 0.0;  // t
    // SHEARCELL
    double maxVisc = 0.0;  // t
    double maxPlasticVisc = 0.0;  // t
    double maxYieldStress = 0.0;  // t
    double maxShearRate = 0.0;  // t
    unsigned int viscSteps = 0;  // t
    unsigned int shearRateSteps = 0;  // t
    // NET, BARRIER: avalanches (net-like)
    double avalanchePosit = 0.0;  // t
    // HK_LARGE: Usman
    double largeFlumeFlowLevel = 0.0;  // t
    // HOURGLASS (mirrors in DEM)
    double hourglassOutletHeight = 0.0;  // t
    // HEAP: continuum heap (mirrors in DEM)
    double heapBaseLevel = 0.0;  // t
    
    // functions for linearized index management
    __host__ __device__ __forceinline__ unsigned int getIndex(const unsigned int& x, const unsigned int& y, const unsigned int& z) {
        return x + y * lbSize[0] + z * lbSize[0] * lbSize[1];
    }
    __host__ __device__ __forceinline__ tVect getPosition(const unsigned int& index) const {
        unsigned int x, y, z;

        // index is calculated in this fashion:
        // index = x + y*X + z*X*Y
        // where X and Y are sizes of the lattice in x and y direction

        // from this stems that
        // x + y*X = index MOD X*Y
        // x = x + y*X MOD X
        // y = x + y*X DIV X
        // z = index DIV X*Y

#ifdef __CUDA_ARCH__
        // div() does not exist in device code
        z = index / (lbSize[0] * lbSize[1]);
        const int t = index % (lbSize[0] * lbSize[1]);
        y = t / lbSize[0];
        x = t % lbSize[0];
#else
        // see online documentation for class div_t (stdlib.h)
        div_t firstDiv, secondDiv;

        firstDiv = div(int(index), int(lbSize[0] * lbSize[1]));
        secondDiv = div(firstDiv.rem, int(lbSize[0]));

        x = secondDiv.rem;
        y = secondDiv.quot;
        z = firstDiv.quot;
#endif
        return tVect(double(x) - 0.5, double(y) - 0.5, double(z) - 0.5);
    }
};

/**
 * Global scope Lattice Boltzmann parameters structure
 * It is defined (e.g. sans 'extern') inside LBParams.cu
 * CUDA constant's can only be defined in global (or file) scope
 * Therefore, the host replacement is also defined this way
 *
 * In order for them to share a name within shared host/device code
 * we selectively define the macro PARAMS, which maps to h_PARAMS
 * or d_PARAMS subject to whether __CUDA_ARCH__ has been defined.
 * __CUDA_ARCH__ is automatically defined by the CUDA compiler whilst
 * it compiles device code, to specify the target architecture,
 * it is not defined whilst the CUDA compiler builds host code.
 * However d_params is behind USE_CUDA, as it must exist in both
 * host and device code if compiled by the CUDA compiler.
 *
 * If you wish to specifically access the host or device versions
 * h_PARAMS and d_PARAMS can be used respectively.
 *
 * It may be possible to replace PARAMS with a reference,
 * however this would require careful CUDA testing.
 */
// Declare the global storage variables
#ifdef USE_CUDA
extern __constant__ LBParams d_PARAMS;
#endif
extern LBParams h_PARAMS;
// Define the PARAMS macro which maps to the corresponding storage variable
#ifdef __CUDA_ARCH__
#define PARAMS d_PARAMS
#else
#define PARAMS h_PARAMS
#endif

#endif  // LBPARAMS_H

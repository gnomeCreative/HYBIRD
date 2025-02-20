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
    // numerical viscosity to avoid energy problems
    double numVisc;

    __host__device__ constexpr tVect prototypes[5][] = { {        
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
        tVect(-sqrt(3) / 2, -1 / 2, 0.0),
        tVect(sqrt(3) / 2, -1 / 2, 0.0)
    }, {
        // prototypes[4]
        tVect(0.0, 0.0, 1.0),
        tVect(0.0, 2.0 * sqrt(2) / 3.0, -1.0 / 3.0),
        tVect(2.0 * sqrt(6) / 6.0, -2.0 * sqrt(2) / 6.0, -1.0 / 3.0),
        tVect(-2.0 * sqrt(6) / 6.0, -2.0 * sqrt(2) / 6.0, -1.0 / 3.0)
    } };
    /**
     * Initialise struct from run args/input file
     */
    void discreteElementGet(GetPot& configFile, GetPot& commandLine, Element2& elmts, Object2& objects);
    // domain size: DEM grid is orthogonal
    doubleList demSize = { 1.0, 1.0, 1.0 };
    // switchers for rotating local system
    bool solveCoriolis;
    bool solveCentrifugal;
    // switcher for static friction
    bool staticFrictionSolve;
    // number of dem time steps before activating LB
    double demInitialRepeat = 0;
    // time step
    unsigned int demTimeStep = 0;
    // total time duration
    double demTime = 0.0;
    // force field (if not constant need to be calculated inside cycle) // should not be here!!!!!
    tVect demF = { 0,0,0 };
    // rotation of the local reference system
    tVect demRot = { 0,0,0 };
    // center of rotation of the local reference system
    tVect demRotCenter = { 0,0,0 };
    // material for elements // should not be here!!!
    material sphereMat;
    // relative difference between time step in DEM and time step in LBM
    unsigned int multiStep = 1;
    // ratio between time step and estimated duration of contacts
    double criticalRatio = 0.1;

    // PROBLEM-SPECIFIC
    // stuff for drum
    double drumSpeed = 0.0;
    // stuff for hong kong
    double hongkongSlitSize;
    bool hongkongSmoothWall;
    // stuff for hourglass
    double hourglassOutletSize;
    double hourglassOutletHeight;
    // stuff for continuum heap
    double heapBaseLevel;
    // stuff for triaxial tests
    double triIsopressure;
    double triBeginX, triBeginY, triBeginZ;
    double triDefSpeed;
    double pressureX, pressureY, pressureZ;
    // stuff for Usman
    bool depositArea = false;
};

#endif  // DEMPARAMS_H
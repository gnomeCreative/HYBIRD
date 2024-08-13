/*
 * File:   DEM.h
 * Author: aleonard
 *
 * Created on November 5, 2013, 4:07 PM
 */

#ifndef DEM_H
#define DEM_H



#include <stdio.h>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <limits>
//
#include "myvector.h"
#include "elmt.h"
#include "utils.h"
#include "LB.h"

class Problem;
using namespace std;
extern ProblemName problemName;

// forward declaration of IO
class IO;

class DEM {
    //IO *io;

private:
    // neighbor list (here COMMENTS are needed)
    // neighbor parameters
    double maxDisp;
    double nebrRange;
    double nebrShell;
    double cellWidth[3];
    unsigned int nCells[3];
    // neighbor tables
    intList cellTable;
    // numerical viscosity to avoid energy problems
    double numVisc;
    // prototype shape collection
    std::vector <vecList> prototypes;
public:
    unsIntList neighborTable;
    unsIntList nearWallTable;
    unsIntList nearObjectTable;
    unsIntList nearCylinderTable;
    // domain size: DEM grid is orthogonal
    doubleList demSize;
    // switchers for rotating local system
    bool solveCoriolis;
    bool solveCentrifugal;
    // switcher for static friction
    bool staticFrictionSolve;
    // number of dem time steps before activating LB
    double demInitialRepeat;
    // time step
    unsigned int demTimeStep;
    // energy
    energy particleEnergy;
    // time step duration
    double deltat;
    // total time duration
    double demTime;
    // force field (if not constant need to be calculated inside cycle) // should not be here!!!!!
    tVect demF;
    // rotation of the local reference system
    tVect demRot;
    // center of rotation of the local reference system
    tVect demRotCenter;
    // material for elements // should not be here!!!
    material sphereMat;
    // relative difference between time step in DEM and time step in LBM
    unsigned int multiStep;
    // ratio between time step and estimated duration of contacts
    double criticalRatio;
    // solid walls
    wallList walls;
    // solid walls
    cylinderList cylinders;
    // periodic boundary conditions
    pbcList pbcs;
    // list with elements
    elmtList elmts;
    unsIntList activeElmts;
    // list with particles
    particleList particles;
    unsIntList activeParticles;
    // list with objects
    objectList objects;
    // particles contains both standard particles and ghost particles. The total number of standard particles needs to be saved
    unsigned int stdPartNumber;
    // number of active particles and elements
    unsigned int actvPartNumber;
    unsigned int actvElmtNumber;
    // mean, max and min radius for particles and objects
    double minPartRadius, maxPartRadius, minObjRadius, maxObjRadius, meanPartRadius, meanObjRadius;
    // list with ghost elements
    ghostList ghosts;
    // true if indices for the LB need to be regenerated
    bool newNeighborList;
    //couple elongation list
    elongVector elongTable;
    elongVector elongTableWall;
    elongVector elongTableObject;
    elongVector elongTableCylinder;

    // Total number of standard objects (i.e. not ghosts) which needs to be saved (Devis))
    unsigned int stdObjects;
    // total force on obstacles
    tVect objMaxTotalForce;
    // Total number of active springs
    unsigned int totSprings;
    // PROBLEM-SPECIFIC
    // stuff for drum
    double drumSpeed;
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
    double triBeginX,triBeginY,triBeginZ;
    double triDefSpeed;
    double pressureX,pressureY,pressureZ;
     // stuff for Usman
    bool depositArea;
    
    // Seguin problem
    bool movingCylinder;
    double coordinateCylinderX;
    double coordinateCylinderY;
    double coordinateCylinderZ;
    double velocityCylinderX;
    double velocityCylinderY;
    double velocityCylinderZ;
public:

    DEM() {
        demInitialRepeat = 0.0;
        demTimeStep = 0;
        drumSpeed = 0.0;
        newNeighborList = false;
        objMaxTotalForce=Zero;
        demSize.resize(3);
        demSize[0] = demSize[1] = demSize[2] = 1.0;
        demTime = 0.0;
        deltat = 1.0;
        demF.reset();
        demRot.reset();
        demRotCenter.reset();
        walls.clear();
        pbcs.clear();
        elmts.clear();
        particles.clear();
        objects.clear();
        stdPartNumber = 0;
        ghosts.clear();
        //        ghosts.clear();
        cellTable.clear();
        neighborTable.clear();
        nearWallTable.clear();
        nearObjectTable.clear();
        nearCylinderTable.clear();
        // elongation variables
        elongTable.clear();
        elongTableWall.clear();
        elongTableObject.clear();
        elongTableCylinder.clear();
        // neighbor list variables
        maxDisp = 0.0;
        nebrRange = 0.0;
        cellWidth[0] = cellWidth[1] = cellWidth[2] = 0.0;
        nCells[0] = nCells[1] = nCells[2] = 0;
        prototypes.clear();
        //        energy.reset();
        criticalRatio = 0.1;
        multiStep = 1;
        numVisc = 0.0;
        // stuff for Usman
        depositArea=false;
    }
    void discreteElementStep();
    void discreteElementGet(GetPot& config_file, GetPot& command_line);
    void discreteElementInit(const Problem& problem, const typeList& externalBoundary, const doubleList& externalSize, const vecList& externalBoundaryLocation,
            const tVect externalAccel, const tVect externalRotation, const tVect externalRotationCenter, const bool externalSolveCoriolis, const bool externalSolveCentrifugal, const double& externalTimeStep);
    void evolveBoundaries();
//    void evolveObj();
    void determineTimeStep(const double& externalTimeStep);
    // energy functions
    void updateEnergy(double& totalKineticEnergy);
private:
    // initialization functions
    void compositeProperties();
    void initializeWalls(const Problem& problem, const typeList& boundary, const vecList& boundaryLocation);
    void initializeCylinders(const Problem& problem);
    void initializePbcs(const typeList& boundary, const vecList& boundaryLocation);
    // integration functions
    void predictor();
    void corrector();
    void evaluateForces();
    void updateParticlesPredicted();
    void updateParticlesCorrected();
    double criticalTimeStep() const;
    // neighbor list functions
    void destroyElements();
    void evalNeighborTable(); // should be private: correct!
    void initNeighborParameters();
    void evalMaxDisp();
    void evalCellTable();
    void updateElongationTable(const unsIntList& oldNeigh, unsIntList& newNeigh, const elongVector& oldElong, elongVector& newElong);
    void evalNearWallTable();
    void evalNearObjectTable();
    void evalNearCylinderTable();
    // periodicity functions
    void createGhosts();
    void pbcShift();
    void periodicObjects();
    // force computation functions
    Elongation* findSpring(const unsigned int& t, const unsigned int& indexI, particle* partj);
    void computeApparentForces();
    void particleParticleContacts();
    void wallParticleContacts();
    void cylinderParticelContacts();
    void objectParticleContacts();
    inline void particleParticleCollision(const particle *partI, const particle *partJ, const tVect& vectorDistance, Elongation* elongation_new);
    inline void wallParticleCollision(wall *walli, const particle *partJ, const double& overlap, Elongation* elongation_new);
    inline void cylinderParticleCollision(cylinder *cylinderI, const particle *partJ, const double& overlap, Elongation* elongation_new);
    inline void objectParticleCollision(object *iObject, const particle *partJ, const tVect& vectorDistance, Elongation* elongation_new);
    double normalContact(const double& overlap, const double& vrelnnorm, const double& rEff, const double& massEff) const;
    double tangentialContact(const double& vreltNorm, const double& fn, const double& effRad, const double& effMass, const double& friction) const;
    tVect FRtangentialContact(const tVect& tangRelVelContact, const double& fn, const double& overlap, const double& effRad, const double& effMass, Elongation* elongation_new, const double& friction, const double& tangStiff, const double& viscTang);
    tVect rollingContact(const tVect& wI, const tVect& wJ, const double& effRad, const double& fn, const double& rolling);
    void saveObjectForces();
    void lubrication(int i, int j, double ri, double rj, tVect x0ij);
    /// this are the (new) ones that save the contacts (Devis))
};

#endif /* DEM_H */


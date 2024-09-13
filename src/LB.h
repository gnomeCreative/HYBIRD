/* 
 * File:   LB.h
 * Author: aleonard
 *
 * Created on November 5, 2013, 4:06 PM
 */

#ifndef LB_H
#define LB_H

//class IO;

#include <stdio.h>
//#include <iostream>
//#include <vector>
#include <stdlib.h>
#include <math.h>
//#include <algorithm>
#include <string.h>
#include <valarray>
#include <sys/time.h>
//
#include "macros.h"
#include "myvector.h"
#include "lattice.h"
#include "node.h"
#include "utils.h"
#include "elmt.h"
#include "getpot.h"

using namespace std;
extern ProblemName problemName;

class LB {
    // Characteristics ///////////////////////////////
    // lattice-Boltzmann equation is solved in a dimensionless space
    // with lattice spacing and time step both unitary
public: //private
    
    // link linearized index with nodes
    nodeMap nodes;
    // model for the fluid
    FluidMaterial fluidMaterial;
    // slip coefficient
    double slipCoefficient;
    // hydrodynamic radius (see Kumnar et al., Mechanics of granular column collapse in fluid at varying slope angles)
    double hydrodynamicRadius;
    // initial velocity for the fluid
    tVect initVelocity;
    // force field (if not constant need to be calculated inside cycle)
    tVect lbF;
    // switchers for rotating local system
    bool solveCoriolis;
    bool solveCentrifugal;
    // rotation speed of the local coordinate system
    tVect rotationSpeed;
    // cebnter of rotation of the local coordinate system
    tVect rotationCenter;
    // total number of nodes
    unsigned int totPossibleNodes;
    // total mass (initial)
    double totalMass;
    // standard neighbors shifting
    intList ne;
    unsIntList shift;
    intList domain;
    // node types (boundary are located halfway between nodes)
    //typeList types;
    //vecList positions;
    // curved info container
    //curveList curves; 
    // boundary conditions
    typeList boundary;
    // list with nodes
    //nodeList nodes;
    // indices of active nodes
    nodeList activeNodes;
    // indices of fluid nodes
    nodeList fluidNodes;
    // indices of interface nodes
    nodeList interfaceNodes;
    // indices of wall nodes (only those in potential contact with the fluid)
    nodeList wallNodes;
    // energy
    energy fluidEnergy, fluidImmersedEnergy;
    // for performance check
    struct timeval startLBStep, endLBStep;
    struct timeval startFreeSurfaceStep, endFreeSurfaceStep;
    struct timeval startCouplingStep, endCouplingStep;
    struct timeval startUpdateMassStep, endUpdateMassStep;
    struct timeval startUpdateInterfaceStep, endUpdateInterfaceStep;
    struct timeval startFindMutantsStep, endFindMutantsStep;
    struct timeval startSmoothenInterfaceStep_1, endSmoothenInterfaceStep_1;
    struct timeval startSmoothenInterfaceStep_2, endSmoothenInterfaceStep_2;
    struct timeval startUpdateMutantsStep, endUpdateMutantsStep;
    struct timeval startRemoveIsolatedStep, endRemoveIsolatedStep;
    struct timeval startRedistributeMassStep, endRedistributeMassStep;
    struct timeval startComputeNormalsStep, endComputeNormalsStep;
public:
    // switcher for restart
    bool lbRestart;
    // restart file
    string lbRestartFile;
    // switcher for imposed volume, and imposed volume
    bool imposeFluidVolume;
    double imposedFluidVolume;
    // switcher for increasing volume, and extra volume and time to add it
    bool increaseVolume;
    double deltaVolume;
    double deltaTime;
    // switcher for topography
    bool lbTopography;
    // switcher for topographic initial level
    bool lbTopographySurface;
    // shifts for topography file
    double translateTopographyX;
    double translateTopographyY;
    double translateTopographyZ;
    // topography file
    string lbTopographyFile;
    // topography container
    topography lbTop;
    // switchers for force field, non-Newtonian and everything
    bool freeSurface;
    bool forceField;
    bool TRTsolver;
    double magicNumber;
    // number of LBM step before activating the free surface
    unsigned int lbmInitialRepeat;
    // absolute time
    unsigned int time;
    // lbm size in cell units
    unsIntList lbSize;
    // lbm size in physical units (with boundaries))
    doubleList lbPhysicalSize;
    // lbm size in physical units (without boundaries))
    doubleList lbInnerPhysicalSize;
    // lbm box boundary locations
    vecList lbBoundaryLocation;
    // standard borders for free surface
    int freeSurfaceBorders[6];
    
    // conversion units /////////////////////////////////////////////////////////////////////////
    // fluid and granular matter are solved in different measure units
    // for a reference, check Feng, Han, Owen, 2007
    measureUnits unit;
    // problem-specific stuff ///////////////////////////////
    // stuff for shear cell
    double maxVisc, maxPlasticVisc, maxYieldStress, maxShearRate;
    unsigned int viscSteps, shearRateSteps;
    // stuff for drum
    double fluidMass;
    // stuff for avalanches (net-like)
    double avalanchePosit;
    // stuff for hourglass (mirrors in DEM)
    double hourglassOutletHeight;
    // stuff for continuum heap (mirrors in DEM)
    double heapBaseLevel;
    // stuff for Usman
    double largeFlumeFlowLevel;
    // lists for "mutant" nodes
    nodeList filledNodes, emptiedNodes;
    // interface nodes created due to interface evolution
    nodeList newInterfaceNodes;
    // age coefficients for new nodes
    const double maxAge=200.0;
    const double ageRatio=1.0/maxAge;
    // standard constructor
public:

    LB() {
        //
        freeSurface = false;
        forceField = false;
        lbRestart = false;
        //
        time = 0;
        lbSize.resize(3);
        lbSize[0] = lbSize[1] = lbSize[2] = 1;
        //
        slipCoefficient = 0.0;
        //
        initVelocity.reset();
        lbF.reset();
        //
        totPossibleNodes = 0;
        //
        shift.resize(3);
        shift[0] = shift[1] = shift[2] = 1;
        domain.resize(3);
        domain[0] = domain[1] = domain[2] = 1;
        //
        boundary.resize(6);
        ne.resize(lbmDirec);
        // initialize variable containers
        nodes.clear();
        //types.clear();
        // initialize lists
        activeNodes.clear();
        fluidNodes.clear();
        interfaceNodes.clear();
        wallNodes.clear();
        //
        lbmInitialRepeat = 0;
        // stuff for shear cell
        maxVisc = 0.0;
        maxPlasticVisc = 0.0;
        maxYieldStress = 0.0;
        maxShearRate = 0.0;
        viscSteps = 0;
        shearRateSteps = 0;
        // stuff for drum
        fluidMass = 0.0;
        // stuff for Usman
        largeFlumeFlowLevel=0.0;
        
    }
    // showing Lattice characteristics
    void LBShow() const;
    void latticeDefinition();
    void latticeBoltzmannGet(GetPot& lbmCfgFile, GetPot& command_line);
    void latticeBolzmannInit(cylinderList& cylinders, wallList& walls, particleList& particles, objectList& objects, const bool externalSolveCoriolis, const bool externalSolveCentrifugal);
    void latticeBolzmannStep(elmtList& elmts, particleList& particles, wallList& walls, objectList& objects);
    void latticeBoltzmannCouplingStep(bool& newNeighborList, elmtList& eltms, particleList& particles);
    void latticeBoltzmannFreeSurfaceStep();
    void accelerateDrum(double& drumSpeed, cylinderList& cylinders, wallList& walls);
    //
    void updateEnergy(double& totalKineticEnergy);
public:
    // node creation and deletion
    //void createNode(const unsigned int& it);
    //void deleteNode(const unsigned int& it);
    void generateNode(const unsigned int& coord, const types& type);
    void eraseNode(const unsigned int& coord);
    // initialize neighbors of a node coherently with the lattice and boundary conditions
    void findNeighbors(unsigned int neighborCoord[], const node* nodeHere);
    //void initializeNodeNeighbors(const int& it, node*& nodeHere);
    // initialization functions
    void initializeTypes(wallList& walls, cylinderList& cylinders, objectList& objects);
    void initializeNodes();
    void initializeLatticeBoundaries();
    void initializeParticleBoundaries(particleList& particles);
    void checkNewInterfaceParticles(elmtList& elmts, particleList& particles);
    void initializeWallBoundaries(wallList& walls);
    void initializeObjectBoundaries(objectList& objects);
    void initializeCylinderBoundaries(cylinderList& cylinders);
    void initializeTopography();
    void initializeCurved(cylinderList& cylinders);
    void setTopographySurface();
    void initializeInterface(double totParticles);
    void restartInterface(ifstream& fluidFileID, unsigned int& totNodes);
    void initializeLists();
    void resetLists();
    void initializeVariables();
    void initializeWalls(wallList& walls, cylinderList& cylinders, objectList& objects);
    // integration functions
    void reconstruction();
    void collision();
    void streaming(wallList& walls, objectList& objects);
    double curvedWallReconstruction(const unsigned int& j, const node* nodeHere, const tVect& wallSpeed) const;
    // error dumping functions
    void cleanLists();
    // coupling functions
    void computeHydroForces(elmtList& elmts, particleList& particles);
    void findNewActive(nodeList& newPopUpNodes, elmtList& elmts, particleList& particles);
    void findNewSolid(nodeList& newSolidNodes, elmtList& elmts, particleList& particles);
    void activeToSolid(unsIntList& newSolidNodes, elmtList& elmts, double& massSurplus);
    void solidToActive(unsIntList& newPopUpNodes, elmtList& elmts, double& massSurplus);
    void updateBoundary(unsIntList& newPopUpNodes, unsIntList& newSolidNodes, elmtList& elmts);
    //    void updateBoundaryOld(intList& newPopUpNodes, intList& newSolidNodes, particle& dummyParticle);
    // interface functions
    void updateMass();
    void updateInterface();
    void findInterfaceMutants();
    void smoothenInterface(double& massSurplus);
    void averageFromNeighbors(node* linkNode, double& newNodeN, tVect& newNodeU,double& newNodeVisc, double newNodeDistributions[]);
    void updateMutants(double& massSurplus);
    void removeIsolated(double& massSurplus);
    void redistributeMass(const double& massSurplus);
    void computeSurfaceNormal();
    void enforceMassConservation();
    
public:
    // function for linearized index management
    double maxHeight(const unsigned int& dir) const;
    unsigned int getIndex(const unsigned int& x, const unsigned int& y, const unsigned int& z);
    tVect getPosition(const unsigned int& index) const;
    unsigned int getX(const unsigned int& index) const;
    unsigned int getY(const unsigned int& index) const;
    unsigned int getZ(const unsigned int& index) const;
    double getPositionX(const unsigned int& index) const;
    double getPositionY(const unsigned int& index) const;
    double getPositionZ(const unsigned int& index) const;

};

#endif /* LB_H */
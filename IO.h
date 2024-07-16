/*
 * File:   IO.h
 * Author: aleonard
 *
 * Created on November 5, 2013, 4:06 PM
 */

#ifndef IO_H
#define IO_H


#include <cstdlib>
#include <iostream>
#include <vector>
#include <climits>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <omp.h>
#include <ctime>
#include <csignal>
#include <iomanip>
// only to create a directory
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
//
#include "LB.h"
#include "DEM.h"
#include "myvector.h"
#include "elmt.h"
#include "node.h"

using namespace std;
extern ProblemName problemName;

// define precision outputs for statistic file
#define wSt int(15)

class IO {
    
public:
    // simulation time in real units
    double realTime;
    // simulation iterations so far
    unsigned int currentTimeStep;
    // switchers for solvers
    bool demSolver, lbmSolver;
    // scwitchers for rotation of the reference frame
    bool coriolisSolver, centrifugalSolver;
    // locations
    string workDirectory, partDirectory, fluidDirectory;
    // time at invocation
    struct tm *now;
    // config file
    string configFileName;   
    // time integration: max number of iterations
    unsigned int maximumTimeSteps;
    // time integration: maximum time
    double maxTime;
    // intervals for output
    double screenExpTime, stateExpTime, fluidExpTime, fluidLagrangianExpTime, fluid2DExpTime, partExpTime, objectExpTime, singleObjectExpTime, cylinderExpTime;
    // recycle intervals
    double fluidRecycleExpTime, partRecycleExpTime, outputExpTime;
    // file for storing indices of single objects to export
    unsIntList singleObjects;
    // FOR HK_SMALL AND HK_LARGE
    // file for storing beginning and end of observation windows
    doubleList flowLevelBegin, flowLevelEnd;
    // file for storing beginning and end of observation windows
    unsIntList objectGroupBegin, objectGroupEnd;
    // general functions
    // 
    double energyStopThreshold;
    bool energyStop;
    double totalKineticEnergy;
    const unsigned int minimumIterations=1000;
    bool energyExit=false;
    
    void initialize();
    void outputStep(LB& lb, DEM& dem); // DEM used to be passed as const
    void outputFinal();
    
private:
    
    // generic file output streams
    string exportFileName;
    ofstream exportFile;
    unsigned int lastScreenExp;
    
    // general export stuff
    void exportMeanViscosity(const LB& lb);
    double meanViscosity(const LB& lb) const;
    double totParticleMass(const elmtList& elmts) const;
    
    // function that groups file creations
    void createFiles(const LB& lb, const DEM& dem);
    
    // RECYCLE /////////////////////////////////////////////////////////////////////////////////////////
    
    // particle recycle file
    unsigned int lastPartRecycleExp;
    string partRecycleFileFormat;
    void exportRecycleParticles(const elmtList& elmts, const pbcList& pbcs, const string& partRecycleFile);
    
    // fluid recycle file    
    unsigned int lastFluidRecycleExp;
    string fluidRecycleFileFormat;
    void exportRecycleFluid(const LB& lb, const string& fluidRecycleFile);
    
    // PARAVIEW /////////////////////////////////////////////////////////////////////////////////////////
    
    // particle paraview file
    unsigned int lastPartExp;
    string partFileFormat;
    void exportParaviewParticles(const elmtList& elmts, const particleList& particles, const string& particleFile);
    
    // object paraview file
    unsigned int lastObjectExp;
    string objectFileFormat;
    void exportParaviewObjects(const objectList& particles, const string& objectFile);
    
    // Cylinder paraview file
    unsigned int lastCylinderExp;
    string cylinderFileFormat;
    void exportParaviewCylinders(const cylinderList& cylinders, const string& cylinderFile);

    // Eulerian fluid paraview file
    unsigned int lastFluidExp;
    string fluidFileFormat;
    void exportEulerianParaviewFluid(const LB& lb, const string& fluidFile);
    
    // Lagrangian fluid paraview file
    unsigned int lastFluidLagrangianExp;
    string fluidLagrangianFileFormat;
    void exportLagrangianParaviewFluid(const LB& lb, const string& fluidFile);
    
    // 2D files for GIS/////////////////////////////////////////////////////////////////////////////////////////
    
    unsigned int lastFluid2DExp;
    string fluid2DFileFormat;
    void initialize2DFile(const LB& lb);
    void update2DFile(const LB& lb);
    void create2DFile(const LB& lb, const string& planarFile);
    
    // 2D planar files - istantaneous values
    doubleSet planarHeight, planarVel, planarPointVel;
    // topography
    doubleSet planarLevel;
    // global maxima quantities
    doubleSet maxPlanarHeight, maxPlanarVel, maxPlanarPointVel;
    // switcher for initialization
    double initializedPlanarFile;
    
    
    // GENERAL OUTPUT FILES ////////////////////////////////////////////////////////////////////////////////////
    
    // single-object data file
    string singleObjectDirectory;
    void exportSingleObjects(const objectList& objects);
    
    // max speed files
    string maxFluidSpeedFileName, maxParticleSpeedFileName;
    ofstream maxFluidSpeedFile, maxParticleSpeedFile;
    void exportMaxSpeedFluid(const LB& lb);
    void exportMaxSpeedParticles(const DEM& dem);
    
    // fluid flow rate file
    string freeSurfaceExtentFileName;
    ofstream freeSurfaceExtentFile;
    void exportFreeSurfaceExtent(const LB& lb);
    
    // fluid flow rate file
    string fluidFlowRateFileName;
    ofstream fluidFlowRateFile;
    void exportFluidFlowRate(const LB& lb);
    
   // particle flow rate file
    string particleFlowRateFileName;
    ofstream particleFlowRateFile;
    void exportParticleFlowRate(const DEM& dem);
    
    // fluid plasticity file
    string plasticityFileName;
    ofstream plasticityFile;
    void exportPlasticity(const LB& lb);
    double totPlastic(const LB& lb) const;
    
    // fluid mass file
    string fluidMassFileName;
    ofstream fluidMassFile;
    void exportFluidMass(const LB& lb);
    double totFluidMass(const LB& lb) const;
    
    // particle center of mass file
    string particleCenterOfMassFileName;
    ofstream particleCenterOfMassFile;
    void exportParticleCenterOfMass(const DEM& dem);
    tVect particleCenterOfMass(const elmtList& elmts) const;
    
    // particle coordination
    string particleCoordinationFileName;
    ofstream particleCoordinationFile;
    void exportParticleCoordination(const DEM& dem);
    
    // fluid center of mass file
    string fluidCenterOfMassFileName;
    ofstream fluidCenterOfMassFile;
    void exportFluidCenterOfMass(const LB& lb);
    tVect fluidCenterOfMass(const LB& lb) const;

    // particle overlap files
    string overlapFileName, dtOverlapFileName;
    ofstream overlapFile, dtOverlapFile;
    void exportParticleOverlap(DEM& dem); // DEM used to be passed as const
        
    // particle force files
    string forceFileName, wallForceFileName, maxWallForceFileName;
    ofstream forceFile, wallForceFile, maxWallForceFile;
    double hydraulicForceTot(const elmtList& elmts) const;
    double collisionForceTot(const elmtList& elmts) const;
    void exportForces(const DEM& dem);
    void exportWallForce(const DEM& dem);
    
    // particle energy file
    string energyFileName;
    void exportEnergy(const DEM& dem, const LB& lb);
    
    // PROBLEM-SPECIFIC OUTPUT FILES ////////////////////////////////////////////////////////////////////////////////////
    
    // SHEARCELL /////
    string  viscFileName;
    ofstream viscFile;
    // time-averaged quantites
    double msdSum, msdSum2;
    unsigned int msdSumN;
    tVect msdVecSum, msdVecSum2;
    void initViscFile();
    void writeViscFile(const LB& lb, const wallList& walls, const elmtList& elmts);
    void exportShearCell(const LB& lb, const DEM& dem);
    void apparentViscosity(const LB& lb, const wallList& walls, double& externalShear, double& wallStress, double& appVisc) const;

    // AVALANCHE /////
    string obstacleFileName;
    void exportForceObstacle(const objectList& objects);
    void exportGroupForce(const objectList& objects);
    tVect totForceObject(const objectList& objects, const unsigned int& groupBegin, const unsigned int& groupEnd) const;

    // OPENBARRIER /////
    string obstacleMomentFileName, obstacleHeightFileName, obstacleKmFileName;
    void exportForceObstacleElement(const DEM& dem);
    void exportMomentObstacleElement(const DEM& dem);
    void exportHPartObstacle(const DEM& dem);
    void exportKmPartObstacle(const DEM& dem);
    tVect totForceObjectElement(const objectList& objects, const unsigned int& obstacleIndex) const;
    tVect totMomentObjectElement(const wallList& walls, const objectList& objects, const unsigned int& obstacleIndex)const;
    
    // HONG KONG /////
    string hongKongFlowFileName, hongKongForceFileName;
    void exportHongKongFlow(DEM& dem);
    void exportHongKongBarrier(DEM& dem);
    double getSurface(doubleList& surfaceParticles);
    
    // USMAN /////
    void exportFlowLevel(const LB& lb);
    string frontFileName;
    void exportFront(const LB& lb);
    
    // MANGENEY /////
    string topFileName;
    void exportTop(const LB& lb);
    
    // INCLINE FLOW /////
    string inclineFlowFileName;
    void exportInclineFlow(const LB& lb);
    
    // HEAP /////
    string heapFileName;
    void exportHeapHeight(const LB& lb);

    // TRIAXIAL /////
    string triaxialFileName;
    void exportTriaxial(DEM& dem);
   

};

#endif /* IO_H */
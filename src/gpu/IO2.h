#ifndef IO2_H
#define IO2_H

#include "IO.h"

class LB2;

class IO2 : public IO {
    
 public:

     // simulation time
     double realTime;
     unsigned int currentTimeStep;

     // solver switches
     bool demSolver, lbmSolver;
     bool coriolisSolver, centrifugalSolver;

     // directories and configuration
     string workDirectory, partDirectory, fluidDirectory;
     struct tm* now;
     string configFileName;

     // simulation duration
     unsigned int maximumTimeSteps;
     double maxTime;

     // stopping criteria
     double energyStopThreshold;
     bool energyStop;
     double totalKineticEnergy;
     const unsigned int minimumIterations = 1000;
     bool energyExit = false;

     double fluidLagrangianExpTime;
     double partExpTime;
     double objectExpTime;
     double partRecycleExpTime;
     double fluidRecycleExpTime;

     // output formats
     enum ParaviewFormat {
         Binary,
         BinaryLowMem,
         Ascii
     };

     ParaviewFormat fluidLagrangianFormat = Binary;
     ParaviewFormat partExpFormat = Binary;

     // selected objects and elements
     unsIntList singleObjects;
     unsIntList singleElements;

     // object observation groups
     unsIntList objectGroupBegin, objectGroupEnd;

     void initialize();
     void outputStep(LB2& lb, DEM& dem);
     void outputFinal();

 private:

     // console and log output
     string exportFileName;
     ofstream exportFile;
     unsigned int lastScreenExp;

    // function that groups file creations
    void createFiles(LB2& lb, const DEM& dem);

    // PARAVIEW /////////////////////////////////////////////////////////////////////////////////////////
    
    // particle paraview file
    unsigned int lastPartExp;
    string partFileFormat;
    void exportParaviewParticles(const elmtList& elmts, const particleList& particles, const string& particleFile);
    void exportParaviewParticles_binaryv3(const elmtList& elmts, const particleList& particles, const string& particleFile);

    // object paraview file
    unsigned int lastObjectExp;
    string objectFileFormat;
    void exportParaviewObjects(const objectList& particles, const string& objectFile);

    // Eulerian fluid paraview file
    unsigned int lastFluidExp;
    string fluidFileFormat;
    // void exportEulerianParaviewFluid(LB2& lb, const string& fluidFile);
    // void exportEulerianParaviewFluid_binary(LB2& lb, const string& fluidFile);
    // void exportEulerianParaviewFluid_binaryv2(LB2& lb, const string& fluidFile);
    void exportEulerianParaviewFluid_binaryv3(LB2& lb, const string& fluidFile); //

    // Lagrangian fluid paraview file
    unsigned int lastFluidLagrangianExp;
    string fluidLagrangianFileFormat;
    void exportLagrangianParaviewFluid(LB2& lb, const string& fluidFile);
    void exportLagrangianParaviewFluid_binaryv3(LB2& lb, const string& fluidFile);
    
    // RESTART & RECYCLE /////////////////////////////////////////////////////////////////////////////////////////
    
    // particle recycle file
    unsigned int lastPartRecycleExp;
    string partRecycleFileFormat;
    void exportRecycleParticles(const elmtList& elmts, const pbcList& pbcs, const string& partRecycleFile);

    // fluid recycle file
    unsigned int lastFluidRecycleExp;
    string fluidRecycleFileFormat;

    // OTHER FILES - LBM /////////////////////////////////////////////////////////////////////////////////////////
    
    // maximum fluid speed
    string maxFluidSpeedFileName;
    ofstream maxFluidSpeedFile;
    void exportMaxSpeedFluid(LB2& lb);

    // free surface extent
    string freeSurfaceExtentFileName;
    ofstream freeSurfaceExtentFile;
    void exportFreeSurfaceExtent(LB2& lb);

    // fluid flow rate
    string fluidFlowRateFileName;
    ofstream fluidFlowRateFile;
    void exportFluidFlowRate(LB2& lb);

    // fluid mass
    string fluidMassFileName;
    ofstream fluidMassFile;
    void exportFluidMass(LB2& lb);

    // fluid centre of mass
    string fluidCenterOfMassFileName;
    ofstream fluidCenterOfMassFile;
    void exportFluidCenterOfMass(LB2& lb);

    // plasticity
    string plasticityFileName;
    ofstream plasticityFile;
    void exportPlasticity(LB2& lb);
    void exportMeanViscosity(LB2& lb);

    // OTHER FILES - DEM /////////////////////////////////////////////////////////////////////////////////////////
    
    // particle flow rate file
    string particleFlowRateFileName;
    ofstream particleFlowRateFile;
    void exportParticleFlowRate(const DEM& dem);

    // maximum particle speed
    string maxParticleSpeedFileName;
    ofstream maxParticleSpeedFile;
    void exportMaxSpeedParticles(const DEM& dem);

    // particle centre of mass
    string particleCenterOfMassFileName;
    ofstream particleCenterOfMassFile;
    void exportParticleCenterOfMass(const DEM& dem);
    tVect particleCenterOfMass(const elmtList& elmts) const;

    // particle coordination
    string particleCoordinationFileName;
    ofstream particleCoordinationFile;
    void exportParticleCoordination(const DEM& dem);

    // particle overlap
    string overlapFileName, dtOverlapFileName;
    ofstream overlapFile, dtOverlapFile;
    void exportParticleOverlap(DEM& dem); // DEM used to be passed as const

    // particle and wall forces
    string forceFileName, wallForceFileName, maxWallForceFileName;
    ofstream forceFile, wallForceFile, maxWallForceFile;
    void exportWallForce(const DEM& dem);
    void exportForces(const DEM& dem);
    double hydraulicForceTot(const elmtList& elmts) const;
    double collisionForceTot(const elmtList& elmts) const;

    // object forces
    string obstacleFileName;
    void exportForceObstacle(const objectList& objects);
    tVect totForceObject(const objectList& objects, const unsigned int& groupBegin, const unsigned int& groupEnd) const;
    
    // selected DEM outputs
    string singleObjectDirectory;
    void exportSingleObjects(const objectList& objects);

    string singleElementDirectory;
    void exportSingleElements(const elmtList& elmts);
    void exportGroupForce(const objectList& objects);

    // OTHER FILES - LBM-DEM /////////////////////////////////////////////////////////////////////////////////////////
    // coupled energy output
    string energyFileName; 
    void exportEnergy(DEM& dem, LB2& lb);

};

#endif /* IO_H */
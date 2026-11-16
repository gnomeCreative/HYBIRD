#ifndef IO2_H
#define IO2_H

#include "IO.h"

class LB2;

class IO2 : public IO {
    
 public:    
    void initialize();
    void outputStep(LB2& lb, DEM& dem);
    void outputFinal();

 private:
    // function that groups file creations
    void createFiles(LB2& lb, const DEM& dem);

    // PARAVIEW /////////////////////////////////////////////////////////////////////////////////////////
    
    // particle paraview file
    void exportParaviewParticles(const elmtList& elmts, const particleList& particles, const string& particleFile);
    void exportParaviewParticles_binaryv3(const elmtList& elmts, const particleList& particles, const string& particleFile);

    
    // Eulerian fluid paraview file
    // void exportEulerianParaviewFluid(LB2& lb, const string& fluidFile);
    // void exportEulerianParaviewFluid_binary(LB2& lb, const string& fluidFile);
    // void exportEulerianParaviewFluid_binaryv2(LB2& lb, const string& fluidFile);
    void exportEulerianParaviewFluid_binaryv3(LB2& lb, const string& fluidFile); //

    // Lagrangian fluid paraview file
    void exportLagrangianParaviewFluid(LB2& lb, const string& fluidFile);
    void exportLagrangianParaviewFluid_binaryv3(LB2& lb, const string& fluidFile);

    // OTHER FILES - LBM /////////////////////////////////////////////////////////////////////////////////////////
    
    // fluid export functions (.dat files)
    void exportMaxSpeedFluid(LB2& lb);
    void exportFreeSurfaceExtent(LB2& lb);
    void exportFluidFlowRate(LB2& lb);
    void exportFluidMass(LB2& lb);
    void exportFluidCenterOfMass(LB2& lb);

    // viscosity and plasticity
    void exportPlasticity(LB2& lb);
    void exportMeanViscosity(LB2& lb);

    // OTHER FILES - DEM /////////////////////////////////////////////////////////////////////////////////////////
    // particle flow rate file
    void exportParticleFlowRate(const DEM& dem);
    void exportMaxSpeedParticles(const DEM& dem);
    void exportParticleCenterOfMass(const DEM& dem);
    tVect particleCenterOfMass(const elmtList& elmts) const;
    void exportParticleCoordination(const DEM& dem);
    void exportParticleOverlap(DEM& dem); // DEM used to be passed as const
    void exportWallForce(const DEM& dem);
    void exportForces(const DEM& dem);
    double hydraulicForceTot(const elmtList& elmts) const;
    double collisionForceTot(const elmtList& elmts) const;
    void exportForceObstacle(const objectList& objects);
    tVect totForceObject(const objectList& objects, const unsigned int& groupBegin, const unsigned int& groupEnd) const;

    // OTHER FILES - LBM-DEM /////////////////////////////////////////////////////////////////////////////////////////
    void exportEnergy(DEM& dem, LB2& lb);

};

#endif /* IO_H */
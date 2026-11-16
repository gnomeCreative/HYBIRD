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
    void exportParticleCenterOfMass(const DEM& dem);
    tVect particleCenterOfMass(const elmtList& elmts) const;
    void exportParticleCoordination(const DEM& dem);

    // OTHER FILES - LBM-DEM /////////////////////////////////////////////////////////////////////////////////////////
    void exportEnergy(DEM& dem, LB2& lb);
    void exportMaxSpeedParticles(const DEM& dem);
};

#endif /* IO_H */
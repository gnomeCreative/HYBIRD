#ifndef IO2_H
#define IO2_H

#include "IO.h"

class LB2;
class DEM2;

class IO2 : public IO {
    
 public:    
    void outputStep(LB2& lb, DEM2& dem);

 private:
    // function that groups file creations
    void createFiles(LB2& lb, const DEM2& dem);

    // PARAVIEW /////////////////////////////////////////////////////////////////////////////////////////
    
    // particle paraview file
    // void exportParaviewParticles(const elmtList& elmts, const particleList& particles, const string& particleFile);
    void exportParaviewParticles_binaryv3(DEM2& dem, const string& particleFile);

    
    // Eulerian fluid paraview file
    // void exportEulerianParaviewFluid(LB2& lb, const string& fluidFile);
    // void exportEulerianParaviewFluid_binary(LB2& lb, const string& fluidFile);
    // void exportEulerianParaviewFluid_binaryv2(LB2& lb, const string& fluidFile);
    void exportEulerianParaviewFluid_binaryv3(LB2& lb, const string& fluidFile);

    // Lagrangian fluid paraview file
    void exportLagrangianParaviewFluid(LB2& lb, const string& fluidFile);
    void exportLagrangianParaviewFluid_binaryv3(LB2& lb, const string& fluidFile);
};

#endif /* IO_H */
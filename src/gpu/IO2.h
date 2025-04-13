#ifndef IO2_H
#define IO2_H

#include "IO.h"

class LBOpenMP;

class IO2 : public IO {
    
 public:    
    void outputStep(LBOpenMP& lb, DEM& dem);

 private:
    // function that groups file creations
    void createFiles(LBOpenMP& lb, const DEM& dem);

    // PARAVIEW /////////////////////////////////////////////////////////////////////////////////////////
    
    // particle paraview file
    void exportParaviewParticles(const elmtList& elmts, const particleList& particles, const string& particleFile);
    void exportParaviewParticles_binaryv3(const elmtList& elmts, const particleList& particles, const string& particleFile);

    
    // Eulerian fluid paraview file
    // void exportEulerianParaviewFluid(LB2& lb, const string& fluidFile);
    // void exportEulerianParaviewFluid_binary(LB2& lb, const string& fluidFile);
    // void exportEulerianParaviewFluid_binaryv2(LB2& lb, const string& fluidFile);
    void exportEulerianParaviewFluid_binaryv3(LBOpenMP& lb, const string& fluidFile);

    // Lagrangian fluid paraview file
    void exportLagrangianParaviewFluid(LBOpenMP& lb, const string& fluidFile);
    void exportLagrangianParaviewFluid_binaryv3(LBOpenMP& lb, const string& fluidFile);

    void exportMaxSpeedFluid(LBOpenMP& lb);
};

#endif /* IO_H */
#ifndef DEM2_H
#define DEM2_H

#include "Cylinder2.h"
#include "DEMParams.h"
#include "Element2.h"
#include "Object2.h"
#include "Particle2.h"
#include "Wall2.h"

/**
 * V2 of DEM.cpp
 * An experimental attempt to produce a combined CPU-CUDA codebase
 * which can be compiled as either CPU or GPU code with as much shared code as possible
 */
class DEM2 {
    // particles contains both standard particles and ghost particles. The total number of standard particles needs to be saved
    // @todo this might move into Particle2
    unsigned int stdPartNumber = 0;
    // number of active particles and elements
    unsigned int actvPartNumber = 0;
    unsigned int actvElmtNumber = 0;

    // Total number of standard objects (i.e. not ghosts) which needs to be saved (Devis))
    unsigned int stdObjects;
 public:
    void discreteElementInit(DEMParams& dem_p, const std::array<types, 6>& externalBoundary, const std::array<double, 3>& externalSize, const std::array<tVect, 6>& externalBoundaryLocation,
        const tVect &externalAccel, const tVect &externalRotation, const tVect &externalRotationCenter, bool externalSolveCoriolis, bool externalSolveCentrifugal, double externalTimeStep);

 public: // Public for convenience
    Particle2 h_particles, hd_particles, * d_particles = nullptr;
    Wall2 h_walls, hd_walls, * d_walls = nullptr;
    Cylinder2 h_cylinders, hd_cylinders, * d_cylinders = nullptr;
    Element2 h_elements, hd_elements, *d_elements = nullptr;
    Object2 h_objects, hd_objects, *d_objects = nullptr;
};
#endif  // DEM2_H
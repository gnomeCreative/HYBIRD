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
    void discreteElementInit(const std::array<types, 6>& externalBoundary, const std::array<double, 3>& externalSize, const std::array<tVect, 6>& externalBoundaryLocation,
        const tVect &externalAccel, const tVect &externalRotation, const tVect &externalRotationCenter, bool externalSolveCoriolis, bool externalSolveCentrifugal, double externalTimeStep);
    /**
     * Initializes all parameters useful for neighbor list algorithm
     * cellWidth, nCells, nebrRange, maxDisp
     */
    void initNeighborParameters();
    void determineTimeStep(double externalTimeStep);
    double criticalTimeStep() const;
    /**
     * (re)construct the neighbour table
     */
    template<int impl>
    void evalNeighborTable();

    void evolveBoundaries();
    void discreteElementStep();

    ///
    /// Discrete element step functions
    ///
    template<int impl>
    void buildActiveLists();
    template<int impl>
    void predictor();
    template<int impl>
    void updateParticlesPredicted();
    template<int impl>
    void evaluateForces();
    template<int impl>
    void particleParticleContacts();
    template<int impl>
    void wallParticleContacts();
    template<int impl>
    void objectParticleContacts();
    template<int impl>
    void cylinderParticleContacts();
    template<int impl>
    void computeApparentForces();
    template<int impl>
    void saveObjectForces();
    template<int impl>
    void newtonEquationsSolution();
    template<int impl>
    void corrector();
    template<int impl>
    void updateParticlesCorrected();

    Particle2 &getParticles();
    Element2 &getElements();
    DEMParams &getParams();
    void syncParamsToDevice();
    /**
     * Copy h_particles to d_particles (and hd_particles)
     * This also allocates memory for d_particles if required
     */
    void syncParticlesToDevice();
    void syncParticlesFromDevice();
    void syncWallsToDevice();
    void syncCylindersToDevice();
    void syncElementsToDevice();
    void syncElementsFromDevice();
    void syncObjectsToDevice();
 public: // Public for convenience
    Particle2 h_particles, hd_particles, * d_particles = nullptr;
    Wall2 h_walls, hd_walls, * d_walls = nullptr;
    Cylinder2 h_cylinders, hd_cylinders, * d_cylinders = nullptr;
    Element2 h_elements, hd_elements, *d_elements = nullptr;
    Object2 h_objects, hd_objects, *d_objects = nullptr;
};
#endif  // DEM2_H
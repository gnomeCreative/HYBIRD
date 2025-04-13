#ifndef LBOPENMP_H
#define LBOPENMP_H

#include "Node2.h"
#include "Particle2.h"
#include "Element2.h"
#include "Object2.h"
#include "Cylinder2.h"
#include "Wall2.h"

class Problem;
class DEM;

/**
 * V3 of LB.cpp
 * An improved architecture to V2, which should hopefully make the CUDA subclass fully hybrid
 * and this version able to be compiled sans CUDA
 */
class LBOpenMP {
 public:
    LBOpenMP() = default;

    void init(Problem& problem, cylinderList& cylinders, wallList& walls, particleList& particles, objectList& objects, bool externalSolveCoriolis, bool externalSolveCentrifugal);
    /**
     * Execute all LB methods from hybird.cpp::goCycle()
     * 1. latticeBoltzmannCouplingStep() if io.demSolver
     * 2. latticeBolzmannStep()
     * 3. latticeBoltzmannFreeSurfaceStep() if lb.freeSurface
     * or if dem.demTime <= dem.demInitialRepeat
     * 1. latticeBoltzmannCouplingStep() if io.demSolver
     * @param dem The DEM instance coupled with LB model
     * @param io_demSolver The result of io.demSolver in the calling method
     * @note io_demSolver may be redundant, surely dem can be probed to detect if DEM is active
     */
    void step(DEM& dem, bool io_demSolver);

    Node2& getNodes() { return h_nodes; }

 private:
    //
    // Init chain
    //
    void allocateNodes(unsigned int count);
    void initializeLatticeBoundaries();
    void initializeTypes(const wallList& walls, const cylinderList& cylinders, const objectList& objects);
    void initializeWallBoundaries(const wallList& walls);
    void initializeObjectBoundaries(const objectList& objects);
    void initializeCylinderBoundaries(const cylinderList& cylinders);
    void initializeTopography();
    void initializeInterface(const Problem& problem);
    void initializeVariables();
    void generateNode(unsigned int coord, types typeHere);
    void initializeWalls();
    void initializeLists();

    //
    // Data ingress from legacy structures
    //
    /**
     * @brief Sync DEM data to structure of arrays format (and device memory)
     * @param elmts Objects, such as walls within the DEM
     * @param particles Spherical particles, which represent a decomposition of the elmts
     * @param walls Wall elmts??
     * @param objects Objects that are not walls, e.g. cylinders??
     * @note Could directly load to h_ structures in future, enabling removal of legacy code
     */
    void syncDEMIn(const elmtList& elmts, const particleList& particles, const wallList& walls, const objectList& objects);
    /**
     * @brief Sync DEM data to structure from arrays format (and device memory)
     * @param elmts Objects, such as walls within the DEM
     * @param particles Spherical particles, which represent a decomposition of the elmts
     * @param walls Wall elmts??
     * @param objects Objects that are not walls, e.g. cylinders??
     * @note Could directly use h_ structures in future, enabling removal of legacy code
     */
    void syncDEMOut(elmtList& elmts, particleList& particles, wallList& walls, objectList& objects) const;
    /**
     * @brief Sync DEM elements list to h_elements
     * @return true if h_elements.componentsData has grown (this may be redundant)
     * @note Could directly load to h_particles in future, enabling removal of legacy code
     */
    bool syncElementsIn(const elmtList& elements);
    /**
     * @brief Sync DEM particle list to h_particles
     * @note Could directly load to h_particles in future, enabling removal of legacy code
     */
    void syncParticlesIn(const particleList& particles);
    /**
     * @brief Sync DEM cylinder list to h_cylinders
     * @note Could directly load to h_cylinders in future, enabling removal of legacy code
     */
    void syncCylindersIn(const cylinderList& cylinders);
    /**
     * @brief Sync DEM wall list to h_walls
     * @note Could directly load to h_walls in future, enabling removal of legacy code
     */
    void syncWallsIn(const wallList& walls);
    /**
     * @brief Sync DEM object list to h_objects
     * @note Could directly load to h_objects in future, enabling removal of legacy code
     */
    void syncObjectsIn(const objectList& objects);
    /**
     * @brief Sync DEM elements list from h_elements
     * @note Could directly use h_elements in future, enabling removal of legacy code
     */
    void syncElementsOut(elmtList& elements) const;
    /**
     * @brief Sync DEM particle list from h_particles
     * @note Could directly use h_particles in future, enabling removal of legacy code
     */
    void syncParticlesOut(particleList& particles) const;
    /**
     * @brief Sync DEM cylinder list from h_cylinders
     * @note Could directly use h_cylinders in future, enabling removal of legacy code
     */
    void syncCylindersOut(cylinderList& cylinders) const;
    /**
     * @brief Sync DEM wall list from h_walls
     * @note Could directly use h_walls in future, enabling removal of legacy code
     */
    void syncWallsOut(wallList& walls) const;
    /**
     * @brief Sync DEM object list from h_objects
     * @note Could directly use h_objects in future, enabling removal of legacy code
     */
    void syncObjectsOut(objectList& objects) const;

    ///
    /// Step stages
    ///
    /**
     * @brief Following the DEM model being stepped, this updates impacted LBM nodes
     * @param newNeighbourList If a new neighbour table has been defined, the indexing will be reinitialised
     */
    void latticeBoltzmannCouplingStep(bool& newNeighbourList);
    /**
     * @brief The main LBM step
     */
    void latticeBoltzmannStep();
    /**
     * @brief free-surface  management functions
     */
    void latticeBoltzmannFreeSurfaceStep();

    //
    // latticeBoltzmannCouplingStep() subroutines
    //
    /**
     * @brief Update all Node2::solid_index to the contained particle
     * @return The total mass of nodes that are inside particles
     *
     * @note Called from latticeBoltzmannCouplingStep() when a new neighbour table has been defined
     * @see latticeBoltzmannCouplingStep()
     */
    double initializeParticleBoundaries();
    /**
     * @brief Check and update whether active nodes are still inside particles
     * @note Called from latticeBoltzmannCouplingStep() when a new neighbour table has not been defined
     * @see latticeBoltzmannCouplingStep()
     */
    void findNewActive();
    /**
     * @brief Find new nodes that are inside particles
     * Checks whether nodes that are neighbouring active nodes marked as inside particle
     * are also inside any particles within that same particles cluster.
     * If so they are marked as inside particle with their solid index updated
     * @see latticeBoltzmannCouplingStep()
     */
    void findNewSolid();
    /**
     * @brief ??
     * @see latticeBoltzmannCouplingStep()
     */
    void checkNewInterfaceParticles();

    //
    // latticeBoltzmannStep() subroutines
    //
    /**
     * @brief combined reconstruct(), computeHydroForces(), collision()
     * Reconstruct macroscopic variables from microscopic distribution
     * Compute interaction forces with DEM elmts
     * Collision step
     * @see latticeBoltzmannStep()
     */
    void reconstructHydroCollide();
    /**
     * @brief Streaming operator
     * @see latticeBoltzmannStep()
     */
    void streaming();
    /**
     * @brief Shift element/wall/object forces and torques to physical units
     * @see latticeBoltzmannStep()
     */
    void shiftToPhysical();

    //
    // latticeBoltzmannFreeSurfaceStep() subroutines
    //
    void enforceMassConservation();
    void redistributeMass(const double& massSurplus);
    void updateMass();
    void updateInterface();
    /**
     * Build a temporary list of h_nodes.interfaceI
     * @todo Profile and see if it's worth optimising this
     * @todo Could optimise to save memory realloc
     */
    std::vector<unsigned int> buildTempNewList(const unsigned int &_max_len);
    /**
     * Rebuild h_nodes.interfaceI
     * @todo Profile and see if it's worth optimising this
     * @todo Could optimise to save memory realloc
     */
    void buildInterfaceList(const unsigned int &_max_len);
    /**
     * Rebuild h_nodes.activeI, h_nodes.interfaceI & h_nodes.fluidI
     * @todo Profile and see if it's worth optimising this
     * @todo Could optimise to save memory realloc
     */
    void buildAllLists(const unsigned int &_max_interface_len, const unsigned int &_max_fluid_len);

 protected:
    // The actual node storage
    // Host copy of node buffers, may not always be current whilst in CUDA mode
    Node2 h_nodes;

    // The temporary host, and device particle storage
    // Until DEM model is moved to CUDA, host copy only acts as a location to build data before copying to device
    Particle2 h_particles;
    Element2 h_elements;
    Wall2 h_walls;
    Cylinder2 h_cylinders;
    Object2 h_objects;

    // Host only parameters used during initialisation
    LBInitParams init_params; 
    // topography container (used during initializeTopography()/initializeInterface())
    topography lbTop = {};
};

#endif // LBOPENMP_H

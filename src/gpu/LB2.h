#ifndef LB2_H
#define LB2_H

#include <array>

#include "LBParams.h"
#include "myvector.h"
#include "Node2.h"
#include "Particle2.h"
#include "Element2.h"
#include "Object2.h"
#include "Cylinder2.h"
#include "Wall2.h"

class DEM;

/**
 * V2 if LB.cpp
 * An experimental attempt to produce a combined CPU-CUDA codebase
 * which can be compiled as either CPU or GPU code with as much shared code as possible
 */
class LB2 {
    struct NewNode {
        types type;
        unsigned int solidIndex = std::numeric_limits<unsigned int>::max();
    };
    public:

    // @todo how do we init params from config?
    LB2() = default;

    void init(cylinderList &cylinders, wallList &walls, particleList &particles, objectList &objects, bool externalSolveCoriolis, bool externalSolveCentrifugal);
    void allocateHostNodes(unsigned int count);
    void initializeLatticeBoundaries();
    void initializeTypes(const wallList &walls, const cylinderList &cylinders, const objectList &objects);
    void initializeWallBoundaries(const wallList &walls);
    void initializeObjectBoundaries(const objectList &objects);
    void initializeCylinderBoundaries(const cylinderList &cylinders);
    void initializeTopography();
    void initializeInterface();
    void initializeVariables();
    void generateNode(unsigned int coord, types typeHere);
    void initializeWalls();
    void initializeLists();

    template<int impl>
    void buildInterfaceList(unsigned int max_len = std::numeric_limits<unsigned int>::max(), bool update_device_struct = true);
    template<int impl>
    void buildFluidList(unsigned int max_len = std::numeric_limits<unsigned int>::max(), bool update_device_struct = true);
    template<int impl>
    unsigned int *buildTempNewList(unsigned int max_len = std::numeric_limits<unsigned int>::max(), bool update_device_struct = true);
    template<int impl>
    void buildActiveList();
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
    void step(DEM &dem, bool io_demSolver);
    /**
     * @brief Sync DEM data to structure of arrays format (and device memory)
     * @param elmts Objects, such as walls within the DEM
     * @param particles Spherical particles, which represent a decomposition of the elmts
     * @param walls Wall elmts??
     * @param objects Objects that are not walls, e.g. cylinders??
     * @note Most of this will be removed once DEM is also moved to CUDA
     */
    void syncDEMIn(const elmtList &elmts, const particleList &particles, const wallList &walls, const objectList &objects);
    /**
     * @brief Sync DEM data to structure from arrays format (and device memory)
     * @param elmts Objects, such as walls within the DEM
     * @param particles Spherical particles, which represent a decomposition of the elmts
     * @param walls Wall elmts??
     * @param objects Objects that are not walls, e.g. cylinders??
     * @note Most of this will be removed once DEM is also moved to CUDA
     * @note Much of this synchronisation is redundant, as the LBM coupling does not impact most DEM properties
     */
    void syncDEMOut(elmtList &elmts, particleList &particles, wallList &walls, objectList &objects);
    /**
     * @brief Following the DEM model being stepped, this updates impacted LBM nodes
     * @param newNeighbourList If a new neighbour table has been defined, the indexing will be reinitialised
     */
    void latticeBoltzmannCouplingStep(bool &newNeighbourList);
    /**
     * @brief The main LBM step
     */
    void latticeBoltzmannStep();
    /**
     * @brief free-surface  management functions
     */
    void latticeBoltzmannFreeSurfaceStep();
    /**
     * @brief Sync DEM elements list to h_elements/d_elements
     * @note Should be redundant once DEM is also ported to CUDA
     */
    template<int impl>
    bool syncElementsIn(const elmtList &elements);
    /**
     * @brief Sync DEM particle list to h_particles/d_particles
     * @note Should be redundant once DEM is also ported to CUDA
     */
    template<int impl>
    void syncParticlesIn(const particleList &particles);
    /**
     * @brief Sync DEM cylinder list to h_cylinders/d_cylinders
     * @note Should be redundant once DEM is also ported to CUDA
     */
    template<int impl>
    void syncCylindersIn(const cylinderList &cylinders);
    /**
     * @brief Sync DEM wall list to h_walls/d_walls
     * @note Should be redundant once DEM is also ported to CUDA
     */
    template<int impl>
    void syncWallsIn(const wallList &walls);
    /**
     * @brief Sync DEM object list to h_objects/d_objects
     * @note Should be redundant once DEM is also ported to CUDA
     */
    template<int impl>
    void syncObjectsIn(const objectList &objects);
    /**
     * @brief Sync DEM elements list from h_elements/d_elements
     * @note Should be redundant once DEM is also ported to CUDA
     */
    template<int impl>
    void syncElementsOut(elmtList &elements);
    /**
     * @brief Sync DEM particle list from h_particles/d_particles
     * @note Should be redundant once DEM is also ported to CUDA
     */
    template<int impl>
    void syncParticlesOut(particleList &particles);
    /**
     * @brief Sync DEM cylinder list from h_cylinders/d_cylinders
     * @note Should be redundant once DEM is also ported to CUDA
     */
    template<int impl>
    void syncCylindersOut(cylinderList &cylinders);
    /**
     * @brief Sync DEM wall list from h_walls/d_walls
     * @note Should be redundant once DEM is also ported to CUDA
     */
    template<int impl>
    void syncWallsOut(wallList &walls);
    /**
     * @brief Sync DEM object list from h_objects/d_objects
     * @note Should be redundant once DEM is also ported to CUDA
     */
    template<int impl>
    void syncObjectsOut(objectList &objects);
    /**
     * Set the LBM params
     * This copies them to GPU storage if required
     * @params params The LBParams structure to overwrite PARAMS with
     *
     * @param skip_sync If true, will not be copied to device
     */
    void setParams(const LBParams& params, const LBInitParams& initParams, bool skip_sync = false);
    void syncParams();

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
    template<int impl>
    double initializeParticleBoundaries();
    /**
     * @brief Check and update whether active nodes are still inside particles
     * @note Called from latticeBoltzmannCouplingStep() when a new neighbour table has not been defined
     * @see latticeBoltzmannCouplingStep()
     */
    template<int impl>
    void findNewActive();
    /**
     * @brief Find new nodes that are inside particles
     * Checks whether nodes that are neighbouring active nodes marked as inside particle
     * are also inside any particles within that same particles cluster.
     * If so they are marked as inside particle with their solid index updated
     * @see latticeBoltzmannCouplingStep()
     */
    template<int impl>
    void findNewSolid();
    /**
     * @brief ??
     * @see latticeBoltzmannCouplingStep()
     */
    template<int impl>
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
    template<int impl>
    void reconstructHydroCollide();
    /**
     * @brief Streaming operator
     * @see latticeBoltzmannStep()
     */
    template<int impl>
    void streaming();
    /**
     * @brief Shift element/wall/object forces and torques to physical units
     * @see latticeBoltzmannStep()
     */
    template<int impl>
    void shiftToPhysical();


    ///
    /// latticeBoltzmannFreeSurfaceStep() subroutines
    ///
    template<int impl>
    void enforceMassConservation();
    template<int impl>
    void redistributeMass(const double &massSurplus);
    template<int impl>
    void updateMass();
    template<int impl>
    void updateInterface();

    Node2 &getNodes();
    void initDeviceNodes();

    /**
     * Member variable philosophy
     * h_<name> are "host" copies of the same information pointed to by <name>
     * In CPU builds, <name> should always point to h_<name>
     * In CUDA builds, <name> will point to device memory, and h_<name> may not have current data at all times
     */
    private:
    // The actual node storage
    // Host copy of node buffers, may not always be current whilst in CUDA mode
    Node2 h_nodes;
    // Host copy of device node buffer pointers
    Node2 hd_nodes;
    // Pointer to device copy of device node buffer pointers
    // In CPU, this is a pointer to h_nodes
    Node2 *d_nodes = nullptr;


    // The temporary host, and device particle storage
    // Until DEM model is moved to CUDA, host copy only acts as a location to build data before copying to device
    Particle2 h_particles, hd_particles, *d_particles = nullptr;
    Element2 h_elements, hd_elements, *d_elements = nullptr;
    Wall2 h_walls, hd_walls, *d_walls = nullptr;
    Cylinder2 h_cylinders, hd_cylinders, * d_cylinders = nullptr;
    Object2 h_objects, hd_objects, *d_objects = nullptr;
    
    LBInitParams init_params; // Host only parameters used during initialisation
    // topography container
    topography lbTop = {};  // t
};

#endif // LB2_H

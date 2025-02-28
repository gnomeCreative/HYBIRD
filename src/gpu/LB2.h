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

class DEM2;

/**
 * V2 of LB.cpp
 * An experimental attempt to produce a combined CPU-CUDA codebase
 * which can be compiled as either CPU or GPU code with as much shared code as possible
 */
class LB2 {
    struct NewNode {
        types type;
        unsigned int solidIndex = std::numeric_limits<unsigned int>::max();
    };
    public:

    LB2(DEM2& _dem) : dem(_dem) { };

    void init(DEM2 &dem, bool externalSolveCoriolis, bool externalSolveCentrifugal);
    void allocateHostNodes(unsigned int count);
    void initializeLatticeBoundaries();
    void initializeTypes(const Wall2 &h_walls, const Cylinder2 &h_cylinders, const Object2 &h_objects);
    void initializeWallBoundaries(const Wall2 &h_walls);
    void initializeObjectBoundaries(const Object2 &h_objects);
    void initializeCylinderBoundaries(const Cylinder2 &h_cylinders);
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
    unsigned int *buildTempNewList(unsigned int max_len = std::numeric_limits<unsigned int>::max());
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
    void step(bool io_demSolver);
    /**
     * @brief Following the DEM model being stepped, this updates impacted LBM nodes
     * @param newNeighbourList If a new neighbour table has been defined, the indexing will be reinitialised
     */
    void latticeBoltzmannCouplingStep(bool newNeighbourList);
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
     * @params params The LBParams structure to overwrite LB_P with
     *
     * @param skip_sync If true, will not be copied to device
     */
    void setParams(const LBParams& params, const LBInitParams& initParams, bool skip_sync = false);
    void syncParamsToDevice();

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
    // Reference to DEM, for accessing particles, elements, walls etc
    DEM2& dem;
    // The actual node storage
    // Host copy of node buffers, may not always be current whilst in CUDA mode
    Node2 h_nodes = {};
    // Host copy of device node buffer pointers
    Node2 hd_nodes = {};
    // Pointer to device copy of device node buffer pointers
    // In CPU, this is a pointer to h_nodes
    Node2 *d_nodes = nullptr;


    // The temporary host, and device particle storage
    // Until DEM model is moved to CUDA, host copy only acts as a location to build data before copying to device
    
    LBInitParams init_params; // Host only parameters used during initialisation
    // topography container
    topography lbTop = {};  // t
};

#endif // LB2_H

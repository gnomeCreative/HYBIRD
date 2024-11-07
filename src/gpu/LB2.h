#ifndef LB2_H
#define LB2_H

#include <array>

#include "LBParams.h"
#include "myvector.h"
#include "Node2.h"
#include "Particle2.h"
#include "Element2.h"

class DEM;

/**
 * V2 if LB.cpp
 * An experimental attempt to produce a combined CPU-CUDA codebase
 * which can be compiled as either CPU or GPU code with as much shared code as possible
 */
enum Implementation {CPU, CUDA};
class LB2 {
    public:

    // @todo how do we init params from config?
    LB2() = default;

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
    void step(const DEM &dem, bool io_demSolver);
    /**
     * @brief Identifies which nodes need to have an update due to particle movement
     * @param newNeighbourList If a new neighbour table has been defined, the indexing will be reinitialised
     * @param eltms ??????? Are these the walls?
     * @param particles The DEM particles
     */
    void latticeBoltzmannCouplingStep(bool &newNeighbourList, const elmtList& eltms, const particleList& particles);

    /**
     * @brief Sync DEM particle list to h_particles/d_particles
     * @note Should be redundant once DEM is also ported to CUDA
     */
    template<int impl>
    void syncParticles(const particleList& particles);
    /**
     * @brief Sync DEM elements list to h_elements/d_elements
     * @note Should be redundant once DEM is also ported to CUDA
     */
    template<int impl>
    bool syncElements(const elmtList& elements);
    /**
     * @brief Update all Node2::solid_index to the contained particle
     * @note Called from latticeBoltzmannCouplingStep() when a new neighbour table has been defined
     */
    template<int impl>
    void initializeParticleBoundaries();
    /**
     * @brief Check and update whether active nodes are still inside particles
     * @note Called from latticeBoltzmannCouplingStep() when a new neighbour table has not been defined
     */
    template<int impl>
    void findNewActive();
    /**
     * @brief Find new nodes that are inside particles
     * Checks whether nodes that are neighbouring active nodes marked as inside particle
     * are also inside any particles within that same particles cluster.
     * If so they are marked as inside particle with their solid index updated
     */
    template<int impl>
    void findNewSolid();
    /**
     * @brief ??
     */
    template<int impl>
    void checkNewInterfaceParticles();
    /**
     * Initialise dynamic lattice params that scale with the model configuration
     * @see lattice.h for the static lattice params
     */
    void latticeDefinition(); // @todo Call Params version and copy to PARAMS?

    void LBShow() { };


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
    Node2 *d_nodes;


    // The temporary host, and device particle storage
    // Until DEM model is moved to CUDA, host copy only acts as a location to build data before copying to device
    Particle2 h_particles, hd_particles, *d_particles;
    Element2 h_elements, hd_elements, *d_elements;
    
    LBParams h_params;  // Host copy of PARAMS

    // switchers for force field, non-Newtonian and everything
    bool freeSurface = false;
};

#endif // LB2_H

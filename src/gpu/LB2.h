#ifndef LB2_H
#define LB2_H

#include <array>

#include "lattice.h"

class DEM;

/**
 * V2 if LB.cpp
 * An experimental attempt to produce a combined CPU-CUDA codebase
 * which can be compiled as either CPU or GPU code with as much shared code as possible
 */
enum Implementation {CPU, CUDA};
class LB2 {
    public:
    /**
     * Lattice Boltzmann model parameters
     * These should not change during after initialisation
     */
    struct Params {
        // LB::lbSize
        std::array<unsigned int, 3> lbSize = {1,1,1};
        /**
         * standard neighbors shifting
         */
        // LB::ne
        std::array<int, lbmDirec> ne = {};
        // LB::shift
        std::array<unsigned int, 3> shift = {1,1,1};
        // LB::domain
        std::array<int, 3> domain = {1,1,1};
        /**
         * Initialise dynamic lattice params that scale with the model configuration
         * @see lattice.h for the static lattice params
         */
        void latticeDefinition();
    };

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
        Params h_params;  // @see PARAMS
};

#endif // LB2_H

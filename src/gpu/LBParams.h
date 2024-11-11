#ifndef LBPARAMS_H
#define LBPARAMS_H

#include <array>

#include "lattice.h"
#include "MeasureUnits.h"

/**
 * Lattice Boltzmann model parameters
 * These should not change during after initialisation
 * @note This struct has been moved into a seperate file to avoid a circular dependency
 *       whilst maximising how much code can be inlined
 */
struct LBParams {
    // hydrodynamic radius (see Kumnar et al., Mechanics of granular column collapse in fluid at varying slope angles)
    double hydrodynamicRadius = 1;
    // LB::lbSize
    std::array<unsigned int, 3> lbSize = { 1,1,1 };
    /**
     * standard neighbors shifting
     */
     // LB::ne
    std::array<int, lbmDirec> ne = {};
    // LB::shift
    std::array<unsigned int, 3> shift = { 1,1,1 };
    // LB::domain
    std::array<int, 3> domain = { 1,1,1 };
    /**
     * Initialise dynamic lattice params that scale with the model configuration
     * @see lattice.h for the static lattice params
     */
    void latticeDefinition();
    // conversion units /////////////////////////////////////////////////////////////////////////
    // fluid and granular matter are solved in different measure units
    // for a reference, check Feng, Han, Owen, 2007
    MeasureUnits unit;
};

/**
 * Global scope Lattice Boltzmann parameters structure
 * It is defined (e.g. sans 'extern') inside LBParams.cu
 * CUDA constant's can only be defined in global (or file) scope
 * Therefore, the host replacement is also defined this way
 * @todo Current implementation relies on undefined behaviour or the host copy being accessible (see compiler warning)
 */
#ifdef USE_CUDA
extern __constant__ LBParams PARAMS;
#else
extern LBParams PARAMS;
#endif

#endif  // LBPARAMS_H

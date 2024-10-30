#include "LB2.h"

#ifdef __CUDACC__
#include <cuda.h>
#endif  // __CUDACC__


#ifdef __CUDACC__
__constant__ LB2::Params PARAMS;
#else
LB2::Params PARAMS;
#endif  // __CUDACC__


void LB2::step(const DEM &dem, bool io_demSolver) {

}

void LB2::Params::latticeDefinition() {
    // LATTICE PARAMETERS  ////////////////
    //size-dependent; the others are in the lattice.h

    // domain movement variables
    shift[2] = lbSize[0] * lbSize[1];
    shift[1] = lbSize[0];
    shift[0] = 1;
    //
    domain[2] = lbSize[2] * lbSize[1] * lbSize[0] - 2 * shift[2];
    domain[1] = lbSize[1] * lbSize[0] - 2 * shift[1];
    domain[0] = lbSize[0] - 2 * shift[0];

    // standard neighbors shifting
    // order is O,x,y,z,xy,yz,zx.
    ne[0] = 0;
    //
    ne[1] = shift[0];
    ne[2] = -shift[0];
    //
    ne[3] = shift[1];
    ne[4] = -shift[1];
    //
    ne[5] = shift[2];
    ne[6] = -shift[2];
    //
    ne[7] = shift[0] + shift[1];
    ne[8] = -shift[0] - shift[1];
    ne[9] = -shift[0] + shift[1];
    ne[10] = shift[0] - shift[1];
    //
    ne[11] = shift[1] + shift[2];
    ne[12] = -shift[1] - shift[2];
    ne[13] = -shift[1] + shift[2];
    ne[14] = shift[1] - shift[2];
    //
    ne[15] = shift[2] + shift[0];
    ne[16] = -shift[2] - shift[0];
    ne[17] = -shift[2] + shift[0];
    ne[18] = shift[2] - shift[0];

}
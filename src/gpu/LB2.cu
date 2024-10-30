#include "LB2.h"

#ifdef __CUDACC__
#include <cuda.h>
#endif  // __CUDACC__

#include "DEM.h"

#ifdef __CUDACC__
__constant__ LB2::Params PARAMS;
#else
LB2::Params PARAMS;
#endif  // __CUDACC__

__host_ __device__ tVect Node2::getPosition(const unsigned int& index) const {
    unsigned int x, y, z;

    // index is calculated in this fashion:
    // index = x + y*X + z*X*Y
    // where X and Y are sizes of the lattice in x and y direction

    // from this stems that
    // x + y*X = index MOD X*Y
    // x = x + y*X MOD X
    // y = x + y*X DIV X
    // z = index DIV X*Y

    // see online documentation for class div_t (stdlib.h)
    div_t firstDiv, secondDiv;

    firstDiv = div(int(index), int(PARAMS.lbSize[0] * PARAMS.lbSize[1]));
    secondDiv = div(firstDiv.rem, int(PARAMS.lbSize[0]));

    x = secondDiv.rem;
    y = secondDiv.quot;
    z = firstDiv.quot;

    return tVect(double(x) - 0.5, double(y) - 0.5, double(z) - 0.5);
}

void LB2::step(const DEM &dem, bool io_demSolver) {
    if (io_demSolver) {
        this->latticeBoltzmannCouplingStep(dem.newNeighborList, dem.elmts, dem.particles);
    }

    if (dem.demTime >= dem.demInitialRepeat) {
        this->latticeBolzmannStep(dem.elmts, dem.particles, dem.walls, dem.objects);

        // Lattice Boltzmann core steps
        if (this->freeSurface) {
            this->latticeBoltzmannFreeSurfaceStep();
        }
    }
}
void LB2::latticeBoltzmannCouplingStep(bool &newNeighbourList, const elmtList& elmts, const particleList& particles) {
    // identifies which nodes need to have an update due to particle movement
    // the complexity arises from trying to keep the scaling as close to linear as possible
    // maybe the best idea is to do this in two step:
    // 1) the first is to check for new active nodes and initialise them
    // 2) the second is checking for new solid nodes.
    // this automatically check also for the hideous case of particle to particle double transition

    // first we check if a new neighbour table has been defined. In that case, the indexing needs to be reinitialised
    if (newNeighbourList) {
        cout << endl << "New neighbour list" << endl;
        this->initializeParticleBoundaries();
        newNeighbourList = false;
    }


    // declaring list for mutant nodes
    nodeList newPopUpNodes, newSolidNodes;
    // emptying list of mutant nodes
    newPopUpNodes.clear();
    newSolidNodes.clear();
    // double massSurplus=0.0;

    // SOLID TO ACTIVE CHECK
    findNewActive(newPopUpNodes, elmts, particles);

    //    solidToActive(newPopUpNodes, elmts, massSurplus);

    // ACTIVE TO SOLID CHECK
    findNewSolid(newSolidNodes, elmts, particles);

    if (freeSurface && time % 1 == 0) {
        checkNewInterfaceParticles(elmts, particles);
    }

    //activeToSolid(newSolidNodes, elmts, massSurplus);

    // redistributing extra mass due to popping of nodes
    // redistributeMass(massSurplus);
}

template<>
void LB2::initializeParticleBoundaries<CPU>(const particleList& particles) {
    // Reset all nodes to outside (std::numeric_limits<unsigned short>::max())
    // (The type short is 2 bytes long, but memset requires 4 bytes, so we pass 0xff..)
    memset(d_nodes.solidIndex, 0xffffffff, d_nodes.nodeCount * sizeof(unsigned short));
    
    for (const auto &particle : particles) {
        const tVect convertedPosition = particle.x0 / PARAMS.unit.Length;
        const double convertedRadius = PARAMS.hydrodynamicRadius * particle.r / PARAMS.unit.Length;
        const unsigned int particleIndexHere = particle.particleIndex;
#pragma omp parallel for
        for (size_t i = 0; i < d_nodes.activeNodeCount; ++i) {
            unsigned int nodeIndexHere = d_nodes.coord[i];
            // checking if node is inside a particle
            if (d_nodes.getPosition(i).insideSphere(convertedPosition, convertedRadius)) { //-0.5?
                d_nodes.solidIndex[i] = particleIndexHere;
            }
        }
    }
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
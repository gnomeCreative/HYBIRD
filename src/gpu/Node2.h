#ifndef NODE2_H
#define NODE2_H

#include <cuda.h>

#include <limits>


struct Node2 {
    // The total number of active nodes
    size_t activeNodeCount;

    unsigned int *activeNodesI;
    // The total number of nodes
    size_t nodeCount;
    /**
     * @brief linear coordinate of node
     */
    unsigned int *coord;
    /**
     * @brief The index of the particle the node is inside
     * @note If val = std::numeric_limits<unsigned short>::max() node is outside a particle
     */
    unsigned short *solidIndex;

    /**
     * @brief Return true if the node's soldIndex denotes that it's inside a DEM particle
     */
    __host__ __device__ bool isInsideParticle(unsigned int i) const {
        return solidIndex[i] != std::numeric_limits<unsigned short>::max();
    }

    /**
     * \brief Return the node
     * \param index 
     * \return 
     */
    __host_ __device__ tVect getPosition(const unsigned int& index) const {
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

        firstDiv = div((int)coord[index], (int)(PARAMS.lbSize[0] * PARAMS.lbSize[1]));
        secondDiv = div(firstDiv.rem, (int)PARAMS.lbSize[0]);

        return tVect(
            secondDiv.rem - 0.5,   // x
            secondDiv.quot - 0.5,  // y
            firstDiv.quot - 0.5); // z
    }
};

#endif  // NODE2_H
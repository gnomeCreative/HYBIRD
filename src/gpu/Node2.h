#ifndef NODE2_H
#define NODE2_H

#include <limits>

#include "cuda_helper.h"
#include "myvector.h"

/**
 * @brief A node is a discrete cell within the LBM model's environment.
 *
 * As the LBM environment is sparse, with many inactive cells that contain
 * gas. Nodes are stored compactly, rather than in a dense array, to reduce
 * memory requirements.
 */
struct Node2 {
    // The total number of active nodes
    // A node is active whilst it contains either fluid or interface
    // A node is inactive whilst it contains gas.
    size_t activeCount = 0;
    // Index of nodes marked as active (@todo when is this generated?)
    unsigned int *activeI = nullptr;
    // The total number of interface nodes
    size_t interfaceCount = 0;
    // Index of nodes marked as interface (@todo when is this generated?)
    unsigned int* interfaceI = nullptr;
    // The total number of nodes
    size_t count = 0;
    /**
     * @brief linear coordinate of node
     */
    unsigned int *coord = nullptr;
    /**
     * @brief The index of the particle the node is inside
     * @note If val = std::numeric_limits<unsigned int>::max() node is outside a particle
     */
    unsigned int *solidIndex = nullptr;
    /**
     * @brief neighbour node indexes
     *       Length == count*lbmDirec
     * @note Each node has lbmDirec (19) neighbour nodes, which are stored in a strided pattern
     *       such that indices 0 to N-1, hold every node's first neighbour
     *       indices N to 2N-1, hold every node's second neighbour and so forth
     * @note If nodes are compacted/sorted, as this indexes other nodes it will need regenerating
     */
    unsigned int *d = nullptr;
    /**
     * @brief Identification of node type
     * @see types
     */
    types *type;
    /**
     * @brief particle flag. If true node is inside particle, else false
     */
    bool *p;

    /**
     * @brief Return true if the node's soldIndex denotes that it's inside a DEM particle
     */
    __host__ __device__ bool isInsideParticle(unsigned int i) const {
        return p[i];
    }
    __host__ __device__ void setInsideParticle(unsigned int i, bool t) {
        p[i] = t;
    }
    __host__ __device__ bool isActive(unsigned int i) const {
        return type[i] == LIQUID || type[i] == INTERFACE;
    }

    /**
     * \brief Return the node
     * \param index 
     * \return 
     */
    __host__ __device__ inline tVect getPosition(const unsigned int& index) const;
};

__host__ __device__ inline tVect Node2::getPosition(const unsigned int& index) const {
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

    firstDiv = div(int(coord[index]), int(PARAMS.lbSize[0] * PARAMS.lbSize[1]));
    secondDiv = div(firstDiv.rem, int(PARAMS.lbSize[0]));

    x = secondDiv.rem;
    y = secondDiv.quot;
    z = firstDiv.quot;

    return tVect(double(x) - 0.5, double(y) - 0.5, double(z) - 0.5);
}
#endif  // NODE2_H
#ifndef PARTICLE2_H
#define PARTICLE2_H

/**
 * \brief Elements within the DEM are decomposed into spherical particles.
 *
 * This leads to simpler maths for elements that would contain hard edges.
 * @todo/note original member vars coord and d are fixed at runtime
 */
struct Particle2 {
    // The total number of particles
    // @note This is variable at runtime
    unsigned int count = 0;
    // belonging element index
    unsigned int *clusterIndex = nullptr;
    // particle radius
    double *r = nullptr;
    // position of the particle
    tVect *x0 = nullptr;
    // vector connecting center of element to center of particle
    tVect *radiusVec = nullptr;
};
#endif  // PARTICLE2_H

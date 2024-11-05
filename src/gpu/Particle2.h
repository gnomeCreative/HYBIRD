#ifndef PARTICLE2_H
#define PARTICLE2_H

/**
 * \brief ???
 */
struct Particle2 {
    // The total number of particles
    size_t count = 0;
    // belonging element index
    unsigned int *clusterIndex = nullptr;
    // particle radius
    double *r = nullptr;
    // position of the particle
    tVect *x0 = nullptr;
};
#endif  // PARTICLE2_H

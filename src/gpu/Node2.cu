#include "Node2.h"

__host__ void Node2::initialize(const unsigned int i, const double density, const tVect &velocity, const double massFunction, const double viscosity, const tVect &F, const double ageHere, const tVect &rotationSpeed) {

    //macroscopic properties of the initial condition
    n[i] = density;
    mass[i] = massFunction;
    u[i] = velocity; //FORCE PART  lbF / 2.0 / n+
    visc[i] = viscosity;
    friction[i] = 0.0;
    age[i] = static_cast<float>(ageHere);


    const tVect coriolisAcceleration = computeCoriolis(u[i], rotationSpeed);

    const tVect force = F + centrifugalForce[i] + coriolisAcceleration;

    setEquilibrium(i, density, velocity);
    shiftVelocity(i, Zero);
    reconstruct(i);

    // add force term to new distributions
    addForce(i, force);

}
__host__ void Node2::setEquilibrium(const unsigned int i, const double nHere, const tVect& velHere) {
    const double usq = velHere.norm2();

    for (int j = 0; j < lbmDirec; ++j) {
        constexpr double C1 = 3.0;
        constexpr double C2 = 4.5;
        constexpr double C3 = 1.5;
        // the following lines initialize f to be the local equilibrium values
        const double vu = velHere.dot(v[j]);
        f[i*lbmDirec + j] = fs[i * lbmDirec + j] = coeff[j] * nHere * (1.0 + C1 * vu + C2 * vu * vu - C3 * usq);
    }
}

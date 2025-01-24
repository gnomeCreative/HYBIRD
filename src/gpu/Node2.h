#ifndef NODE2_H
#define NODE2_H

#include <algorithm>
#include <array>
#include <limits>
#include <cub/cub.cuh>

#include "cuda_helper.h"
#include "cub_temp_mem.h"
#include "LBParams.h"
#include "myvector.h"
#include "node.h"

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
    unsigned int activeCount = 0;
    // The allocated size of activeI
    // We don't shrink the buffer if the number of active nodes decreases
    unsigned int activeAlloc = 0;
    // Index of nodes marked as active
    unsigned int *activeI = nullptr;

    // The total number of interface nodes
    unsigned int interfaceCount = 0;
    // Index of nodes marked as interface
    unsigned int* interfaceI = nullptr;

    // The total number of fluid nodes
    unsigned int fluidCount = 0;
    // Index of nodes marked as fluid
    unsigned int *fluidI = nullptr;

    // The total number of wall nodes
    unsigned int wallCount = 0;
    // Index of nodes marked as walls
    unsigned int *wallI = nullptr;

    // The total number of nodes
    unsigned int count = 0;
    /**
     * @brief linear coordinate of node
     */
    // unsigned int *coord = nullptr; // TODO nolonger required
    /**
     * @brief probability density functions (one per lattice direction)
     * @note lbmDirec elements per node, not stored interleaved
     */
    double *f = nullptr;
    /**
     * @brief support variables for streaming
     * @note lbmDirec elements per node, not stored interleaved
     */
    double *fs = nullptr;
    /**
     * @brief density in the node
     */
    double *n = nullptr;
    /**
     * @brief velocity of fluid in node
     */
    tVect *u = nullptr;
    /**
     * @brief DEM-coupling forces on the node
     */
    tVect *hydroForce = nullptr;
    /**
     * @brief centrifugal force on the node
     */
    tVect *centrifugalForce = nullptr;
    // mass functions
    double *mass = nullptr;
    double *newMass = nullptr;
    // viscosity functions
    double *visc = nullptr;
    bool *basal = nullptr;
    // friction (for frictional models)
    double *friction = nullptr;
    // for smooth surface transitions
    float *age = nullptr;
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
     * @note Unallocated/gas nodes are represented by std::numeric_limits<unsigned int>::max()
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

    // functions working on distribution functions
    __host__ void initialize(unsigned int i, double density, const tVect &velocity, double massFunction, double viscosity, const tVect &F, double ageHere, const tVect &rotationSpeed);
    __host__ void setEquilibrium(unsigned int i, double nHere, const tVect& velHere);
    
    __host__ __device__ double liquidFraction(const unsigned int i) const {
        return mass[i] / n[i];
    }
    /**
     * @brief Return true if the node's soldIndex denotes that it's inside a DEM particle
     */
    __host__ __device__ bool isInsideParticle(const unsigned int i) const {
        return p[i];
    }
    __host__ __device__ bool isWall(const unsigned int i) const {
        switch (type[i]) {
        case SLIP_STAT_WALL:
        case SLIP_DYN_WALL:
        case STAT_WALL:
        case DYN_WALL:
        case CYL:
        case OBJ:
        case TOPO:
            return true;
        }
        return false;
    }
    __host__ __device__ void setInsideParticle(const unsigned int i, const bool t) {
        p[i] = t;
    }
    __host__ __device__ bool isActive(const unsigned int i) const {
        return type[i] == LIQUID || type[i] == INTERFACE;
    }
    __host__ __device__ void copy(unsigned int dest_index, unsigned int src_index);
    __host__ __device__ void scatterMass(unsigned int index, double& extraMass);
    __host__ __device__ void shiftVelocity(unsigned int index, const tVect& F);
    __host__ __device__ void computeEquilibrium(unsigned int index, std::array<double, lbmDirec> &feq) const;
    __host__ __device__ void computeEquilibriumTRT(unsigned int index, std::array<double, lbmDirec>& feqp, std::array<double, lbmDirec>& feqm) const;
    __host__ __device__ void computeApparentViscosity(unsigned int index, const std::array<double, lbmDirec> &feq);
    __host__ __device__ void solveCollision(unsigned int index, const std::array<double, lbmDirec> &feq);
    __host__ __device__ void solveCollisionTRT(unsigned int index, const std::array<double, lbmDirec> &feqp, const std::array<double, lbmDirec> &feqm);
    __host__ __device__ void addForce(unsigned int index, const tVect &F);
    __host__ __device__ void addForce(unsigned int index, const std::array<double, lbmDirec>& feq, const tVect& F);
    __host__ __device__ void addForceTRT(unsigned int index, const tVect &F);
    __host__ __device__ tVect bounceBackForce(unsigned int index, unsigned int j, const std::array<double, lbmDirec> &staticPres, double BBi) const;
    __host__ __device__ double massStream(unsigned int index, unsigned int sdir) const;

    /**
     * @brief Copy the force variable (f) to the streaming support variable (fs) for all nodes 
     */
    template<int impl>
    __host__ __device__ void store();
    /**
     * @brief Return the specified node's position within the grid
     * @param index The storage index of the specified node
     */
    __host__ __device__ tVect getPosition(unsigned int index) const;
    /**
     * @brief Returns the grid cell coordinate
     * Replacement for getX(), getY(), getZ()
     */
    __host__ __device__ std::array<int, 3> getGridPosition(unsigned int index) const;
    /**
     * reconstruction of macroscopic variables from microscopic distribution
     * this step is necessary to proceed to the collision step
     */
    __host__ __device__ void reconstruct(unsigned int index);
    /**
     * collision operator
     */
    __host__ __device__ void collision(unsigned int index);

    /**
     * Regenerate activeI by combining and sorting fluidI and interfaceI
     */
    template<int impl>
    void cleanLists();

    __host__ __device__ std::array<unsigned int, lbmDirec> findNeighbors(unsigned int index) const;
    __host__ __device__ void generateNode(unsigned int index, types typeHere);
    __host__ __device__ void eraseNode(unsigned int index);
};

__host__ __device__ __forceinline__ tVect Node2::getPosition(const unsigned int index) const {
    unsigned int x, y, z;

    // index is calculated in this fashion:
    // index = x + y*X + z*X*Y
    // where X and Y are sizes of the lattice in x and y direction

    // from this stems that
    // x + y*X = index MOD X*Y
    // x = x + y*X MOD X
    // y = x + y*X DIV X
    // z = index DIV X*Y

#ifdef __CUDA_ARCH__
    // div() does not exist in device code
    z = index / (PARAMS.lbSize[0] * PARAMS.lbSize[1]);
    const int t = index % (PARAMS.lbSize[0] * PARAMS.lbSize[1]);
    y = t / PARAMS.lbSize[0];
    x = t % PARAMS.lbSize[0];
#else
    // see online documentation for class div_t (stdlib.h)
    div_t firstDiv, secondDiv;

    firstDiv = div(int(index), int(PARAMS.lbSize[0] * PARAMS.lbSize[1]));
    secondDiv = div(firstDiv.rem, int(PARAMS.lbSize[0]));

    x = secondDiv.rem;
    y = secondDiv.quot;
    z = firstDiv.quot;
#endif

    return { double(x) - 0.5, double(y) - 0.5, double(z) - 0.5 };
}
__host__ __device__ __forceinline__ std::array<int, 3> Node2::getGridPosition(const unsigned int index) const {
    unsigned int x, y, z;

    // index is calculated in this fashion:
    // index = x + y*X + z*X*Y
    // where X and Y are sizes of the lattice in x and y direction

    // from this stems that
    // x + y*X = index MOD X*Y
    // x = x + y*X MOD X
    // y = x + y*X DIV X
    // z = index DIV X*Y
#ifdef __CUDA_ARCH__
    // div() does not exist in device code
    z = index / (PARAMS.lbSize[0] * PARAMS.lbSize[1]);
    const int t = index % (PARAMS.lbSize[0] * PARAMS.lbSize[1]);
    y = t / PARAMS.lbSize[0];
    x = t % PARAMS.lbSize[0];
#else
    // see online documentation for class div_t (stdlib.h)
    div_t firstDiv, secondDiv;

    firstDiv = div(int(index), int(PARAMS.lbSize[0] * PARAMS.lbSize[1]));
    secondDiv = div(firstDiv.rem, int(PARAMS.lbSize[0]));

    x = secondDiv.rem;
    y = secondDiv.quot;
    z = firstDiv.quot;
#endif

    return { static_cast<int>(x), static_cast<int>(y), static_cast<int>(z) };
}
__host__ __device__ __forceinline__ void Node2::copy(unsigned int dest_index, unsigned int src_index) {
    //macroscopic properties of the initial condition
    this->n[dest_index] = this->n[src_index];
    this->mass[dest_index] = 0.1 * this->mass[src_index];
    this->u[dest_index] = this->u[src_index];
    this->visc[dest_index] = this->visc[src_index];
    this->friction[dest_index] = this->friction[src_index];
    this->centrifugalForce[dest_index] = this->centrifugalForce[src_index];

    //for (int j = 0; j < lbmDirec; ++j) {
    //    this->f[dest_index * lbmDirec + j] = this->f[src_index * lbmDirec + j];
    //    this->fs[dest_index * lbmDirec + j] = this->fs[src_index * lbmDirec + j];
    //}
    memcpy(this->f + (dest_index * lbmDirec), this->f + (src_index * lbmDirec), sizeof(double) * lbmDirec);
    memcpy(this->fs + (dest_index * lbmDirec), this->fs + (src_index * lbmDirec), sizeof(double) * lbmDirec);
}
__host__ __device__ __forceinline__ void Node2::scatterMass(unsigned int index, double& extraMass) {

    unsigned int totNodes = 0;
    for (int j = 1; j < lbmDirec; ++j) {
        unsigned int ln_i = d[j];
        if (ln_i != std::numeric_limits<unsigned int>::max()) {
            if (this->type[ln_i] == INTERFACE) {
                ++totNodes;
            }
        }
    }
    if (totNodes) {
        for (int j = 1; j < lbmDirec; ++j) {
            unsigned int ln_i = d[j];
            if (ln_i != std::numeric_limits<unsigned int>::max()) {
                if (this->type[ln_i] == INTERFACE) {
                    this->extraMass[ln_i] += extraMass / double(totNodes);
                }
            }
        }
        extraMass = 0.0;
    }
    // otherwise mass is not redistributed, but rather given to the global redistribution function
}
__host__ __device__ __forceinline__ void Node2::shiftVelocity(const unsigned int index, const tVect& F) {
    const tVect totalForce = F + this->hydroForce[index];
    this->u[index] += 0.5 * this->mass[index] * totalForce;
}
__host__ __device__ __forceinline__ void Node2::computeEquilibrium(const unsigned int index, std::array<double, lbmDirec> &feq) const {
    constexpr double C1 = 3.0;
    constexpr double C2 = 4.5;
    constexpr double C3 = 1.5;

    const double usq = this->u[index].norm2();

    for (int j = 0; j < lbmDirec; ++j) {
        const double vu = this->u[index].dot(v[j]);
        feq[j] = coeff[j] * this->n[index] * (1.0 + C1 * vu + C2 * vu * vu - C3 * usq);
    }    
}
__host__ __device__ __forceinline__ void Node2::computeEquilibriumTRT(const unsigned int index, std::array<double, lbmDirec> &feqp, std::array<double, lbmDirec> &feqm) const {
    constexpr double C1 = 3.0;
    constexpr double C2 = 4.5;
    constexpr double C3 = 1.5;

    const double usq = this->u[index].norm2();
    double posSum = 0.0;
    for (int j = 1; j < lbmDirec; ++j) {
        const double vu = this->u[index].dot(v[j]);
        feqp[j] = coeff[j] * this->n[index] * (1.0 + C2 * vu * vu - C3 * usq);
        feqm[j] = coeff[j] * this->n[index] * (C1 * vu);
        posSum += feqp[j];
    }
    feqp[0] = this->n[index] - posSum;
    feqm[0] = 0.0;
}
__host__ __device__ __forceinline__ void Node2::computeApparentViscosity(const unsigned int index, const std::array<double, lbmDirec> &feq) {
    // minimum and maximum viscosity
    const double tau = 0.5 + 3.0 * this->visc[index];

    tMat gamma(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    for (int j = 1; j < lbmDirec; ++j) {
        gamma += (this->f[index * lbmDirec + j] - feq[j]) * vv[j];
    }
    gamma *= -1.5 / (tau * this->n[index]);

    // shear rate (second invariant)
    const double shearRate = 2.0 * gamma.magnitude();

    // Bingham model
    double nuApp = 0.0;
    switch (PARAMS.fluidMaterial.rheologyModel) {
        case NEWTONIAN:
        {
            nuApp = PARAMS.fluidMaterial.initDynVisc;
            break;
        }
        case BINGHAM:
        {
            nuApp = PARAMS.fluidMaterial.plasticVisc + PARAMS.fluidMaterial.yieldStress / shearRate;
            break;
        }
        case FRICTIONAL:
        {
            // p=c_s^2 * n, scaled by atmospheric pressure (n_atm=1.0)
            const double pressure = std::max(PARAMS.fluidMaterial.minimumPressure, 0.33333333 * (this->n[index] - 1.0));
            if (this->basal[index]) {
                this->friction[index] = PARAMS.fluidMaterial.basalFrictionCoefFluid;
            } else {
                this->friction[index] = PARAMS.fluidMaterial.frictionCoefFluid;
            }
            nuApp = PARAMS.fluidMaterial.initDynVisc + this->friction[index] * pressure / shearRate;
            break;
        }
        case VOELLMY:
        {
            // p=c_s^2 * n, scaled by atmospheric pressure (n_atm=1.0)
            const double pressure = std::max(PARAMS.fluidMaterial.minimumPressure, 0.33333333 * (this->n[index] - 1.0));
            if (this->basal[index]) {
                this->friction[index] = PARAMS.fluidMaterial.basalFrictionCoefFluid;
            } else {
                this->friction[index] = PARAMS.fluidMaterial.frictionCoefFluid;
            }
            nuApp = this->friction[index] * pressure / shearRate + PARAMS.fluidMaterial.rhod2 * shearRate;
            break;
        }
        case BAGNOLD:
        {
            // p=c_s^2 * n, scaled by atmospheric pressure (n_atm=1.0)
            nuApp = PARAMS.fluidMaterial.rhod2 * shearRate;
            break;
        }
        case MUI:
        {
            // p=c_s^2 * n, scaled by atmospheric pressure (n_atm=1.0)
            const double pressure = std::max(PARAMS.fluidMaterial.minimumPressure, 0.33333333 * (this->n[index] - 1.0)); //smoothedPressure
            const double inertialNumber = PARAMS.fluidMaterial.particleDiameter * shearRate / sqrt(pressure / PARAMS.fluidMaterial.particleDensity);
            const double regularizationFactor = 1.0;//-exp(-shearRate/0.00005);
            if (this->basal[index]) {
                this->friction[index] = PARAMS.fluidMaterial.basalFrictionCoefFluid * regularizationFactor + PARAMS.fluidMaterial.deltaFriction / (PARAMS.fluidMaterial.baseInertial / inertialNumber + 1.0);
            } else {
                this->friction[index] = PARAMS.fluidMaterial.frictionCoefFluid * regularizationFactor + PARAMS.fluidMaterial.deltaFriction / (PARAMS.fluidMaterial.baseInertial / inertialNumber + 1.0);
            }
            nuApp = this->friction[index] * pressure / shearRate;
            if (this->type[index] == INTERFACE) { //(pressure<1.5*fluidMaterial.minimumPressure) {
                nuApp = PARAMS.fluidMaterial.lbMinVisc;
            }
            break;
        }
    }

    // Smagorinsky turbulence model
    if (PARAMS.fluidMaterial.turbulenceOn) {
        const double nuTurb = PARAMS.fluidMaterial.turbConst * shearRate;
        nuApp += nuTurb;
    }

    // limiting for stability
    this->visc[index] = std::max(PARAMS.fluidMaterial.lbMinVisc, std::min(PARAMS.fluidMaterial.lbMaxVisc, nuApp));
}
__host__ __device__ __forceinline__ void Node2::solveCollision(const unsigned int index, const std::array<double, lbmDirec> &feq) {
    // relaxation frequency
    const double omega = 1.0 / (0.5 + 3.0 * this->visc[index]);

    // temporary pointer to start of this node's array
    const auto _f = this->f + (index * lbmDirec);  // @todo Is _f the most suitable name?
    for (unsigned int j = 0; j < lbmDirec; ++j) {
        _f[j] += omega * (feq[j] - _f[j]);
    }
}
__host__ __device__ __forceinline__ void Node2::solveCollisionTRT(const unsigned int index, const std::array<double, lbmDirec> &feqp, const std::array<double, lbmDirec> &feqm) {
    // relaxation frequency
    const double omegap = 1.0 / (0.5 + 3.0 * this->visc[index]);
    const double omegam = 6.0 * this->visc[index] / (2.0 * PARAMS.magicNumber + 3.0 * this->visc[index]);

    // symmetric probability density functions (one per lattice direction)
    std::array<double, lbmDirec> fp;
    // antisymmetric probability density functions (one per lattice direction)
    std::array<double, lbmDirec> fm;

    // temporary pointer to start of this node's array
    const auto _f = this->f + (index * lbmDirec);  // @todo Is _f the most suitable name?
    fp[0] = _f[0];
    fm[0] = 0;
    for (unsigned int j = 1; j < lbmDirec; ++j) {
        // symmetric part
        fp[j] = 0.5 * (_f[j] + _f[opp[j]]);
        // antisymmetric part
        fm[j] = 0.5 * (_f[j] - _f[opp[j]]);
    }

    for (unsigned int j = 0; j < lbmDirec; ++j) {
        _f[j] += omegam * (feqm[j] - fm[j]) + omegap * (feqp[j] - fp[j]);
    }
}
__host__ __device__ __forceinline__ void Node2::addForce(const unsigned int index, const tVect &F) {
    constexpr double F1 = 3.0;
    constexpr double F2 = 9.0;

    const double omegaf = 1.0 - 1.0 / (1.0 + 6.0 * this->visc[index]);
    const tVect totalForce = F + this->hydroForce[index];

    // temporary pointer to start of this node's array
    const auto _f = this->f + (index * lbmDirec);  // @todo Is _f the most suitable name?
    for (int j = 0; j < lbmDirec; j++) {
        const double vu = this->u[index].dot(v[j]);
        const tVect vmu = v[j] - this->u[index];
        const tVect forcePart = F2 * vu * v[j] + F1*vmu;
        _f[j] +=  omegaf * coeff[j] * forcePart.dot(this->mass[index] *totalForce);
    }
}
__host__ __device__ __forceinline__ void Node2::addForce(const unsigned int index, const std::array<double, lbmDirec> &feq, const tVect &F) {
    constexpr double F1 = 3.0;

    const double omegaf = 1.0 - 1.0 / (1.0 + 6.0 * this->visc[index]);
    const tVect totalForce = F + this->hydroForce[index];

    // temporary pointer to start of this node's array
    const auto _f = this->f + (index * lbmDirec);  // @todo Is _f the most suitable name?
    for (int j = 0; j < lbmDirec; j++) {
        const tVect vmu = v[j] - this->u[index];
        const tVect forcePart = feq[j]*F1*vmu / this->n[index];
        _f[j] +=  omegaf * forcePart.dot(this->mass[index]*totalForce);
    }
}
__host__ __device__ __forceinline__ void Node2::addForceTRT(const unsigned int index, const tVect& F) {
    constexpr double F1 = 3.0;
    constexpr double F2 = 9.0;

    constexpr double Lambda = 0.25;//9*visc*visc;
    const double omegam = 6.0 * this->visc[index] / (2.0 * Lambda + 3.0 * this->visc[index]);
    const double omegaf = 1.0 - 0.5 * omegam;

    const tVect totalForce = F + this->hydroForce[index];

    // temporary pointer to start of this node's array
    const auto _f = this->f + (index * lbmDirec);  // @todo Is _f the most suitable name?
    for (int j = 0; j < lbmDirec; j++) {
        const double vu = this->u[index].dot(v[j]);
        const tVect vmu = v[j] - this->u[index];
        const tVect forcePart = F2 * vu * v[j] + F1 * vmu;
        _f[j] += omegaf * coeff[j] * forcePart.dot(totalForce);
    }
}
__host__ __device__ __forceinline__ tVect Node2::bounceBackForce(const unsigned int index, const unsigned int j, const std::array<double, lbmDirec>& staticPres, const double BBi) const {
    if (visc[index] == PARAMS.fluidMaterial.lbMaxVisc) {
        return (2.0 * (fs[index * lbmDirec + j] + coeff[j] * (n[index] - 1.0) * (PARAMS.fluidMaterial.earthPressureCoeff - 1.0) - staticPres[j]) - BBi) * v[j];
    } else {
        return (2.0 * (fs[index * lbmDirec + j] - staticPres[j]) - BBi) * v[j];
    }
}

__host__ __device__ __forceinline__ double Node2::massStream(const unsigned int index, const unsigned int sdir) const {
    return (f[index * lbmDirec + opp[sdir]] - fs[index * lbmDirec + sdir]);
}
template <>
inline void Node2::store<CPU>() {
    // Saving in streaming support variables f->fs
    memcpy(fs, f, sizeof(double) * lbmDirec * count);
}
#ifdef USE_CUDA
template <>
inline void Node2::store<CUDA>() {
#ifdef __CUDA_ARCH__
    // memcpy(fs, f, sizeof(double) * lbmDirec * count);  // This won't work, would need to be lbmDirec per thread
#else
    // Saving in streaming support variables f->fs
    cudaMemcpy(fs, f, sizeof(double) * lbmDirec * count, cudaMemcpyDeviceToDevice);
#endif
}
#endif
__host__ __device__ __forceinline__ void Node2::reconstruct(const unsigned int index) {
    // reconstruction of macroscopical physical variables
    this->n[index] = 0.0;
    this->u[index].reset();
    for (unsigned int j = 0; j < lbmDirec; ++j) {
        const double _f = this->f[index * lbmDirec + j];
        // density
        this->n[index] += _f;
        // momentum
        this->u[index] += _f * v[j];
    }

    // velocity
    this->u[index] /= this->n[index];
}
__host__ __device__ __forceinline__ void Node2::collision(const unsigned int index) {
    if (!PARAMS.TRTsolver) {
        // equilibrium distributions
        std::array<double, lbmDirec> feq;
        // force field
        tVect force = PARAMS.lbF;
        if (PARAMS.solveCentrifugal) {
            force += this->centrifugalForce[index];
        }
        if (PARAMS.solveCoriolis) {
            force += computeCoriolis(this->u[index], PARAMS.rotationSpeed);
        }
        // shift velocity field to F/2
        this->shiftVelocity(index, force);
        this->computeEquilibrium(index, feq);
        if (PARAMS.fluidMaterial.rheologyModel != NEWTONIAN || PARAMS.fluidMaterial.turbulenceOn) {
            // compute shear rate tensor, find invariant calculate viscosity (Bingham)
            this->computeApparentViscosity(index, feq);
        }
        // compute new distributions
        this->solveCollision(index, feq);
        // add force term to new distributions
        this->addForce(index, feq, force);
    } else {
        // equilibrium distributions
        std::array<double, lbmDirec> feqp;
        std::array<double, lbmDirec> feqm;
        // force field
        tVect force = PARAMS.lbF;
        if (PARAMS.solveCentrifugal) {
            force += this->centrifugalForce[index];
        }
        if (PARAMS.solveCoriolis) {
            force += computeCoriolis(this->u[index], PARAMS.rotationSpeed);
        }
        // shift velocity field to F/2
        this->shiftVelocity(index, force);
        this->computeEquilibriumTRT(index, feqp, feqm);
        if (PARAMS.fluidMaterial.rheologyModel != NEWTONIAN || PARAMS.fluidMaterial.turbulenceOn) {
            // compute shear rate tensor, find invariant calculate viscosity (Bingham)
            this->computeApparentViscosity(index, feqp);
        }
        else this->visc[index] = PARAMS.fluidMaterial.initDynVisc;
        // compute new distributions
        this->solveCollisionTRT(index, feqp, feqm);
        // add force term to new distributions
        this->addForceTRT(index, force);
    }
}

/**
 * cleanLists()
 */
template<>
inline void Node2::cleanLists<CPU>() {
    // If the buffer needs to be grown
    if (interfaceCount + fluidCount > activeAlloc) {
        activeAlloc = interfaceCount + fluidCount;
        activeCount = activeAlloc;
        if (activeI) {
            free(activeI);
        }
        activeI = static_cast<unsigned int*>(malloc(sizeof(unsigned int) * activeAlloc));
    } else {
        activeCount = interfaceCount + fluidCount;
    }
    // Copy interface IDs
    memcpy(activeI, interfaceI, sizeof(unsigned int) * interfaceCount);
    // Copy fluid IDs (there shouldn't be duplicates)
    memcpy(activeI + interfaceCount, fluidI, sizeof(unsigned int) * fluidCount);
    // Sort (@todo it is assumed, not known, that nodes are in coord order)
    std::sort(activeI, activeI + activeCount);
}
#ifdef USE_CUDA
template<>
inline void Node2::cleanLists<CUDA>() {
    // If the buffer needs to be grown
    if (interfaceCount + fluidCount > activeAlloc) {
        activeAlloc = interfaceCount + fluidCount;
        activeCount = activeAlloc;
        if (activeI) {
            CUDA_CALL(cudaFree(activeI));
        }
        CUDA_CALL(cudaMalloc(&activeI, sizeof(unsigned int) * activeAlloc));
    } else {
        activeCount = interfaceCount + fluidCount;
    }
    // Copy interface IDs
    CUDA_CALL(cudaMemcpy(activeI, interfaceI, sizeof(unsigned int)* interfaceCount, cudaMemcpyDeviceToDevice));
    // Copy fluid IDs (there shouldn't be duplicates)
    CUDA_CALL(cudaMemcpy(activeI + interfaceCount, fluidI, sizeof(unsigned int)* fluidCount, cudaMemcpyDeviceToDevice));
    /**
     * Sort (@todo it is assumed, not known, that nodes are in coord order)
     */
    // Prepare the sort, calculate how much and allocate temporary memory for cub
    const int max_bit = static_cast<int>(floor(log2(count))) + 1;
    size_t temp_storage_bytes = 0;
    auto& ctb = CubTempMem::GetBufferSingleton();
    ctb.resize(activeCount * sizeof(unsigned int));
    cub::DeviceRadixSort::SortKeys(nullptr, temp_storage_bytes, activeI, static_cast<unsigned int *>(ctb.getPtr()), activeCount, 0, max_bit);
    auto &ctm = CubTempMem::GetTempSingleton();
    ctm.resize(temp_storage_bytes);
    // Actually perform the sort
    cub::DeviceRadixSort::SortKeys(ctm.getPtr(), ctm.getSize(), activeI, static_cast<unsigned int*>(ctb.getPtr()), activeCount, 0, max_bit);
    // Swap the input and output buffers
    activeAlloc = static_cast<unsigned int>(ctb.swapPtr(activeI, activeAlloc * sizeof(unsigned int)) / sizeof(unsigned int));
}
#endif


__host__ __device__ __forceinline__ std::array<unsigned int, lbmDirec> Node2::findNeighbors(unsigned int index) const {
    std::array<unsigned int, lbmDirec> neighborCoord;
    // assign boundary characteristic to nodes (see class)
    // if not differently defined, type is 0 (fluid)

    for (int j = 1; j < lbmDirec; ++j) {
        neighborCoord[j] = index + PARAMS.ne[j];
    }

    // BOUNDARY CONDITIONS ///////////////////////////
    // nodes on the boundary have no neighbors
    if (this->isWall(index)) {
        const std::array<int, 3> pos = this->getGridPosition(index);// getPosition() adds + 0.5x
        if (pos[0] == 0) {
            for (unsigned int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Xp) < 0.0) {
                    neighborCoord[j] = index;
                }
            }
        }
        if (pos[0] == PARAMS.lbSize[0] - 1) {
            for (unsigned int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Xp) > 0.0) {
                    neighborCoord[j] = index;
                }
            }
        }
        if (pos[1] == 0) {
            for (unsigned int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Yp) < 0.0) {
                    neighborCoord[j] = index;
                }
            }
        }
        if (pos[1] == PARAMS.lbSize[1] - 1) {
            for (unsigned int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Yp) > 0.0) {
                    neighborCoord[j] = index;
                }
            }
        }
        if (pos[2] == 0) {
            for (unsigned int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Zp) < 0.0) {
                    neighborCoord[j] = index;
                }
            }
        }
        if (pos[2] == PARAMS.lbSize[2] - 1) {
            for (unsigned int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Zp) > 0.0) {
                    neighborCoord[j] = index;
                }
            }
        }
    }// PERIODICITY ////////////////////////////////////////////////
        // assigning periodicity conditions (this needs to be done after applying boundary conditions)
        // runs through free cells and identifies neighboring cells. If neighbor cell is
        // a special cell (periodic) then the proper neighboring condition is applied
        // calculates the effect of periodicity
    else if (this->isActive(index)) {
        // neighboring and periodicity vector for boundary update
        std::array<unsigned int, lbmDirec> pbc = {};
        if (this->type[neighborCoord[1]] == PERIODIC) {
            for (unsigned int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Xp) > 0.0) {
                    pbc[j] -= PARAMS.domain[0];
                }
            }
        }
        if (this->type[neighborCoord[2]] == PERIODIC) {
            for (unsigned int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Xp) < 0.0) {
                    pbc[j] += PARAMS.domain[0];
                }
            }
        }
        if (this->type[neighborCoord[3]] == PERIODIC) {
            for (unsigned int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Yp) > 0.0) {
                    pbc[j] -= PARAMS.domain[1];
                }
            }
        }
        if (this->type[neighborCoord[4]] == PERIODIC) {
            for (unsigned int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Yp) < 0.0) {
                    pbc[j] += PARAMS.domain[1];
                }
            }
        }
        if (this->type[neighborCoord[5]] == PERIODIC) {
            for (unsigned int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Zp) > 0.0) {
                    pbc[j] -= PARAMS.domain[2];
                }
            }
        }
        if (this->type[neighborCoord[6]] == PERIODIC) {
            for (unsigned int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Zp) < 0.0) {
                    pbc[j] += PARAMS.domain[2];
                }
            }
        }

        // apply periodicity
        for (unsigned int j = 1; j < lbmDirec; ++j) {
            neighborCoord[j] += pbc[j];
        }

    }

    return neighborCoord;
}
/**
 * CUDA capable version of generateNode()
 */
__host__ __device__ __forceinline__ void Node2::generateNode(unsigned int index, types typeHere) {
    // set type
    this->type[index] = typeHere;
    this->p[index] = false;  // setOutsideParticle()
    this->age[index] = 0.0;

    // find neighbor indices
    const std::array<unsigned int, lbmDirec> neighborindex = this->findNeighbors(index);

    this->basal[index] = false;

    // set centrifugal acceleration
    this->centrifugalForce[index] = computeCentrifugal(this->getPosition(index), PARAMS.rotationCenter, PARAMS.rotationSpeed);

    // assign neighbor nodes
    for (unsigned int j = 1; j < lbmDirec; ++j) {
        // linearized indexinate of neighbor nodes
        const unsigned int link = neighborindex[j];
        // check if node at that location exists
        if (link < this->count && this->type[link] != GAS) {
            // assign neighbor for local node
            this->d[j * this->count + index] = link;
            // if neighbor node is also active, link it to local node
            if (this->isActive(index)) {
                /// TODO potential race condition
                this->d[opp[j] * this->count + link] = index;
                if (this->isWall(link)) {
                    this->basal[index] = true;
                }
            }
        } else {
            this->d[j * this->count + index] = std::numeric_limits<unsigned int>::max();
        }
    }
}
__host__ __device__ __forceinline__ void Node2::eraseNode(const unsigned int index) {
    // find neighbor indices
    std::array<unsigned int, lbmDirec> neighborCoord = findNeighbors(index);

    //remove from map;
    this->type[index] = GAS;

    // assign neighbor nodes
    for (int j = 1; j < lbmDirec; ++j) {
        // linearized coordinate of neighbor nodes
        const unsigned int ln_i = neighborCoord[j];
        // check if node at that location exists
        if (this->type[ln_i] != std::numeric_limits<unsigned int>::max()) {
            // if neighbor node is active, remove link to local node
            if (this->isActive(ln_i)) {
                this->d[this->count * opp[j] + ln_i] = std::numeric_limits<unsigned int>::max();
            }
        }
    }
}
#endif  // NODE2_H
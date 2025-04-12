#ifndef LBCOMMON_H
#define LBCOMMON_H

/**
 * This header contains common implementation usable by both CPU/OPENMP and CUDA implementations
 * All functions must be marked __hdi__, which is a macro that should avoid linker issues
 */

#include "Node2.h"
#include "Particle2.h"

/**
 * __hdi__ should map to a different token dependent on whether it's seen by C/C++ or CUDA compiler
 * C/C++: Just mark the method inline to avoid linker errors
 * CUDA: Force inline and compile for both host and device code
 */
#ifdef __CUDACC__
#define __hdi__ __host__ __device__ __forceinline__
#else
#define __hdi__ inline
#endif

//
// latticeBoltzmannCouplingStep() subroutines
//
__hdi__ double common_initializeParticleBoundaries(const unsigned int i, Node2* nodes, Particle2* particles) {
    // Fetch the index of the (active) node being processed
    const unsigned int an_i = nodes->activeI[i];
    const tVect node_position = nodes->getPosition(an_i);
    for (unsigned int p_i = 0; p_i < particles->count; ++p_i) {
        const tVect convertedPosition = particles->x0[p_i] / PARAMS.unit.Length;
        // @todo pre-compute PARAMS.hydrodynamicRadius / PARAMS.unit.Length ?
        const double convertedRadius = particles->r[p_i] * PARAMS.hydrodynamicRadius / PARAMS.unit.Length;
        if (node_position.insideSphere(convertedPosition, convertedRadius)) { //-0.5?
            nodes->setInsideParticle(an_i, true);
            nodes->solidIndex[an_i] = p_i;
            return nodes->mass[an_i];  // @todo in original code it doesn't break after setting
        }
    }
    return 0.0;
}
__hdi__ void common_findNewActive(const unsigned int i, Node2* nodes, Particle2* particles, Element2* elements) {
    // Fetch the index of the (active) node being processed
    const unsigned int an_i = nodes->activeI[i];
    if (nodes->p[an_i]) {
        const tVect nodePosition = nodes->getPosition(an_i);
        // solid index to identify cluster
        const unsigned int particleIndex = nodes->solidIndex[an_i];
        const unsigned int clusterIndex = particles->clusterIndex[particleIndex];
        // in this case check if it has been uncovered (must be out of all particles of the cluster) - we start with a true hypothesis
        // cycling through component particles
        const unsigned int first_component = elements->componentsIndex[clusterIndex];
        const unsigned int last_component = elements->componentsIndex[clusterIndex + 1];
        for (unsigned int j = first_component; j < last_component; ++j) {
            // getting indexes from particle composing the cluster
            const unsigned int componentIndex = elements->componentsData[j];
            // checking if it has been uncovered in component j of the cluster
            // radius need to be increased by half a lattice unit
            // this is because solid boundaries are located halfway between solid and fluid nodes
            const tVect convertedPosition = particles->x0[componentIndex] / PARAMS.unit.Length;
            // @todo pre-compute PARAMS.hydrodynamicRadius / PARAMS.unit.Length ?
            const double convertedRadius = particles->r[componentIndex] * PARAMS.hydrodynamicRadius / PARAMS.unit.Length;
            if (nodePosition.insideSphere(convertedPosition, convertedRadius)) { //-0.5?
                // if the node is still inside the element, the hypothesis of new active is not true anymore
                // and we can get out of the cycle
                return;
            }
        }
        // turning up the cell as we didn't exit early
        nodes->setInsideParticle(an_i, false);
    }
}
__hdi__ void common_findNewSolid(const unsigned int i, Node2* nodes, Particle2* particles, Element2* elements) {
    const unsigned int an_i = nodes->activeI[i];
    if (nodes->isInsideParticle(an_i)) {  // If node is inside particle
        // solid index to identify cluster
        const unsigned int particleIndex = nodes->solidIndex[an_i];
        const unsigned int clusterIndex = particles->clusterIndex[particleIndex];
        // cycle through first neighbors
        const unsigned int nodeCount = nodes->count;
        for (int k = 1; k < lbmMainDirec; ++k) {
            const unsigned int l_i = nodes->d[nodeCount * k + an_i];
            if (l_i != std::numeric_limits<unsigned int>::max()) {
                // checking if solid particle is close to an active one -> we have an active node to check
                if (!nodes->isInsideParticle(l_i) && nodes->isActive(l_i)) {
                    const tVect linkPosition = nodes->getPosition(l_i);
                    // check if neighbors has been covered (by any of the particles of the cluster) - we start with a false hypothesis
                    // cycling through all components of the cluster
                    const unsigned int first_component = elements->componentsIndex[clusterIndex];
                    const unsigned int last_component = elements->componentsIndex[clusterIndex + 1];
                    for (unsigned int j = first_component; j < last_component; ++j) {
                        // getting component particle index
                        const unsigned int componentIndex = elements->componentsData[j];
                        // check if it getting inside
                        // radius need to be increased by half a lattice unit
                        // this is because solid boundaries are located halfway between solid and fluid nodes
                        // @todo pre-compute PARAMS.hydrodynamicRadius / PARAMS.unit.Length ?
                        if (linkPosition.insideSphere(particles->x0[componentIndex] / PARAMS.unit.Length, particles->r[componentIndex] * PARAMS.hydrodynamicRadius / PARAMS.unit.Length)) { //-0.5?
                            // if so, then the false hypothesis does not hold true anymore
                            nodes->solidIndex[l_i] = componentIndex;
                            // By setting particle to inside, it won't be checked again, newSolidNodes hence becomes redundant
                            nodes->setInsideParticle(l_i, true);  // @todo Is this a race condition? Multiple nodes may share a link node?
                            // and we exit the cycle
                            break;
                        }
                    }
                }
            }
        }
    }
}
__hdi__ void common_checkNewInterfaceParticles(const unsigned int e_i, Node2* nodes, Particle2* particles, Element2* elements) {
    // INITIAL PARTICLE POSITION ////////////////////////
    if (elements->FHydro[e_i].norm2() == 0.0) {
        const unsigned int first_component = elements->componentsIndex[e_i];
        const unsigned int last_component = elements->componentsIndex[e_i + 1];
        for (unsigned int n = first_component; n < last_component; ++n) {
            const unsigned int componentIndex = elements->componentsData[n];
            const tVect convertedPosition = particles->x0[componentIndex] / PARAMS.unit.Length;
            // @todo pre-compute PARAMS.hydrodynamicRadius / PARAMS.unit.Length ?
            const double convertedRadius = particles->r[componentIndex] * PARAMS.hydrodynamicRadius / PARAMS.unit.Length;
            for (unsigned int i_i = 0; i_i < nodes->interfaceCount; ++i_i) {
                const unsigned int nodeHere = nodes->interfaceI[i_i];
                if (!nodes->isInsideParticle(nodeHere)) {
                    // checking if node is inside a particle
                    const tVect nodePosition = nodes->getPosition(nodeHere);
                    if (nodePosition.insideSphere(convertedPosition, convertedRadius)) { //-0.5?
                        nodes->setInsideParticle(nodeHere, true);
                        nodes->solidIndex[nodeHere] = componentIndex;
                    }
                }
            }
        }
    }
}

//
// latticeBoltzmannStep() subroutines
//
__hdi__ void common_computeHydroForces(const unsigned int an_i, Node2* nodes, Particle2* particles, Element2* elements) {
    // resetting hydrodynamic forces on nodes
    nodes->hydroForce[an_i].reset();
    if (nodes->isInsideParticle(an_i)) {
        // getting the index of the particle to compute force in the right object
        const unsigned int index = an_i;
        const unsigned int particleIndex = nodes->solidIndex[an_i];
        const unsigned int clusterIndex = particles->clusterIndex[particleIndex];
        // calculating velocity of the solid boundary at the node (due to rotation of particles)
        // vectorized radius (real units)
        const tVect radius = nodes->getPosition(index) - particles->x0[particleIndex] / PARAMS.unit.Length + particles->radiusVec[particleIndex] / PARAMS.unit.Length;
        // update velocity of the particle node (u=v_center+omega x radius) (real units)
        const tVect localVel = elements->x1[clusterIndex] / PARAMS.unit.Speed + (elements->wGlobal[clusterIndex].cross(radius)) / PARAMS.unit.AngVel;

        // calculate differential velocity
        const tVect diffVel = nodes->age[an_i] * nodes->age[an_i] * nodes->liquidFraction(an_i) * (nodes->u[an_i] - localVel);

        // force on fluid
        nodes->hydroForce[an_i] += -1.0 * diffVel;

        // force on particle
#ifdef __CUDA_ARCH__
        // CUDA atomics
        atomicAdd(&elements->fluidVolume[clusterIndex], nodes->mass[an_i]);
        atomicAdd(&elements->FHydro[clusterIndex].x, 1.0 * diffVel.x);
        atomicAdd(&elements->FHydro[clusterIndex].y, 1.0 * diffVel.y);
        atomicAdd(&elements->FHydro[clusterIndex].z, 1.0 * diffVel.z);
        const tVect t = 1.0 * radius.cross(diffVel);
        atomicAdd(&elements->MHydro[clusterIndex].x, t.x);
        atomicAdd(&elements->MHydro[clusterIndex].y, t.y);
        atomicAdd(&elements->MHydro[clusterIndex].z, t.z);
#else
        // CPU atomics
#pragma omp atomic update
        elements->fluidVolume[clusterIndex] += nodes->mass[an_i];
#pragma omp atomic update
        elements->FHydro[clusterIndex] += 1.0 * diffVel;
#pragma omp atomic update
        elements->MHydro[clusterIndex] += 1.0 * radius.cross(diffVel);
#endif
    }
}
__hdi__ void common_streaming(const unsigned int i, Node2* nodes, Wall2* walls) {
    // Convert index to active node index
    const unsigned int an_i = nodes->activeI[i];

    // coefficient for free-surface
    constexpr double C2x2 = 9.0;
    constexpr double C3x2 = 3.0;
    // coefficient for slip conditions
    const double S1 = PARAMS.slipCoefficient;
    const double S2 = (1.0 - PARAMS.slipCoefficient);
    // creating list for collision function @todo can this be precomputed, rather than once per node?
    std::array<double, lbmDirec> staticPres;
    for (int j = 0; j < lbmDirec; j++) {
        staticPres[j] = PARAMS.fluidMaterial.initDensity * coeff[j];
    }

    // coefficient for bounce-back
    constexpr double BBCoeff = 2.0 * 3.0;

    const unsigned int A_OFFSET = an_i * lbmDirec;
    // cycling through neighbours
    for (unsigned int j = 1; j < lbmDirec; ++j) {
        // getting neighbour index
        const unsigned int ln_i = nodes->d[nodes->count * j + an_i];
        // if neighbour is normal fluid cell what follows is true

        if (ln_i == std::numeric_limits<unsigned int>::max()) { // is gas
            // additional variables for equilibrium f computation
            const double usq = nodes->u[an_i].norm2();
            const double vuj = nodes->u[an_i].dot(v[j]);
            // streaming with constant pressure interface
            nodes->f[A_OFFSET + opp[j]] = -nodes->fs[A_OFFSET + j] + coeff[j] * PARAMS.fluidMaterial.initDensity * (2.0 + C2x2 * (vuj * vuj) - C3x2 * usq);
        } else {
            const unsigned int L_OFFSET = ln_i * lbmDirec;
            // @todo this could be improved by stacking matching cases to reduce divergence
            switch (nodes->type[ln_i]) {
            case LIQUID:
            {
                nodes->f[A_OFFSET + opp[j]] = nodes->fs[L_OFFSET + opp[j]];
                break;
            }
            case INTERFACE:
            {
#ifdef DEBUG
                // TEST USING AGE //////////////////////////////////////
                const double usq = nodes->u[an_i].norm2();
                const double vuj = nodes->u[an_i].dot(v[j]);
                nodes->f[A_OFFSET + opp[j]] = nodes->age[ln_i] * nodes->fs[L_OFFSET + opp[j]] +
                    (1.0 - nodes->age[ln_i]) * (-nodes->fs[A_OFFSET + j] + coeff[j] * PARAMS.fluidMaterial.initDensity * (2.0 + C2x2 * (vuj * vuj) - C3x2 * usq));
#else

                nodes->f[A_OFFSET + opp[j]] = nodes->fs[L_OFFSET + opp[j]];
#endif
                break;

            }
            // for walls there is simple bounce-back
            case STAT_WALL:
            {
#ifndef DEBUG 
                if (nodes->type[an_i] == INTERFACE) {
                    // additional variables for equilibrium f computation
                    const double usq = nodes->u[an_i].norm2();
                    const double vuj = nodes->u[an_i].dot(v[j]);
                    //streaming with constant pressure interface
                    nodes->f[A_OFFSET + opp[j]] = -nodes->fs[A_OFFSET + j] + coeff[j] * PARAMS.fluidMaterial.initDensity * (2.0 + C2x2 * (vuj * vuj) - C3x2 * usq);
                    break;
                }
#endif      
                // getting the index of the wall to compute force in the right object
                const unsigned int solidIndex = nodes->solidIndex[ln_i];

                // static pressure is subtracted in order to correctly compute buoyancy for floating objects
                const tVect BBforce = nodes->bounceBackForce(an_i, j, staticPres, 0.0);
                // updating force and torque on the object (lattice units). This point is critical since many nodes update the force on the same object (lattice units)
#ifdef __CUDA_ARCH__
                    // CUDA atomics
                atomicAdd(&walls->FHydro[solidIndex].x, BBforce.x);
                atomicAdd(&walls->FHydro[solidIndex].y, BBforce.y);
                atomicAdd(&walls->FHydro[solidIndex].z, BBforce.z);
#else
                    // CPU atomics
#pragma omp atomic update
                walls->FHydro[solidIndex] += BBforce;
#endif
                nodes->f[A_OFFSET + opp[j]] = nodes->fs[A_OFFSET + j];
                break;
            }
            // for curved walls there is the rule of Mei-Luo-Shyy
            case TOPO:
            {
                nodes->f[A_OFFSET + opp[j]] = nodes->fs[A_OFFSET + j];
                break;
            }
            case OUTLET:
            {
                nodes->f[A_OFFSET + opp[j]] = std::min(nodes->fs[A_OFFSET + opp[j]], nodes->fs[A_OFFSET + j]);
                break;
            }
            // for moving walls there is simple bounce-back with velocity correction
            case DYN_WALL:
            {
                // getting the index of the wall to compute force in the right object
                const unsigned int solidIndex = nodes->solidIndex[ln_i];
                // velocity of the wall
                const tVect vel = nodes->u[ln_i];
                // variation in Bounce-Back due to moving object
                const double BBi = BBCoeff * nodes->n[an_i] * coeff[j] * vel.dot(v[j]); // mass!!!!!

                // static pressure is subtracted in order to correctly compute buoyancy for floating objects
                const tVect BBforce = nodes->bounceBackForce(an_i, j, staticPres, BBi);
                // updating force and torque on the object (lattice units). This point is critical since many nodes update the force on the same object (lattice units)
#ifdef __CUDA_ARCH__
                    // CUDA atomics
                atomicAdd(&walls->FHydro[solidIndex].x, BBforce.x);
                atomicAdd(&walls->FHydro[solidIndex].y, BBforce.y);
                atomicAdd(&walls->FHydro[solidIndex].z, BBforce.z);
#else
                    // CPU atomics
#pragma omp atomic update
                walls->FHydro[solidIndex] += BBforce;
#endif
                nodes->f[A_OFFSET + opp[j]] = nodes->fs[A_OFFSET + j] - BBi;
                // adding the extra mass to the surplus //@todo extraMass required for parity
                // extraMass = BBi * nodes->mass[an_i];  // redistributeMass() currently not used, so this isn't implemented properly
                break;
            }// for walls there is simple bounce-back
            case OBJ:
            {
                // getting the index of the wall to compute force in the right object
                const unsigned int solidIndex = nodes->solidIndex[ln_i];
                // static pressure is subtracted in order to correctly compute buoyancy for floating objects
                const tVect BBforce = nodes->bounceBackForce(an_i, j, staticPres, 0.0);
                // updating force and torque on the object (lattice units). This point is critical since many nodes update the force on the same object (lattice units)
#ifdef __CUDA_ARCH__
                    // CUDA atomics
                atomicAdd(&walls->FHydro[solidIndex].x, BBforce.x);
                atomicAdd(&walls->FHydro[solidIndex].y, BBforce.y);
                atomicAdd(&walls->FHydro[solidIndex].z, BBforce.z);
#else
                    // CPU atomics
#pragma omp atomic update
                walls->FHydro[solidIndex] += BBforce;
#endif
                nodes->f[A_OFFSET + opp[j]] = nodes->fs[A_OFFSET + j];
                break;
            }
            case SLIP_STAT_WALL:
            {
                if (j > 6) {
                    const unsigned int nodeCheck1 = nodes->d[slip1Check[j] * nodes->count + an_i];
                    const unsigned int nodeCheck2 = nodes->d[slip2Check[j] * nodes->count + an_i];
                    // check for the environment
                    const bool active1 = nodeCheck1 != std::numeric_limits<unsigned int>::max() && nodes->isActive(nodeCheck1);
                    const bool active2 = nodeCheck2 != std::numeric_limits<unsigned int>::max() && nodes->isActive(nodeCheck2);
                    // given the environment, perform the right operation
                    if (active1 && !active2) {
                        // first
                        nodes->f[A_OFFSET + opp[j]] = S1 * nodes->fs[nodeCheck1 * lbmDirec + slip1[j]] + S2 * nodes->fs[A_OFFSET + j];
                    }
                    else if (!active1 && active2) {
                        // second
                        nodes->f[A_OFFSET + opp[j]] = S1 * nodes->fs[nodeCheck2 * lbmDirec + slip2[j]] + S2 * nodes->fs[A_OFFSET + j];
                    }
                    else {
                        // standard BB
                        nodes->f[A_OFFSET + opp[j]] = nodes->fs[A_OFFSET + j];
                    }
                }
                else {
                    // standard BB
                    nodes->f[A_OFFSET + opp[j]] = nodes->fs[A_OFFSET + j];
                }
                break;
            }
            case SLIP_DYN_WALL:
            {
                // velocity of the wall
                const tVect vel = nodes->u[ln_i];
                // variation in Bounce-Back due to moving object
                const double BBi = BBCoeff * nodes->n[an_i] * coeff[j] * vel.dot(v[j]);
                if (j > 6) {
                    const unsigned int nodeCheck1 = nodes->d[slip1Check[j] * nodes->count + an_i];
                    const unsigned int nodeCheck2 = nodes->d[slip2Check[j] * nodes->count + an_i];
                    // check for the environment
                    const bool active1 = nodeCheck1 != std::numeric_limits<unsigned int>::max() && nodes->isActive(nodeCheck1);
                    const bool active2 = nodeCheck2 != std::numeric_limits<unsigned int>::max() && nodes->isActive(nodeCheck2);
                    // given the environment, perform the right operation
                    if (active1 && !active2) {
                        // first
                        nodes->f[A_OFFSET + opp[j]] = S1 * nodes->fs[nodeCheck1 * lbmDirec + slip1[j]] + S2 * (nodes->fs[A_OFFSET + j] - BBi);
                        // adding the extra mass to the surplus //@todo extraMass required for parity
                        // extraMass += S2 * nodes->mass[an_i] * BBi;  // redistributeMass() currently not used, so this isn't implemented properly
                    }
                    else if (!active1 && active2) {
                        // second
                        nodes->f[A_OFFSET + opp[j]] = S1 * nodes->fs[nodeCheck2 * lbmDirec + slip2[j]] + S2 * (nodes->fs[A_OFFSET + j] - BBi);
                        // adding the extra mass to the surplus //@todo extraMass required for parity
                        // extraMass += S2 * nodes->mass[an_i] * BBi;  // redistributeMass() currently not used, so this isn't implemented properly
                    }
                    else {
                        // standard BB
                        nodes->f[A_OFFSET + opp[j]] = nodes->fs[A_OFFSET + j] - BBi;
                        // adding the extra mass to the surplus //@todo extraMass required for parity
                        // extraMass += nodes->mass[an_i] * BBi;  // redistributeMass() currently not used, so this isn't implemented properly
                    }
                }
                else {
                    // standard BB
                    nodes->f[A_OFFSET + opp[j]] = nodes->fs[A_OFFSET + j] - BBi;
                    // adding the extra mass to the surplus //@todo extraMass required for parity
                    // extraMass += nodes->mass[an_i] * BBi;  // redistributeMass() currently not used, so this isn't implemented properly
                }
                break;
            }
            case UNUSED:
            case GAS:
            case PERIODIC:
            case CYL:
            default:
            {
                {
                    // @todo This may print out of order if multiple threads break in parallel
                    tVect pos = nodes->getPosition(an_i);
                    printf("%u(%f, %f, %f) %s TYPE ERROR:\n", an_i, pos.x, pos.y, pos.z, typeString(nodes->type[an_i]));
                    for (unsigned int k = 1; k < lbmDirec; ++k) {
                        printf("before error: j=%u link=%u\n", k, nodes->d[k * nodes->count + an_i]);
                    }
                    pos = nodes->getPosition(ln_i);
                    printf("(%f, %f, %f) %s TYPE ERROR\n", pos.x, pos.y, pos.z, typeString(nodes->type[ln_i]));
                    // @todo aborting from CUDA is harder, especially if the printf() is to be saved
#ifndef __CUDA_ARCH__
                    std::abort();
#endif
                    return;
                }
                break;

            }
            }
        }
    }
}

//
// latticeBoltzmannFreeSurfaceStep() subroutines
//
__hdi__ void common_updateMassInterface(const unsigned int in_i, Node2 *nodes) {
    // mass for interface nodes is regulated by the evolution equation
    nodes->newMass[in_i] = nodes->mass[in_i];
    // additional mass streaming to/from interface
    double deltaMass = 0.0;
    const unsigned int nodeCount = nodes->count;
    // cycling through neighbors
    for (unsigned int j = 1; j < lbmDirec; ++j) {
        // getting neighbor index
        const unsigned int ln_i = nodes->d[nodeCount * j + in_i];
        // average liquid fraction
        if (ln_i == std::numeric_limits<unsigned int>::max()) {
            // do nothing
        } else if (nodes->type[ln_i] == INTERFACE) {
            // average liquid fraction
            const double averageMass = 0.5 * (nodes->mass[ln_i] / nodes->n[ln_i] + nodes->mass[in_i] / nodes->n[in_i]);
            deltaMass += averageMass * nodes->massStream(in_i, j);
        } else if (nodes->type[ln_i] == LIQUID) {
            const double averageMass = 1.0;
            deltaMass += averageMass * nodes->massStream(in_i, j);
        } else if (nodes->type[ln_i] == DYN_WALL) {
            const double averageMass = 1.0 * nodes->mass[in_i];
            deltaMass += averageMass * nodes->massStream(in_i, j);
        } else if (nodes->type[ln_i] == CYL) {
            const double averageMass = 1.0 * nodes->mass[in_i];
            deltaMass += averageMass * nodes->massStream(in_i, j);
        } else if (nodes->type[ln_i] == SLIP_DYN_WALL) {
            if (j > 6) {
                bool active1 = false;
                bool active2 = false;
                const unsigned int c1_i = nodes->d[nodeCount * slip1Check[j] + in_i];
                const unsigned int c2_i = nodes->d[nodeCount * slip2Check[j] + in_i];
                // check for the environment
                if (c1_i != std::numeric_limits<unsigned int>::max()) {
                    if (nodes->isActive(c1_i)) {
                        active1 = true;
                    }
                }
                if (c2_i != std::numeric_limits<unsigned int>::max()) {
                    if (nodes->isActive(c2_i)) {
                        active2 = true;
                    }
                }
                // given the environment, perform the right operation
                double averageMass = 0.0;
                if (active1 && !active2) {
                    // adding the extra mass to the surplus
                    averageMass += 1.0 * (1.0 - PARAMS.slipCoefficient) * nodes->mass[in_i];
                } else if (!active1 && active2) {
                    // adding the extra mass to the surplus
                    averageMass += 1.0 * (1.0 - PARAMS.slipCoefficient) * nodes->mass[in_i];
                } else {
                    // adding the extra mass to the surplus
                    averageMass += 1.0 * nodes->mass[in_i];
                }
                deltaMass += averageMass * nodes->massStream(in_i, j);
            } else {
                // adding the extra mass to the surplus
                const double averageMass = 1.0 * nodes->mass[in_i];
                deltaMass += averageMass * nodes->massStream(in_i, j);
            }
        }
    }
    nodes->newMass[in_i] += deltaMass;

    nodes->mass[in_i] = nodes->newMass[in_i];
    nodes->age[in_i] = min(nodes->age[in_i] + PARAMS.ageRatio, 1.0f);
}
__hdi__ void common_updateMassFluid(const unsigned int fn_i, Node2 *nodes) {
    nodes->mass[fn_i] = nodes->n[fn_i];
    nodes->age[fn_i] = min(nodes->age[fn_i] + PARAMS.ageRatio, 1.0f);
}
__hdi__ void common_findInterfaceMutants(const unsigned int in_i, Node2* nodes) {
    // CHECKING FOR NEW FLUID NODES from filling
    if (nodes->mass[in_i] > nodes->n[in_i]) {
        // updating type
        nodes->type[in_i] = INTERFACE_FILLED;
    }// CHECKING FOR NEW GAS NODES from emptying
    else if (nodes->mass[in_i] < 0.0) {
        // updating type
        nodes->type[in_i] = INTERFACE_EMPTY;
    }
}
__hdi__ void common_smoothenInterface_find(const unsigned int in_i, Node2* nodes) {
    // CHECKING FOR NEW INTERFACE NODES from neighboring a new fluid node
    if (nodes->type[in_i] == INTERFACE_FILLED) {
        // neighor indices
        const std::array<unsigned int, lbmDirec> neighborCoord = nodes->findNeighbors(in_i);
        // cycling through neighbors
        for (int j = 1; j < lbmDirec; ++j) {
            // neighbor index
            const unsigned int ln_i = neighborCoord[j];
            // checking if node is gas (so to be transformed into interface)
            if (ln_i < nodes->count && nodes->type[ln_i] == GAS) { // @todo this should probably include INTERFACE_EMPTY (see issue #5)
                nodes->type[ln_i] = GAS_TO_INTERFACE;
                // Track source node for the copy in part 2
                nodes->d[ln_i] = in_i; // race condition, but eh
            }
        }
    }

    // CHECKING FOR NEW INTERFACE NODES from neighboring a new gas node
    // tested unordered_set, was slower
    if (nodes->type[in_i] == INTERFACE_EMPTY) {
        // neighor indices
        const std::array<unsigned int, lbmDirec> neighborCoord = nodes->findNeighbors(in_i);
        // cycling through neighbors
        for (int j = 1; j < lbmDirec; ++j) {
            // neighbor node
            const unsigned int ln_i = neighborCoord[j];
            if (ln_i < nodes->count && (nodes->type[ln_i] == LIQUID || nodes->type[ln_i] == INTERFACE_FILLED)) {
                nodes->type[ln_i] = FLUID_TO_INTERFACE;
            }
        }
    }
}
__hdi__ void common_smoothenInterface_update(const unsigned int in_i, Node2* nodes) {
    constexpr double marginalMass = 1.0e-2;
    // CHECKING FOR NEW INTERFACE NODES from neighboring a new fluid node
    if (nodes->type[in_i] == GAS_TO_INTERFACE) {
        // create new interface node
        nodes->generateNode(in_i, INTERFACE);
        // node is becoming active and needs to be initialized
        double massSurplusHere = -marginalMass * PARAMS.fluidMaterial.initDensity;
        // same density and velocity; 1% of the mass
        nodes->copy(in_i, nodes->d[in_i]); // d[0] contains src node
        nodes->mass[in_i] = -massSurplusHere;
        // the 1% of the mass is taken form the surplus
        nodes->scatterMass(in_i, massSurplusHere);  // @TODO race condition on extraMass (not currently enabled as redundant)?
        // massSurplus += massSurplusHere;
    }


    // CHECKING FOR NEW INTERFACE NODES from neighboring a new gas node
    // tested unordered_set, was slower
    else if (nodes->type[in_i] == FLUID_TO_INTERFACE) {
        // ln_i should equal nodes->d[in_i * nodes.count + j];
        nodes->type[in_i] = INTERFACE;
        double massSurplusHere = marginalMass * nodes->n[in_i];
        // characteristics are inherited by previous fluid cell. Only mass must be updated to 99% of initial mass
        nodes->mass[in_i] = nodes->n[in_i] - massSurplusHere;
        // the remaining 1% of the mass is added to the surplus
        nodes->scatterMass(in_i, massSurplusHere);
        //massSurplus += massSurplusHere;
    }
}
__hdi__ void common_updateMutants(const unsigned int in_i, Node2* nodes, double *massSurplus) {
    // resetting new gas macroscopic quantities
    if (nodes->type[in_i] == INTERFACE_EMPTY) {
        // updating mass surplus
#ifdef __CUDA_ARCH__
        // CUDA atomics
        atomicAdd(massSurplus, nodes->mass[in_i]);
#else
        // CPU atomics
        #pragma omp atomic update
        *massSurplus += nodes->mass[in_i];
#endif
        // deleting node
        nodes->eraseNode(in_i);
    }

    // resetting new fluid macroscopic quantities
    if (nodes->type[in_i] == INTERFACE_FILLED) {
        // updating mass surplus
#ifdef __CUDA_ARCH__
        // CUDA atomics
        atomicAdd(massSurplus, nodes->mass[in_i] - nodes->n[in_i]);
#else
        // CPU atomics
        #pragma omp atomic update
        *massSurplus += (nodes->mass[in_i] - nodes->n[in_i]);
#endif
        // setting liquid fraction for new fluid cell (other macroscopic characteristics stay the same)
        nodes->mass[in_i] = nodes->n[in_i];
        // Complete it's conversion to type LIQUID
        nodes->type[in_i] = LIQUID;
    }
}
__hdi__ void common_removeIsolated(const unsigned int in_i, Node2* nodes, double *massSurplus) {
    // remove isolated interface cells (surrounded either by only fluid or only solid cells)

    // checking if it is surrounded by fluid (in that case is converted to fluid). Solid is an exception
    // reverse cycle is needed because of deletion function
    {
        bool surroundedFluid = true;
        for (int j = 1; j < lbmDirec; ++j) {
            const unsigned int ln_i = nodes->d[nodes->count * j + in_i];
            if (ln_i == std::numeric_limits<unsigned int>::max() || nodes->type[ln_i] == GAS) {
                surroundedFluid = false;
                break;
            }
        }
        if (surroundedFluid) {
            // update mass storage for balance
#ifdef __CUDA_ARCH__
            // CUDA atomics
            atomicAdd(massSurplus, nodes->mass[in_i] - nodes->n[in_i]);
#else
            // CPU atomics
            #pragma omp atomic update
            *massSurplus += (nodes->mass[in_i] - nodes->n[in_i]);
#endif
            // update characteristics (inherited from the gas node)
            nodes->mass[in_i] = nodes->n[in_i];
            nodes->type[in_i] = LIQUID;
        }
    }

    // checking if it is surrounded by gas (in that case is converted to gas)
    // or, better, if it is not connected to fluid (could be connected to walls or particles)
    {
        bool surroundedGas = true;
        for (int j = 1; j < lbmDirec; ++j) {
            const unsigned int ln_i = nodes->d[nodes->count * j + in_i];
            if (ln_i != std::numeric_limits<unsigned int>::max()) {
                if (nodes->type[ln_i] == LIQUID) {
                    surroundedGas = false;
                    break;
                }
            }
        }
        // updating mass surplus
        if (surroundedGas) {
            // update mass
#ifdef __CUDA_ARCH__
            // CUDA atomics
            atomicAdd(massSurplus, nodes->mass[in_i]);
#else
            // CPU atomics
            #pragma omp atomic update
            *massSurplus += nodes->mass[in_i];
#endif
            nodes->eraseNode(in_i);
        }
    }
}

#endif // LBCOMMON_H
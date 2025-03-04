#include "Element2.h"

#include <vector>

#include "myvector.h"
#include "DEMParams.h"
#include "Particle2.h"


void Element2::allocComponentsData() {
    if (count) {
        if (componentsIndex && componentsIndex[count]) {
            if (componentsData) {
                free(componentsData);
            }
            componentsData = (unsigned int*)malloc(componentsIndex[count] * sizeof(unsigned int));
        }
    }
}
void Element2::initialize(const unsigned int index, const double partDensity) {
    // translational degrees of freedom
    xp0[index] = x0[index];
    xp1[index] = x1[index];
    xp2[index].reset();
    x2[index].reset();
    xp3[index].reset();
    x3[index].reset();
    xp4[index].reset();
    x4[index].reset();
    xp5[index].reset();
    x5[index].reset();

    // rotational degrees of freedom
    qp0[index] = q0[index];
    qp1[index] = q1[index];
    qp2[index] = q2[index] =tQuat(0.0,0.0,0.0,0.0);
    qp3[index] = q3[index] =tQuat(0.0,0.0,0.0,0.0);
    qp4[index] = q4[index] =tQuat(0.0,0.0,0.0,0.0);
    qp5[index] = q5[index] =tQuat(0.0,0.0,0.0,0.0);
    
    if (wSolver[index]) {
        wpGlobal[index] = wGlobal[index] = wpLocal[index] = wLocal[index] = w0[index];
    } else {
        const tQuat q0adj=q0[index].adjoint();
        wp0[index] = w0[index] = wpGlobal[index] = wGlobal[index] = 2.0*quat2vec( q1[index].multiply(q0adj) );
        wpLocal[index] = wLocal[index] = 2.0*quat2vec( q0adj.multiply(q1[index]));
    }
    
    wp1[index].reset();
    w1[index].reset();
    wp2[index].reset();
    w2[index].reset();
    wp3[index].reset();
    w3[index].reset();
    wp4[index].reset();
    w4[index].reset();
    wp5[index].reset();
    w5[index].reset();
    
    // tQuat wQuat=vec2quat(wGlobal[index]);
    // qp1[index] = q1[index] = 0.5*q0[index].multiply(wQuat);

    // calculated variables (the element is supposed to be a sphere for the moment)
    // mass
    const double singleMass=4.0/3.0*partDensity*M_PI*radius[index]*radius[index]*radius[index];
    m[index] = size[index] * singleMass;
    // inertia moment (diagonal) - Huygens-Steiner theorem
    // inertia of single spheres
    I[index] = size[index] * 2.0/5.0 * singleMass * radius[index] * radius[index] * tVect(1.0,1.0,1.0);
    // transport components
    for (unsigned int n=0; n<size[index]; ++n) {
        I[index] += singleMass*radius[index] *radius[index] * prototypes[size[index]][n].transport();
    }

    active[index] = true;
    
    // initialize forces
    FHydro[index].reset();
    FParticle[index].reset();
    FWall[index].reset();
    FGrav[index] = DEM_P.demF*m[index];
    MHydro[index].reset();
    MParticle[index].reset();
    MWall[index].reset();    
}

void Element2::generateParticles(const unsigned int e_i, unsigned int& p_i, Particle2 &particles) {
    for (unsigned int i = 0; i < size[e_i]; ++i) {
        // add standard particle
        particles.particleIndex[p_i] = p_i;
        componentsData[componentsIndex[e_i] + i] = p_i;
        particles.clusterIndex[p_i] = index[e_i];
        particles.r[p_i] = radius[e_i];
        particles.protoIndex[p_i] = i;
        particles.isGhost[p_i] = false;
        particles.active[p_i] = true;
        //@todo how do springs look? How will they exist on gpu
        /*
        particles.springs[p_i].resize(4);
        for (int t = 1; t < 4; ++t) {
            particles.springs[t].clear();
        }
        */
        particles.updateCorrected(p_i, this); // was predicted
        ++p_i;
    }
}

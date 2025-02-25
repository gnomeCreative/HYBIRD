#include "Element2.h"

#include <vector>

#include "myvector.h"
#include "DEMParams.h"
#include "Particle2.h"

template<>
inline void Element2::memoryAlloc<CPU>(unsigned int num) {
    alloc = num;

    if(wSolver) {
        free(wSolver);
        free(index);
        free(componentsIndex);
        free(size);
        free(m);
        free(I);
        free(x0);
        free(x1);
        free(x2);
        free(x3);
        free(x4);
        free(x5);
        free(xp5);
        free(xp0);
        free(xp1);
        free(xp2);
        free(xp3);
        free(xp4);
        free(xp5);
        free(x0history);
        free(q0);
        free(q1);
        free(q2);
        free(q3);
        free(q4);
        free(q5);
        free(qp5);
        free(qp0);
        free(qp1);
        free(qp2);
        free(qp3);
        free(qp4);
        free(qp5);
        free(wGlobal);
        free(wLocal);
        free(wpLocal);
        free(wpGlobal);
        free(w0);
        free(w1);
        free(w2);
        free(w3);
        free(w4);
        free(w5);
        free(wp5);
        free(wp0);
        free(wp1);
        free(wp2);
        free(wp3);
        free(wp4);
        free(wp5);
        free(FHydro);
        free(FParticle);
        free(FWall);
        free(FGrav);
        free(FSpringP);
        free(FSpringW);
        free(MHydro);
        free(MParticle);
        free(MWall);
        free(MRolling);
        free(solidIntensity);
        free(coordination);
        free(ACoriolis);
        free(ACentrifugal);
        free(fluidVolume);
        free(maxOverlap);
        free(maxDtOverlap);
        free(slippingCase);
    }

    wSolver = (bool*)malloc(alloc * sizeof(bool));
    index = (unsigned int*)malloc(alloc * sizeof(unsigned int));
    active = (bool*)malloc(alloc * sizeof(bool));
    componentsIndex = (unsigned int*)malloc((alloc + 1) * sizeof(unsigned int));    
    // componentsData = (unsigned int*)malloc(? * sizeof(unsigned int)); // This is handled later once index is full
    size = (unsigned int*)malloc(alloc * sizeof(unsigned int));
    radius = (double*)malloc(alloc * sizeof(double));
    m = (double*)malloc(alloc * sizeof(double));
    I = (tVect*)malloc(alloc * sizeof(tVect));
    x0 = (tVect*)malloc(alloc * sizeof(tVect));
    x1 = (tVect*)malloc(alloc * sizeof(tVect));
    x2 = (tVect*)malloc(alloc * sizeof(tVect));
    x3 = (tVect*)malloc(alloc * sizeof(tVect));
    x4 = (tVect*)malloc(alloc * sizeof(tVect));
    x5 = (tVect*)malloc(alloc * sizeof(tVect));
    xp0 = (tVect*)malloc(alloc * sizeof(tVect));
    xp1 = (tVect*)malloc(alloc * sizeof(tVect));
    xp2 = (tVect*)malloc(alloc * sizeof(tVect));
    xp3 = (tVect*)malloc(alloc * sizeof(tVect));
    xp4 = (tVect*)malloc(alloc * sizeof(tVect));
    xp5 = (tVect*)malloc(alloc * sizeof(tVect));
    x0history = (tVect*)malloc(alloc * sizeof(tVect));
    q0 = (tQuat*)malloc(alloc * sizeof(tQuat));
    q1 = (tQuat*)malloc(alloc * sizeof(tQuat));
    q2 = (tQuat*)malloc(alloc * sizeof(tQuat));
    q3 = (tQuat*)malloc(alloc * sizeof(tQuat));
    q4 = (tQuat*)malloc(alloc * sizeof(tQuat));
    q5 = (tQuat*)malloc(alloc * sizeof(tQuat));
    qp0 = (tQuat*)malloc(alloc * sizeof(tQuat));
    qp1 = (tQuat*)malloc(alloc * sizeof(tQuat));
    qp2 = (tQuat*)malloc(alloc * sizeof(tQuat));
    qp3 = (tQuat*)malloc(alloc * sizeof(tQuat));
    qp4 = (tQuat*)malloc(alloc * sizeof(tQuat));
    qp5 = (tQuat*)malloc(alloc * sizeof(tQuat));
    wGlobal = (tVect*)malloc(alloc * sizeof(tVect));
    wLocal = (tVect*)malloc(alloc * sizeof(tVect));
    wpLocal = (tVect*)malloc(alloc * sizeof(tVect));
    wpGlobal = (tVect*)malloc(alloc * sizeof(tVect));
    w0 = (tVect*)malloc(alloc * sizeof(tVect));
    w1 = (tVect*)malloc(alloc * sizeof(tVect));
    w2 = (tVect*)malloc(alloc * sizeof(tVect));
    w3 = (tVect*)malloc(alloc * sizeof(tVect));
    w4 = (tVect*)malloc(alloc * sizeof(tVect));
    w5 = (tVect*)malloc(alloc * sizeof(tVect));
    wp0 = (tVect*)malloc(alloc * sizeof(tVect));
    wp1 = (tVect*)malloc(alloc * sizeof(tVect));
    wp2 = (tVect*)malloc(alloc * sizeof(tVect));
    wp3 = (tVect*)malloc(alloc * sizeof(tVect));
    wp4 = (tVect*)malloc(alloc * sizeof(tVect));
    wp5 = (tVect*)malloc(alloc * sizeof(tVect));
    FHydro = (tVect*)malloc(alloc * sizeof(tVect));
    FParticle = (tVect*)malloc(alloc * sizeof(tVect));
    FWall = (tVect*)malloc(alloc * sizeof(tVect));
    FGrav = (tVect*)malloc(alloc * sizeof(tVect));
    FSpringP = (tVect*)malloc(alloc * sizeof(tVect));
    FSpringW = (tVect*)malloc(alloc * sizeof(tVect));
    MHydro = (tVect*)malloc(alloc * sizeof(tVect));
    MParticle = (tVect*)malloc(alloc * sizeof(tVect));
    MWall = (tVect*)malloc(alloc * sizeof(tVect));
    MRolling = (tVect*)malloc(alloc * sizeof(tVect));
    solidIntensity = (tVect*)malloc(alloc * sizeof(tVect));
    coordination = (unsigned int*)malloc(alloc * sizeof(unsigned int));
    ACoriolis = (tVect*)malloc(alloc * sizeof(tVect));
    ACentrifugal = (tVect*)malloc(alloc * sizeof(tVect));
    fluidVolume = (double*)malloc(alloc * sizeof(double));
    maxOverlap = (double*)malloc(alloc * sizeof(double));
    maxDtOverlap = (double*)malloc(alloc * sizeof(double));
    slippingCase = (int*)malloc(alloc * sizeof(int));
}
template <>
inline void Element2::initElements<CPU>() {
    //initializing this time step hydrodynamic force
    memset(FHydro, 0, sizeof(tVect) * count);
    memset(MHydro, 0, sizeof(tVect) * count);
    // initializing the fluid mass for buoyancy
    memset(fluidVolume, 0, sizeof(double) * count);
}
#ifdef USE_CUDA
template <>
inline void Element2::initElements<CUDA>() {
    //initializing this time step hydrodynamic force
    CUDA_CALL(cudaMemset(FHydro, 0, sizeof(tVect) * count));
    CUDA_CALL(cudaMemset(MHydro, 0, sizeof(tVect) * count));
    // initializing the fluid mass for buoyancy
    CUDA_CALL(cudaMemset(fluidVolume, 0, sizeof(double) * count));
}
#endif

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
        I[index] += singleMass*radius[index] *radius[index] * DEM_P.prototypes[size[index]][n].transport();
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
        particles.clusterIndex = index;
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
        particles.updateCorrected(p_i, *this, e_i); // was predicted
        ++p_i;
    }
}

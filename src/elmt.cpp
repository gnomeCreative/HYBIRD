
#include "myvector.h"
#include "macros.h"
#include "elmt.h"
//#include "LB.h"

using namespace std;

void elmt::elmtShow() const {
    cout<<"Element number "<<index<<" with "<<size<<" particles of radius "<<radius<<"\n";
    cout<<"Position: "; x0.show();
    cout<<"; Velocity: "; x1.show(); cout<<";\n";
    cout<<"Orientation: "; q0.show(); cout<<";\n";
    cout<<"Inertia: "; I.show(); cout<<";\n";
    cout<<"Mass: "<<m<<";\n";
}

void elmt::initialize(const double& partDensity, std::vector <vecList>& prototypes, tVect& demF, doubleList& demSize) {

    // translational degrees of freedom
    xp0=x0;
    xp1=x1;
    xp2.reset();
    x2.reset();
    xp3.reset();
    x3.reset();
    xp4.reset();
    x4.reset();
    xp5.reset();
    x5.reset();

    // rotational degrees of freedom
    qp0=q0;
    qp1=q1;
    qp2=q2=tQuat(0.0,0.0,0.0,0.0);
    qp3=q3=tQuat(0.0,0.0,0.0,0.0);
    qp4=q4=tQuat(0.0,0.0,0.0,0.0);
    qp5=q5=tQuat(0.0,0.0,0.0,0.0);
    
    if (wSolver) {
        wpGlobal=wGlobal=wpLocal=wLocal=w0;
    }
    else {
        const tQuat q0adj=q0.adjoint();
        wp0=w0=wpGlobal=wGlobal=2.0*quat2vec( q1.multiply(q0adj) );
        wpLocal=wLocal=2.0*quat2vec( q0adj.multiply(q1) );
    }
    
    wp1.reset();
    w1.reset();
    wp2.reset();
    w2.reset();
    wp3.reset();
    w3.reset();
    wp4.reset();
    w4.reset();
    wp5.reset();
    w5.reset();
    
    // tQuat wQuat=vec2quat(wGlobal);
    // qp1=q1=0.5*q0.multiply(wQuat);

    // calculated variables (the element is supposed to be a sphere for the moment)
    // mass
    // mass correction particles must behave like cylinders or whatever shape
    double singleMass = 0.0;
    switch (problemName){
        case (TBAR): { // cylinders
            singleMass = partDensity*radius*radius*M_PI* (demSize[1]);
            break;
        }
        case (SEGUIN):{ // cylinders
            singleMass = partDensity*radius*radius*M_PI* (demSize[1]);
            break;
        }
        default:{ // spheres
            singleMass=4.0/3.0*partDensity*M_PI*radius*radius*radius;
            break;
        }
    }
//    const double singleMass=4.0/3.0*partDensity*M_PI*radius*radius*radius;
    m=size*singleMass;
    // inertia moment (diagonal) - Huygens-Steiner theorem
    // inertia of single spheres
    // Inertia correction particles must behave like cylinders or whatever shape
    switch (problemName){
        case (TBAR):{ // cylinders
            I=size*0.5*singleMass*radius*radius*tVect(1.0,1.0,1.0);
            break;
            
        }
        case (SEGUIN): { // cylinders
            I=size*0.5*singleMass*radius*radius*tVect(1.0,1.0,1.0);
            break;
        } 
        default: { // spheres
            I=size*2.0/5.0*singleMass*radius*radius*tVect(1.0,1.0,1.0);
            break;
        }
    }
//    I=size*2.0/5.0*singleMass*radius*radius*tVect(1.0,1.0,1.0);
    // transport components
    for (int n=0; n<size; ++n) {
        I+=singleMass*radius*radius*prototypes[size][n].transport();
    }

    active=true;
    
    // initialize forces
    FHydro.reset();
    FParticle.reset();
    FWall.reset();
    FGrav=demF*m;
    MHydro.reset();
    MParticle.reset();
    MWall.reset();
    
}

void elmt::resetVelocity() {
    // resets translational and rotational velocity to zero

    // translational
    xp1.reset();
    x1.reset();
    xp2.reset();
    x2.reset();
    xp3.reset();
    x3.reset();
    xp4.reset();
    x4.reset();
    xp5.reset();
    x5.reset();
    // rotational
    wGlobal.reset();
    wLocal.reset();
    wpGlobal.reset();
    wpLocal.reset();
    wp1.reset();
    w1.reset();
    wp2.reset();
    w2.reset();
    wp3.reset();
    w3.reset();
    wp4.reset();
    w4.reset();
    wp5.reset();
    w5.reset();
    tQuat wQuat=vec2quat(wGlobal);
    qp1=q1=0.5*q0.multiply(wQuat);
    qp2=q2=tQuat(0.0,0.0,0.0,0.0);
    qp3=q3=tQuat(0.0,0.0,0.0,0.0);
    qp4=q4=tQuat(0.0,0.0,0.0,0.0);
    qp5=q5=tQuat(0.0,0.0,0.0,0.0);
}

void elmt::generateParticles(unsigned int& globalIndex, particleList& particles, const std::vector<vecList>& prototypes) {
    components.resize(size);

    for (int i = 0; i < size; ++i) {
        // add standard particle
        particle dummyPart;
        dummyPart.particleIndex = globalIndex;
        components[i] = dummyPart.particleIndex;
        dummyPart.clusterIndex = index;
        dummyPart.r = radius;
        dummyPart.protoIndex = i;
        dummyPart.isGhost = false;
        dummyPart.active=true;
        dummyPart.springs.resize(4);
        for (int t=1; t<4; ++t) {
            dummyPart.springs[t].clear();
        }
        dummyPart.updateCorrected(*this, prototypes); // was predicted
        particles.push_back(dummyPart);
        ++globalIndex;
    }
}

void elmt::predict(const double c1[], const double c2[]) {

    xp0 = x0 + x1 * c2[0] + x2 * c2[1] + x3 * c2[2] + x4 * c2[3] + x5 * c2[4];
    xp1 = x1 + x2 * c2[0] + x3 * c2[1] + x4 * c2[2] + x5 * c2[3];
    xp2 = x2 + x3 * c2[0] + x4 * c2[1] + x5 * c2[2];
    xp3 = x3 + x4 * c2[0] + x5 * c2[1];
    xp4 = x4 + x5 * c2[0];
    xp5 = x5;

    qp0 = q0 + q1 * c2[0] + q2 * c2[1] + q3 * c2[2] + q4 * c2[3] + q5 * c2[4];
    qp1 = q1 + q2 * c2[0] + q3 * c2[1] + q4 * c2[2] + q5 * c2[3];
    qp2 = q2 + q3 * c2[0] + q4 * c2[1] + q5 * c2[2];
    qp3 = q3 + q4 * c2[0] + q5 * c2[1];
    qp4 = q4 + q5 * c2[0];
    qp5 = q5;

    qp0.normalize();

    wp0 = w0 + w1 * c1[0] + w2 * c1[1] + w3 * c1[2] + w4 * c1[3] + w5 * c1[4];
    wp1 = w1 + w2 * c1[0] + w3 * c1[1] + w4 * c1[2] + w5 * c1[3];
    wp2 = w2 + w3 * c1[0] + w4 * c1[1] + w5 * c1[2];
    wp3 = w3 + w4 * c1[0] + w5 * c1[1];
    wp4 = w4 + w5 * c1[0];
    wp5 = w5;

    if (wSolver) {
        wpGlobal = wp0;
        wpLocal = project(wpGlobal, qp0.adjoint());
    } else {
        //q1.forceStability(q0);
        const tQuat qp0adj = qp0.adjoint();
        wpGlobal = 2.0 * quat2vec(qp1.multiply(qp0adj));
        wpLocal = 2.0 * quat2vec(qp0adj.multiply(qp1));
    }

}

void elmt::correct(const double coeff1ord[], const double coeff2ord[]) {

    const tVect x2Corr = x2 - xp2;

    x0 = xp0 + x2Corr * coeff2ord[0];
    x1 = xp1 + x2Corr * coeff2ord[1];
    // x2 calculated directly at the end of force routine
    x3 = xp3 + x2Corr * coeff2ord[3];
    x4 = xp4 + x2Corr * coeff2ord[4];
    x5 = xp5 + x2Corr * coeff2ord[5];

    xp0 = x0;
    xp1 = x1;
    xp2 = x2;
    xp3 = x3;
    xp4 = x4;
    xp5 = x5;

    const tVect w1Corr = w1 - wp1;

    w0 = wp0 + w1Corr * coeff1ord[0];
    // w1 calculated directly at the end of force routine
    w2 = wp2 + w1Corr * coeff1ord[2];
    w3 = wp3 + w1Corr * coeff1ord[3];
    w4 = wp4 + w1Corr * coeff1ord[4];
    w5 = wp5 + w1Corr * coeff1ord[5];

    wp0 = w0;
    wp1 = w1;
    wp2 = w2;
    wp3 = w3;
    wp4 = w4;
    wp5 = w5;

    const tQuat q2Corr = q2 - qp2;

    q0 = qp0 + q2Corr * coeff2ord[0];
    q1 = qp1 + q2Corr * coeff2ord[1];
    // q2 calculated directly at the end of force routine
    q3 = qp3 + q2Corr * coeff2ord[3];
    q4 = qp4 + q2Corr * coeff2ord[4];
    q5 = qp5 + q2Corr * coeff2ord[5];
    //normalization of q0
    q0.normalize();

    qp0 = q0;
    qp1 = q1;
    qp2 = q2;
    qp3 = q3;
    qp4 = q4;
    qp5 = q5;


    if (wSolver) {
        wGlobal = w0;
        wLocal = project(wGlobal, q0.adjoint());
    } else {
        //q1.forceStability(q0);
        const tQuat q0adj = q0.adjoint();
        wGlobal = 2.0 * quat2vec(q1.multiply(q0adj));
        wLocal = 2.0 * quat2vec(q0adj.multiply(q1));
    }
}

void elmt::translate(const tVect& transVec) {
    x0+=transVec;
    xp0+=transVec;
}

void particle::updatePredicted(const elmt& motherElmt, const std::vector <vecList>& prototypes) {

    // updating position and velocity for simple case
    x0=motherElmt.xp0;
    radiusVec=Zero;
    x1=motherElmt.xp1;
            
    if (motherElmt.size>1) {
        x0=x0+r*project(prototypes[motherElmt.size][protoIndex],motherElmt.qp0);
        // updating radius (distance of particle center of mass to element center of mass)
        radiusVec=x0-motherElmt.xp0;
        x1=x1+motherElmt.wpGlobal.cross(radiusVec);
    }
}

void particle::updateCorrected(const elmt& motherElmt, const std::vector <vecList>& prototypes) {
    
    // updating position and velocity for simple case
    x0=motherElmt.x0;
    radiusVec=Zero;
    x1=motherElmt.x1;
            
    if (motherElmt.size>1) {
        x0=x0+r*project(prototypes[motherElmt.size][protoIndex],motherElmt.q0);
        // updating radius (distance of particle center of mass to element center of mass)
        radiusVec=x0-motherElmt.x0;
        x1=x1+motherElmt.wGlobal.cross(radiusVec);
    }
    
}

void particle::ghostUpdate(particle& originParticle, tVect& pbcVector) {
        // updating position
        x0=originParticle.x0+pbcVector;
        // updating radius (distance of particle barycentre from element barycentre)
        radiusVec=originParticle.radiusVec;
        // updating particle speed
        x1=originParticle.x1;
}

void object::updateMax(const tVect& direction, const double& time) {
    if (maxFParticle.dot(direction) < FParticle.dot(direction)) {
        maxFParticle = FParticle;
        timeMaxFParticle = time;
    }
}

void object::saveForces() {
        savedFParticle = FParticle;
}

void Elongation::reset(){
     e.reset();
     p.reset();
     slipping=false;
}
void Elongation::copy(const Elongation& e2){
     e=e2.e;
     p=e2.p;
     slipping=e2.slipping;
}
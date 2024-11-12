/*
 * File:   elmt.h
 * Author: aleonard
 *
 * Created on November 5, 2013, 4:11 PM
 */

#ifndef ELMT_H
#define	ELMT_H

#include "myvector.h"

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// CLASS  DEFINITIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/


class particle{
public:
    // particle index
    unsigned int particleIndex;
    // belonging element index
    unsigned int clusterIndex;
    // index of the component of the prototype
    unsigned int protoIndex;
    // is it active? (false=destroyed)
    bool active;
    // is it a ghost?
    bool isGhost;
    // belonging cell for neighbor list
    unsigned int tableCell;
    // particle radius
    double r;
    // position of the particle
    tVect  x0;
    // velocity of the particle
    tVect x1;
    // vector connecting center of element to center of particle
    tVect radiusVec;
    // list with springs (0=particle-particle; 1=particle-wall; 2=particle-cylinder; 3=particle-object)
    std::vector<elongVector> springs;
    particle() {
        particleIndex=0;
        clusterIndex=0;
        protoIndex=0;
        tableCell=0;
        r=0.0;
        x0.reset();
        x1.reset();
        radiusVec.reset();
        //springs.clear();
        isGhost=false;
    }
    ~particle (){
        for (int i=0; i<springs.size(); ++i) {
            springs[i].clear();
        }
        
    }
    void updatePredicted(const elmt& motherElmt, const std::vector <vecList>& prototypes);
    void updateCorrected(const elmt& motherElmt, const std::vector <vecList>& prototypes);
    void ghostUpdate(particle& originParticle, tVect& pbcVector);
};

class elmt{
public:
    bool wSolver;
    // element index
    unsigned int index;
    // is it active? (false=destroyed)
    bool active;
    // constitutive particles indexes
    intList components;
    // number of constitutive particles
    unsigned int size;
    // radius of constitutive particles (supposed, for the moment, of being constant)
    double radius;
    // mass stored in the element
    double m;
    // inertia tensor on principal axis (diagonal)
    tVect I;
    // Position and derivatives
    // position of the center of mass of the element
    tVect x0;
    // velocity of the center of the element
    tVect x1;
    // acceleration of the center of the element
    tVect x2;
    // other velocity derivative for Gear scheme
    tVect x3,x4,x5;
    // predicted quantities for Gear scheme
    tVect xp0,xp1,xp2,xp3,xp4,xp5;
    // total displacement
    tVect x0history;
    // Orientation and derivatives
    // orientation of the center of mass of the element
    tQuat q0;
    // quaternion rates (related to angular velocity, acceleration...)
    tQuat q1,q2,q3,q4,q5;
    // predicted quantities for Gear scheme
    tQuat qp0,qp1,qp2,qp3,qp4,qp5;
    // angular velocity in global and local reference frame
    tVect wGlobal,wLocal;
    // predicted angular velocities
    tVect wpGlobal,wpLocal;
    // angular velocity rates (global)
    tVect w0,w1,w2,w3,w4,w5;
    // predicted angular velocity rates (global)
    tVect wp0,wp1,wp2,wp3,wp4,wp5;
    // forces and moments
    tVect FHydro,FParticle,FWall,FGrav,FSpringP,FSpringW;
    tVect MHydro,MParticle,MWall,MRolling;
    // force intensities
    tVect solidIntensity;
    // connectivity (coordination number)
    unsigned int coordination;
    // apparent accelerations
    tVect ACoriolis,ACentrifugal;
    // fluid mass entrained by the particle
    double fluidVolume;
    // container for maximum overlap, used for output
    double maxOverlap;
    // container for maximum overlap between time steps, used for output
    double maxDtOverlap;
    // whether the element is slipping on a wall
    int slippingCase;
    // default constructor
    elmt() {
        //
        wSolver=true;
        //
        active=true;
        size=1;
        radius=1.0;
        x0=x1=x2=x3=x4=x5=tVect(0.0,0.0,0.0);
        q0=tQuat(1.0,0.0,0.0,0.0);
        q1=q2=q3=q4=q5=tQuat(0.0,0.0,0.0,0.0);
        wGlobal=tVect(0.0,0.0,0.0);
        wLocal=tVect(0.0,0.0,0.0);
        wpGlobal=tVect(0.0,0.0,0.0);
        wpLocal=tVect(0.0,0.0,0.0);
        w0=w1=w2=w3=w4=w5=tVect(0.0,0.0,0.0);
        xp0=xp1=xp2=xp3=xp4=xp5=tVect(0.0,0.0,0.0);
        qp0=tQuat(1.0,0.0,0.0,0.0);
        qp1=qp2=qp3=qp4=qp5=tQuat(0.0,0.0,0.0,0.0);
        wp0=wp1=wp2=wp3=wp4=wp5=tVect(0.0,0.0,0.0);
        index=0;
        m=1.0;
        I=tVect(1.0,1.0,1.0);
        FHydro.reset();
        FWall.reset();
        FParticle.reset();
        FSpringP.reset();
        FSpringW.reset();
        MHydro.reset();
        MParticle.reset();
        MWall.reset();
        MRolling.reset();
        ACoriolis.reset();
        ACentrifugal.reset();
        components.resize(size);
        fluidVolume=0.0;
        maxOverlap=0.0;
    }
    void elmtShow()const;
    void initialize(const double& partDensity, std::vector <vecList>& prototypes, tVect& demF);
    void resetVelocity();
    void generateParticles(unsigned int& globalIndex, particleList& particles,  const std::vector <vecList>& prototypes);
    void predict(const double c1[],const double c2[]);
    void correct(const double coeff1ord[],const double coeff2ord[]);
    void translate(const tVect& transVec);
};

enum ContactModel {LINEAR, HERTZIAN};
   
class material{
public:
    // density [ mass / length² ]
    double density;

    // contact model
    ContactModel contactModel;

    // linear model ////////////////////
    // stiffness
    double linearStiff;

    // Hertzian model /////////////////////////
    // Young modulus [ force/length² ]
    double youngMod;
    // Poisson ratio [ / ]
    double poisson;
    // constant part of normal stiffness (dependent variable)
    double knConst;
    // constant part of shear stiffness (dependent variable) ->this is not physically sound
    double ksConst;

    // normal damping ///////////////////////////
    // normal viscous coefficient
    double restitution;
    double dampCoeff;

    // tangential model ////////////////////////
    // tangential viscous coefficient
    double viscTang;
    // particle-particle friction
    double frictionCoefPart;
    // particle-wall friction
    double frictionCoefWall;
    // particle-object friction
    double frictionCoefObj;
    // rolling model ////////////////////////  
    // particle-particle rolling friction
    double rollingCoefPart;

    // default constructor
    material(){
        contactModel=LINEAR;
        density=1.0;
        youngMod=1.0;
        poisson=0.3;
        knConst=1.0;
        ksConst=1.0;
        linearStiff=1.0;
        restitution=0.5;
        dampCoeff=0.5;
        viscTang=0.5;
        frictionCoefPart=0.3;
        frictionCoefObj=0.3;
        frictionCoefWall=0.3;
        rollingCoefPart=0.02;
    }
};

class ghost{
public:
    unsigned int ghostIndex;
    tVect pbcVector;
    ghost() {
        ghostIndex=0;
        pbcVector=tVect(0.0,0.0,0.0);
    }
};

class object{
public:
    // object index
    unsigned int index;
    // original object index in case of copies
    unsigned int originalIndex;
    // element index
    unsigned int ElID;
    // particle radius
    double r;
    // position of the object
    tVect x0;
    // velocity of the object
    tVect x1;
    // force on the object
    tVect FParticle,FHydro;
    // force on the object, maximum over simulation time, and time of occurrence
    tVect maxFParticle;
    double timeMaxFParticle;
    // force on the object, saved when necessary
    tVect savedFParticle;
//     is it translating?
    bool translating;
//     translation vector
    tVect trans;
    object() {
        index=ElID=0;
        originalIndex=0;
        r=0.0;
        x0=Zero;
        x1=Zero;
        FParticle=FHydro=Zero;
        maxFParticle=Zero;
        savedFParticle=Zero;
        timeMaxFParticle=0.0;
        translating = true;
//        trans = tVect(0.0,0.0,0.0);
    }
    void updateMax(const tVect& direction, const double& time);
    void saveForces();
};

class Elongation{
public:
    tVect e; //elongation vector
    tVect p; //contact point between particles
    bool slipping;
    unsigned int slippingCase;
    unsigned int indexI;
    unsigned int indexJ;
    bool active;
    
    Elongation() {
        e=tVect (0.0,0.0,0.0);
        p=tVect (0.0,0.0,0.0);
        slipping=false;
        slippingCase=0;
        active=false;
    }
    
    void reset();
    void copy (const Elongation& e2);
};

#endif	/* ELMT_H */


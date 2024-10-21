/*
 * File:   node.h
 * Author: aleonard
 *
 * Created on November 5, 2013, 4:10 PM
 */

#ifndef NODE_H
#define	NODE_H

#include <stdlib.h>

#include "lattice.h"
#include "myvector.h"

enum Rheology {NEWTONIAN, BINGHAM, FRICTIONAL, BAGNOLD, MUI, VOELLMY};

class FluidMaterial{
public:
    // maximum and minimum viscosity due to limitations to tau
    double lbMinVisc;
    double lbMaxVisc;
    // initial density
    double initDensity;
    // initial viscosity
    double initDynVisc;
    // rheology model
    Rheology rheologyModel;
    // parameters for Bingham fluid (rheologyModel=BINGHAM)
    double yieldStress, plasticVisc;    
    // parameter for friction (rheologyModel=FRICTIONAL,VOELLMY,MUI)
    double frictionCoefFluid, basalFrictionCoefFluid;
    // range of frictional values (rheologyModel=MUI)
    double deltaFriction;
    // basic inertial number (rheologyModel=MUI)
    double baseInertial;
    // particle diameter for computing inertial number or bagnold viscosity (rheologyModel=MUI, BAGNOLD)
    double particleDiameter;
    // square of particle diameter for computing bagnold viscosity (rheologyModel=BAGNOLD)
    double rhod2;
    // particle density for computing inertial number (rheologyModel=MUI)
    double particleDensity;
    // minimum pressure for stabilizing frictional models (rheologyModel=MUI,FRICTIONAL)
    double minimumPressure;
    // activation of Smagorinsky turbulence model
    bool turbulenceOn;
    // parameter for Smagorinsky turbulence model
    double turbConst;
    // earth pressure coefficient for static load on structures 
    double earthPressureCoeff;
            // default constructor
    FluidMaterial(){
        initDensity=1.0;
        initDynVisc=(1.0-0.5)/3.0;
        rheologyModel=NEWTONIAN;
        yieldStress=0.0;
        plasticVisc=0.0;
        frictionCoefFluid=0.0;
        basalFrictionCoefFluid=0.0;
        turbulenceOn=false;
        turbConst=0.0;
        particleDiameter=0.0;
        baseInertial=0.0;
        deltaFriction=0.0;
        earthPressureCoeff=1.0;
    }
};

class node{
    // lattice node declaration
public:
    // linear coordinate of node
    unsigned int coord;
    // probability density functions (one per lattice direction)
    double f[lbmDirec];
    // support variables for streaming
    double fs[lbmDirec];
    // density in the node
    double n;
    // smoothed pressure in the node
    double smoothedPressure;
    // velocity of fluid in node
    tVect u;
    // DEM-coupling forces on the node
    tVect hydroForce;
    // surface normal
    tVect surfaceNormal;
    // centrifugal force on the node
    tVect centrifugalForce;
    // mass functions
    double mass, newMass, extraMass;
    // viscosity functions
    double visc;
    bool basal;
    // friction (for frictional models)
    double friction;
    // for smooth surface transitions
    float age;
    // neighbor nodes
    node* d[lbmDirec];
    // informations for curved boundaries
    curve* curved;
    // identification of node type
    // 0 = liquid cells
    // 3 = interface
    // 4 = periodic
    // 5 = slip stationary planar walls
    // 6 = slip moving planar walls
    // 7 = no-slip stationary planar walls
    // 8 = no-slip moving planar walls
    // 9 = cylinder walls
    // 10 = object walls
    // 11 = topographic walls
    types type;
    // particle flag
    bool p;
    // solid index
    unsigned short int solidIndex;
    // default constructor
    node(){
        n=0.0;
        smoothedPressure=0.0;
        u.reset();
        hydroForce.reset();
        centrifugalForce.reset();
        newMass=mass=extraMass=0.0;
        visc=1.0;
        friction=0.0;
        age=0.0;
        f[0]=f[1]=f[2]=f[3]=f[4]=f[5]=f[6]=f[7]=f[8]=f[9]=0.0;
        f[10]=f[11]=f[12]=f[13]=f[14]=f[15]=f[16]=f[17]=f[18]=0.0;
//        fp[0]=fp[1]=fp[2]=fp[3]=fp[4]=fp[5]=fp[6]=fp[7]=fp[8]=fp[9]=0.0;
//        fp[10]=fp[11]=fp[12]=fp[13]=fp[14]=fp[15]=fp[16]=fp[17]=fp[18]=0.0;
//        fm[0]=fm[1]=fm[2]=fm[3]=fm[4]=fm[5]=fm[6]=fm[7]=fm[8]=fm[9]=0.0;
//        fm[10]=fm[11]=fm[12]=fm[13]=fm[14]=fm[15]=fm[16]=fm[17]=fm[18]=0.0;
        fs[0]=fs[1]=fs[2]=fs[3]=fs[4]=fs[5]=fs[6]=fs[7]=fs[8]=fs[9]=0.0;
        fs[10]=fs[11]=fs[12]=fs[13]=fs[14]=fs[15]=fs[16]=fs[17]=fs[18]=0.0;
        d[0]=d[1]=d[2]=d[3]=d[4]=d[5]=d[6]=d[7]=d[8]=d[9]=0;
        d[10]=d[11]=d[12]=d[13]=d[14]=d[15]=d[16]=d[17]=d[18]=0;
        curved=0;
        type=LIQUID;
        solidIndex=0;
        basal=false;
    }
    // getting liquid fraction
    double liquidFraction() const;
    // shows the population in case debugging requires it
    void showPopulations() const;
    // functions working on distribution functions
    void initialize(const double& density, const tVect& velocity, const double& massFunction, const double& viscosity, const tVect& F, const double& ageHere, const tVect& rotationSpeed);
    void copy(const node*& copyNode);
    void scatterMass(double& extraMass);
    void restart(const double& restDensity, const tVect& restVelocity, const double& restMassFunction, const double& restViscosity, const double restF[]) ;
    void setEquilibrium(const tVect& F, const double& nHere, const tVect& velHere);
    void reconstruct();
    void shiftVelocity(const tVect& F);
    void computeEquilibrium(double feq[]);
    void computeEquilibriumTRT(double feqp[], double feqm[]);
    void computeApparentViscosity(const double feq[], const FluidMaterial& fluidMaterial);
    void resetViscosity(const double& initViscosity);
    void solveCollision(const double feq[]);
    void solveCollisionTRT(const double feqp[], const double feqm[], const double& magicNumber);
    void addForce(const tVect& F);
    void addForce(double feq[], const tVect& F);
    void addForceTRT(const tVect& F);
    void collideFNN(const double& turbConst, const tVect& F, const double& plasticVisc, const double& yieldStress);
    void collideNN(const double& turbConst, const double& plasticVisc, const double& yieldStress);
    void collideF(const tVect& F, const double& initVisc);
    void collide(const double& initVisc);
    void store();
    tVect bounceBackForce(const unsigned int& j, const double staticPres[], const double& BBi, const double& earthPressureCoeff, const double& maxVisc) const;
    double massStream(const unsigned int& sdir) const;
    // type
            // functions for identification of type
    bool isActive() const;
    bool isWall() const;
    bool isCurvedWall() const;
    inline bool isSlipStatWall() const { return type == SLIP_STAT_WALL; }
    inline bool isSlipDynWall() const { return type == SLIP_DYN_WALL; }
    inline bool isStatWall() const { return type == STAT_WALL; }
    inline bool isDynWall() const { return type == DYN_WALL; }
    inline bool isCylinderWall() const { return type == CYL; }
    inline bool isObjectWall() const { return type == OBJ; }
    inline bool isTopography() const { return type == TOPO; }
    inline bool isOutlet() const { return type == OUTLET; }
    inline bool isInterface() const { return type == INTERFACE; }
    inline bool isFluid() const { return type == LIQUID; }
    inline bool isGas() const { return type == GAS; }
    inline bool isPeriodic() const { return type == PERIODIC; }
    // functions for change of type
    void setSlipStatWall();
    void setSlipDynWall();
    void setStatWall();
    void setDynWall();
    void setCylinderWall();
    void setObjectWall();
    void setTopography();
    void setOutlet();
    void setInterface();
    void setFluid();
    void setGas();
    void setPeriodic();
    void setType(const types& typ);
    //unsigned int getType() const;
    // particle flag
    inline bool isInsideParticle() const { return p; };
    void setInsideParticle();
    void setOutsideParticle();
    // functions to get and assign index
    unsigned short int getSolidIndex() const;
    void setSolidIndex(const unsigned short int& ind);
};

class curve{
public:
    tVect wallNormal;
    double delta[lbmDirec];
    double m1[lbmDirec], m2[lbmDirec], m3[lbmDirec];
    curve(){
        for (int j=1; j<lbmDirec; ++j) {
            delta[j]=0.0;
            m1[j]=0.0;
            m2[j]=0.0;
            m3[j]=0.0;
        }
    }
    void computeCoefficients();
    double getChi(const unsigned int& j, const double& tau) const;
};

class measureUnits{
public:
    // conversion units /////////////////////////////////////////////////////////////////////////
    // fluid and granular matter are solved in different measure units
    // for a reference, check Feng, Han, Owen, 2007
    // unit length
    double Length;
    //unit time
    double Time;
    //unit density
    double Density;
    // secondary units: speed, acceleration, mass, force, flow rate
    double Volume, Speed, Accel, AngVel, Force, Torque, Mass, KinVisc, DynVisc, Stress, Pressure, FlowRate, Energy;
    // inverted unit length for efficiency
    double invLength;
    // constructor
    measureUnits(){
        Volume=1.0;
        Length=1.0;
        Time=1.0;
        Density=1.0;
        Speed=1.0;
        Accel=1.0;
        AngVel=1.0;
        Force=1.0;
        Torque=1.0;
        Mass=1.0;
        KinVisc=1.0;
        DynVisc=1.0;
        Stress=1.0;
        Pressure=1.0;
    }
    void setComposite();
};


types boundary2Type(const string& typeString);

//types int2Type(const unsigned int& typeID);

#endif /* NODE_H */


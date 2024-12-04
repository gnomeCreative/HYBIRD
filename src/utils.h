/* 
 * File:   utils.h
 * Author: aleonard
 *
 * Created on November 5, 2013, 4:11 PM
 */

#ifndef UTILS_H
#define	UTILS_H

#include "myvector.h"
#include "macros.h"
#include <limits>
#include <algorithm>

class elmt;

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// CLASS  DEFINITIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/

class wall {
public:
    // index
    unsigned int index;
    // normal vector
    tVect n;
    // basic point
    tVect p;
    // center of rotation
    tVect rotCenter;
    // rotational speed
    tVect omega;
    // translational speed
    tVect vel;
    // hydraulic force on wall
    tVect FHydro;
    // collision force on wall
    tVect FParticle;
    // max hydraulic force on wall
    tVect maxFHydro;
    // max collision force on wall
    tVect maxFParticle;
    // is it a moving wall? (true = moving; false = fixed);
    bool moving;
    // is it a slipping wall?
    bool slip;
    // is it translating?
    bool translating;
    // translation vector
    tVect trans;
    // limits for walls
    bool limited;
    double xMin,xMax,yMin,yMax,zMin,zMax;
    // default constructor
    wall(){
        n=tVect(1.0,0.0,0.0);
        p=tVect(0.0,0.0,0.0);
        rotCenter=tVect(0.0,0.0,0.0);
        omega=tVect(0.0,0.0,0.0);
        vel=tVect(0.0,0.0,0.0);
        FHydro=tVect(0.0,0.0,0.0);
        FParticle=tVect(0.0,0.0,0.0);
        maxFHydro=tVect(0.0,0.0,0.0);
        maxFParticle=tVect(0.0,0.0,0.0);
        moving=false;
        slip=false;
        translating=false;
        trans=tVect(0.0,0.0,0.0);
        limited=false;
        xMin=yMin=zMin=-numeric_limits<double>::max();
        xMax=yMax=zMax=numeric_limits<double>::max();
        index=0;
    }
    // better constructor
    wall(tVect ip, tVect in){
        n=in;
        p=ip;
        rotCenter=tVect(0.0,0.0,0.0);
        omega=tVect(0.0,0.0,0.0);
        vel=tVect(0.0,0.0,0.0);
        moving=false;
        slip=false;
        translating=false;
        trans=tVect(0.0,0.0,0.0);
        limited=false;
        xMin=yMin=zMin=0.0;
        xMax=yMax=zMax=HUGE_VAL;
        index=0;
    }
    // distance point-wall
    double dist(tVect x) const;
    // show wall characteristics
    void wallShow() const;
    // computes the speed of a point of the wall
    tVect getSpeed(tVect pt) const;
};

class topography {
public:
    // index
    unsigned int sizeX = 0, sizeY = 0;
    // coordinate vectors
    doubleList coordX = {};
    doubleList coordY = {};
    vecSet points = {};
    vecSet surface = {};
    double spacingX = 0,spacingY = 0;
    double corner1X = 0,corner1Y = 0;
    double corner2X = 0,corner2Y = 0;
    constexpr topography() = default;
    void readFromFile(string& fileName, double& transX, double& transY, double& transZ);
    void show();
    // computes the distance sfrom the topographical surface
    void getReferenceTriangle(const tVect& point, tVect& point1, tVect& point2, tVect& point3, tVect& planeNormal)  const;
    void getSurfaceTriangle(const tVect& point, tVect& point1, tVect& point2, tVect& point3, tVect& planeNormal)  const;
    double distance(const tVect& point) const;
    tVect surfaceNormal(const tVect& point) const;
    double surfaceIsoparameter(const tVect& point) const;
    double directionalDistance(const tVect& point, const tVect& dir) const;
};

enum cylinderType {EMPTY, FULL};

class cylinder {
    // Cylinder class, defined by two points and a radius. Cylinder is of infinite length
public:
    //index
    unsigned int index;
    // two axes point
    tVect p1,p2;
    // radius
    double R;
    // axes vector and unit vector (defined in init_axes for speed)
    tVect axes, naxes;
    // rotational speed
    tVect omega;
    // is it a rotating wall (true = rotating; false = fixed);
    bool moving;
    // is it a slip cylinder?
    bool slip;
    // inner or outer
    cylinderType type;
    // function for distance calculation (norm)
    double dist(tVect pt) const;
    // function for distance calculation (scalar)
    tVect vecDist(tVect pt) const;
    // limits for cylinders
    bool limited;
    double xMin,xMax,yMin,yMax,zMin,zMax;
    // is it translating?
    bool translating;
    // translation vector
    tVect trans;
    
    // default constructor
    cylinder(){
        index=0;
        p1 = tVect(0.0,0.0,0.0);
        p2 = tVect(1.0,0.0,0.0);
        R=1.0;
        axes =p2-p1;
        naxes=axes/axes.norm();
        omega=tVect(0.0,0.0,0.0);
        moving=false;
        slip=false;
        type=EMPTY;
        translating=false;
        trans=tVect(0.0,0.0,0.0);
        }
    // constructor
    cylinder(tVect ip1, tVect ip2, double radius){
        p1=ip1;
        p2=ip2;
        R=radius;
        axes =p2-p1;
        naxes=axes/axes.norm();
        omega=tVect(0.0,0.0,0.0);
        moving=false;
        slip=false;
    }
    // get tangential speed of a point rotating with the cylinder
    tVect getSpeed(tVect pt) const;
    // get angle and distance from axes in the cylindrical coordinates defined by the cylinder
    double getTheta(const tVect& pt, const tVect& zeroDir, const tVect& thetaDir) const;
    double getR(const tVect& pt) const;
    // axes variables initialization
    void initAxes();
    void cylinderShow() const;
    double segmentIntercept(const tVect& start, const tVect& dir) const;
};

class pbc {
 // periodic boundary class, defined by a point and a distance vector between the two planes
public:
    // index
    int index;
    // point of a plane
    tVect p;
    // distance between two planes
    tVect v;
    // two periodic planes (to be defined)
    wall pl1, pl2;
    void pbcShow() const;
    //function for the definition of two planes
    void setPlanes();

    pbc(){
        v=tVect(1.0,0.0,0.0);
        p=tVect(0.0,0.0,0.0);
    }
};

class energy {
public:
    // total mass
    double mass;
    // total kinetic energy
    double kin;
    // rotational kinetic energy
    double rotKin;
    // translational kinetic energy
    double trKin;
    // gravitational potential energy
    double grav;
    // elastic potential energy
    double elastic;
    // total energy
    double total;
    energy(){
        kin=0.0;
        rotKin=0.0;
        trKin=0.0;
        grav=0.0;
        elastic=0.0;
        total=0.0;
    }
    void reset();
//    void updateKinetic(elmtList& elmts);
//    void updatePotential(elmtList& elmts, tVect& gravity, doubleList& dim);
//    void updateElastic();
    void updateTotal();
    void show();
};


/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// TYPE DEFINITIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/





#endif	/* UTILS_H */


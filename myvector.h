/* 
 * File:   vector.h
 * Author: aleonard
 *
 * Created on November 5, 2013, 4:11 PM
 */

//#define DEBUG
#undef DEBUG

#ifndef VECTOR_H
#define	VECTOR_H

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>

using namespace std;
enum ProblemName {NONE, SHEARCELL, AVALANCHE, DRUM, NET, BARRIER, 
    ZHOU, OPENBARRIER, HONGKONG, STVINCENT, STAVA, NIGRO, CAROLINE, DAMBREAK, GRAY_DAMBREAK, GRAY_DAMBREAK_2D, INCLINEFLOW,
    HOURGLASS, IERVOLINO,IERVOLINO_2D, IERVOLINO_CYLINDERTEST, HEAP, TRIAXIAL, JOP, WILL, WILL_SETTLING, MANGENEY, GRAY, ESERCITAZIONE, FILIPPO_SILOS, HK_SMALL, HK_LARGE, KELVIN, SHEARCELL2022};
    
enum types : unsigned char {LIQUID, UNUSED, GAS, INTERFACE, PERIODIC, SLIP_STAT_WALL, SLIP_DYN_WALL, STAT_WALL, DYN_WALL, CYL, OBJ, TOPO, OUTLET};

inline const char* typeString(types v)
{
    switch (v)
    {
        case LIQUID:   return "LIQUID";
        case UNUSED:   return "UNUSED";
        case INTERFACE:      return "INTERFACE";
        case PERIODIC:      return "PERIODIC";
        case SLIP_STAT_WALL:      return "SLIP_STAT_WALL";
        case SLIP_DYN_WALL:      return "SLIP_DYN_WALL";
        case STAT_WALL:      return "STAT_WALL";
        case DYN_WALL:      return "DYN_WALL";
        case CYL:      return "CYL";
        case OBJ:      return "OBJ";
        case TOPO:      return "TOPO";
        case OUTLET:      return "OUTLET";

        default:       return "[Unknown type]";
    }
}


/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// BASIC TYPE DEFINITIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/

// lists of  basic variables
typedef std::vector<unsigned int> unsIntList;
typedef std::vector<int> intList;
typedef std::vector<double> doubleList;
// double lists of basic variables
typedef std::vector< std::vector<unsigned int> > unsIntSet;
typedef std::vector< std::vector<int> > intSet;
typedef std::vector< std::vector<double> > doubleSet;
typedef std::vector< types > typeList;

// vector variable type
// defined with 3 double numbers
//typedef class tinyVector tVect;
// quaternion variable type
// defined with 4 double numbers (should be normalized, ||q||=1)
typedef class quaternion tQuat;
// matrix variable type
// defined with 9 double numbers
typedef class tinyMatrix tMat;

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// CLASS  DEFINITIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/

// VECTOR CLASS
class tVect {
public:
    double x,y,z;
public:
    tVect(){
        x=0.0;
        y=0.0;
        z=0.0;
    }
    //tinyVector(double a, double b, double c){
    tVect(const double& a, const double& b, const double& c){
        x=a;
        y=b;
        z=c;
    }
    void get(std::ifstream& inputFile);
    void show() const;
    void printLine(std::ofstream& outputFile) const;
    void print(std::ofstream& outputFile) const;
    void printFixedLine(std::ofstream& outputFile) const;
    void printFixed(std::ofstream& outputFile) const;
    double max() const;
    void reset();
    // overloaded operators
    tVect operator+(const tVect& vec) const;
    tVect& operator+=(const tVect& vec);
    tVect operator-(const tVect& vec) const;
    tVect& operator-=(const tVect& vec);
    tVect operator*(const double& scalar) const;
    tVect& operator*=(const double& scalar);
    tVect operator/(const double& scalar) const;
    tVect& operator/=(const double& scalar);
    //friend tVect operator*(const double& scalar, const tVect& vec);
    // mathematical operations
    tVect transport() const;
    tVect abs() const;
    double dot(const tVect& vec) const;
    double dot2(const tVect& vec) const;
    tVect cross(const tVect& vec) const;
    tVect compProd(const tVect& vec) const;
    tMat outer(const tVect& vec) const;
    double norm() const;
    double norm2() const;
    int linearizePosition(double cellWidth[], unsigned int nCells[]) const;
    //friend tVect newtonAcc(tVect moment, tVect I, tVect wBf);
    //friend tQuat vec2quat(tVect vec);
    //friend tVect project(tVect& vec, tQuat quat);
    // geometric position functions
    bool insideSphere(const tVect& center, const double& radius) const;
    bool insideCylinder(const tVect& p1, const tVect& naxes, const double& R1, const double& R2) const;
    bool insidePlane(const tVect& p, const tVect& n) const;
    double distance2Plane(const tVect& p, const tVect& n) const;
    bool close2Plane(const tVect& p, const tVect& n, const double& distance) const;
    double isoparameterSphere(const tVect& center, const double& radius) const;
};


// direction vectors
#define Xp tVect(1.0,0.0,0.0)
#define Yp tVect(0.0,1.0,0.0)
#define Zp tVect(0.0,0.0,1.0)
#define Xm tVect(-1.0,0.0,0.0)
#define Ym tVect(0.0,-1.0,0.0)
#define Zm tVect(0.0,0.0,-1.0)
const tVect Zero=tVect(0.0,0.0,0.0);

// QUATERNION CLASS
class quaternion{
public:
    double q0,q1,q2,q3;
public:
    quaternion(){
        q0=1.0;
        q1=q2=q3=0.0;
    }
    quaternion(double a, double b, double c, double d){
        q0=a;
        q1=b;
        q2=c;
        q3=d;
    }
    void show() const;
    void printLine(std::ofstream& outputFile) const;
    void print(std::ofstream& outputFile) const;
    // overloaded operators
    tQuat operator+(const tQuat& quat) const;
    tQuat operator-(const tQuat& quat) const;
    tQuat operator*(const double& scalar) const;
    tQuat operator/(const double& scalar) const;
    friend tQuat operator*(const double& scalar, const tQuat& quat);
    // mathematical operations
    void normalize();
    void forceStability(tQuat& q0);
    void resetHard();
    void resetSoft();
    tQuat adjoint() const;
    double norm() const;
    double norm2() const;
    tQuat inverse() const;
    tQuat multiply(tQuat r) const;
    double dot(tQuat r) const;
    //friend tVect quat2vec(tQuat quat);
    //friend tVect project(tVect vec, tQuat quat);
    //friend tQuat quatAcc(tVect waBf, tQuat Q1);
};

// MATRIX CLASS
class tinyMatrix {
private:
    double m00,m01,m02,m10,m11,m12,m20,m21,m22;
public:
    tinyMatrix(){
        m00=m01=m02=0.0;
        m10=m11=m12=0.0;
        m20=m21=m22=0.0;
    }
    tinyMatrix(double a00, double a01, double a02,
                         double a10, double a11, double a12,
                         double a20, double a21, double a22){
        m00=a00;
        m01=a01;
        m02=a02;
        m10=a10;
        m11=a11;
        m12=a12;
        m20=a20;
        m21=a21;
        m22=a22;
    }
    tinyMatrix(tVect v1, tVect v2){
        tMat matrix=v1.outer(v2);
        m00=matrix.m00;
        m01=matrix.m01;
        m02=matrix.m02;
        m10=matrix.m10;
        m11=matrix.m11;
        m12=matrix.m12;
        m20=matrix.m20;
        m21=matrix.m21;
        m22=matrix.m22;
    }
    void show() const;
    // overloaded operators
    tMat operator+(const tMat& mat) const;
    tMat& operator+=(const tMat& mat);
    tMat operator-(const tMat& mat) const;
    tMat& operator-=(const tMat& mat);
    tMat operator*(const double& scalar) const;
    tMat& operator*=(const double& scalar);
    tMat operator/(const double& scalar) const;
    tMat& operator/=(const double& scalar);
    friend tMat operator*(const double& scalar, const tMat& mat);
    // mathematical operations
    double magnitude() const;
};

// INTER-CLASS FUNCTIONS
// overloaded operators
tVect operator*(const double& scalar, const tVect& vec);
tQuat operator*(const double& scalar, const tQuat& quat);
tMat operator*(const double& scalar, const tMat& mat);
// mathematical operations
tVect quat2vec(tQuat quat);
tQuat vec2quat(tVect vec);
tVect project(tVect vec, tQuat quat);
tVect newtonAcc(tVect moment, tVect I, tVect wBf);
tQuat quatAcc(tVect waBf, tQuat Q1);

tVect computeCentrifugal(const tVect& position, const tVect& rotationCenter, const tVect& rotationSpeed);
tVect computeCoriolis(const tVect& velocity, const tVect& rotationSpeed);

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// COMPOSITE TYPE DEFINITIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/

// lists of vectors and list of lists of vectors
typedef std::vector<tVect> vecList;
typedef std::vector< std::vector<tVect> > vecSet;

// list of matrices
typedef std::vector<tMat> matList;

// list of particles, elements and ghosts
class particle;
typedef std::vector<particle> particleList;
typedef std::vector<particle*> particlePointerList;
class elmt;
typedef std::vector<elmt> elmtList;
class ghost;
typedef std::vector<ghost> ghostList;
class object;
typedef std::vector<object> objectList;
class Elongation;
typedef std::vector <Elongation> elongVector;

// list of nodes
class node;
typedef std::vector<node*> nodeList;
typedef std::map<const unsigned int,node> nodeMap;

// list of curved
class curve;
typedef std::vector<curve*> curveList;

// list of geometric entities
class wall;
typedef std::vector<wall> wallList;
class cylinder;
typedef std::vector<cylinder> cylinderList;
class pbc;
typedef std::vector<pbc> pbcList;

#endif	/* VECTOR_H */


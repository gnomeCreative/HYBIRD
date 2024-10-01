

#include "myvector.h"

using namespace std;

/*/////////////////////////////////////////////////////////////////////////////////////////
// VECTOR
/////////////////////////////////////////////////////////////////////////////////////////*/

// basic functions

inline void tVect::show() const {
    cout<<"("<<x<<", "<<y<<", "<<z<<")";
}

inline void tVect::printLine(std::ofstream& outputFile) const {
    outputFile<<x<<"\t"<<y<<"\t"<<z<<"\n";
}

inline void tVect::print(std::ofstream& outputFile) const {
    outputFile<<x<<"\t"<<y<<"\t"<<z<<"\t";
}

inline void tVect::printFixedLine(std::ofstream& outputFile) const {
    outputFile<<std::setprecision(8)<<std::scientific<<x<<" "<<y<<" "<<z<<"\n";
}

inline void tVect::printFixed(std::ofstream& outputFile) const {
    outputFile<<std::setprecision(8)<<std::fixed<<x<<" "<<y<<" "<<z<<" ";
}


inline void tVect::get(std::ifstream& inputFile) {
        inputFile>>x;
        inputFile>>y;
        inputFile>>z;
    }

inline double tVect::max() const {
    return std::max(std::abs(x),std::max(std::abs(y),std::abs(z)));
}

inline void tVect::reset() {
    x=0.0;
    y=0.0;
    z=0.0;
}

// overloaded operators

inline tVect tVect::operator+(const tVect& vec) const {
    return tVect(
        x+vec.x,
        y+vec.y,
        z+vec.z);
}

inline tVect& tVect::operator+=(const tVect& vec) {
    x+=vec.x;
    y+=vec.y;
    z+=vec.z;
    return *this;
}

inline tVect tVect::operator-(const tVect& vec) const {
    return tVect(
        x-vec.x,
        y-vec.y,
        z-vec.z);
}

inline tVect& tVect::operator-=(const tVect& vec) {
    x-=vec.x;
    y-=vec.y;
    z-=vec.z;
    return *this;
}

inline tVect tVect::operator*(const double& scalar) const {
    return tVect(
        x*scalar,
        y*scalar,
        z*scalar);
}

inline tVect& tVect::operator*=(const double& scalar) {
    x*=scalar;
    y*=scalar;
    z*=scalar;
    return *this;
}

inline tVect tVect::operator/(const double& scalar) const {
    return tVect(
        x/scalar,
        y/scalar,
        z/scalar);
}

inline tVect& tVect::operator/=(const double& scalar) {
    x/=scalar;
    y/=scalar;
    z/=scalar;
    return *this;
}

inline tVect operator *(const double& scalar, const tVect& vec) {
    return tVect(
        vec.x*scalar,
        vec.y*scalar,
        vec.z*scalar);
}

// mathematical functions

inline tVect tVect::transport() const {
    return tVect(y*y+z*z,z*z+x*x,x*x+y*y);
}

inline tVect tVect::abs() const {
    return tVect(std::abs(x),std::abs(y),std::abs(z));
}

inline double tVect::dot(const tVect& vec) const {
    return x*vec.x+y*vec.y+z*vec.z;
}

inline double tVect::dot2(const tVect& vec) const {
    return pow(this->dot(vec),2.0);
}

inline tVect tVect::cross(const tVect& vec) const {
    return tVect(
            y*vec.z-z*vec.y,
            z*vec.x-x*vec.z,
            x*vec.y-y*vec.x);
}

inline tVect tVect::compProd(const tVect& vec) const {
    return tVect(
            x*vec.x,
            y*vec.y,
            z*vec.z);
}

inline tMat tVect::outer(const tVect& vec) const {
    return tMat(x*vec.x,x*vec.y,x*vec.z,y*vec.x,y*vec.y,y*vec.z,z*vec.x,z*vec.y,z*vec.z);
}

inline double tVect::norm() const {
    return sqrt(x*x+y*y+z*z);
}

inline double tVect::norm2() const {
    return x*x+y*y+z*z;
}

inline int tVect::linearizePosition(double cellWidth[], unsigned int nCells[]) const  {
    const int xc=floor(x/cellWidth[0])+1;
    const int yc=floor(y/cellWidth[1])+1;
    const int zc=floor(z/cellWidth[2])+1;
    const int index=xc+nCells[0]*(yc+ nCells[1]*zc );
    return index;
    cout<<"x="<<x<<endl;
    cout<<"y="<<y<<endl;
    cout<<"z="<<z<<endl;
    cout<<"floor(x/cellWidth[0])+1="<<floor(x/cellWidth[0])+1<<endl;
    cout<<"floor(y/cellWidth[1])+1="<<floor(y/cellWidth[1])+1<<endl;
    cout<<"floor(z/cellWidth[2]))+1="<<floor(z/cellWidth[2])+1<<endl;
    cout<<"index="<<index<<endl;
    

}

// geometric position functions

inline bool tVect::insideSphere(const tVect& center, const double& radius) const {
    //distance between center of the sphere and point
    tVect dist=*this -center;
//    return ((z-center.z)*(z-center.z)+(y-center.y)*(y-center.y)+(x-center.x)*(x-center.x)<radius*radius);
    return (dist.norm2()<radius*radius);
}

inline bool tVect::insideCylinder(const tVect& p1, const tVect& naxes, const double& R1, const double& R2) const {
    // distance to point 1 of axes
    const tVect p1dist=*this-p1;
    // same but projected on the axes
    const tVect p1distax=(p1dist.dot(naxes))*naxes;
    // distance center to cylinder
    const tVect p1distcylinder=p1dist-p1distax;
    // square of distance to center
    const double p1distcylinderSquare=p1distcylinder.norm2();
    // condition for being inside
    return (p1distcylinderSquare>R1*R1 && p1distcylinderSquare<R2*R2);
}



inline bool tVect::insidePlane(const tVect& p, const tVect& n) const {
    tVect pdist=*this-p;
    return n.dot(pdist)<0.0;
}

inline bool tVect::close2Plane(const tVect& p, const tVect& n, const double& distance) const {
    const double distanceHere=this->distance2Plane(p,n);
    return distanceHere<distance;
}

inline double tVect::distance2Plane(const tVect& p, const tVect& n) const {
    tVect pdist=*this-p;
    return n.dot(pdist);
}

inline double tVect::isoparameterSphere(const tVect& center, const double& radius) const {
    //distance between center of the sphere and point
    tVect dist=*this -center;
//    return ((z-center.z)*(z-center.z)+(y-center.y)*(y-center.y)+(x-center.x)*(x-center.x)<radius*radius);
    return (radius*radius)/dist.norm2();
}

/*/////////////////////////////////////////////////////////////////////////////////////////
// QUATERNION
/////////////////////////////////////////////////////////////////////////////////////////*/

// basic functions

inline void quaternion::show() const {
    cout<<" ("<<q0<<", "<<q1<<", "<<q2<<", "<<q3<<")";
}

inline void quaternion::printLine(std::ofstream& outputFile) const {
    cout.precision(17);
    outputFile<<q0<<"\t"<<q1<<"\t"<<q2<<"\t"<<q3<<"\n";
}

inline void quaternion::print(std::ofstream& outputFile) const {
    cout.precision(17);
    outputFile<<q0<<"\t"<<q1<<"\t"<<q2<<"\t"<<q3<<"\t";
}

inline void quaternion::forceStability(tQuat& quat){
    const static double relax=0.0;
    q0=q0;//-1.0*this->dot(quat)/quat.q0;
    q1=q1;//-relax*this->dot(quat)/quat.q1;
    q2=q2;//-relax*this->dot(quat)/quat.q2;
    q3=q3;//-relax*this->dot(quat)/quat.q3;

}

inline void quaternion::resetHard() {
    q0=0.0;
    q1=0.0;
    q2=0.0;
    q3=0.0;
}

inline void quaternion::resetSoft() {
    q0=1.0;
    q1=0.0;
    q2=0.0;
    q3=0.0;
}

// overloaded operators

inline tQuat quaternion::operator+(const tQuat& quat) const {
    return tQuat(
        q0+quat.q0,
        q1+quat.q1,
        q2+quat.q2,
        q3+quat.q3);
}

inline tQuat quaternion::operator-(const tQuat& quat) const {
    return tQuat(
        q0-quat.q0,
        q1-quat.q1,
        q2-quat.q2,
        q3-quat.q3);
}

inline tQuat quaternion::operator*(const double& scalar) const {
    return tQuat(
        q0*scalar,
        q1*scalar,
        q2*scalar,
        q3*scalar);
}

inline tQuat quaternion::operator/(const double& scalar) const {
    return tQuat(
        q0/scalar,
        q1/scalar,
        q2/scalar,
        q3/scalar);
}

inline tQuat operator *(const double& scalar, const tQuat& quat){
    return tQuat(
        quat.q0*scalar,
        quat.q1*scalar,
        quat.q2*scalar,
        quat.q3*scalar);
}

// mathematical functions

inline void quaternion::normalize(){
    double norm=this->norm();
    q0/=norm;
    q1/=norm;
    q2/=norm;
    q3/=norm;
}

inline double quaternion::norm() const {
    return sqrt(q0*q0+q1*q1+q2*q2+q3*q3);
}

inline double quaternion::norm2() const {
    return q0*q0+q1*q1+q2*q2+q3*q3;
}

inline tQuat quaternion::adjoint() const {
    return tQuat(
        q0,
        -q1,
        -q2,
        -q3);
}

inline tQuat quaternion::inverse() const {
    // hopefully this is not needed, since for a unit quaternion inverse and adjoint are equal
    tQuat dummyQ;
    dummyQ=dummyQ.adjoint();
    dummyQ.normalize();
    return dummyQ;
}

inline tQuat quaternion::multiply(quaternion quat) const {
    // note that quaternion multiplication is non-commutative
    // here the product q*r is solved, which in general is different from r*q
    // refer to the internet (mathworks page on quaternion multiplication) for documentation
    return tQuat(
        quat.q0*q0-quat.q1*q1-quat.q2*q2-quat.q3*q3,
        quat.q0*q1+quat.q1*q0-quat.q2*q3+quat.q3*q2,
        quat.q0*q2+quat.q1*q3+quat.q2*q0-quat.q3*q1,
        quat.q0*q3-quat.q1*q2+quat.q2*q1+quat.q3*q0);
}

inline double quaternion::dot(quaternion quat) const {
    // note that quaternion multiplication is non-commutative
    // here the product q*r is solved, which in general is different from r*q
    // refer to the internet (mathworks page on quaternion multiplication) for documentation
    return quat.q0*q0+quat.q1*q1+quat.q2*q2+quat.q3*q3;
}

/*/////////////////////////////////////////////////////////////////////////////////////////
// MATRIX
/////////////////////////////////////////////////////////////////////////////////////////*/

// basic functions

inline void tinyMatrix::show() const {
    cout<<" (";
    cout<<m00<<", "<<m01<<", "<<m02<<"; ";
    cout<<m10<<", "<<m11<<", "<<m12<<"; ";
    cout<<m20<<", "<<m21<<", "<<m22<<") ";
}

// overloaded operators

inline tMat tinyMatrix::operator+(const tMat& mat) const {
    return tMat(
        m00+mat.m00,
        m01+mat.m01,
        m02+mat.m02,
        m10+mat.m10,
        m11+mat.m11,
        m12+mat.m12,
        m20+mat.m20,
        m21+mat.m21,
        m22+mat.m22);
}

inline tMat& tinyMatrix::operator+=(const tMat& mat){
    m00+=mat.m00;
    m01+=mat.m01;
    m02+=mat.m02;
    m10+=mat.m10;
    m11+=mat.m11;
    m12+=mat.m12;
    m20+=mat.m20;
    m21+=mat.m21;
    m22+=mat.m22;
    return *this;
}

inline tMat tinyMatrix::operator-(const tMat& mat) const {
    return tMat(
        m00-mat.m00,
        m01-mat.m01,
        m02-mat.m02,
        m10-mat.m10,
        m11-mat.m11,
        m12-mat.m12,
        m20-mat.m20,
        m21-mat.m21,
        m22-mat.m22);
}

inline tMat& tinyMatrix::operator-=(const tMat& mat){
    m00-=mat.m00;
    m01-=mat.m01;
    m02-=mat.m02;
    m10-=mat.m10;
    m11-=mat.m11;
    m12-=mat.m12;
    m20-=mat.m20;
    m21-=mat.m21;
    m22-=mat.m22;
    return *this;
}

inline tMat tinyMatrix::operator*(const double& scalar) const {
    return tMat(
        m00*scalar,
        m01*scalar,
        m02*scalar,
        m10*scalar,
        m11*scalar,
        m12*scalar,
        m20*scalar,
        m21*scalar,
        m22*scalar);
}

inline tMat& tinyMatrix::operator*=(const double& scalar){
    m00*=scalar;
    m01*=scalar;
    m02*=scalar;
    m10*=scalar;
    m11*=scalar;
    m12*=scalar;
    m20*=scalar;
    m21*=scalar;
    m22*=scalar;
    return *this;
}

inline tMat tinyMatrix::operator/(const double& scalar) const {
    return tMat(
        m00/scalar,
        m01/scalar,
        m02/scalar,
        m10/scalar,
        m11/scalar,
        m12/scalar,
        m20/scalar,
        m21/scalar,
        m22/scalar);
}

inline tMat& tinyMatrix::operator/=(const double& scalar){
    m00/=scalar;
    m01/=scalar;
    m02/=scalar;
    m10/=scalar;
    m11/=scalar;
    m12/=scalar;
    m20/=scalar;
    m21/=scalar;
    m22/=scalar;
    return *this;
}

inline tMat operator *(const double& scalar, const tMat& mat){
    return tMat(
        mat.m00*scalar,
        mat.m01*scalar,
        mat.m02*scalar,
        mat.m10*scalar,
        mat.m11*scalar,
        mat.m12*scalar,
        mat.m20*scalar,
        mat.m21*scalar,
        mat.m22*scalar);
}

// mathematical functions

inline double tinyMatrix::magnitude() const {
//    return 2.0*sqrt((m01*m10+m20*m02+m12*m21-(m00*m11+m11*m22+m22*m00))); // this gives nan
    return sqrt(0.5*(m00*m00+m11*m11+m22*m22+2.0*(m01*m10+m20*m02+m12*m21))); // this WORKS
//    return 2.0*sqrt(m00*m00+m11*m11+m22*m22+2.0*(m01*m10+m20*m02+m12*m21)); // this does not work
}

/*/////////////////////////////////////////////////////////////////////////////////////////
// INTER-CLASS
/////////////////////////////////////////////////////////////////////////////////////////*/

inline tVect quat2vec(tQuat quat) {
    //if (abs(quat.q0)>0.001*abs(quat.q1)){
//    cout<<"transform error="<<quat.q0<<"\n";
    //}
    return tVect(
            quat.q1,
            quat.q2,
            quat.q3);
}

inline tQuat vec2quat(tVect vec){
    return tQuat(
            0.0,
            vec.x,
            vec.y,
            vec.z);
}

inline tVect project(tVect vec, tQuat quat){
    // projection of a vector  v from a reference frame to another.
    // the reference frame is identified by a quaternion q
    // v'=qvq*
    // v=q*v'q
    const tQuat qAdjoint=quat.adjoint();
    const tQuat vQuat=vec2quat(vec);
    tQuat rotQuat=quat.multiply(vQuat);
    rotQuat=rotQuat.multiply(qAdjoint);
    return quat2vec(rotQuat);
}

inline tVect newtonAcc(tVect moment, tVect I, tVect wBf){
    const double xA=(moment.x+(I.y-I.z)*wBf.y*wBf.z)/I.x;
    const double yA=(moment.y+(I.z-I.x)*wBf.z*wBf.x)/I.y;
    const double zA=(moment.z+(I.x-I.y)*wBf.x*wBf.y)/I.z;
    return tVect(xA,yA,zA);
}

inline tQuat quatAcc(tVect waBf, tQuat Q1){
    tQuat waQuat=vec2quat(waBf);
    waQuat.q0=-2.0*Q1.norm2();
    return waQuat;
}

inline tVect computeCentrifugal(const tVect& position, const tVect& rotationCenter, const tVect& rotationSpeed) {
    
    // Calculate the centrifugal acceleration
    // vector from rotating center to the particle position
    const tVect dist2Center = position - rotationCenter;
    // epurate from component aligned with the rotation axes
    //const tVect part2Axis = dist2Center - unitRot * (dist2Center.dot(unitRot));
    // vector for the first part of the centrifugal acceleration
    const tVect VecCenPar = rotationSpeed.cross(dist2Center);
    // centrifugal acceleration
    return -1.0 * rotationSpeed.cross(VecCenPar);

}

inline tVect computeCoriolis(const tVect& velocity, const tVect& rotationSpeed) {

    // Calculate Coriolis acceleration
    // Get Coriolis acceleration: -2 w X v
    return -2.0 * rotationSpeed.cross(velocity);
}
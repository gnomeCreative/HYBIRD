

#include <vector>
#include <iostream>
//
#include "utils.h"

#include <string>

using namespace std;

/*/////////////////////////////////////////////////////////////////////////////////////////
// WALL
/////////////////////////////////////////////////////////////////////////////////////////*/

double wall::dist(tVect pt) const{ // ACHTUNG!! This assume vector n is unitary
    const tVect d=pt-p;
    return n.dot(d);
}

void wall::wallShow() const {
    cout<<"Wall number "<<index<<": ";
    cout<<"with base point:";
    p.show();
    cout<<" and normal vector:";
    n.show();
    cout<<" moving:"<<moving<<"\n";
    if (moving) {
        cout<<"with translational velocity:";
        vel.show();
        cout<<"\n";
        cout<<"and rotational velocity:";
        omega.show();
        cout<<" around:";
        rotCenter.show();
        cout<<"\n";
    }
}

tVect wall::getSpeed(tVect pt) const {
    if (moving) {
        // distance rotation center
        const tVect distFromCenter=pt-rotCenter;
        // distance from axes
        const tVect distFromAxes=distFromCenter-(distFromCenter.dot(n))*n;
        // tangential velocity
        return vel+omega.cross(distFromAxes);
    }
    else {
        return tVect(0.0,0.0,0.0);
    }
}

/*/////////////////////////////////////////////////////////////////////////////////////////
// CYLINDER
/////////////////////////////////////////////////////////////////////////////////////////*/

double cylinder::dist(tVect pt) const {
    // calculates the distance between a given point and the inner surface of the cylinder
    
    // distance to point 1 of axis
    const tVect p1dist=pt-p1;
    // same but projected on the axis
    const tVect p1distax=(p1dist.dot(naxes))*naxes;
    // distance of point from cylinder axis
    const tVect p1distcylinder=p1distax-p1dist;
    // norm of the distance
    if (type==EMPTY) return (R-p1distcylinder.norm());
    else if (type==FULL) return (p1distcylinder.norm()-R);
    else {
        cout<<"Problem with cylinder type"<<endl;
        exit(0);
    }
}

tVect cylinder::vecDist(tVect pt) const {
    // calculates the distance between a given point and the surface of the cylinder

    // distance to point one
    const tVect pt1=pt-p1;
    // projected point on the cylinder axes
    const tVect projectedPoint=p1+(naxes.dot(pt1))*naxes;
    // distance of point from axes
    const tVect distFromAxes=pt-projectedPoint;
    // vectorized distance of point from the cylinder surface
    const tVect d = pt-(projectedPoint+distFromAxes/distFromAxes.norm()*R);
    return d;

}

tVect cylinder::getSpeed(tVect pt) const {
    if (moving) {
        // distance to point one
        const tVect pt1=pt-p1;
        // projected point on the cylinder axes
        const tVect projectedPoint=p1+(naxes.dot(pt1))*naxes;
        // distance of point from axes
        const tVect distFromAxes=pt-projectedPoint;
        // tangential velocity
        return omega.cross(distFromAxes);
    }
    else {
        return tVect(0.0,0.0,0.0);
    }
}

double cylinder::getTheta(const tVect& pt, const tVect& zeroDir, const tVect& thetaDir) const {
    // calculates angular coordinate, the zero being a point located directly under the axes (coord. Z)

    // distance to axis point 1
    const tVect p1dist=pt-p1;
    // same but projected on the axis
    const tVect p1distax=(p1dist.dot(naxes))*naxes;
    // distance of point from cylinder axis
    const tVect p1distcylinder=p1distax-p1dist;

    // sin of the angle
    const double cos=p1distcylinder.dot(zeroDir);
    //cos of the angle
    const double sin=p1distcylinder.dot(thetaDir);
    // angular coordinate
    return atan(sin/cos);

}

double cylinder::getR(const tVect& pt) const {
    // calculates distance from cylinder axes

    // distance to point 1 of axes
    const tVect p1dist=pt-p1;
    // same but projected on the axes
    const tVect p1distax=(p1dist.dot(naxes))*naxes;
    // distance center to cylinder
    const tVect p1distcylinder=p1distax-p1dist;
    // norm of the distance
    return (p1distcylinder.norm());
}

void cylinder::initAxes() {
    // initializes the axes variables (outside dist function for improved performance)

    // axes vector
    axes =p2-p1;
    // axes unit vector
    naxes=axes/axes.norm();
}

void cylinder::cylinderShow() const {
    if (type==EMPTY)
        cout<<"Inner ";
    else if (type==FULL) {
        cout<<"Outer ";
    }
    cout<<"cylinder, number "<<index;
    cout<<" with base point: p1";
    p1.show();
    cout<<", p2";
    p2.show();
    cout<<", and radius: "<<R<<endl;
    if (moving) {
        cout<<"rotating with speed:";
        omega.show();
        cout<<"\n";
    }
    
}

double cylinder::segmentIntercept(const tVect& start, const tVect& dir) const {
    /// Line segment VS cylinder
    // - cylinder (A, B, r) (start point, end point, radius) -> in our case (p1, p2, R)
    // - line has starting point (x0, y0, z0) and ending point (x0+ux, y0+uy, z0+uz) ((ux, uy, uz) is "direction")
    // for our purposes, the starting point is always OUTSIDE, the final point INSIDE
    //   optimize? (= don't care for t > 1)
    // <= t  = "time" of intersection
    //   norm = surface normal of intersection point
    //  t = NaN;

  // Solution : http://www.gamedev.net/community/forums/topic.asp?topic_id=467789 <- Thanks!
//  double cxmin, cymin, czmin, cxmax, cymax, czmax;
//  if (p1.z < p2.z) { czmin = p1.z - r; czmax = p2.z + r; } else { czmin = p2.z - r; czmax = p1.z + r; }
//  if (p1.y < p2.y) { cymin = p1.y - r; cymax = p2.y + r; } else { cymin = p2.y - r; cymax = p1.y + r; }
//  if (p1.x < p2.x) { cxmin = p1.x - r; cxmax = p2.x + r; } else { cxmin = p2.x - r; cxmax = p1.x + r; }
//  if (optimize) {
//   if (start.z >= czmax && (start.z + dir.z) > czmax) return;
//   if (start.z <= czmin && (start.z + dir.z) < czmin) return;
//   if (start.y >= cymax && (start.y + dir.y) > cymax) return;
//   if (start.y <= cymin && (start.y + dir.y) < cymin) return;
//   if (start.x >= cxmax && (start.x + dir.x) > cxmax) return;
//   if (start.x <= cxmin && (start.x + dir.x) < cxmin) return;
//  }

    // in my case A=p1, B=p2
    const tVect AB=p2-p1;
    const tVect AO=start-p1;
    const tVect AOxAB=AO.cross(AB);
    const tVect VxAB=dir.cross(AB);
    const double ab2=AB.dot(AB);
    const double a=VxAB.dot(VxAB);
    const double b=2*VxAB.dot(AOxAB);
    const double c=AOxAB.dot(AOxAB)-(R*R*ab2);
    const double delta = b * b - 4 * a * c;
    if (delta < 0) {
        cout<<"Error in cylindrical boundary, delta<0"<<endl;
        return 0.0;
    }
    const double time1 = (-b -sqrt(delta)) / (2 * a);
    const double time2 = (-b -sqrt(delta)) / (2 * a);
    if (time1<0 || time2<0) {
        cout<<"Error in cylindrical boundary, t1<0 or t2<0"<<endl;
        return 0.0;
    }

    const double time=std::min(time1,time2);
//    const tVect intersection = start + dir * time;        /// intersection point
    return time;
//    tVect projection = p1 + (AB.dot(intersection - p1) / ab2) * AB; /// intersection projected onto cylinder axis
//    if ((projection - p1).norm() + (p2 - projection).norm() > AB.norm()) return; /// THIS IS THE SLOW SAFE WAY
    //if (projection.z > czmax - r || projection.z < czmin + r ||
    // projection.y > cymax - r || projection.y < cymin + r ||
    // projection.x > cxmax - r || projection.x < cxmin + r ) return; /// THIS IS THE FASTER BUGGY WAY

//    normal = (intersection - projection);
//    normal.normalise();
//    t = time; /// at last
}

/*/////////////////////////////////////////////////////////////////////////////////////////
// TOPOGRAPHY
/////////////////////////////////////////////////////////////////////////////////////////*/

void topography::readFromFile(string& topographyFileName, double& transX, double& transY, double& transZ) {
    ifstream topographyFileID;
    cout << "Reading " << topographyFileName.c_str() << endl;
    // open topography file
    topographyFileID.open(topographyFileName.c_str(), ios::in);
    // check if it exists
    ASSERT(topographyFileID.is_open());
    // skip header
    topographyFileID.ignore(10000,'\n');
    // read first two lines
    double dummy,firstX,firstY,secondX,secondY;
    topographyFileID>>firstX;
    topographyFileID>>firstY;
    topographyFileID>>dummy;
    topographyFileID>>dummy;
    topographyFileID>>secondX;
    topographyFileID>>secondY;
    topographyFileID>>dummy;
    topographyFileID>>dummy;
    // translate
    firstX+=transX;
    secondX+=transX;
    firstY+=transY;
    secondY+=transY;
    cout<<"firstX="<<firstX<<" firstY="<<firstY<<" secondX="<<secondX<<" secondY="<<secondY<<endl;
    bool progressIsX=false;
    if (firstX!=secondX) {
        progressIsX=true;
    } else if (firstY!=secondY) {
        progressIsX=false;
    }
    // advance in Y first is not supported
    ASSERT(progressIsX=true);
    spacingX=secondX-firstX;
    topographyFileID.close();
    
    // find number of topography points
    int totPoints = 0;
    string line;
    topographyFileID.open(topographyFileName.c_str(), ios::in);
    topographyFileID.ignore(10000,'\n');
    while (std::getline(topographyFileID, line)) {
        ++totPoints;
    }
    std::cout << "Number of lines in text file: " << totPoints<<endl;
    topographyFileID.close();
    
    // find rows and columns
    topographyFileID.open(topographyFileName.c_str(), ios::in);
    topographyFileID.ignore(10000,'\n');
    sizeX=0;
    bool foundEndX=false;
    while (foundEndX==false) {
        double coordXHere;
        double coordYHere;
        topographyFileID>>coordXHere;
        topographyFileID>>coordYHere;
        topographyFileID>>dummy;
        topographyFileID>>dummy;
        // translate
        coordXHere += transX;
        coordYHere += transY;
        if (coordYHere == firstY) {
            coordX.push_back(coordXHere);
            sizeX++;
        }
        else {
            secondY=coordYHere;
            foundEndX=true;
        }
    }
    spacingY=secondY-firstY;
    sizeY=totPoints/sizeX;
    cout<<"spacingX="<<spacingX<<" spacingY="<<spacingY<<" sizeX="<<sizeX<<" sizeY="<<sizeY<<endl;    
    
    topographyFileID.close();
    
    // read data
    topographyFileID.open(topographyFileName.c_str(), ios::in);
    topographyFileID.ignore(10000,'\n');
    for (int iy = 0; iy < sizeY; iy++) {
        vecList pointRow;
        vecList surfaceRow;
        for (int ix = 0; ix < sizeX; ix++) {
            double x, y, z, h;
            topographyFileID>>x;
            topographyFileID>>y;
            topographyFileID>>z;
            topographyFileID>>h;
            // add to row
            x+=transX;
            y+=transY;
            z+=transZ;
            tVect coordHere(x, y, z);
            tVect surfCoordHere(x, y, z+h);
            pointRow.push_back(coordHere);
            surfaceRow.push_back(surfCoordHere);
            // save coordinates
            if (ix == 0) {
                coordY.push_back(y);
            }
        }
        points.push_back(pointRow);
        surface.push_back(surfaceRow);
    }
    corner1X=coordX[0];
    corner1Y=coordY[0];
    corner2X=coordX[sizeX-1];
    corner2Y=coordY[sizeY-1];
    
    ASSERT(spacingX != 0);
    ASSERT(spacingY != 0);
    
    // fix case with negative spacing
    if (spacingY < 0.0) {
        std::reverse(points.begin(), points.end());
        std::reverse(surface.begin(), surface.end());
        std::reverse(coordY.begin(), coordY.end());
        corner1Y = coordY[0];
        corner2Y = coordY[sizeY - 1];
        spacingY = coordY[1] - coordY[0];
    }

    if (spacingX < 0.0) {
        for (int ix = 0; ix < sizeX; ix++) {
            std::reverse(points[ix].begin(), points[ix].end());
            std::reverse(surface[ix].begin(), surface[ix].end());
        }
        std::reverse(coordX.begin(), coordX.end());
        corner1X = coordX[0];
        corner2X = coordX[sizeX - 1];
        spacingX = coordX[1] - coordX[0];
    }
    
    ASSERT(spacingX > 0);
    ASSERT(spacingY > 0);
    

}

void topography::show() {

    cout<<"Topography with size "<<sizeX<<" x "<<sizeY<<endl;
    cout<<"Points with spacing "<<spacingX<<" x "<<spacingY<<":"<<endl;
    cout<<"X range: from "<<coordX[0]<<" to "<<coordX[sizeX - 1]<<endl;
    cout<<"Y range: from "<<coordY[0]<<" to "<<coordY[sizeY - 1]<<endl;
//    for (int iy = 0; iy < sizeY; iy++) {
//        for (int ix = 0; ix < sizeX; ix++) {
//            cout<<"("<<ix<<","<<iy<<")"<<" "<<points[iy][ix].dot(Xp)<<" "<<points[iy][ix].dot(Yp)<<" "<<points[iy][ix].dot(Zp)<<" "<<fluidLevel[iy][ix]<<endl;
//        }
//    }

}


void topography::getReferenceTriangle(const tVect& point, tVect& point1, tVect& point2, tVect& point3, tVect& planeNormal)  const{
    
        // coordinates of the point
    const double pointX = point.dot(Xp);
    const double pointY = point.dot(Yp);

    // approximate indices of the point, projected on the xy plane
    const double indexX = (pointX - corner1X) / spacingX;
    //cout<<"pointX="<<pointX<<" indexX="<<indexX<<endl;
    const double indexY = (pointY - corner1Y) / spacingY;
    //cout<<"pointY="<<pointY<<" indexY="<<indexY<<endl;

    // indices of the inscribing rectangle
    const int negX = floor(indexX);
    int posX = ceil(indexX);
    if (negX==posX) posX++;
    const int negY = floor(indexY);
    int posY = ceil(indexY);
    if (negY==posY) posY++;
    //cout<<negX<<" "<<posX<<" "<<negY<<" "<<posY<<endl;
    // use the approximate indices to decide which triangle to use
    bool useTop = false;
    if ((indexX - double(negX))<(indexY - double(negY))) {
        useTop = true;
    }

    // setup triangle to use
    point1 = points[negY][negX];
    if (useTop) {
        point2 = points[posY][posX];
        point3 = points[posY][negX];
    } else {
        point2 = points[negY][posX];
        point3 = points[posY][posX];
    }

    // compute triangle normal
    const tVect side1=point2-point1;
    const tVect side2=point3-point1;
    
    planeNormal=side1.cross(side2);
    planeNormal=planeNormal/planeNormal.norm();

}

void topography::getSurfaceTriangle(const tVect& point, tVect& point1, tVect& point2, tVect& point3, tVect& planeNormal)  const{
    
        // coordinates of the point
    const double pointX = point.dot(Xp);
    const double pointY = point.dot(Yp);

    // approximate indices of the point, projected on the xy plane
    const double indexX = (pointX - corner1X) / spacingX;
    const double indexY = (pointY - corner1Y) / spacingY;

    // indices of the inscribing rectangle
    const int negX = floor(indexX);
    int posX = ceil(indexX);
    if (negX==posX) posX++;
    const int negY = floor(indexY);
    int posY = ceil(indexY);
    if (negY==posY) posY++;
    //cout<<negX<<" "<<posX<<" "<<negY<<" "<<posY<<endl;
    // use the approximate indices to decide which triangle to use
    bool useTop = false;
    if ((indexX - double(negX))<(indexY - double(negY))) {
        useTop = true;
    }

    // setup triangle to use
    point1 = surface[negY][negX];
    if (useTop) {
        point2 = surface[posY][posX];
        point3 = surface[posY][negX];
    } else {
        point2 = surface[negY][posX];
        point3 = surface[posY][posX];
    }

    // compute triangle normal
    const tVect side1=point2-point1;
    const tVect side2=point3-point1;
    
    planeNormal=side1.cross(side2);
    planeNormal=planeNormal/planeNormal.norm();

}

double topography::distance(const tVect& point)  const{

    tVect point1,point2,point3;
    tVect planeNormal;

    getReferenceTriangle(point,point1,point2,point3,planeNormal);

    // vector connecting the interrogation point and one point of the triangle
    const tVect connection=point-point1;
    // distance is the dot product of these two
    return connection.dot(planeNormal);
    
}

tVect topography::surfaceNormal(const tVect& point)  const{

    tVect point1,point2,point3;
    tVect planeNormal;
    
    getReferenceTriangle(point,point1,point2,point3,planeNormal);

    // distance is the dot product of these two
    return planeNormal;
    
}

double topography::surfaceIsoparameter(const tVect& point)  const{

    tVect topPoint1,topPoint2,topPoint3;
    tVect topPlaneNormal;
    getReferenceTriangle(point,topPoint1,topPoint2,topPoint3,topPlaneNormal);
    
    tVect surfPoint1,surfPoint2,surfPoint3;
    tVect surfPlaneNormal;
    getSurfaceTriangle(point,surfPoint1,surfPoint2,surfPoint3,surfPlaneNormal);

    // vector connecting the interrogation point and one point of the surface triangle
    const tVect surfaceConnection=point-surfPoint1;
    // vector connecting the interrogation point and one point of the topographical surface triangle
    const tVect topographyConnection=point-topPoint1;
    
    const double surfDistance=surfaceConnection.dot(surfPlaneNormal);
    const double topDistance=topographyConnection.dot(topPlaneNormal);
    
    // distance is the dot product of these two
    return -1.0*topDistance/(surfDistance-topDistance);
    
}

double topography::directionalDistance(const tVect& point, const tVect& dir)  const{
    // distance point to topography, measured along a line defined by dir (normalized vector)
    // see https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
    
    tVect point1,point2,point3;
    tVect planeNormal;
    getReferenceTriangle(point,point1,point2,point3,planeNormal);

    // vector connecting the interrogation point and one point of the triangle
    const tVect connection=point-point1;
    // distance is the dot product of these two
    return -1.0*connection.dot(planeNormal)/dir.dot(planeNormal);
    
}


/*/////////////////////////////////////////////////////////////////////////////////////////
// PERIODIC BOUNDARY CONDITION (PBC)
/////////////////////////////////////////////////////////////////////////////////////////*/

void pbc::setPlanes(){
    pl1.p=p;
    pl2.p=p+v;
    pl1.n=v/v.norm();
    pl2.n=-1.0*pl1.n;
}

void pbc::pbcShow() const {
    cout<<"Periodic condition number "<<index<<"\n";
    cout<<"with base point: ";
    p.show();
    cout<<" and translation vector:";
    v.show();
    cout<<"\n";

    cout<<"plane 1: ";
    pl1.p.show();
    cout<<" and normal vector:";
    pl1.n.show();
    cout<<"\n";

    cout<<"plane 2: ";
    pl2.p.show();
    cout<<" and normal vector:";
    pl2.n.show();
    cout<<"\n";

}

/*/////////////////////////////////////////////////////////////////////////////////////////
// ENERGY
/////////////////////////////////////////////////////////////////////////////////////////*/

void energy::reset(){
    // reset all energy types for
    kin=0.0;
    rotKin=0.0;
    trKin=0.0;
    grav=0.0;
    elastic=0.0;
    total=0.0;
}

//void energy::updateKinetic(elmtList& elmts){
//
//    double tKin, rKin;
//    tKin=rKin=0.0;
//    for (int i=0; i<elmts.size(); i++){
//        tKin+=0.5*elmts[i].m*elmts[i].x1.norm2();
//        // adjoint of orientation quaternion
//        const tQuat q0adj=elmts[i].q0.adjoint();
//        // rotational velocity (body-fixed reference frame)
//        const tVect w=2.0*quat2vec( q0adj.multiply( elmts[i].q1 ) );
////        Tvect w=2.0*quat2vec( elmts[i].q1.multiply( elmts[i].q0.adjoint() ) );
//        const tVect wSquare=w.compProd(w);
//        rKin+=elmts[i].I.dot(wSquare);
//    }
//
//    rotationalKinetic=rKin;
//    translationalKinetic=tKin;
//}

//void energy::updatePotential(elmtList& elmts, tVect& gravity, doubleList& dim){
//    int i;
//    wall zeroWall;
//    double grav,gravityNorm,heigth;
//    gravityNorm=gravity.norm();
//    grav=0.0;
//    if (gravityNorm!=0.0) {
//        zeroWall.n=(-1.0)*gravity/gravityNorm;
//        if (zeroWall.n.x>0) zeroWall.p.x=0.0;
//        else zeroWall.p.x=dim[0]-1.0;
//        if (zeroWall.n.y>0) zeroWall.p.y=0.0;
//        else zeroWall.p.y=dim[1]-1.0;
//        if (zeroWall.n.z>0) zeroWall.p.z=0.0;
//        else zeroWall.p.z=dim[2]-1.0;
//        for (i=0; i<elmts.size(); i++){
//            heigth=zeroWall.dist(elmts[i].x0);
//            grav+=elmts[i].m*heigth*gravityNorm;
//        }
//    }
//    gravitational=grav;
//}

void energy::updateTotal(){
    kin=rotKin+trKin;
    total=rotKin+trKin+grav+elastic;
}

void energy::show(){
    cout<<"Total energy ="<<total<<"; (Kinetic= "<<kin<<"; Gravitational= "<<grav<<"; Elastic= "<<elastic<<");\t";
}

#ifndef CYLINDER2_H
#define CYLINDER2_H
#include "myvector.h"
#include "utils.h"

struct DEMParams;

struct Cylinder2 {
    unsigned int alloc = 0;
    unsigned int count = 0;

    //index
    unsigned int *index = nullptr;
    // two axes point
    tVect *p1 = nullptr;
    tVect *p2 = nullptr;
    // radius
    double *R = nullptr;
    // axes vector and unit vector (defined in init_axes for speed)
    tVect *axes = nullptr;
    tVect *naxes = nullptr;
    // rotational speed
    tVect *omega = nullptr;
    // is it a rotating wall (true = rotating; false = fixed);
    bool *moving = nullptr;
    // is it a slip cylinder?
    bool *slip = nullptr;
    // inner or outer
    cylinderType *type = nullptr;
    // limits for cylinders
    bool *limited = nullptr;
    double *xMin = nullptr, *xMax = nullptr;
    double *yMin = nullptr, *yMax = nullptr;
    double *zMin = nullptr, *zMax = nullptr;
    // is it translating?
    bool *translating = nullptr;
    // translation vector
    tVect *trans = nullptr;

    /**
     * Does not update count
     * Initialises based on original constructor
     */
    template<int impl>
    void memoryAlloc(unsigned int num);
    void initialize(const DEMParams& dem_p);
    // axes variables initialization
    void initAxes(unsigned int i);
    void cylinderShow(unsigned int i) const;

    // get tangential speed of a point rotating with the cylinder
    __host__ __device__ tVect getSpeed(unsigned int i, tVect pt) const {
        if (moving[i]) {
            // distance to point one
            const tVect pt1 = pt - p1[i];
            // projected point on the cylinder axes
            const tVect projectedPoint = p1[i] + (naxes[i].dot(pt1)) * naxes[i];
            // distance of point from axes
            const tVect distFromAxes = pt - projectedPoint;
            // tangential velocity
            return omega[i].cross(distFromAxes);
        }
        return tVect(0.0, 0.0, 0.0);
    }
    double segmentIntercept(const unsigned int i, const tVect& start, const tVect& dir) const {
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
        const tVect AB=p2[i] - p1[i];
        const tVect AO=start-p1[i];
        const tVect AOxAB=AO.cross(AB);
        const tVect VxAB=dir.cross(AB);
        const double ab2=AB.dot(AB);
        const double a=VxAB.dot(VxAB);
        const double b=2*VxAB.dot(AOxAB);
        const double c=AOxAB.dot(AOxAB)-(R[i] *R[i] *ab2);
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
};

template<>
inline void Cylinder2::memoryAlloc<CPU>(unsigned int num) {
    alloc = num;
    if (!num) return;

    index = (unsigned int*)malloc(alloc * sizeof(unsigned int));
    p1 = (tVect*)malloc(alloc * sizeof(tVect));
    p2 = (tVect*)malloc(alloc * sizeof(tVect));
    R = (double*)malloc(alloc * sizeof(double));
    axes = (tVect*)malloc(alloc * sizeof(tVect));
    naxes = (tVect*)malloc(alloc * sizeof(tVect));
    omega = (tVect*)malloc(alloc * sizeof(tVect));
    moving = (bool*)malloc(alloc * sizeof(bool));
    slip = (bool*)malloc(alloc * sizeof(bool));
    type = (cylinderType*)malloc(alloc * sizeof(cylinderType));
    limited = (bool*)malloc(alloc * sizeof(bool));
    xMin = (double*)malloc(alloc * sizeof(double));
    xMax = (double*)malloc(alloc * sizeof(double));
    yMin = (double*)malloc(alloc * sizeof(double));
    yMax = (double*)malloc(alloc * sizeof(double));
    zMin = (double*)malloc(alloc * sizeof(double));
    zMax = (double*)malloc(alloc * sizeof(double));
    translating = (bool*)malloc(alloc * sizeof(bool));
    trans = (tVect*)malloc(alloc * sizeof(tVect));
    // Init
    memset(index, 0, alloc * sizeof(unsigned int));
    memset(p1, 0, alloc * sizeof(tVect));
    std::fill(p2, p2 + alloc, tVect(1.0, 0.0, 0.0));
    std::fill(R, R + alloc, 1.0);
    std::fill(axes, axes + alloc, p2[0]-p1[0]);
    std::fill(naxes, naxes + alloc, axes[0]/axes[0].norm());
    memset(omega, 0, alloc * sizeof(tVect));
    std::fill(moving, moving + alloc, false);
    std::fill(slip, slip + alloc, false);
    std::fill(type, type + alloc, EMPTY);
    // limited?
    // x-z min/max?
    std::fill(translating, translating + alloc, false);
    memset(trans, 0, alloc * sizeof(tVect));
}
#endif  // CYLINDER2_H

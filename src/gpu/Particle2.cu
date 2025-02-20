#include "Particle2.h"

#include "Element2.h"

void Particle2::updateCorrected(const unsigned int p_i, const Element2 &elements, const unsigned int e_i, const DEMParams& dem_p) {
    // updating position and velocity for simple case
    x0[p_i] = elements.x0[e_i];
    radiusVec[p_i] = Zero;
    x1[p_i] = elements.x1[e_i];

    if (elements.size[e_i] > 1) {
        x0[p_i] = x0[p_i] + r[p_i] * project(dem_p.prototypes[elements.size[e_i]][protoIndex[p_i]], elements.q0[e_i]);
        // updating radius (distance of particle center of mass to element center of mass)
        radiusVec[p_i] = x0[p_i] - elements.x0[e_i];
        x1[p_i] = x1[p_i] + elements.wGlobal[e_i].cross(radiusVec[p_i]);
    }

}
#include "Wall2.h"

#include "LBParams.h"
#include "DEMParams.h"

extern ProblemName problemName;

void Wall2::initialize(const std::array<types, 6> &externalBoundary, const std::array<tVect, 6> &boundaryLocation) {
    // boundary directors
    vecList boundaryDirec;
    boundaryDirec.resize(6);
    boundaryDirec[0] = tVect(1.0, 0.0, 0.0);
    boundaryDirec[1] = tVect(-1.0, 0.0, 0.0);
    boundaryDirec[2] = tVect(0.0, 1.0, 0.0);
    boundaryDirec[3] = tVect(0, 0. - 1.0, 0.0);
    boundaryDirec[4] = tVect(0.0, 0.0, 1.0);
    boundaryDirec[5] = tVect(0.0, 0.0, -1.0);

    // Count basic walls
    unsigned int basic_wall_count = 0;
    for (int i = 0; i < externalBoundary.size(); ++i) {
        if ((externalBoundary[i] == SLIP_STAT_WALL) || (externalBoundary[i] == SLIP_DYN_WALL) || (externalBoundary[i] == STAT_WALL) || (externalBoundary[i] == DYN_WALL)) {
            ++basic_wall_count;
        }
    }
    unsigned int index = basic_wall_count;

    // additional walls
    switch (problemName) {
        case HK_LARGE:
        {
            // 1-2 extra walls
            memoryAlloc<CPU>(basic_wall_count + (DEM_P.depositArea ? 2 : 1));
            // storage container
            this->p[index] = tVect(10.0, 0.0, 0.0);
            this->n[index] = tVect(0.17365, 0.98481, 0.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating[index] = false;
            this->limited[index] = false;
            ++index;
            if (DEM_P.depositArea) {
                // deposit area
                this->p[index] = tVect(25.2795, 0.5008, 0.0);
                this->n[index] = tVect(-0.342020, 0.939693, 0.0);
                this->index[index] = index;
                this->moving[index] = false;
                this->rotCenter[index].reset();
                this->omega[index].reset();
                this->vel[index].reset();
                this->translating[index] = false;
                this->limited[index] = false;
                ++index;
            }
            break;
        }
        case KELVIN:
        {
            // 2 extra walls
            memoryAlloc<CPU>(basic_wall_count + 2);
            // ramp
            this->p[index] = tVect(1.42535, 0.04048, 0.0);
            this->n[index] = tVect(-0.43837, 0.89879, 0.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating[index] = false;
            this->limited[index] = true;
            this->xMin[index] = 1.42535;
            this->xMax[index] = 1.65005;
            this->yMin[index] = -1.0;
            this->yMax[index] = 0.15008;
            ++index;
            // ramp
            this->p[index] = tVect(1.65005, 0.15008, 0.0);
            this->n[index] = tVect(0.89879, 0.43837, 0.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating[index] = false;
            this->limited[index] = true;
            this->xMin[index] = 1.65005;
            this->xMax[index] = 1.72324;
            this->yMin[index] = -1.0;
            this->yMax[index] = 0.15008;
            ++index;
            break;
        }
        case AVALANCHE:
        {
            // 1 extra walls
            memoryAlloc<CPU>(basic_wall_count + 1);
            this->p[index] = tVect(60, 15, 15 - 14.86);
            this->n[index] = tVect(0.42262, 0.0, 0.9063);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating[index] = false;
            this->limited[index] = false;
            ++index;
            break;
        }
        case DRUM:
        {
            unsigned int extra_count = 0;
            // left wall, just next to basic wall 2
            if ((externalBoundary[2] == SLIP_STAT_WALL) || (externalBoundary[2] == SLIP_DYN_WALL) || (externalBoundary[2] == STAT_WALL) || (externalBoundary[2] == DYN_WALL)) {
                ++extra_count;
            }
            // right wall, just next to basic wall 3
            if ((externalBoundary[3] == SLIP_STAT_WALL) || (externalBoundary[3] == SLIP_DYN_WALL) || (externalBoundary[3] == STAT_WALL) || (externalBoundary[3] == DYN_WALL)) {
                ++extra_count;
            }
            // extra walls
            memoryAlloc<CPU>(basic_wall_count + extra_count);
            break;
        }
        case ZHOU:
        {
            // 1 extra walls
            memoryAlloc<CPU>(basic_wall_count + 1);
            this->p[index] = tVect(0.053, 0.01, 0.21);
            this->n[index] = tVect(0.0, 0.0, 1.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating[index] = false;
            this->limited[index] = true;
            this->xMin[index] = 0.0525;
            this->xMax[index] = 0.3675;
            ++index;
            break;
        }
        case OPENBARRIER:
        {
            // 2 extra walls
            memoryAlloc<CPU>(basic_wall_count + 2);
            this->p[index] = tVect(0.0, 0.0375, 0.0);
            this->n[index] = tVect(0.0, 1.0, 0.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating[index] = false;
            this->limited[index] = false;
            ++index;

            this->p[index] = tVect(0.0, 2.0625, 0.0);
            this->n[index] = tVect(0.0, -1.0, 0.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating[index] = false;
            this->limited[index] = false;
            ++index;

            break;
        }
        case HONGKONG:
        {
            const double width = 0.2;

            if (DEM_P.hongkongSlitSize < 0.9 * width && DEM_P.hongkongSmoothWall) {
                // 4 extra walls
                memoryAlloc<CPU>(basic_wall_count + 4);

                const double edgeRadius = 0.005;
                const double xBarrierPosition = 1.4;
                const double yBarrierSize = (width - DEM_P.hongkongSlitSize) / 2.0;
                const double barrierLeftWing = yBarrierSize;
                const double barrierRightWing = yBarrierSize + DEM_P.hongkongSlitSize;
                
                this->p[index] = tVect(xBarrierPosition - edgeRadius * 0.95, 0.0, 0.0);
                this->n[index] = tVect(-1.0, 0.0, 0.0);
                this->index[index] = index;
                this->moving[index] = false;
                this->rotCenter[index].reset();
                this->omega[index].reset();
                this->vel[index].reset();
                this->translating[index] = false;
                this->limited[index] = true;
                this->yMin[index] = -100.0;
                this->yMax[index] = barrierLeftWing - edgeRadius;
                ++index;

                this->p[index] = tVect(xBarrierPosition - edgeRadius * 0.95, 0.0, 0.0);
                this->n[index] = tVect(-1.0, 0.0, 0.0);
                this->index[index] = index - 1;
                this->moving[index] = false;
                this->rotCenter[index].reset();
                this->omega[index].reset();
                this->vel[index].reset();
                this->translating[index] = false;
                this->limited[index] = true;
                this->yMin[index] = barrierRightWing + edgeRadius;
                this->yMax[index] = 100.0;
                ++index;

                this->p[index] = tVect(xBarrierPosition + edgeRadius * 0.95, 0.0, 0.0);
                this->n[index] = tVect(1.0, 0.0, 0.0);
                this->index[index] = index-2;
                this->moving[index] = false;
                this->rotCenter[index].reset();
                this->omega[index].reset();
                this->vel[index].reset();
                this->translating[index] = false;
                this->limited[index] = true;
                this->yMin[index] = -100.0;
                this->yMax[index] = barrierLeftWing - edgeRadius;
                ++index;

                this->p[index] = tVect(xBarrierPosition - edgeRadius * 0.95, 0.0, 0.0);
                this->n[index] = tVect(-1.0, 0.0, 0.0);
                this->index[index] = index - 3;
                this->moving[index] = false;
                this->rotCenter[index].reset();
                this->omega[index].reset();
                this->vel[index].reset();
                this->translating[index] = false;
                this->limited[index] = true;
                this->yMin[index] = barrierRightWing + edgeRadius;
                this->yMax[index] = 100.0;
                ++index;
            } else {
                memoryAlloc<CPU>(basic_wall_count);
            }

            break;
        }
        case HOURGLASS:
        {
            // 2 extra walls
            memoryAlloc<CPU>(basic_wall_count + 2);

            const double outletCenter = DEM_P.demSize[0] / 2.0;
            
            this->p[index] = tVect(0.0, 0.0, DEM_P.hourglassOutletHeight);
            this->n[index] = tVect(0.0, 0.0, 1.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating[index] = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = true;
            this->xMin[index] = -numeric_limits<double>::max();
            this->xMax[index] = outletCenter - DEM_P.hourglassOutletSize;
            ++index;

            this->p[index] = tVect(0.0, 0.0, DEM_P.hourglassOutletHeight);
            this->n[index] = tVect(0.0, 0.0, 1.0);
            this->index[index] = index - 1;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating[index] = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = true;
            this->xMin[index] = outletCenter + DEM_P.hourglassOutletSize;
            this->xMax[index] = numeric_limits<double>::max();
            ++index;

            break;
        }
        case IERVOLINO_2D:
        {
            // 3 extra walls
            memoryAlloc<CPU>(basic_wall_count + 3);

            const double reservoirX = 0.825;
            const double erodibleSizeX = 0.51;
            const double erodibleHeight = 0.02;
            const double edgeRadius = 0.005;
            
            // obstacle
            this->p[index] = tVect(reservoirX + erodibleSizeX, 0.0, 0);
            this->n[index] = tVect(-1.0, 0.0, 0.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating[index] = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = true;
            this->xMin[index] = -numeric_limits<double>::max();
            this->xMax[index] = +numeric_limits<double>::max();
            this->yMin[index] = -numeric_limits<double>::max();
            this->yMax[index] = +numeric_limits<double>::max();
            this->zMin[index] = -numeric_limits<double>::max();
            this->zMax[index] = +numeric_limits<double>::max();
            ++index;

            // bottom of reservoir
            this->p[index] = tVect(0.0, 0.0, erodibleHeight);
            this->n[index] = tVect(0.0, 0.0, 1.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating[index] = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = true;
            this->xMin[index] = -numeric_limits<double>::max();
            this->xMax[index] = reservoirX - edgeRadius;
            this->yMin[index] = -numeric_limits<double>::max();
            this->yMax[index] = +numeric_limits<double>::max();
            this->zMin[index] = -numeric_limits<double>::max();
            this->zMax[index] = +numeric_limits<double>::max();
            ++index;

            this->p[index] = tVect(reservoirX, 0.0, 0);
            this->n[index] = tVect(1.0, 0.0, 0.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating[index] = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = true;
            this->xMin[index] = -numeric_limits<double>::max();
            this->xMax[index] = +numeric_limits<double>::max();
            this->yMin[index] = -numeric_limits<double>::max();
            this->yMax[index] = +numeric_limits<double>::max();
            this->zMin[index] = -numeric_limits<double>::max();
            this->zMax[index] = erodibleHeight - edgeRadius;
            ++index;
            break;
        }
        case IERVOLINO:
        {
            // 14 extra walls
            memoryAlloc<CPU>(basic_wall_count + 14);

            const double centerY = DEM_P.demSize[1] / 2.0;
            const double reservoirX = 0.8;
            const double wallSizeX = 0.025;
            const double outletSizeY = 0.3;
            const double obstacleSizeX = 0.155;
            const double obstacleSizeY = 0.3;
            const double erodibleSizeX = 0.51;
            const double erodibleHeight = 0.02;
            const double edgeRadius = 0.005;

            // reservoir left wall
            this->p[index] = tVect(reservoirX, 0.0, 0);
            this->n[index] = tVect(-1.0, 0.0, 0.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = true;
            this->xMin[index] = -numeric_limits<double>::max();
            this->xMax[index] = +numeric_limits<double>::max();
            this->yMin[index] = -numeric_limits<double>::max();
            this->yMax[index] = centerY - outletSizeY / 2.0 - edgeRadius;
            this->zMin[index] = -numeric_limits<double>::max();
            this->zMax[index] = +numeric_limits<double>::max();
            ++index;

            this->p[index] = tVect(reservoirX + wallSizeX, 0.0, 0);
            this->n[index] = tVect(+1.0, 0.0, 0.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = true;
            this->xMin[index] = -numeric_limits<double>::max();
            this->xMax[index] = +numeric_limits<double>::max();
            this->yMin[index] = -numeric_limits<double>::max();
            this->yMax[index] = centerY - outletSizeY / 2.0 - edgeRadius;
            this->zMin[index] = -numeric_limits<double>::max();
            this->zMax[index] = +numeric_limits<double>::max();
            ++index;

            this->p[index] = tVect(reservoirX + wallSizeX / 2, centerY - outletSizeY / 2.0, 0);
            this->n[index] = tVect(0.0, +1.0, 0.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = true;
            this->xMin[index] = reservoirX + edgeRadius;
            this->xMax[index] = reservoirX + wallSizeX - edgeRadius;
            this->yMin[index] = -numeric_limits<double>::max();
            this->yMax[index] = +numeric_limits<double>::max();
            this->zMin[index] = -numeric_limits<double>::max();
            this->zMax[index] = +numeric_limits<double>::max();
            ++index;


            // reservoir right wall
            this->p[index] = tVect(reservoirX, 0.0, 0);
            this->n[index] = tVect(-1.0, 0.0, 0.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = true;
            this->xMin[index] = -numeric_limits<double>::max();
            this->xMax[index] = +numeric_limits<double>::max();
            this->yMin[index] = centerY + outletSizeY / 2.0 + edgeRadius;
            this->yMax[index] = +numeric_limits<double>::max();
            this->zMin[index] = -numeric_limits<double>::max();
            this->zMax[index] = +numeric_limits<double>::max();
            ++index;

            this->p[index] = tVect(reservoirX + wallSizeX, 0.0, 0);
            this->n[index] = tVect(+1.0, 0.0, 0.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = true;
            this->xMin[index] = -numeric_limits<double>::max();
            this->xMax[index] = +numeric_limits<double>::max();
            this->yMin[index] = centerY + outletSizeY / 2.0 + edgeRadius;
            this->yMax[index] = +numeric_limits<double>::max();
            this->zMin[index] = -numeric_limits<double>::max();
            this->zMax[index] = +numeric_limits<double>::max();
            ++index;

            this->p[index] = tVect(reservoirX + wallSizeX / 2, centerY + outletSizeY / 2.0, 0);
            this->n[index] = tVect(0.0, -1.0, 0.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = true;
            this->xMin[index] = reservoirX + edgeRadius;
            this->xMax[index] = reservoirX + wallSizeX - edgeRadius;
            this->yMin[index] = -numeric_limits<double>::max();
            this->yMax[index] = +numeric_limits<double>::max();
            this->zMin[index] = -numeric_limits<double>::max();
            this->zMax[index] = +numeric_limits<double>::max();
            ++index;


            // obstacle
            this->p[index] = tVect(reservoirX + wallSizeX + erodibleSizeX, 0.0, 0);
            this->n[index] = tVect(-1.0, 0.0, 0.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = true;
            this->xMin[index] = -numeric_limits<double>::max();
            this->xMax[index] = +numeric_limits<double>::max();
            this->yMin[index] = centerY - obstacleSizeY / 2.0 + edgeRadius;
            this->yMax[index] = centerY + obstacleSizeY / 2.0 - edgeRadius;
            this->zMin[index] = -numeric_limits<double>::max();
            this->zMax[index] = +numeric_limits<double>::max();
            ++index;

            this->p[index] = tVect(reservoirX + wallSizeX + erodibleSizeX + obstacleSizeX, 0.0, 0);
            this->n[index] = tVect(1.0, 0.0, 0.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = true;
            this->xMin[index] = -numeric_limits<double>::max();
            this->xMax[index] = +numeric_limits<double>::max();
            this->yMin[index] = centerY - obstacleSizeY / 2.0 + edgeRadius;
            this->yMax[index] = centerY + obstacleSizeY / 2.0 - edgeRadius;
            this->zMin[index] = -numeric_limits<double>::max();
            this->zMax[index] = +numeric_limits<double>::max();
            ++index;

            this->p[index] = tVect(0.0, centerY + obstacleSizeY / 2.0, 0);
            this->n[index] = tVect(0.0, 1.0, 0.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = true;
            this->xMin[index] = reservoirX + wallSizeX + erodibleSizeX + edgeRadius;
            this->xMax[index] = reservoirX + wallSizeX + erodibleSizeX + obstacleSizeX - edgeRadius;
            this->yMin[index] = -numeric_limits<double>::max();
            this->yMax[index] = +numeric_limits<double>::max();
            this->zMin[index] = -numeric_limits<double>::max();
            this->zMax[index] = +numeric_limits<double>::max();
            ++index;

            this->p[index] = tVect(0.0, centerY - obstacleSizeY / 2.0, 0);
            this->n[index] = tVect(0.0, -1.0, 0.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = true;
            this->xMin[index] = reservoirX + wallSizeX + erodibleSizeX + edgeRadius;
            this->xMax[index] = reservoirX + wallSizeX + erodibleSizeX + obstacleSizeX - edgeRadius;
            this->yMin[index] = -numeric_limits<double>::max();
            this->yMax[index] = +numeric_limits<double>::max();
            this->zMin[index] = -numeric_limits<double>::max();
            this->zMax[index] = +numeric_limits<double>::max();
            ++index;

            // bottom of reservoir
            this->p[index] = tVect(0.0, 0.0, erodibleHeight);
            this->n[index] = tVect(0.0, 0.0, 1.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = true;
            this->xMin[index] = -numeric_limits<double>::max();
            this->xMax[index] = reservoirX + wallSizeX - edgeRadius;
            this->yMin[index] = -numeric_limits<double>::max();
            this->yMax[index] = +numeric_limits<double>::max();
            this->zMin[index] = -numeric_limits<double>::max();
            this->zMax[index] = +numeric_limits<double>::max();
            ++index;

            this->p[index] = tVect(reservoirX + wallSizeX, 0.0, 0);
            this->n[index] = tVect(1.0, 0.0, 0.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = true;
            this->xMin[index] = -numeric_limits<double>::max();
            this->xMax[index] = +numeric_limits<double>::max();
            this->yMin[index] = -numeric_limits<double>::max();
            this->yMax[index] = +numeric_limits<double>::max();
            this->zMin[index] = -numeric_limits<double>::max();
            this->zMax[index] = erodibleHeight - edgeRadius;
            ++index;

            // bottom of non-erodible floodplain
            this->p[index] = tVect(0.0, 0.0, erodibleHeight);
            this->n[index] = tVect(0.0, 0.0, 1.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = true;
            this->xMin[index] = reservoirX + wallSizeX + erodibleSizeX + edgeRadius;
            this->xMax[index] = numeric_limits<double>::max();
            this->yMin[index] = -numeric_limits<double>::max();
            this->yMax[index] = +numeric_limits<double>::max();
            this->zMin[index] = -numeric_limits<double>::max();
            this->zMax[index] = +numeric_limits<double>::max();
            ++index;

            this->p[index] = tVect(reservoirX + wallSizeX + erodibleSizeX, 0.0, 0);
            this->n[index] = tVect(-1.0, 0.0, 0.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = true;
            this->xMin[index] = -numeric_limits<double>::max();
            this->xMax[index] = +numeric_limits<double>::max();
            this->yMin[index] = -numeric_limits<double>::max();
            this->yMax[index] = +numeric_limits<double>::max();
            this->zMin[index] = -numeric_limits<double>::max();
            this->zMax[index] = erodibleHeight - edgeRadius;
            ++index;

            break;
        }
        case HEAP:
        {
            // 1 extra walls
            memoryAlloc<CPU>(basic_wall_count + 1);

            const double outletSize = DEM_P.demSize[0] / 4.0;
            
            this->p[index] = tVect(0.0, 0.0, DEM_P.heapBaseLevel);
            this->n[index] = tVect(0.0, 0.0, 1.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = true;
            this->xMin[index] = outletSize;
            this->xMax[index] = DEM_P.demSize[0] - outletSize;
            ++index;

            break;
        }
        case TRIAXIAL:
        {
            // 3 extra walls
            memoryAlloc<CPU>(basic_wall_count + 3);
            
            this->p[index] = tVect(DEM_P.triBeginX, 0.0, 0.0);
            this->n[index] = tVect(-1.0, 0.0, 0.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = false;
            ++index;

            this->p[index] = tVect(0.0, DEM_P.triBeginY, 0.0);
            this->n[index] = tVect(0.0, -1.0, 0.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = false;
            ++index;

            this->p[index] = tVect(0.0, 0.0, DEM_P.triBeginZ);
            this->n[index] = tVect(0.0, 0.0, -1.0);
            this->index[index] = index;
            this->moving[index] = false;
            this->rotCenter[index].reset();
            this->omega[index].reset();
            this->vel[index].reset();
            this->translating = false;
            this->slip[index] = false;
            this->moving[index] = false;
            this->limited[index] = false;
            ++index;

            break;
        }
        case default: 
        {
            // no extra walls
            memoryAlloc<CPU>(basic_wall_count);
        }
    }
    this->count = index;
    index = 0;

    // Create basic walls
    for (unsigned int i = 0; i < externalBoundary.size(); ++i) {
        if ((externalBoundary[i] == SLIP_STAT_WALL) || (externalBoundary[i] == SLIP_DYN_WALL) || (externalBoundary[i] == STAT_WALL) || (externalBoundary[i] == DYN_WALL)) {
            this->p[index] = boundaryLocation[i]; //tVect(0.5*unit.Length,0.0,0.0);
            this->n[index] = boundaryDirec[i];
            this->index[index] = index;
            this->translating = false;
            this->trans[index].reset();
            this->limited[index] = false;
            if (externalBoundary[i] == SLIP_STAT_WALL) {
                this->moving[index] = false;
                this->slip[index] = true;
            } else if (externalBoundary[i] == SLIP_DYN_WALL) {
                this->moving[index] = true;
                this->slip[index] = true;
                this->rotCenter[index].reset();
                this->omega[index].reset();
                this->vel[index].reset();
            } else if (externalBoundary[i] == STAT_WALL) {
                this->moving[index] = false;
                this->slip[index] = false;
            } else if (externalBoundary[i] == DYN_WALL) {
                this->moving[index] = true;
                this->slip[index] = false;
                this->rotCenter[index].reset();
                this->omega[index].reset();
                this->vel[index] = tVect(0.0, 0.0, 0.0);
            }
            ++index;
        }
    }

    // Drum edge case, init extra walls after, as they depend on basic
    if (problemName == DRUM) {
        // left wall, just next to basic wall 2
        if ((externalBoundary[2] == SLIP_STAT_WALL) || (externalBoundary[2] == SLIP_DYN_WALL) || (externalBoundary[2] == STAT_WALL) || (externalBoundary[2] == DYN_WALL)) {
            // Copy wall 2
            this->n[index] = this->n[2];
            this->translating[index] = this->translating[2];
            this->trans[index] = this->trans[2];
            this->limited[index] = this->limited[2];
            this->translating[index] = this->translating[2];
            this->moving[index] = this->moving[2];
            this->slip[index] = this->slip[2];
            this->rotCenter[index] = this->rotCenter[2];
            this->omega[index] = this->omega[2];
            this->vel[index] = this->vel[2];
            this->p[index] = this->p[2] + tVect(0.6 + 0.35, 0.0001, 1.260);
            this->index[index] = index;
            if (externalBoundary[2] == SLIP_STAT_WALL) {
                this->moving[index] = false;
                this->slip[index] = true;
                this->rotCenter[index].reset();
                this->omega[index].reset();
            } else if (externalBoundary[2] == SLIP_DYN_WALL) {
                this->moving[index] = true;
                this->slip[index] = true;
                this->rotCenter[index] = this->p[index];
                this->omega[index] = tVect(0.0, DEM_P.drumSpeed, 0.0);
            } else if (externalBoundary[2] == STAT_WALL) {
                this->moving[index] = false;
                this->slip[index] = false;
                this->rotCenter[index].reset();
                this->omega[index].reset();
            } else if (externalBoundary[2] == DYN_WALL) {
                this->moving[index] = true;
                this->slip[index] = false;
                this->rotCenter[index] = this->p[index];
                this->omega[index] = tVect(0.0, DEM_P.drumSpeed, 0.0);
            }
            ++index;
        }
        // right wall, just next to basic wall 3
        if ((externalBoundary[3] == SLIP_STAT_WALL) || (externalBoundary[3] == SLIP_DYN_WALL) || (externalBoundary[3] == STAT_WALL) || (externalBoundary[3] == DYN_WALL)) {
            // Copy wall 3
            this->n[index] = this->n[3];
            this->translating[index] = this->translating[3];
            this->trans[index] = this->trans[3];
            this->limited[index] = this->limited[3];
            this->translating[index] = this->translating[3];
            this->moving[index] = this->moving[3];
            this->slip[index] = this->slip[3];
            this->rotCenter[index] = this->rotCenter[3];
            this->omega[index] = this->omega[3];
            this->vel[index] = this->vel[3];
            this->p[index] = this->p[3] + tVect(0.6 + 0.35, -0.0001, 1.260);
            this->index[index] = index;
            this->limited[index] = false;
            if (externalBoundary[3] == SLIP_STAT_WALL) {
                this->moving[index] = false;
                this->slip[index] = true;
                this->rotCenter[index].reset();
                this->omega[index].reset();
            } else if (externalBoundary[3] == SLIP_DYN_WALL) {
                this->moving[index] = true;
                this->slip[index] = true;
                this->rotCenter[index] = this->p[index];
                this->omega[index] = tVect(0.0, DEM_P.drumSpeed, 0.0);
            } else if (externalBoundary[3] == STAT_WALL) {
                this->moving[index] = false;
                this->slip[index] = false;
                this->rotCenter[index].reset();
                this->omega[index].reset();
            } else if (externalBoundary[3] == DYN_WALL) {
                this->moving[index] = true;
                this->slip[index] = false;
                this->rotCenter[index] = this->p[index];
                this->omega[index] = tVect(0.0, DEM_P.drumSpeed, 0.0);
            }
            ++index;
        }
    }

    for (unsigned int i = 0; i < this->count; ++i) {
        this->wallShow(i);
    }
}
void Wall2::wallShow(const unsigned int i) const {
    cout << "Wall number " << index[i] << ": ";
    cout << "with base point:";
    p[i].show();
    cout << " and normal vector:";
    n[i].show();
    cout << " moving:" << moving[i] << "\n";
    if (moving[i]) {
        cout << "with translational velocity:";
        vel[i].show();
        cout << "\n";
        cout << "and rotational velocity:";
        omega[i].show();
        cout << " around:";
        rotCenter[i].show();
        cout << "\n";
    }
}

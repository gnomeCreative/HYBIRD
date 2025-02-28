#include "Cylinder2.h"

#include "DEMParams.h"

extern ProblemName problemName;

void Cylinder2::initialize() {
    
    unsigned int index = 0;

    switch (problemName) {
        case HK_LARGE:
        {
            if (DEM_P.depositArea) {
                this->memoryAlloc<CPU>(1);
                this->p1[index] = tVect(21.0858, 13.3466, 0.0);
                this->p2[index] = tVect(21.0858, 13.3466, 100.0);
                this->R[index] = 13.513015;
                this->omega[index] = tVect(0.0, 0.0, 0.0);
                this->initAxes(index);
                this->type[index] = cylinderType::EMPTY;
                this->moving[index] = false;
                this->slip[index] = false;
                this->limited[index] = true;
                this->xMin[index] = 23.2;
                this->xMax[index] = 25.2795;
                this->yMin[index] = 0.0;
                this->yMax[index] = 0.5008;
                this->zMin[index] = -1.0;
                this->zMax[index] = 2.0 * DEM_P.demSize[2];
                ++index;
            }
            break;
        }
        case IERVOLINO_2D:
        {

            const double reservoirX = 0.8;
            const double wallSizeX = 0.025;
            // const double erodibleSizeX = 0.51;
            const double erodibleHeight = 0.02;
            const double edgeRadius = 0.005;

            this->memoryAlloc<CPU>(1);
            this->R[index] = edgeRadius;
            this->omega[index] = tVect(0.0, 0.0, 0.0);
            this->type[index] = FULL;
            this->moving[index] = false;
            this->slip[index] = false;
            this->limited[index] = false;


            // erodible part borders
            this->p1[index] = tVect(reservoirX + wallSizeX - edgeRadius, 0.0, erodibleHeight - edgeRadius);
            this->p2[index] = tVect(reservoirX + wallSizeX - edgeRadius, 1.0, erodibleHeight - edgeRadius);
            this->initAxes(index);
            ++index;

            break;
        }
        case IERVOLINO:
        {

            const double centerY = DEM_P.demSize[1] / 2.0;
            const double reservoirX = 0.8;
            const double wallSizeX = 0.025;
            const double outletSizeY = 0.3;
            const double obstacleSizeX = 0.155;
            const double obstacleSizeY = 0.3;
            const double erodibleSizeX = 0.51;
            const double erodibleHeight = 0.02;
            const double edgeRadius = 0.005;

            this->memoryAlloc<CPU>(10);

            // Common features
            for (unsigned int i = 0; i < 10; ++i) {
                this->R[index] = edgeRadius;
                this->omega[index] = tVect(0.0, 0.0, 0.0);
                this->type[index] = FULL;
                this->moving[index] = false;
                this->slip[index] = false;
                this->limited[index] = false;
            }

            // left wall borders
            this->p1[index] = tVect(reservoirX + edgeRadius, centerY - outletSizeY / 2.0 - edgeRadius, 0.0);
            this->p2[index] = tVect(reservoirX + edgeRadius, centerY - outletSizeY / 2.0 - edgeRadius, 1.0);
            this->initAxes(index);
            ++index;
            
            this->p1[index] = tVect(reservoirX + wallSizeX - edgeRadius, centerY - outletSizeY / 2.0 - edgeRadius, 0.0);
            this->p2[index] = tVect(reservoirX + wallSizeX - edgeRadius, centerY - outletSizeY / 2.0 - edgeRadius, 1.0);
            this->initAxes(index);
            ++index;

            // right wall borders
            this->p1[index] = tVect(reservoirX + edgeRadius, centerY + outletSizeY / 2.0 + edgeRadius, 0.0);
            this->p2[index] = tVect(reservoirX + edgeRadius, centerY + outletSizeY / 2.0 + edgeRadius, 1.0);
            this->initAxes(index);
            ++index;
            
            this->p1[index] = tVect(reservoirX + wallSizeX - edgeRadius, centerY + outletSizeY / 2.0 + edgeRadius, 0.0);
            this->p2[index] = tVect(reservoirX + wallSizeX - edgeRadius, centerY + outletSizeY / 2.0 + edgeRadius, 1.0);
            this->initAxes(index);
            ++index;

            // obstacle corners
            this->p1[index] = tVect(reservoirX + wallSizeX + erodibleSizeX + edgeRadius, centerY - obstacleSizeY / 2.0 + edgeRadius, 0.0);
            this->p2[index] = tVect(reservoirX + wallSizeX + erodibleSizeX + edgeRadius, centerY - obstacleSizeY / 2.0 + edgeRadius, 1.0);
            this->initAxes(index);
            ++index;
            
            this->p1[index] = tVect(reservoirX + wallSizeX + erodibleSizeX + obstacleSizeX - edgeRadius, centerY - obstacleSizeY / 2.0 + edgeRadius, 0.0);
            this->p2[index] = tVect(reservoirX + wallSizeX + erodibleSizeX + obstacleSizeX - edgeRadius, centerY - obstacleSizeY / 2.0 + edgeRadius, 1.0);
            this->initAxes(index);
            ++index;
            
            this->p1[index] = tVect(reservoirX + wallSizeX + erodibleSizeX + edgeRadius, centerY + obstacleSizeY / 2.0 - edgeRadius, 0.0);
            this->p2[index] = tVect(reservoirX + wallSizeX + erodibleSizeX + edgeRadius, centerY + obstacleSizeY / 2.0 - edgeRadius, 1.0);
            this->initAxes(index);
            ++index;
            
            this->p1[index] = tVect(reservoirX + wallSizeX + erodibleSizeX + obstacleSizeX - edgeRadius, centerY + obstacleSizeY / 2.0 - edgeRadius, 0.0);
            this->p2[index] = tVect(reservoirX + wallSizeX + erodibleSizeX + obstacleSizeX - edgeRadius, centerY + obstacleSizeY / 2.0 - edgeRadius, 1.0);
            this->initAxes(index);
            ++index;

            // erodible part borders
            this->p1[index] = tVect(reservoirX + wallSizeX - edgeRadius, 0.0, erodibleHeight - edgeRadius);
            this->p2[index] = tVect(reservoirX + wallSizeX - edgeRadius, 1.0, erodibleHeight - edgeRadius);
            this->initAxes(index);
            ++index;
            
            this->p1[index] = tVect(reservoirX + wallSizeX + erodibleSizeX + edgeRadius, 0.0, erodibleHeight - edgeRadius);
            this->p2[index] = tVect(reservoirX + wallSizeX + erodibleSizeX + edgeRadius, 1.0, erodibleHeight - edgeRadius);
            this->initAxes(index);
            ++index;

            break;
        }
        case IERVOLINO_CYLINDERTEST:
        {
            this->memoryAlloc<CPU>(2);
            
            this->p1[index] = tVect(1.5, 0, 1.5);
            this->p2[index] = tVect(1.5, 1.0, 1.5);
            this->R[index] = 1.0;
            this->omega[index] = tVect(0.0, 0.0, 0.0);
            this->initAxes(index);
            this->type[index] = EMPTY;
            this->moving[index] = false;
            this->slip[index] = false;
            this->limited[index] = false;
            ++index;
            
            this->p1[index] = tVect(1.5, 0, 1.5);
            this->p2[index] = tVect(1.5, 1.0, 1.5);
            this->R[index] = 0.5;
            this->omega[index] = tVect(0.0, 0.0, 0.0);
            this->initAxes(index);
            this->type[index] = FULL;
            this->moving[index] = false;
            this->slip[index] = false;
            this->limited[index] = false;
            ++index;

            break;

        }
        case WILL_SETTLING:
        {
            this->memoryAlloc<CPU>(1);
            this->p1[index] = tVect(0, 0, 0);
            this->p2[index] = tVect(0, 0, 1);
            this->R[index] = 0.108 / 2.0;
            this->omega[index] = tVect(0.0, 0.0, 0.0);
            this->initAxes(index);
            this->type[index] = EMPTY;
            this->moving[index] = false;
            this->slip[index] = false;
            this->limited[index] = false;
            ++index;

            break;

        }
        case KELVIN:
        {
            this->memoryAlloc<CPU>(1);
            this->p1[index] = tVect(1.2504, 0.4, 0.0);
            this->p2[index] = tVect(1.2504, 0.4, 100.0);
            this->R[index] = 0.400;
            this->omega[index] = tVect(0.0, 0.0, 0.0);
            this->initAxes(index);
            this->type[index] = EMPTY;
            this->moving[index] = false;
            this->slip[index] = false;
            this->limited[index] = true;
            this->xMin[index] = 1.250;
            this->xMax[index] = 1.42535;
            this->yMin[index] = -10.0;
            this->yMax[index] = 0.04048;
            this->zMin[index] = -1.0;
            this->zMax[index] = 2.0 * DEM_P.demSize[2];
            ++index;
            break;
        }
        case DRUM:
        {
            this->memoryAlloc<CPU>(1);
            this->p1[index] = tVect(0.0 + 0.6 + 0.35, 0.0, 1.26);
            this->p2[index] = tVect(0.0 + 0.6 + 0.35, 1.0, 1.26);
            this->R[index] = 1.243;
            this->omega[index] = tVect(0.0, DEM_P.drumSpeed, 0.0);
            this->initAxes(index);
            this->type[index] = EMPTY;
            this->moving[index] = true;
            this->slip[index] = false;
            this->limited[index] = false;
            ++index;
            break;
        }
        case AVALANCHE:
        {
            this->memoryAlloc<CPU>(1);
            this->p1[index] = tVect(0.0, 15, 15); //tVect(0.0,14.86,14.86);
            this->p2[index] = tVect(1.0, 15, 15);
            this->R[index] = 14.86;
            this->omega[index].reset();
            this->initAxes(index);
            this->type[index] = EMPTY;
            this->moving[index] = false;
            this->slip[index] = false;
            this->limited[index] = false;
            ++index;
            break;
        }
        case NET:
        {
            this->memoryAlloc<CPU>(1);
            this->p1[index] = tVect(0.0, 4.0, 4.0);
            this->p2[index] = tVect(1.0, 4.0, 4.0);
            this->R[index] = 4.0;
            this->omega[index].reset();
            this->initAxes(index);
            this->type[index] = EMPTY;
            this->moving[index] = false;
            this->slip[index] = false;
            this->limited[index] = false;
            ++index;
            break;
        }
        case INTRUDER:
        {
            const double intruderX = 0.02; //0.0675;   // X and Z coordinates of the intruder (in the muddle of the box X and Z)
            const double intruderZ = 0.0675;
            const double intruderd = 16e-3; // diameter of the intruder
            this->memoryAlloc<CPU>(1);
            this->p1[index] = tVect(intruderX, 0.0, intruderZ);
            this->p2[index] = tVect(intruderX, 100.0, intruderZ);
            this->R[index] = intruderd / 2.0;
            this->omega[index] = tVect(0.0, 0.0, 0.0);
            this->initAxes(index);
            this->moving[index] = false;
            this->slip[index] = false;
            this->limited[index] = false;
            this->type[index] = FULL;
            this->translating[index] = true;
            this->trans[index] = tVect(7.5, 0.0, 0.0);
            ++index;
            break;
        }
    }
    count = alloc;
    for (unsigned int i = 0; i < this->count; ++i) {
        this->cylinderShow(i);
    }
}

void Cylinder2::initAxes(const unsigned int i) {
    // initializes the axes variables (outside dist function for improved performance)

    // axes vector
    axes[i] = p2[i] - p1[i];
    // axes unit vector
    naxes[i] = axes[i] / axes[i].norm();
}

void Cylinder2::cylinderShow(const unsigned int i) const {
    if (type[i] == EMPTY) {
        cout << "Inner ";
    } else if (type[i] == FULL) {
        cout << "Outer ";
    }
    cout << "cylinder, number " << i;
    cout << " with base point: p1";
    p1[i].show();
    cout << ", p2";
    p2[i].show();
    cout << ", and radius: " << R << endl;
    if (moving[i]) {
        cout << "rotating with speed:";
        omega[i].show();
        cout << "\n";
    }
}
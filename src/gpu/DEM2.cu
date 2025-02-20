#include "DEM2.h"

#include <numeric>

void DEM2::discreteElementInit(DEMParams &dem_p, const std::array<types, 6>& externalBoundary, const std::array<double, 3>& externalSize, const std::array<tVect, 6>& externalBoundaryLocation,
                               const tVect& externalAccel, const tVect& externalRotation, const tVect& externalRotationCenter, const bool externalSolveCoriolis, const bool externalSolveCentrifugal, const double externalTimeStep) {

    // initializing DEM parameters from external parameters (LES or LBM))

    // domain size is the same of LBM. This is expressed in lattice coordinates
    dem_p.demSize[0] = externalSize[0];
    dem_p.demSize[1] = externalSize[1];
    dem_p.demSize[2] = externalSize[2];

    // switchers for apparent accelerations
    dem_p.solveCoriolis = externalSolveCoriolis;
    dem_p.solveCentrifugal = externalSolveCentrifugal;

    // acceleration field
    dem_p.demF = externalAccel;

    // roation of the reference frame
    dem_p.demRotCenter = externalRotationCenter;
    dem_p.demRot = externalRotation;

    // initializing particles
    const double partDensity = dem_p.sphereMat.density;
    
    unsigned int globalIndex = 0;
    unsigned int particleCount = 0;
    for (unsigned int i = 0; i < h_elements.count; ++i) {
        // initialize element
        h_elements.initialize(i, partDensity, dem_p);
        // Count required particles
        particleCount += h_elements.size[i];

    }
    // Allocate memory for particles
    h_particles.memoryAlloc<CPU>(particleCount);
    // Allocate memory for element components
    h_elements.allocComponentsData();

    for (unsigned int i = 0; i < h_elements.count; ++i) {
        h_elements.generateParticles(i, globalIndex, h_particles, dem_p);
    }

    // the number of standard particles (=not ghosts) is now fixed, and WILL NOT BE CHANGED
    stdPartNumber = h_particles.count;
    actvPartNumber = stdPartNumber;

    // initializing wall for DEM
    h_walls.initialize(dem_p, externalBoundary, externalBoundaryLocation);
    // initializing cylinders for DEM
    h_cylinders.initialize(dem_p);
    // initializing periodic boundary conditions for DEM
    initializePbcs(externalBoundary, externalBoundaryLocation);
    // initializing destroy planes for DEM
    //initializeDestroy(externalBoundary, externalBoundaryLocation);

    // used to get tot number of objects for CG (coarse graining)
    stdObjects = h_objects.count;

    // initialize neighbor list parameters (also generate particles and set ghosts)
    if (h_elements.count) {
        initNeighborParameters();
        periodicObjects();
        evalNeighborTable();
    }


    const double totMass = std::reduce(h_elements.m, h_elements.m + h_elements.count, 0.0);
    const double minPartRadius = *std::min_element(h_elements.radius, h_elements.radius + h_objects.count);
    const double maxPartRadius = *std::max_element(h_elements.radius, h_elements.radius + h_objects.count);
    const double meanPartRadius = h_elements.count ? std::reduce(h_elements.radius, h_elements.radius + h_elements.count) / h_elements.count : 0;

    const double minObjRadius = *std::min_element(h_objects.r, h_objects.r + h_objects.count);
    const double maxObjRadius = *std::max_element(h_objects.r, h_objects.r + h_objects.count);
    const double meanObjRadius = h_objects.count ? std::reduce(h_objects.r, h_objects.r + h_objects.count) / h_objects.count : 0;

    // DEM time step
    // if multistep is 0, it should be calculated by the program here
    determineTimeStep(externalTimeStep);

    cout << "DEM parameters\n";
    cout << "domain size: xdim =" << dem_p.demSize[0] << "; ydim= " << dem_p.demSize[1] << "; zdim= " << dem_p.demSize[2] << ";" << endl;
    cout << "Tot elements: " << h_elements.count << ";\t";
    cout << "Tot objects: " << h_objects.count << ";\t";
    cout << "Object radius: Mean=" << meanObjRadius << " Min=" << minObjRadius << " Max=" << maxObjRadius << ";" << endl;
    cout << "Tot standard particles: " << stdPartNumber << endl;
    cout << "Particle radius: Mean=" << meanPartRadius << " Min=" << minPartRadius << " Max=" << maxPartRadius << ";" << endl;
    cout << "Total particle mass=" << totMass << endl;
    if (dem_p.multiStep > 0) {
        cout << "Deltat =" << deltat << ", therefore " << multiStep << " multiple steps, as a " << criticalRatio << ":1 ratio to estimated collision time" << endl;
    } else {
        cout << "Deltat =" << deltat << ", by imposing " << multiStep << " substeps" << endl;
    }

    switch (dem_p.sphereMat.contactModel) {
        case LINEAR:
        {
            cout << "Contact model: linear dashpot" << endl;
            cout << "Normal stiffness = " << dem_p.sphereMat.linearStiff << endl;
            break;
        }
        case HERTZIAN:
        {
            cout << "Contact model: damped Hertzian contact" << endl;
            break;
        }
    }
    cout << "Damping ratio = " << dem_p.sphereMat.dampCoeff << ", equivalent to a coefficient of restitution of c=" << sphereMat.restitution << endl;
    cout << "Tangential viscosity = " << dem_p.sphereMat.viscTang << endl;
    cout << "Particle-particle friction = " << dem_p.sphereMat.frictionCoefPart
            << ", wall-particle friction = " << dem_p.sphereMat.frictionCoefWall
            << ", object-particle friction = " << dem_p.sphereMat.frictionCoefObj << endl;
    cout << "Rolling coefficient = " << dem_p.sphereMat.rollingCoefPart << endl;
    cout << "Numerical viscosity =" << dem_p.numVisc << endl;

}

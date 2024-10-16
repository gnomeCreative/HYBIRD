
#include "DEM.h"
#include "macros.h"


#include "Problem.h"

#define USE_MATH_DEFINES

using namespace std;

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
 // PUBLIC FUNCTIONS
 //////////////////////////////////////////////////////////////////////////////////////////////////////*/

void DEM::discreteElementGet(GetPot& configFile, GetPot& commandLine) {
    // getting material properties
    PARSE_CLASS_MEMBER(configFile, sphereMat.density, "particleDensity", 0.0);
    ASSERT(sphereMat.density > 0.0);

    string contactModelString;
    PARSE_CLASS_MEMBER(configFile, contactModelString, "contactModel", "none");
    if (contactModelString == "HERTZIAN") sphereMat.contactModel = HERTZIAN;
    else if (contactModelString == "LINEAR") sphereMat.contactModel = LINEAR;


    // Hertzian contact model /////////////
    PARSE_CLASS_MEMBER(configFile, sphereMat.youngMod, "youngMod", 1.0);
    ASSERT(sphereMat.youngMod > 0.0);
    PARSE_CLASS_MEMBER(configFile, sphereMat.poisson, "poisson", 0.3);
    ASSERT(sphereMat.poisson >= 0.0);
    // normal stiffness (constant part, calculated here to avoid repetition)
    sphereMat.knConst = 2.0 / 3.0 * sphereMat.youngMod / (1.0 - sphereMat.poisson * sphereMat.poisson);
    // there was a bracket mistake here
    sphereMat.ksConst = 2.0 * sphereMat.youngMod / (2.0 - sphereMat.poisson) / (1.0 + sphereMat.poisson);

    // linear contact model /////////////
    PARSE_CLASS_MEMBER(configFile, sphereMat.linearStiff, "linearStiff", 1.0);
    ASSERT(sphereMat.linearStiff >= 0.0);

    // normal damping ///////////////////////
    PARSE_CLASS_MEMBER(configFile, sphereMat.restitution, "restitution", 1.0);
    ASSERT(sphereMat.restitution > 0.0);
    ASSERT(sphereMat.restitution <= 1.0);
    // calculating coefficient for normal damping
    switch (sphereMat.contactModel) {
        case HERTZIAN:
        {
            // see "On the Determination of the Damping Coefficient of Non-linear Spring-dashpot System to Model Hertz Contact for Simulation by Discrete Element Method"
            // Hu, Hu, Jian, Liu, Wan, Journal of Computers, 6 (2011) OR BETTER Antypov & Elliott
            sphereMat.dampCoeff = -1.0 * sqrt(5) * log(sphereMat.restitution) / sqrt((log(sphereMat.restitution) * log(sphereMat.restitution) + M_PI * M_PI));
            break;
        }
        case LINEAR:
        {
            //sphereMat.dampCoeff=-1.0*sqrt(2.0)*log(sphereMat.restitution)/sqrt((log(sphereMat.restitution)*log(sphereMat.restitution)+M_PI));
            sphereMat.dampCoeff = -1.0 * log(sphereMat.restitution) / sqrt((log(sphereMat.restitution) * log(sphereMat.restitution) + M_PI * M_PI));
            ASSERT(sphereMat.dampCoeff < 1.0);
            break;
        }
    }

    // tangential model //////////////////////
    PARSE_CLASS_MEMBER(configFile, sphereMat.viscTang, "viscTang", 0.0);
    ASSERT(sphereMat.viscTang >= 0.0);
    PARSE_CLASS_MEMBER(configFile, sphereMat.frictionCoefPart, "frictionCoefPart", 0.0);
    ASSERT(sphereMat.frictionCoefPart >= 0.0);
    PARSE_CLASS_MEMBER(configFile, sphereMat.frictionCoefWall, "frictionCoefWall", 0.0);
    ASSERT(sphereMat.frictionCoefWall >= 0.0);
    PARSE_CLASS_MEMBER(configFile, sphereMat.frictionCoefObj, "frictionCoefObj", 0.0);
    if (sphereMat.frictionCoefObj == 0.0) {
        sphereMat.frictionCoefObj = sphereMat.frictionCoefWall;
    }
    // rolling model //////////////////////
    PARSE_CLASS_MEMBER(configFile, sphereMat.rollingCoefPart, "rollingCoefPart", 0.0);
    ASSERT(sphereMat.rollingCoefPart >= 0.0);

    // particle initial state //////////////////////
    string particleFile;
    PARSE_CLASS_MEMBER(configFile, particleFile, "particleFile", "particles.dat");
    double translateX(0.0), translateY(0.0), translateZ(0.0);
    PARSE_CLASS_MEMBER(configFile, translateX, "particleTranslateX", 0.0);
    PARSE_CLASS_MEMBER(configFile, translateY, "particleTranslateY", 0.0);
    PARSE_CLASS_MEMBER(configFile, translateZ, "particleTranslateZ", 0.0);
    tVect translate(translateX, translateY, translateZ);
    double scale = 1.0;
    PARSE_CLASS_MEMBER(configFile, scale, "particleScale", 1.0);

    ifstream particleFileID;
    particleFileID.open(particleFile.c_str(), ios::in);
    cout << "Reading " << particleFile.c_str() << "...";
    ASSERT(particleFileID.is_open());
    unsigned int totElmt;
    particleFileID>>totElmt;
    for (int n = 0; n < totElmt; ++n) {
        elmt dummyElmt;

        // import variables
        particleFileID >> dummyElmt.index;
        particleFileID >> dummyElmt.size;
        particleFileID >> dummyElmt.radius;
        dummyElmt.radius = dummyElmt.radius*scale;
        // position
        double x0, y0, z0;
        particleFileID>>x0;
        particleFileID>>y0;
        particleFileID>>z0;
        dummyElmt.x0 = tVect(x0*scale, y0*scale, z0 * scale) + translate;
        // translational velocity
        double x1, y1, z1;
        particleFileID>>x1;
        particleFileID>>y1;
        particleFileID>>z1;
        dummyElmt.x1 = tVect(x1, y1, z1); //.reset();
        // rotational velocity
        double w01, w02, w03;
        particleFileID>>w01;
        particleFileID>>w02;
        particleFileID>>w03;
        dummyElmt.w0 = tVect(w01, w02, w03); //.reset();
        //dummyElmt.x1.reset();
        // orientation
        double p0, q0, r0, s0;
        particleFileID>>p0;
        particleFileID>>q0;
        particleFileID>>r0;
        particleFileID>>s0;
        dummyElmt.q0 = tQuat(p0, q0, r0, s0);
        //dummyElmt.q0.resetSoft();
        // translational velocity (in quaternion rates))
        double p1, q1, r1, s1;
        particleFileID>>p1;
        particleFileID>>q1;
        particleFileID>>r1;
        particleFileID>>s1;
        dummyElmt.q1 = tQuat(p1, q1, r1, s1);
        //dummyElmt.q1.resetHard();
        dummyElmt.active = true;
        // add to list
        elmts.push_back(dummyElmt);

    }
    cout << " done" << endl;

    // objects initial state //////////////////////
    string objectFile;
    PARSE_CLASS_MEMBER(configFile, objectFile, "objectFile", "objects.dat");
    ifstream objectFileID;
    objectFileID.open(objectFile.c_str(), ios::in);
    ASSERT(objectFileID.is_open());
    cout << "Reading " << objectFile.c_str() << "...";
    unsigned int totObjects;
    objectFileID>>totObjects;

    for (int n = 0; n < totObjects; ++n) {
        object dummyObject;

        // import variables
        objectFileID >> dummyObject.index;
        // this is used to identify objects belonging to different groups
        objectFileID >> dummyObject.ElID; // must be one
        objectFileID >> dummyObject.r;
        double x0, y0, z0;
        objectFileID>>x0;
        objectFileID>>y0;
        objectFileID>>z0;
        dummyObject.x0 = tVect(x0, y0, z0);
        double x1, y1, z1;
        objectFileID>>x1;
        objectFileID>>y1;
        objectFileID>>z1;
        dummyObject.x1 = tVect(x1, y1, z1);
        // the next eight values are for rotation, and are not used
        double trash;
        objectFileID>>trash;
        objectFileID>>trash;
        objectFileID>>trash;
        objectFileID>>trash;
        objectFileID>>trash;
        objectFileID>>trash;
        objectFileID>>trash;
        objectFileID>>trash;
        objectFileID>>trash;
        objectFileID>>trash;
        objectFileID>>trash;
        dummyObject.originalIndex = dummyObject.index;
        objects.push_back(dummyObject);
    }
    cout << " done" << endl;

    // numerical viscosity for stability
    PARSE_CLASS_MEMBER(configFile, numVisc, "numVisc", 0.0);
    // set multiplication number (sets the number of DEM steps between two fluid steps)
    PARSE_CLASS_MEMBER(configFile, multiStep, "multiStep", 1);
    // set ratio between time step and estimated duration of contacts (only if multiStep=0)
    PARSE_CLASS_MEMBER(configFile, criticalRatio, "criticalRatio", 0.1);
    
    // problem specific variables
    // restart
    PARSE_CLASS_MEMBER(configFile, movingCylinder, "movingCylinder", false);
    PARSE_CLASS_MEMBER(configFile, velocityCylinderX, "velocityCylinderX", 0.0);
    PARSE_CLASS_MEMBER(configFile, velocityCylinderY, "velocityCylinderY", 0.0);
    PARSE_CLASS_MEMBER(configFile, velocityCylinderZ, "velocityCylinderZ", 0.0);

}

void DEM::discreteElementInit(const Problem &problem, const typeList& externalBoundary, const doubleList& externalSize, const vecList& externalBoundaryLocation,
        const tVect externalAccel, const tVect externalRotation, const tVect externalRotationCenter, const bool externalSolveCoriolis, const bool externalSolveCentrifugal, const double& externalTimeStep) {

    // initializing DEM parameters from external parameters (LES or LBM))

    // domain size is the same of LBM. This is expressed in lattice coordinates
    demSize[0] = externalSize[0];
    demSize[1] = externalSize[1];
    demSize[2] = externalSize[2];

    // switchers for apparent accelerations
    solveCoriolis = externalSolveCoriolis;
    solveCentrifugal = externalSolveCentrifugal;

    // acceleration field
    demF = externalAccel;

    // roation of the reference frame
    demRotCenter = externalRotationCenter;
    demRot = externalRotation;

    // initializing particles
    const double partDensity = sphereMat.density;

    // initializing composite particle properties
    compositeProperties();
    // clear particle list
    particles.clear();

    unsigned int globalIndex = 0;
    for (int n = 0; n < elmts.size(); ++n) {
        // initialize element
        elmts[n].initialize(partDensity, prototypes, demF, demSize);
        // generate particles
        elmts[n].generateParticles(globalIndex, particles, prototypes);
    }

    // the number of standard particles (=not ghosts) is now fixed, and WILL NOT BE CHANGED
    stdPartNumber = particles.size();
    actvPartNumber = stdPartNumber;

    // initializing wall for DEM
    initializeWalls(problem, externalBoundary, externalBoundaryLocation);
    // initializing cylinders for DEM
    initializeCylinders(problem);
    // initializing periodic boundary conditions for DEM
    initializePbcs(externalBoundary, externalBoundaryLocation);
    // initializing destroy planes for DEM
    //initializeDestroy(externalBoundary, externalBoundaryLocation);

    // used to get tot number of objects for CG (coarse graining)
    stdObjects = objects.size();

    // initialize neighbor list parameters (also generate particles and set ghosts)
    if (elmts.size()) {
        initNeighborParameters();
        periodicObjects();
        evalNeighborTable();
    }


    double totMass(0.0);
    minPartRadius = 1.0e99;
    maxPartRadius = 0.0;
    meanPartRadius = 0.0;
    for (int n = 0; n < elmts.size(); ++n) {
        const double radiusHere = elmts[n].radius;
        // calculate mass
        totMass += elmts[n].m;
        meanPartRadius += radiusHere;
        if (radiusHere > maxPartRadius) {
            maxPartRadius = radiusHere;
        }
        if (radiusHere < minPartRadius) {
            minPartRadius = radiusHere;
        }
    }
    if (elmts.size() > 0.0) {
        meanPartRadius /= double(elmts.size());
    }

    minObjRadius = 1.0e99;
    maxObjRadius = 0.0;
    meanObjRadius = 0.0;
    for (int o = 0; o < objects.size(); ++o) {
        const double radiusHere = objects[o].r;
        // calculate mass
        meanObjRadius += radiusHere;
        if (radiusHere > maxObjRadius) {
            maxObjRadius = radiusHere;
        }
        if (radiusHere < minObjRadius) {
            minObjRadius = radiusHere;
        }
    }
    if (objects.size() > 0) {
        meanObjRadius /= double(objects.size());
    }

    // DEM time step
    // if multistep is 0, it should be calculated by the program here
    determineTimeStep(externalTimeStep);

    cout << "DEM parameters\n";
    cout << "domain size: xdim =" << demSize[0] << "; ydim= " << demSize[1] << "; zdim= " << demSize[2] << ";" << endl;
    cout << "Tot elements: " << elmts.size() << ";\t";
    cout << "Tot objects: " << objects.size() << ";\t";
    cout << "Object radius: Mean=" << meanObjRadius << " Min=" << minObjRadius << " Max=" << maxObjRadius << ";" << endl;
    cout << "Tot standard particles: " << stdPartNumber << endl;
    cout << "Particle radius: Mean=" << meanPartRadius << " Min=" << minPartRadius << " Max=" << maxPartRadius << ";" << endl;
    cout << "Total particle mass=" << totMass << endl;
    if (multiStep > 0) {
        cout << "Deltat =" << deltat << ", therefore " << multiStep << " multiple steps, as a " << criticalRatio << ":1 ratio to estimated collision time" << endl;
    } else {
        cout << "Deltat =" << deltat << ", by imposing " << multiStep << " substeps" << endl;
    }

    switch (sphereMat.contactModel) {
        case LINEAR:
        {
            cout << "Contact model: linear dashpot" << endl;
            cout << "Normal stiffness = " << sphereMat.linearStiff << endl;
            break;
        }
        case HERTZIAN:
        {
            cout << "Contact model: damped Hertzian contact" << endl;
            break;
        }
    }
    cout << "Damping ratio = " << sphereMat.dampCoeff << ", equivalent to a coefficient of restitution of c=" << sphereMat.restitution << endl;
    cout << "Tangential viscosity = " << sphereMat.viscTang << endl;
    cout << "Particle-particle friction = " << sphereMat.frictionCoefPart
            << ", wall-particle friction = " << sphereMat.frictionCoefWall
            << ", object-particle friction = " << sphereMat.frictionCoefObj << endl;
    cout << "Rolling coefficient = " << sphereMat.rollingCoefPart << endl;
    cout << "Numerical viscosity =" << numVisc << endl;

}

void DEM::determineTimeStep(const double& externalTimeStep) {

    // if there are particles, compute a deltat
    if (elmts.size()) {
        // if multistep is 0, it should be calculated by the program here
        if (multiStep == 0) {
            // find critical deltaT
            const double crit = criticalRatio * criticalTimeStep();
            // if the critical time is bigger than the external time step, then just use the LBM time step
            if (crit >= externalTimeStep) {
                multiStep = 1;
                deltat = externalTimeStep;
            }// if it is lower, calculate the number of substeps
            else {
                const double ratio = externalTimeStep / crit;
                ASSERT(ratio >= 1);
                multiStep = std::floor(ratio) + 1;
                deltat = externalTimeStep / (double(multiStep));
            }
        }// multistep can also be imposed by the user
        else {
            deltat = externalTimeStep / (double(multiStep));
        }
    }// otherwise just take the fluid time step
    else {
        deltat = externalTimeStep;
        multiStep = 1;
    }
}

// solver DEM
void DEM::discreteElementStep() {

    // set trigger for new neighbor list
    static const double neighListTrigger = 0.25 * nebrRange;

    for (int demIter = 0; demIter < multiStep; ++demIter) {

        demTimeStep++;
        demTime += deltat;
        // neighbor management
        //        cout<<"1.0 ("<<demIter<<") "<<"a="<<activeElmts[214]<<endl;
        evalMaxDisp();
        //        cout<<"1.1 ("<<demIter<<") "<<"a="<<activeElmts[214]<<endl;
        if (maxDisp > neighListTrigger) {
            maxDisp = 0.0;
            //cout<<"new neighbor list"<<endl;
            evalNeighborTable();
        }
        //cout<<"1.2"<<endl;
        //        cout<<"1.2 ("<<demIter<<") "<<"a="<<activeElmts[214]<<endl;
        // predictor step
        predictor();
        //        cout<<"1.3 ("<<demIter<<") "<<"a="<<activeElmts[214]<<endl;
        // particles generation
        updateParticlesPredicted();

        // force evaluation
        evaluateForces();

        // corrector step
        corrector();

        // particles re-generation
        updateParticlesCorrected();

    }

}

void DEM::evolveBoundaries() {

    

    switch (problemName) {
        
        case TRIAXIAL:
        {
            double frac = 1.0 / demInitialRepeat / ((double) multiStep);
            const double reduceFactor = 10.0 / (deltat * double(multiStep));

            // wall in X
            {
                const double sideY = walls[7].p.dot(Yp) - walls[2].p.dot(Yp);
                const double sideZ = walls[8].p.dot(Zp) - walls[4].p.dot(Zp);
                const double area = sideY*sideZ;

                const double expectedForce = triIsopressure*area;
                const double measuredForce = walls[6].FParticle.dot(-1.0 * walls[6].n);
                pressureX = measuredForce / area;
                const double deltaForce = expectedForce - measuredForce;

                const double displacement = deltaForce / sphereMat.linearStiff / reduceFactor;

                walls[6].p = walls[6].p + displacement * walls[6].n;
                //cout << "X: deltaForce " << deltaForce << " displacement " << displacement << " position "<< walls[6].p.dot(Xp)<< endl;

            }
            {
                // wall in Y
                const double sideZ = walls[8].p.dot(Zp) - walls[4].p.dot(Zp);
                const double sideX = walls[6].p.dot(Xp) - walls[0].p.dot(Xp);
                const double area = sideX*sideZ;

                const double expectedForce = triIsopressure*area;
                const double measuredForce = walls[7].FParticle.dot(-1.0 * walls[7].n);
                pressureY = measuredForce / area;
                const double deltaForce = expectedForce - measuredForce;

                const double displacement = deltaForce / sphereMat.linearStiff / reduceFactor;

                walls[7].p = walls[7].p + displacement * walls[7].n;
                //cout << "Y: deltaForce " << deltaForce << " displacement " << displacement << endl;
            }
            {
                // wall in Z
                const double sideX = walls[6].p.dot(Xp) - walls[0].p.dot(Xp);
                const double sideY = walls[7].p.dot(Yp) - walls[2].p.dot(Yp);
                const double area = sideX*sideY;

                const double expectedForce = triIsopressure*area;
                const double measuredForce = walls[8].FParticle.dot(-1.0 * walls[8].n);
                pressureZ = measuredForce / area;


                if (triDefSpeed == 0.0) {
                    const double deltaForce = expectedForce - measuredForce;
                    const double displacement = deltaForce / sphereMat.linearStiff / reduceFactor;
                    walls[8].p = walls[8].p + displacement * walls[8].n;
                } else {
                    walls[8].p = walls[8].p + triDefSpeed * deltat * double(multiStep) * walls[8].n;
                }
                //cout << "Z: deltaForce " << deltaForce << " displacement " << displacement << endl;
            }

            break;
        }
    }


    // update the position of the walls (rigid body)
    for (int n = 0; n < walls.size(); ++n) {
        if (walls[n].translating) {
            walls[n].p = walls[n].p + deltat * walls[n].trans;
        }
    }
    // update the position of the cylinders (rigid body)
    for (int c = 0; c < cylinders.size(); c++) {
        if (cylinders[c].translating == true) {
            // update p1 and p2
            cylinders[c].p1 = cylinders[c].p1 + deltat * cylinders[c].trans;
            cylinders[c].p2 = cylinders[c].p2 + deltat * cylinders[c].trans;

        }
    }
    
    // update the position of the object (rigid body)
    for (unsigned int o = 0; o < objects.size(); o++) {
        // Check if the velocity vector is != 0
        double objTrans = objects[o].x1.x + objects[o].x1.y + objects[o].x1.z;
        // update object position
        if (objects[o].translating == true && objTrans != 0.0) {
            objects[o].x0 = objects[o].x0 + deltat * objects[o].x1;
        }

    }

    //evalNeighborTable();
}

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
 // PRIVATE FUNCTIONS
 //////////////////////////////////////////////////////////////////////////////////////////////////////*/

// initialization functions

void DEM::compositeProperties() {

    vecList prototype1, prototype2, prototype3, prototype4, prototype5;


    // prototypes for shapes
    // every vector defines the position of a particle in the object reference frame
    // unit is radius

    prototypes.resize(5);
    prototype1.resize(1);
    prototype1[0].reset();
    prototypes[1] = prototype1;
    prototype2.resize(2);
    prototype2[0] = tVect(0.5, 0.0, 0.0);
    prototype2[1] = tVect(-0.5, 0.0, 0.0);
    prototypes[2] = prototype2;
    prototype3.resize(3);
    prototype3[0] = tVect(0.0, 1.0, 0.0);
    prototype3[1] = tVect(-sqrt(3) / 2, -1 / 2, 0.0);
    prototype3[2] = tVect(sqrt(3) / 2, -1 / 2, 0.0);
    prototypes[3] = prototype3;
    prototype4.resize(4);
    prototype4[0] = tVect(0.0, 0.0, 1.0);
    prototype4[1] = tVect(0.0, 2.0 * sqrt(2) / 3.0, -1.0 / 3.0);
    prototype4[2] = tVect(2.0 * sqrt(6) / 6.0, -2.0 * sqrt(2) / 6.0, -1.0 / 3.0);
    prototype4[3] = tVect(-2.0 * sqrt(6) / 6.0, -2.0 * sqrt(2) / 6.0, -1.0 / 3.0);
    prototypes[4] = prototype4;

    // prototype for immersed cylinders (submarine)
    prototypes[5] = prototype5;
    prototype5.resize(5);
    prototype5[0] = tVect(0.0, 0.0, 0.0);
    prototype5[1] = tVect(1.0, 0.0, 0.0);
    prototype5[2] = tVect(2.00, 0.0, 0.0);
    prototype5[3] = tVect(3.0, 0.0, 0.0);
    prototype5[4] = tVect(4.0, 0.0, 0.0);
    prototypes[5] = prototype5;

}

void DEM::initializeWalls(const Problem &problem, const typeList& externalBoundary, const vecList& boundaryLocation) {
    walls.clear();

    // boundary directors
    vecList boundaryDirec;
    boundaryDirec.resize(6);
    boundaryDirec[0] = tVect(1.0, 0.0, 0.0);
    boundaryDirec[1] = tVect(-1.0, 0.0, 0.0);
    boundaryDirec[2] = tVect(0.0, 1.0, 0.0);
    boundaryDirec[3] = tVect(0, 0. - 1.0, 0.0);
    boundaryDirec[4] = tVect(0.0, 0.0, 1.0);
    boundaryDirec[5] = tVect(0.0, 0.0, -1.0);

    unsigned int index = 0;
    // basic walls
    for (int i = 0; i < externalBoundary.size(); ++i) {
        if ((externalBoundary[i] == SLIP_STAT_WALL) || (externalBoundary[i] == SLIP_DYN_WALL) || (externalBoundary[i] == STAT_WALL) || (externalBoundary[i] == DYN_WALL)) {
            wall dummyWall;
            dummyWall.p = boundaryLocation[i]; //tVect(0.5*unit.Length,0.0,0.0);
            dummyWall.n = boundaryDirec[i];
            dummyWall.index = index;
            dummyWall.translating = false;
            dummyWall.trans.reset();
            dummyWall.limited = false;
            if (externalBoundary[i] == SLIP_STAT_WALL) {
                dummyWall.moving = false;
                dummyWall.slip = true;
            } else if (externalBoundary[i] == SLIP_DYN_WALL) {
                dummyWall.moving = true;
                dummyWall.slip = true;
                dummyWall.rotCenter.reset();
                dummyWall.omega.reset();
                dummyWall.vel.reset();
            } else if (externalBoundary[i] == STAT_WALL) {
                dummyWall.moving = false;
                dummyWall.slip = false;
            } else if (externalBoundary[i] == DYN_WALL) {
                dummyWall.moving = true;
                dummyWall.slip = false;
                dummyWall.rotCenter.reset();
                dummyWall.omega.reset();
                dummyWall.vel = tVect(0.0, 0.0, 0.0);
            }
            ++index;
            walls.push_back(dummyWall);
        }
    }

    // additional walls
    if (!problem.file.empty()) {  // Problem file was loaded
        // Copy additional walls directly from problem file, and correct their indices
        for (size_t i = 0; i < problem.walls.size(); ++i) {
            wall t = problem.walls[i];
            t.index = static_cast<unsigned int>(walls.size());
            walls.push_back(t);
        }
    } else {
        switch (problemName) {
        case HK_LARGE:
        {
            // storage container
            wall dummyWall;
            dummyWall.p = tVect(10.0, 0.0, 0.0);
            dummyWall.n = tVect(0.17365, 0.98481, 0.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.limited = false;
            ++index;
            walls.push_back(dummyWall);
            if (depositArea) {
                // deposit area
                dummyWall.p = tVect(25.2795, 0.5008, 0.0);
                dummyWall.n = tVect(-0.342020, 0.939693, 0.0);
                dummyWall.index = index;
                dummyWall.moving = false;
                dummyWall.rotCenter.reset();
                dummyWall.omega.reset();
                dummyWall.vel.reset();
                dummyWall.translating = false;
                dummyWall.limited = false;
                ++index;
                walls.push_back(dummyWall);
            }
            break;
        }
        case KELVIN:
        {
            // ramp
            wall dummyWall;
            dummyWall.p = tVect(1.42535, 0.04048, 0.0);
            dummyWall.n = tVect(-0.43837, 0.89879, 0.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.limited = true;
            dummyWall.xMin = 1.42535;
            dummyWall.xMax = 1.65005;
            dummyWall.yMin = -1.0;
            dummyWall.yMax = 0.15008;
            ++index;
            walls.push_back(dummyWall);
            // ramp
            dummyWall.p = tVect(1.65005, 0.15008, 0.0);
            dummyWall.n = tVect(0.89879, 0.43837, 0.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.limited = true;
            dummyWall.xMin = 1.65005;
            dummyWall.xMax = 1.72324;
            dummyWall.yMin = -1.0;
            dummyWall.yMax = 0.15008;
            ++index;
            walls.push_back(dummyWall);
            break;
        }
        case AVALANCHE:
        {
            wall dummyWall;
            dummyWall.p = tVect(60, 15, 15 - 14.86);
            dummyWall.n = tVect(0.42262, 0.0, 0.9063);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.limited = false;
            ++index;
            walls.push_back(dummyWall);
            break;
        }
        case DRUM:
        {
            // left wall, just next to basic wall 2
            if ((externalBoundary[2] == SLIP_STAT_WALL) || (externalBoundary[2] == SLIP_DYN_WALL) || (externalBoundary[2] == STAT_WALL) || (externalBoundary[2] == DYN_WALL)) {
                wall dummyWall = walls[2];
                dummyWall.p = dummyWall.p + tVect(0.6 + 0.35, 0.0001, 1.260);
                dummyWall.index = index;
                if (externalBoundary[2] == SLIP_STAT_WALL) {
                    dummyWall.moving = false;
                    dummyWall.slip = true;
                    dummyWall.rotCenter.reset();
                    dummyWall.omega.reset();
                }
                else if (externalBoundary[2] == SLIP_DYN_WALL) {
                    dummyWall.moving = true;
                    dummyWall.slip = true;
                    dummyWall.rotCenter = dummyWall.p;
                    dummyWall.omega = tVect(0.0, drumSpeed, 0.0);
                }
                else if (externalBoundary[2] == STAT_WALL) {
                    dummyWall.moving = false;
                    dummyWall.slip = false;
                    dummyWall.rotCenter.reset();
                    dummyWall.omega.reset();
                }
                else if (externalBoundary[2] == DYN_WALL) {
                    dummyWall.moving = true;
                    dummyWall.slip = false;
                    dummyWall.rotCenter = dummyWall.p;
                    dummyWall.omega = tVect(0.0, drumSpeed, 0.0);
                }
                ++index;
                walls.push_back(dummyWall);
            }
            // right wall, just next to basic wall 3
            if ((externalBoundary[3] == SLIP_STAT_WALL) || (externalBoundary[3] == SLIP_DYN_WALL) || (externalBoundary[3] == STAT_WALL) || (externalBoundary[3] == DYN_WALL)) {
                wall dummyWall = walls[3];
                dummyWall.p = dummyWall.p + tVect(0.6 + 0.35, -0.0001, 1.260);
                dummyWall.index = index;
                dummyWall.limited = false;
                if (externalBoundary[3] == SLIP_STAT_WALL) {
                    dummyWall.moving = false;
                    dummyWall.slip = true;
                    dummyWall.rotCenter.reset();
                    dummyWall.omega.reset();
                }
                else if (externalBoundary[3] == SLIP_DYN_WALL) {
                    dummyWall.moving = true;
                    dummyWall.slip = true;
                    dummyWall.rotCenter = dummyWall.p;
                    dummyWall.omega = tVect(0.0, drumSpeed, 0.0);
                }
                else if (externalBoundary[3] == STAT_WALL) {
                    dummyWall.moving = false;
                    dummyWall.slip = false;
                    dummyWall.rotCenter.reset();
                    dummyWall.omega.reset();
                }
                else if (externalBoundary[3] == DYN_WALL) {
                    dummyWall.moving = true;
                    dummyWall.slip = false;
                    dummyWall.rotCenter = dummyWall.p;
                    dummyWall.omega = tVect(0.0, drumSpeed, 0.0);
                }
                ++index;
                walls.push_back(dummyWall);
            }
            break;
        }
        case ZHOU:
        {
            wall dummyWall;
            dummyWall.p = tVect(0.053, 0.01, 0.21);
            dummyWall.n = tVect(0.0, 0.0, 1.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.limited = true;
            dummyWall.xMin = 0.0525;
            dummyWall.xMax = 0.3675;
            ++index;
            walls.push_back(dummyWall);
            break;
        }
        case OPENBARRIER:
        {
            wall dummyWall;
            dummyWall.p = tVect(0.0, 0.0375, 0.0);
            dummyWall.n = tVect(0.0, 1.0, 0.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.limited = false;
            ++index;
            walls.push_back(dummyWall);

            dummyWall.p = tVect(0.0, 2.0625, 0.0);
            dummyWall.n = tVect(0.0, -1.0, 0.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.limited = false;
            ++index;
            walls.push_back(dummyWall);

            break;
        }
        case HONGKONG:
        {
            const double width = 0.2;

            if (hongkongSlitSize < 0.9 * width && hongkongSmoothWall) {
                const double edgeRadius = 0.005;
                const double xBarrierPosition = 1.4;
                const double yBarrierSize = (width - hongkongSlitSize) / 2.0;
                const double barrierLeftWing = yBarrierSize;
                const double barrierRightWing = yBarrierSize + hongkongSlitSize;

                wall dummyWall;
                dummyWall.p = tVect(xBarrierPosition - edgeRadius * 0.95, 0.0, 0.0);
                dummyWall.n = tVect(-1.0, 0.0, 0.0);
                dummyWall.index = index;
                dummyWall.moving = false;
                dummyWall.rotCenter.reset();
                dummyWall.omega.reset();
                dummyWall.vel.reset();
                dummyWall.translating = false;
                dummyWall.limited = true;
                dummyWall.yMin = -100.0;
                dummyWall.yMax = barrierLeftWing - edgeRadius;
                ++index;
                walls.push_back(dummyWall);

                dummyWall.yMin = barrierRightWing + edgeRadius;
                dummyWall.yMax = 100.0;
                ++index;
                walls.push_back(dummyWall);

                dummyWall.p = tVect(xBarrierPosition + edgeRadius * 0.95, 0.0, 0.0);
                dummyWall.n = tVect(1.0, 0.0, 0.0);
                dummyWall.yMin = -100.0;
                dummyWall.yMax = barrierLeftWing - edgeRadius;
                ++index;
                walls.push_back(dummyWall);

                dummyWall.yMin = barrierRightWing + edgeRadius;
                dummyWall.yMax = 100.0;
                ++index;
                walls.push_back(dummyWall);
            }

            break;
        }
        case HOURGLASS:
        {
            double a = 0;

            const double outletCenter = demSize[0] / 2.0;

            wall dummyWall;
            dummyWall.p = tVect(0.0, 0.0, hourglassOutletHeight);
            dummyWall.n = tVect(0.0, 0.0, 1.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.slip = false;
            dummyWall.moving = false;
            dummyWall.limited = true;
            dummyWall.xMin = -numeric_limits<double>::max();
            dummyWall.xMax = outletCenter - hourglassOutletSize;
            ++index;
            walls.push_back(dummyWall);

            dummyWall.xMin = outletCenter + hourglassOutletSize;
            dummyWall.xMax = numeric_limits<double>::max();
            ++index;
            walls.push_back(dummyWall);

            break;
        }
        case IERVOLINO_2D:
        {
            const double reservoirX = 0.825;
            const double erodibleSizeX = 0.51;
            const double erodibleHeight = 0.02;
            const double edgeRadius = 0.005;

            wall dummyWall;

            // obstacle
            dummyWall.p = tVect(reservoirX + erodibleSizeX, 0.0, 0);
            dummyWall.n = tVect(-1.0, 0.0, 0.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.slip = false;
            dummyWall.moving = false;
            dummyWall.limited = true;
            dummyWall.xMin = -numeric_limits<double>::max();
            dummyWall.xMax = +numeric_limits<double>::max();
            dummyWall.yMin = -numeric_limits<double>::max();
            dummyWall.yMax = +numeric_limits<double>::max();
            dummyWall.zMin = -numeric_limits<double>::max();
            dummyWall.zMax = +numeric_limits<double>::max();
            ++index;
            walls.push_back(dummyWall);

            // bottom of reservoir
            dummyWall.p = tVect(0.0, 0.0, erodibleHeight);
            dummyWall.n = tVect(0.0, 0.0, 1.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.slip = false;
            dummyWall.moving = false;
            dummyWall.limited = true;
            dummyWall.xMin = -numeric_limits<double>::max();
            dummyWall.xMax = reservoirX - edgeRadius;
            dummyWall.yMin = -numeric_limits<double>::max();
            dummyWall.yMax = +numeric_limits<double>::max();
            dummyWall.zMin = -numeric_limits<double>::max();
            dummyWall.zMax = +numeric_limits<double>::max();
            ++index;
            walls.push_back(dummyWall);

            dummyWall.p = tVect(reservoirX, 0.0, 0);
            dummyWall.n = tVect(1.0, 0.0, 0.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.slip = false;
            dummyWall.moving = false;
            dummyWall.limited = true;
            dummyWall.xMin = -numeric_limits<double>::max();
            dummyWall.xMax = +numeric_limits<double>::max();
            dummyWall.yMin = -numeric_limits<double>::max();
            dummyWall.yMax = +numeric_limits<double>::max();
            dummyWall.zMin = -numeric_limits<double>::max();
            dummyWall.zMax = erodibleHeight - edgeRadius;
            ++index;
            walls.push_back(dummyWall);
            break;
        }
        case IERVOLINO:
        {
            const double centerY = demSize[1] / 2.0;
            const double reservoirX = 0.8;
            const double wallSizeX = 0.025;
            const double outletSizeY = 0.3;
            const double obstacleSizeX = 0.155;
            const double obstacleSizeY = 0.3;
            const double erodibleSizeX = 0.51;
            const double erodibleHeight = 0.02;
            const double edgeRadius = 0.005;

            // reservoir left wall
            wall dummyWall;
            dummyWall.p = tVect(reservoirX, 0.0, 0);
            dummyWall.n = tVect(-1.0, 0.0, 0.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.slip = false;
            dummyWall.moving = false;
            dummyWall.limited = true;
            dummyWall.xMin = -numeric_limits<double>::max();
            dummyWall.xMax = +numeric_limits<double>::max();
            dummyWall.yMin = -numeric_limits<double>::max();
            dummyWall.yMax = centerY - outletSizeY / 2.0 - edgeRadius;
            dummyWall.zMin = -numeric_limits<double>::max();
            dummyWall.zMax = +numeric_limits<double>::max();
            ++index;
            walls.push_back(dummyWall);

            dummyWall.p = tVect(reservoirX + wallSizeX, 0.0, 0);
            dummyWall.n = tVect(+1.0, 0.0, 0.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.slip = false;
            dummyWall.moving = false;
            dummyWall.limited = true;
            dummyWall.xMin = -numeric_limits<double>::max();
            dummyWall.xMax = +numeric_limits<double>::max();
            dummyWall.yMin = -numeric_limits<double>::max();
            dummyWall.yMax = centerY - outletSizeY / 2.0 - edgeRadius;
            dummyWall.zMin = -numeric_limits<double>::max();
            dummyWall.zMax = +numeric_limits<double>::max();
            ++index;
            walls.push_back(dummyWall);

            dummyWall.p = tVect(reservoirX + wallSizeX / 2, centerY - outletSizeY / 2.0, 0);
            dummyWall.n = tVect(0.0, +1.0, 0.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.slip = false;
            dummyWall.moving = false;
            dummyWall.limited = true;
            dummyWall.xMin = reservoirX + edgeRadius;
            dummyWall.xMax = reservoirX + wallSizeX - edgeRadius;
            dummyWall.yMin = -numeric_limits<double>::max();
            dummyWall.yMax = +numeric_limits<double>::max();
            dummyWall.zMin = -numeric_limits<double>::max();
            dummyWall.zMax = +numeric_limits<double>::max();
            ++index;
            walls.push_back(dummyWall);


            // reservoir right wall
            dummyWall.p = tVect(reservoirX, 0.0, 0);
            dummyWall.n = tVect(-1.0, 0.0, 0.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.slip = false;
            dummyWall.moving = false;
            dummyWall.limited = true;
            dummyWall.xMin = -numeric_limits<double>::max();
            dummyWall.xMax = +numeric_limits<double>::max();
            dummyWall.yMin = centerY + outletSizeY / 2.0 + edgeRadius;
            dummyWall.yMax = +numeric_limits<double>::max();
            dummyWall.zMin = -numeric_limits<double>::max();
            dummyWall.zMax = +numeric_limits<double>::max();
            ++index;
            walls.push_back(dummyWall);

            dummyWall.p = tVect(reservoirX + wallSizeX, 0.0, 0);
            dummyWall.n = tVect(+1.0, 0.0, 0.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.slip = false;
            dummyWall.moving = false;
            dummyWall.limited = true;
            dummyWall.xMin = -numeric_limits<double>::max();
            dummyWall.xMax = +numeric_limits<double>::max();
            dummyWall.yMin = centerY + outletSizeY / 2.0 + edgeRadius;
            dummyWall.yMax = +numeric_limits<double>::max();
            dummyWall.zMin = -numeric_limits<double>::max();
            dummyWall.zMax = +numeric_limits<double>::max();
            ++index;
            walls.push_back(dummyWall);

            dummyWall.p = tVect(reservoirX + wallSizeX / 2, centerY + outletSizeY / 2.0, 0);
            dummyWall.n = tVect(0.0, -1.0, 0.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.slip = false;
            dummyWall.moving = false;
            dummyWall.limited = true;
            dummyWall.xMin = reservoirX + edgeRadius;
            dummyWall.xMax = reservoirX + wallSizeX - edgeRadius;
            dummyWall.yMin = -numeric_limits<double>::max();
            dummyWall.yMax = +numeric_limits<double>::max();
            dummyWall.zMin = -numeric_limits<double>::max();
            dummyWall.zMax = +numeric_limits<double>::max();
            ++index;
            walls.push_back(dummyWall);


            // obstacle
            dummyWall.p = tVect(reservoirX + wallSizeX + erodibleSizeX, 0.0, 0);
            dummyWall.n = tVect(-1.0, 0.0, 0.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.slip = false;
            dummyWall.moving = false;
            dummyWall.limited = true;
            dummyWall.xMin = -numeric_limits<double>::max();
            dummyWall.xMax = +numeric_limits<double>::max();
            dummyWall.yMin = centerY - obstacleSizeY / 2.0 + edgeRadius;
            dummyWall.yMax = centerY + obstacleSizeY / 2.0 - edgeRadius;
            dummyWall.zMin = -numeric_limits<double>::max();
            dummyWall.zMax = +numeric_limits<double>::max();
            ++index;
            walls.push_back(dummyWall);

            dummyWall.p = tVect(reservoirX + wallSizeX + erodibleSizeX + obstacleSizeX, 0.0, 0);
            dummyWall.n = tVect(1.0, 0.0, 0.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.slip = false;
            dummyWall.moving = false;
            dummyWall.limited = true;
            dummyWall.xMin = -numeric_limits<double>::max();
            dummyWall.xMax = +numeric_limits<double>::max();
            dummyWall.yMin = centerY - obstacleSizeY / 2.0 + edgeRadius;
            dummyWall.yMax = centerY + obstacleSizeY / 2.0 - edgeRadius;
            dummyWall.zMin = -numeric_limits<double>::max();
            dummyWall.zMax = +numeric_limits<double>::max();
            ++index;
            walls.push_back(dummyWall);

            dummyWall.p = tVect(0.0, centerY + obstacleSizeY / 2.0, 0);
            dummyWall.n = tVect(0.0, 1.0, 0.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.slip = false;
            dummyWall.moving = false;
            dummyWall.limited = true;
            dummyWall.xMin = reservoirX + wallSizeX + erodibleSizeX + edgeRadius;
            dummyWall.xMax = reservoirX + wallSizeX + erodibleSizeX + obstacleSizeX - edgeRadius;
            dummyWall.yMin = -numeric_limits<double>::max();
            dummyWall.yMax = +numeric_limits<double>::max();
            dummyWall.zMin = -numeric_limits<double>::max();
            dummyWall.zMax = +numeric_limits<double>::max();
            ++index;
            walls.push_back(dummyWall);

            dummyWall.p = tVect(0.0, centerY - obstacleSizeY / 2.0, 0);
            dummyWall.n = tVect(0.0, -1.0, 0.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.slip = false;
            dummyWall.moving = false;
            dummyWall.limited = true;
            dummyWall.xMin = reservoirX + wallSizeX + erodibleSizeX + edgeRadius;
            dummyWall.xMax = reservoirX + wallSizeX + erodibleSizeX + obstacleSizeX - edgeRadius;
            dummyWall.yMin = -numeric_limits<double>::max();
            dummyWall.yMax = +numeric_limits<double>::max();
            dummyWall.zMin = -numeric_limits<double>::max();
            dummyWall.zMax = +numeric_limits<double>::max();
            ++index;
            walls.push_back(dummyWall);

            // bottom of reservoir
            dummyWall.p = tVect(0.0, 0.0, erodibleHeight);
            dummyWall.n = tVect(0.0, 0.0, 1.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.slip = false;
            dummyWall.moving = false;
            dummyWall.limited = true;
            dummyWall.xMin = -numeric_limits<double>::max();
            dummyWall.xMax = reservoirX + wallSizeX - edgeRadius;
            dummyWall.yMin = -numeric_limits<double>::max();
            dummyWall.yMax = +numeric_limits<double>::max();
            dummyWall.zMin = -numeric_limits<double>::max();
            dummyWall.zMax = +numeric_limits<double>::max();
            ++index;
            walls.push_back(dummyWall);

            dummyWall.p = tVect(reservoirX + wallSizeX, 0.0, 0);
            dummyWall.n = tVect(1.0, 0.0, 0.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.slip = false;
            dummyWall.moving = false;
            dummyWall.limited = true;
            dummyWall.xMin = -numeric_limits<double>::max();
            dummyWall.xMax = +numeric_limits<double>::max();
            dummyWall.yMin = -numeric_limits<double>::max();
            dummyWall.yMax = +numeric_limits<double>::max();
            dummyWall.zMin = -numeric_limits<double>::max();
            dummyWall.zMax = erodibleHeight - edgeRadius;
            ++index;
            walls.push_back(dummyWall);

            // bottom of non-erodible floodplain
            dummyWall.p = tVect(0.0, 0.0, erodibleHeight);
            dummyWall.n = tVect(0.0, 0.0, 1.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.slip = false;
            dummyWall.moving = false;
            dummyWall.limited = true;
            dummyWall.xMin = reservoirX + wallSizeX + erodibleSizeX + edgeRadius;
            dummyWall.xMax = numeric_limits<double>::max();
            dummyWall.yMin = -numeric_limits<double>::max();
            dummyWall.yMax = +numeric_limits<double>::max();
            dummyWall.zMin = -numeric_limits<double>::max();
            dummyWall.zMax = +numeric_limits<double>::max();
            ++index;
            walls.push_back(dummyWall);

            dummyWall.p = tVect(reservoirX + wallSizeX + erodibleSizeX, 0.0, 0);
            dummyWall.n = tVect(-1.0, 0.0, 0.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.slip = false;
            dummyWall.moving = false;
            dummyWall.limited = true;
            dummyWall.xMin = -numeric_limits<double>::max();
            dummyWall.xMax = +numeric_limits<double>::max();
            dummyWall.yMin = -numeric_limits<double>::max();
            dummyWall.yMax = +numeric_limits<double>::max();
            dummyWall.zMin = -numeric_limits<double>::max();
            dummyWall.zMax = erodibleHeight - edgeRadius;
            ++index;
            walls.push_back(dummyWall);

            break;
        }
        case HEAP:
        {

            const double outletSize = demSize[0] / 4.0;

            wall dummyWall;
            dummyWall.p = tVect(0.0, 0.0, heapBaseLevel);
            dummyWall.n = tVect(0.0, 0.0, 1.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.slip = false;
            dummyWall.moving = false;
            dummyWall.limited = true;
            dummyWall.xMin = outletSize;
            dummyWall.xMax = demSize[0] - outletSize;
            ++index;
            walls.push_back(dummyWall);

            break;
        }
        case TRIAXIAL:
        {

            wall dummyWall;
            dummyWall.p = tVect(triBeginX, 0.0, 0.0);
            dummyWall.n = tVect(-1.0, 0.0, 0.0);
            dummyWall.index = index;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            dummyWall.slip = false;
            dummyWall.moving = false;
            dummyWall.limited = false;
            ++index;
            walls.push_back(dummyWall);

            dummyWall.p = tVect(0.0, triBeginY, 0.0);
            dummyWall.n = tVect(0.0, -1.0, 0.0);
            dummyWall.index = index;
            ++index;
            walls.push_back(dummyWall);

            dummyWall.p = tVect(0.0, 0.0, triBeginZ);
            dummyWall.n = tVect(0.0, 0.0, -1.0);
            dummyWall.index = index;
            ++index;
            walls.push_back(dummyWall);

            break;
        }
        }
    }
    for (int n = 0; n < walls.size(); ++n) {
        walls[n].wallShow();
    }
}

void DEM::initializeCylinders(const Problem &problem) {
    cylinders.clear();
    unsigned int index = 0;
    if (!problem.file.empty()) {  // Problem file was loaded
        // Copy cylinders directly from problem file
        //cylinders.insert(cylinders.end(), problem.cylinders.begin(), problem.cylinders.end());
        for (size_t i = 0; i < problem.cylinders.size(); ++i) {
            cylinder t = problem.cylinders[i];
            t.index = static_cast<unsigned int>(cylinders.size());
            cylinders.push_back(t);
        }
    } else {
        switch (problemName) {
        case HK_LARGE:
        {
            if (depositArea) {
                cylinder dummyCylinder;
                dummyCylinder.index = index;
                dummyCylinder.p1 = tVect(21.0858, 13.3466, 0.0);
                dummyCylinder.p2 = tVect(21.0858, 13.3466, 100.0);
                dummyCylinder.R = 13.513015;
                dummyCylinder.omega = tVect(0.0, 0.0, 0.0);
                dummyCylinder.initAxes();
                dummyCylinder.type = EMPTY;
                dummyCylinder.moving = false;
                dummyCylinder.slip = false;
                dummyCylinder.limited = true;
                dummyCylinder.xMin = 23.2;
                dummyCylinder.xMax = 25.2795;
                dummyCylinder.yMin = 0.0;
                dummyCylinder.yMax = 0.5008;
                dummyCylinder.zMin = -1.0;
                dummyCylinder.zMax = 2.0 * demSize[2];
                cylinders.push_back(dummyCylinder);
                ++index;
            }
            break;
        }
        case IERVOLINO_2D:
        {

            const double reservoirX = 0.8;
            const double wallSizeX = 0.025;
            const double erodibleSizeX = 0.51;
            const double erodibleHeight = 0.02;
            const double edgeRadius = 0.005;

            cylinder dummyCylinder;
            dummyCylinder.R = edgeRadius;
            dummyCylinder.omega = tVect(0.0, 0.0, 0.0);
            dummyCylinder.type = FULL;
            dummyCylinder.moving = false;
            dummyCylinder.slip = false;
            dummyCylinder.limited = false;


            // erodible part borders
            dummyCylinder.index = index;
            dummyCylinder.p1 = tVect(reservoirX + wallSizeX - edgeRadius, 0.0, erodibleHeight - edgeRadius);
            dummyCylinder.p2 = tVect(reservoirX + wallSizeX - edgeRadius, 1.0, erodibleHeight - edgeRadius);
            dummyCylinder.initAxes();
            cylinders.push_back(dummyCylinder);
            ++index;

            break;
        }
        case IERVOLINO:
        {

            const double centerY = demSize[1] / 2.0;
            const double reservoirX = 0.8;
            const double wallSizeX = 0.025;
            const double outletSizeY = 0.3;
            const double obstacleSizeX = 0.155;
            const double obstacleSizeY = 0.3;
            const double erodibleSizeX = 0.51;
            const double erodibleHeight = 0.02;
            const double edgeRadius = 0.005;

            // left wall borders
            cylinder dummyCylinder;
            dummyCylinder.index = index;
            dummyCylinder.p1 = tVect(reservoirX + edgeRadius, centerY - outletSizeY / 2.0 - edgeRadius, 0.0);
            dummyCylinder.p2 = tVect(reservoirX + edgeRadius, centerY - outletSizeY / 2.0 - edgeRadius, 1.0);
            dummyCylinder.R = edgeRadius;
            dummyCylinder.omega = tVect(0.0, 0.0, 0.0);
            dummyCylinder.initAxes();
            dummyCylinder.type = FULL;
            dummyCylinder.moving = false;
            dummyCylinder.slip = false;
            dummyCylinder.limited = false;
            cylinders.push_back(dummyCylinder);
            ++index;

            dummyCylinder.index = index;
            dummyCylinder.p1 = tVect(reservoirX + wallSizeX - edgeRadius, centerY - outletSizeY / 2.0 - edgeRadius, 0.0);
            dummyCylinder.p2 = tVect(reservoirX + wallSizeX - edgeRadius, centerY - outletSizeY / 2.0 - edgeRadius, 1.0);
            dummyCylinder.initAxes();
            cylinders.push_back(dummyCylinder);
            ++index;

            // right wall borders
            dummyCylinder.index = index;
            dummyCylinder.p1 = tVect(reservoirX + edgeRadius, centerY + outletSizeY / 2.0 + edgeRadius, 0.0);
            dummyCylinder.p2 = tVect(reservoirX + edgeRadius, centerY + outletSizeY / 2.0 + edgeRadius, 1.0);
            dummyCylinder.initAxes();
            cylinders.push_back(dummyCylinder);
            ++index;

            dummyCylinder.index = index;
            dummyCylinder.p1 = tVect(reservoirX + wallSizeX - edgeRadius, centerY + outletSizeY / 2.0 + edgeRadius, 0.0);
            dummyCylinder.p2 = tVect(reservoirX + wallSizeX - edgeRadius, centerY + outletSizeY / 2.0 + edgeRadius, 1.0);
            dummyCylinder.initAxes();
            cylinders.push_back(dummyCylinder);
            ++index;

            // obstacle corners
            dummyCylinder.index = index;
            dummyCylinder.p1 = tVect(reservoirX + wallSizeX + erodibleSizeX + edgeRadius, centerY - obstacleSizeY / 2.0 + edgeRadius, 0.0);
            dummyCylinder.p2 = tVect(reservoirX + wallSizeX + erodibleSizeX + edgeRadius, centerY - obstacleSizeY / 2.0 + edgeRadius, 1.0);
            dummyCylinder.initAxes();
            cylinders.push_back(dummyCylinder);
            ++index;

            dummyCylinder.index = index;
            dummyCylinder.p1 = tVect(reservoirX + wallSizeX + erodibleSizeX + obstacleSizeX - edgeRadius, centerY - obstacleSizeY / 2.0 + edgeRadius, 0.0);
            dummyCylinder.p2 = tVect(reservoirX + wallSizeX + erodibleSizeX + obstacleSizeX - edgeRadius, centerY - obstacleSizeY / 2.0 + edgeRadius, 1.0);
            dummyCylinder.initAxes();
            cylinders.push_back(dummyCylinder);
            ++index;

            dummyCylinder.index = index;
            dummyCylinder.p1 = tVect(reservoirX + wallSizeX + erodibleSizeX + edgeRadius, centerY + obstacleSizeY / 2.0 - edgeRadius, 0.0);
            dummyCylinder.p2 = tVect(reservoirX + wallSizeX + erodibleSizeX + edgeRadius, centerY + obstacleSizeY / 2.0 - edgeRadius, 1.0);
            dummyCylinder.initAxes();
            cylinders.push_back(dummyCylinder);
            ++index;

            dummyCylinder.index = index;
            dummyCylinder.p1 = tVect(reservoirX + wallSizeX + erodibleSizeX + obstacleSizeX - edgeRadius, centerY + obstacleSizeY / 2.0 - edgeRadius, 0.0);
            dummyCylinder.p2 = tVect(reservoirX + wallSizeX + erodibleSizeX + obstacleSizeX - edgeRadius, centerY + obstacleSizeY / 2.0 - edgeRadius, 1.0);
            dummyCylinder.initAxes();
            cylinders.push_back(dummyCylinder);
            ++index;

            // erodible part borders
            dummyCylinder.index = index;
            dummyCylinder.p1 = tVect(reservoirX + wallSizeX - edgeRadius, 0.0, erodibleHeight - edgeRadius);
            dummyCylinder.p2 = tVect(reservoirX + wallSizeX - edgeRadius, 1.0, erodibleHeight - edgeRadius);
            dummyCylinder.initAxes();
            cylinders.push_back(dummyCylinder);
            ++index;

            dummyCylinder.index = index;
            dummyCylinder.p1 = tVect(reservoirX + wallSizeX + erodibleSizeX + edgeRadius, 0.0, erodibleHeight - edgeRadius);
            dummyCylinder.p2 = tVect(reservoirX + wallSizeX + erodibleSizeX + edgeRadius, 1.0, erodibleHeight - edgeRadius);
            dummyCylinder.initAxes();
            cylinders.push_back(dummyCylinder);
            ++index;

            break;
        }
        case IERVOLINO_CYLINDERTEST:
        {
            cylinder dummyCylinder;
            dummyCylinder.index = index;
            dummyCylinder.p1 = tVect(1.5, 0, 1.5);
            dummyCylinder.p2 = tVect(1.5, 1.0, 1.5);
            dummyCylinder.R = 1.0;
            dummyCylinder.omega = tVect(0.0, 0.0, 0.0);
            dummyCylinder.initAxes();
            dummyCylinder.type = EMPTY;
            dummyCylinder.moving = false;
            dummyCylinder.slip = false;
            dummyCylinder.limited = false;
            cylinders.push_back(dummyCylinder);
            ++index;

            dummyCylinder.index = index;
            dummyCylinder.p1 = tVect(1.5, 0, 1.5);
            dummyCylinder.p2 = tVect(1.5, 1.0, 1.5);
            dummyCylinder.R = 0.5;
            dummyCylinder.omega = tVect(0.0, 0.0, 0.0);
            dummyCylinder.initAxes();
            dummyCylinder.type = FULL;
            dummyCylinder.moving = false;
            dummyCylinder.slip = false;
            dummyCylinder.limited = false;
            cylinders.push_back(dummyCylinder);
            ++index;

            break;

        }
        case WILL_SETTLING:
        {
            cylinder dummyCylinder;
            dummyCylinder.index = index;
            dummyCylinder.p1 = tVect(0, 0, 0);
            dummyCylinder.p2 = tVect(0, 0, 1);
            dummyCylinder.R = 0.108 / 2.0;
            dummyCylinder.omega = tVect(0.0, 0.0, 0.0);
            dummyCylinder.initAxes();
            dummyCylinder.type = EMPTY;
            dummyCylinder.moving = false;
            dummyCylinder.slip = false;
            dummyCylinder.limited = false;
            cylinders.push_back(dummyCylinder);
            ++index;

            break;

        }
        case KELVIN:
        {
            cylinder dummyCylinder;
            dummyCylinder.index = index;
            dummyCylinder.p1 = tVect(1.2504, 0.4, 0.0);
            dummyCylinder.p2 = tVect(1.2504, 0.4, 100.0);
            dummyCylinder.R = 0.400;
            dummyCylinder.omega = tVect(0.0, 0.0, 0.0);
            dummyCylinder.initAxes();
            dummyCylinder.type = EMPTY;
            dummyCylinder.moving = false;
            dummyCylinder.slip = false;
            dummyCylinder.limited = true;
            dummyCylinder.xMin = 1.250;
            dummyCylinder.xMax = 1.42535;
            dummyCylinder.yMin = -10.0;
            dummyCylinder.yMax = 0.04048;
            dummyCylinder.zMin = -1.0;
            dummyCylinder.zMax = 2.0 * demSize[2];
            cylinders.push_back(dummyCylinder);
            ++index;
            break;
        }
        case DRUM:
        {
            cylinder dummyCylinder;
            dummyCylinder.index = index;
            dummyCylinder.p1 = tVect(0.0 + 0.6 + 0.35, 0.0, 1.26);
            dummyCylinder.p2 = tVect(0.0 + 0.6 + 0.35, 1.0, 1.26);
            dummyCylinder.R = 1.243;
            dummyCylinder.omega = tVect(0.0, drumSpeed, 0.0);
            dummyCylinder.initAxes();
            dummyCylinder.type = EMPTY;
            dummyCylinder.moving = true;
            dummyCylinder.slip = false;
            dummyCylinder.limited = false;
            cylinders.push_back(dummyCylinder);
            ++index;
            break;
        }
        case AVALANCHE:
        {
            cylinder dummyCylinder;
            dummyCylinder.index = index;
            dummyCylinder.p1 = tVect(0.0, 15, 15); //tVect(0.0,14.86,14.86);
            dummyCylinder.p2 = tVect(1.0, 15, 15);
            dummyCylinder.R = 14.86;
            dummyCylinder.omega.reset();
            dummyCylinder.initAxes();
            dummyCylinder.type = EMPTY;
            dummyCylinder.moving = false;
            dummyCylinder.slip = false;
            dummyCylinder.limited = false;
            cylinders.push_back(dummyCylinder);
            ++index;
            break;
        }
        case NET:
        {
            cylinder dummyCylinder;
            dummyCylinder.index = index;
            dummyCylinder.p1 = tVect(0.0, 4.0, 4.0);
            dummyCylinder.p2 = tVect(1.0, 4.0, 4.0);
            dummyCylinder.R = 4.0;
            dummyCylinder.omega.reset();
            dummyCylinder.initAxes();
            dummyCylinder.type = EMPTY;
            dummyCylinder.moving = false;
            dummyCylinder.slip = false;
            dummyCylinder.limited = false;
            cylinders.push_back(dummyCylinder);
            ++index;
            break;
        }
        case INTRUDER:
        {
            const double intruderX = 0.02; //0.0675;   // X and Z coordinates of the intruder (in the muddle of the box X and Z)
            const double intruderZ = 0.0675;
            const double intruderd = 16e-3; // diameter of the intruder
            cylinder intruder;
            intruder.index = index;
            intruder.p1 = tVect(intruderX, 0.0, intruderZ);
            intruder.p2 = tVect(intruderX, 100.0, intruderZ);
            intruder.R = intruderd / 2.0;
            intruder.omega = tVect(0.0, 0.0, 0.0);
            intruder.initAxes();
            intruder.moving = false;
            intruder.slip = false;
            intruder.limited = false;
            intruder.type = FULL;
            intruder.translating = true;
            intruder.trans = tVect(7.5, 0.0, 0.0);
            cylinders.push_back(intruder);
            ++index;
            break;
        }
        case TBAR:
        {
            const double intruderX = 0.07; //0.0675;   // X and Z coordinates of the intruder (in the muddle of the box X and Z)
            const double intruderZ = 0.10; // to be determinated after column collpase (probelm name = none)
            const double intruderd = 16e-3; // diameter of the intruder
//            const double velocityIntruder = -30.0e-3;
            
            cylinder intruder;
            intruder.index = index;
            intruder.p1 = tVect(intruderX, -0.5, intruderZ);
            intruder.p2 = tVect(intruderX, 0.5, intruderZ);
            intruder.R = intruderd / 2.0;
            intruder.omega = tVect(0.0, 0.0, 0.0);
            intruder.initAxes();
            intruder.moving = true;
            intruder.slip = false;
            intruder.limited = false;
            intruder.type = FULL;
            intruder.translating = true;
            intruder.trans = tVect(velocityCylinderX, velocityCylinderZ, velocityCylinderZ);
            cylinders.push_back(intruder);
            ++index;
            break;
        }
        
        case SEGUIN: {
            const double intruderX = 0.00; //0.0675;   // X and Z coordinates of the intruder (in the muddle of the box X and Z)
            const double intruderZ = 0.05; // to be determinated after column collpase (probelm name = none)
            const double intruderd = 16e-3; // diameter of the intruder
//            const double velocityIntruder = -30.0e-3;
            
            cylinder intruder;
            intruder.index = index;
            intruder.p1 = tVect(intruderX, -0.05, intruderZ);
            intruder.p2 = tVect(intruderX, 0.05, intruderZ);
            intruder.R = intruderd / 2.0;
            intruder.omega = tVect(0.0, 0.0, 0.0);
            intruder.initAxes();
            intruder.moving = true;
            intruder.slip = false;
            intruder.limited = false;
            intruder.type = FULL;
            intruder.translating = true;
            intruder.trans = tVect(velocityCylinderX, velocityCylinderZ, velocityCylinderZ);
            cylinders.push_back(intruder);
            ++index;
            break;
        }
    }
    for (int n = 0; n < cylinders.size(); ++n) {
        cylinders[n].cylinderShow();
    }
}

void DEM::initializePbcs(const typeList& externalBoundary, const vecList& boundaryLocation) {

    pbcs.clear();

    unsigned int index = 0;

    for (int i = 0; i < externalBoundary.size(); i = i + 2) {

        if (externalBoundary[i] == PERIODIC) {
            if (externalBoundary[i + 1] != PERIODIC) {
                cout << "Periodic boundaries error [" << i << "]!" << endl;
                exit(0);
            } else {
                pbc dummyPbc;
                dummyPbc.p = boundaryLocation[i];
                dummyPbc.v = boundaryLocation[i + 1] - boundaryLocation[i];
                dummyPbc.setPlanes();
                dummyPbc.index = index;
                ++index;
                pbcs.push_back(dummyPbc);
            }
        }
    }

    for (int b = 0; b < pbcs.size(); ++b) {
        pbcs[b].pbcShow();
    }
}

void DEM::periodicObjects() {

    if (pbcs.size()) {
        // copy objects
        ghostList objectGhosts;
        objectGhosts.clear();

        cout << "Original objects: " << stdObjects << endl;

        cout << "Finding side object ghosts" << endl;
        for (int b = 0; b < pbcs.size(); ++b) {
            const tVect pbcVector = pbcs[b].v;
            // cycle through objects
            for (int o = 0; o < stdObjects; ++o) {
                // distances from the periodic walls
                const double leftDist = pbcs[b].pl1.dist(objects[o].x0);
                const double rightDist = pbcs[b].pl2.dist(objects[o].x0);
                //            ASSERT(leftDist>0);
                ASSERT(rightDist > 0)
                        // first plane of couple (left)
                if (leftDist < 2.0 * nebrRange && leftDist > -nebrRange) {
                    ghost dummyGhost;
                    dummyGhost.ghostIndex = objects[o].index;
                    dummyGhost.pbcVector = pbcVector;
                    objectGhosts.push_back(dummyGhost);
                }// second plane of couple (right), we copy only one time
                else if (rightDist < 2.0 * nebrRange && rightDist > -nebrRange) {
                    ghost dummyGhost;
                    dummyGhost.ghostIndex = objects[o].index;
                    dummyGhost.pbcVector = -1.0 * pbcVector;
                    objectGhosts.push_back(dummyGhost);
                }
            }
        }
        cout << "Found" << objectGhosts.size() << endl;

        // corner objects
        cout << "Finding corner object ghosts" << endl;
        ghostList cornerObjectsGhosts;
        cornerObjectsGhosts.clear();
        for (int og1 = 0; og1 < objectGhosts.size(); ++og1) {
            //cout<<"og1="<<og1<<endl;
            for (int og2 = og1 + 1; og2 < objectGhosts.size(); ++og2) {
                //cout<<"og2="<<og2<<"("<<objectGhosts[og1].ghostIndex<<" "<<objectGhosts[og2].ghostIndex<<")"<<endl;
                if (objectGhosts[og1].ghostIndex == objectGhosts[og2].ghostIndex) {
                    // we need to create a corner ghost
                    ghost dummyGhost;
                    dummyGhost.ghostIndex = objectGhosts[og1].ghostIndex;
                    dummyGhost.pbcVector = objectGhosts[og1].pbcVector + objectGhosts[og2].pbcVector;
                    cornerObjectsGhosts.push_back(dummyGhost);
                }
            }
        }
        cout << "Found" << cornerObjectsGhosts.size() << endl;


        objectGhosts.insert(objectGhosts.end(), cornerObjectsGhosts.begin(), cornerObjectsGhosts.end());


        if (objectGhosts.size() != 0) {
            cout << "Generating object ghosts" << endl;
            // updating ghost particles
            for (int og = 0; og < objectGhosts.size(); ++og) {
                //getting original particle index
                const unsigned int originObjectIndex = objectGhosts[og].ghostIndex;
                const unsigned int ghostObjectIndex = stdObjects + og;
                // reconstructing ghost objects
                object dummyObject = objects[originObjectIndex];
                dummyObject.x0 = dummyObject.x0 + objectGhosts[og].pbcVector;
                dummyObject.originalIndex = dummyObject.index;
                dummyObject.index = ghostObjectIndex;
                objects.push_back(dummyObject);
            }
        }
    }
}

// integration functions

void DEM::predictor() {
    static const double c1[5] = {deltat, deltat * deltat / 2.0, deltat * deltat * deltat / 6.0, deltat * deltat * deltat * deltat / 24.0, deltat * deltat * deltat * deltat * deltat / 120.0};
    static const double c2[5] = {deltat, deltat * deltat / 2.0, deltat * deltat * deltat / 6.0, deltat * deltat * deltat * deltat / 24.0, deltat * deltat * deltat * deltat * deltat / 120.0};
    //    c[0] = deltat;
    //    c[1] = c[0] * deltat / 2.0;
    //    c[2] = c[1] * deltat / 3.0;
    //    c[3] = c[2] * deltat / 4.0;
    //    c[4] = c[3] * deltat / 5.0;

    //#pragma omp parallel for
    for (int a = 0; a < activeElmts.size(); ++a) {
        //cout<<"pr a= "<<a<<endl;
        unsigned int n = activeElmts[a];
        //cout<<"pr n= "<<n<<endl;
        elmts[n].predict(c1, c2);
//        elmts[n].xp0[1]=0.004;
//        elmts[n].xp1[1]=0.0;
        //cout<<"pr end"<<endl;
    }
}

void DEM::corrector() {
    //static double gear[6] = {3.0/16.0, 251.0/360.0, 1.0, 11.0/18.0, 1.0/6.0, 1.0/60.0};
    //static const double c[5]={deltat, deltat*deltat/2.0, deltat*deltat*deltat/6.0, deltat*deltat*deltat*deltat/24.0, deltat*deltat*deltat*deltat*deltat/120.0};
    //static const double coeff[6]={gear[0]*c[1], gear[1]*c[1]/c[0], gear[2]*c[1]/c[1], gear[3]*c[1]/c[2], gear[4]*c[1]/c[3], gear[5]*c[1]/c[4]};


    // see: Computer Simulation of Liquids By M. P. Allen, D. J. Tildesle 
    static const double gear1ord[6] = {95.0 / 288.0, 1.0, 25.0 / 24.0, 35.0 / 72.0, 5.0 / 48.0, 1.0 / 120.0}; // coefficients for a first-order equation
    static const double gear2ord[6] = {3.0 / 16.0, 251.0 / 360.0, 1.0, 11.0 / 18.0, 1.0 / 6.0, 1.0 / 60.0}; // coefficients for a second-order equation

    static const double c[5] = {deltat, deltat * deltat / 2.0, deltat * deltat * deltat / 6.0, deltat * deltat * deltat * deltat / 24.0, deltat * deltat * deltat * deltat * deltat / 120.0};
    // Careful the c's must be changed for eq of 1st order
    static const double coeff1ord[6] = {gear1ord[0] * c[0], gear1ord[1] * c[0] / c[0], gear1ord[2] * c[0] / c[1], gear1ord[3] * c[0] / c[2], gear1ord[4] * c[0] / c[3], gear1ord[5] * c[0] / c[4]};
    // coeff2ord has different c's sequence 
    static const double coeff2ord[6] = {gear2ord[0] * c[1], gear2ord[1] * c[1] / c[0], gear2ord[2] * c[1] / c[1], gear2ord[3] * c[1] / c[2], gear2ord[4] * c[1] / c[3], gear2ord[5] * c[1] / c[4]};


    //    doubleList coeff;
    //    coeff.resize(6);
    //    coeff[0]=gear[0]*c[1];
    //    coeff[1]=gear[1]*c[1]/c[0];
    //    coeff[2]=gear[2]*c[1]/c[1];
    //    coeff[3]=gear[3]*c[1]/c[2];
    //    coeff[4]=gear[4]*c[1]/c[3];
    //    coeff[5]=gear[5]*c[1]/c[4];

#pragma omp parallel for
    for (int a = 0; a < activeElmts.size(); ++a) {
        unsigned int n = activeElmts[a];
        elmts[n].correct(coeff1ord, coeff2ord);
    }
}

void DEM::evaluateForces() {
    for (int a = 0; a < activeElmts.size(); ++a) {
        unsigned int n = activeElmts[a];
        elmts[n].FParticle.reset();
        elmts[n].FWall.reset();
        elmts[n].FCylinder.reset();
        elmts[n].MParticle.reset();
        elmts[n].MWall.reset();
        elmts[n].MCylinder.reset();
        elmts[n].FSpringP.reset();
        elmts[n].FSpringW.reset();
        elmts[n].MRolling.reset();
        // reset also apparent accellerations
        elmts[n].ACoriolis.reset();
        elmts[n].ACentrifugal.reset();
        // reset also maximum overlap (moved to IO function)
        if (elmts[n].maxOverlap > elmts[n].maxDtOverlap) {
            elmts[n].maxDtOverlap = elmts[n].maxOverlap;
        }
        elmts[n].maxOverlap = 0.0;
        elmts[n].slippingCase = 0;
        // information about phase change and connectivity
        elmts[n].solidIntensity.reset();
        elmts[n].coordination = 0;
    }

    for (int w = 0; w < walls.size(); ++w) {
        walls[w].FParticle.reset();
    }

    for (int o = 0; o < objects.size(); ++o) {
        objects[o].FParticle.reset();
    }
    
    for (int c = 0; c < cylinders.size(); ++c) {
        cylinders[c].FParticle.reset();
    }

    if (staticFrictionSolve) {
        // set tentatively all springs to inactive
        for (int a = 0; a < activeParticles.size(); ++a) {
            unsigned int p = activeParticles[a];
            for (int t = 0; t < 4; ++t) {
                for (int s = 0; s < particles[p].springs[t].size(); ++s) {
                    particles[p].springs[t][s].active = false;
                }
            }
        }
    }

    // forces due to particle overlap and lubrication
    particleParticleContacts();

    // forces due to contact with plane walls and lubrication
    wallParticleContacts();

    // forces due to contact with stationary spheres (objects)
    objectParticleContacts();

    // forces due to contact with cylinders
    cylinderParticelContacts();

    totSprings = 0;
    if (staticFrictionSolve) {
        // erase springs that are still inactive
        for (int a = 0; a < activeParticles.size(); ++a) {
            unsigned int p = activeParticles[a];
            for (int t = 0; t < 4; ++t) {
                for (int s = 0; s < particles[p].springs[t].size(); ++s) {
                    if (particles[p].springs[t][s].active == false) {
                        particles[p].springs[t].erase(particles[p].springs[t].begin() + s);
                        --s;
                    } else {
                        ++totSprings;
                    }
                }
            }
        }
    }

    // compute apparent accelerations
    if (solveCentrifugal || solveCoriolis) {
        computeApparentForces();
    }

    // save info on maximum local force on the particles
    saveObjectForces();

    //  Newton equations solution
    for (int a = 0; a < activeElmts.size(); ++a) {
        unsigned int n = activeElmts[a];

        // numerical viscosity for stability
        // see "Viscous torque on a sphere under arbitrary rotation" by Lei,  Yang, and Wu, Applied Physics Letters 89, 181908 (2006)
        const tVect FVisc = -6.0 * M_PI * numVisc * elmts[n].radius * elmts[n].xp1;
        const tVect MVisc = -8.0 * M_PI * numVisc * elmts[n].radius * elmts[n].radius * elmts[n].radius * elmts[n].wpGlobal;

        // translational motion
        // acceleration (sum of forces / mass + sum of accelerations)
        elmts[n].x2 = (FVisc + elmts[n].FHydro + elmts[n].FParticle + elmts[n].FWall + elmts[n].FCylinder) / elmts[n].m + demF + elmts[n].ACoriolis + elmts[n].ACentrifugal;
        // problem name Y correction ( Y needs to be fixed to have spheres that acts like cylinders)
        switch (problemName) {
            case (TBAR): { // fix acceleration = 0.0 (elements in 2D domain)
                elmts[n].x2.y = 0.0;
                break;
            }
            case (SEGUIN): { // fix acceleration = 0.0 (elements in 2D domain)
                elmts[n].x2.y = 0.0;
                if (elmts[n].index == 0){ // immersed cylinder case
                    elmts[n].x2.y = 0.0;
                }
                break;
            }
        }
        // rotational motion
        // adjoint of orientation quaternion
        //const tQuat q0adj=elmts[n].qp0.adjoint();
        // rotational velocity (body-fixed reference frame)
        //const tVect wBf=2.0*quat2vec( q0adj.multiply( elmts[n].qp1 ) );
        // moment in global reference frame
        const tVect moment = MVisc + elmts[n].MHydro + elmts[n].MParticle + elmts[n].MWall + elmts[n].MRolling + elmts[n].MCylinder;

        // moment in body-fixed reference frame
        //if (elmts[n].size
        const tVect momentBf = project(moment, elmts[n].qp0.adjoint());
        // rotational acceleration (body-fixed reference frame) (Newton equation for principal system)
        const tVect waBf = newtonAcc(momentBf, elmts[n].I, elmts[n].wpLocal);
        // rotational acceleration (vector)
        elmts[n].w1 = project(waBf, elmts[n].qp0);
        // problem name XZ correction ( XZ needs to be fixed to have spheres that acts like cylinders)
        switch (problemName){
            // fix acceleration = 0.0 (elements in 2D domain)
            case (TBAR): {
                elmts[n].w1.x = 0.0;
                elmts[n].w1.z = 0.0;
                break;
            }
            case (SEGUIN): {
                // Fix acceleration = 0.0 (elements in 2D domain)
                elmts[n].w1.x = 0.0;
                elmts[n].w1.z = 0.0;

                // Immersed cylinder
                if (elmts[n].index==0){
                    elmts[n].w1.x = 0.0;
                    elmts[n].w1.z = 0.0;
                }
                break;
            }
        }
        // rotational acceleration (quaternion)
        if (elmts[n].size > 1) {
            switch (problemName) {
                // immersed cylinder
                case (SEGUIN):{
                    elmts[n].q2 = tQuat(0.0,0.0,0.0,0.0);
                    break;
                } 
                // others
                case (NONE):{
                    const tQuat waQuat = quatAcc(waBf, elmts[n].qp1);
                    elmts[n].q2 = 0.5 * elmts[n].qp0.multiply(waQuat);
                    break;
                }
            }
        }
    }

}

void DEM::updateParticlesPredicted() {
    // this subroutine assumes that there is no change in the structure of the particles
    // any creation, deletion or reorganization of particles requires a different subroutine

    //    #pragma omp parallel for
    for (int p = 0; p < stdPartNumber; ++p) {
        if (particles[p].active) {
            //getting belonging element index
            const unsigned int clusterIndex = particles[p].clusterIndex;
            particles[p].updatePredicted(elmts[clusterIndex], prototypes);
        }
    }
    if (ghosts.size() != 0) {
        // updating ghost particles
        for (int g = 0; g < ghosts.size(); ++g) {
            //getting original particle index
            const unsigned int originParticleIndex = ghosts[g].ghostIndex;
            const unsigned int ghostParticleIndex = stdPartNumber + g;
            // reconstructing ghost particle
            particles[ghostParticleIndex].ghostUpdate(particles[originParticleIndex], ghosts[g].pbcVector);
        }
    }
}

void DEM::updateParticlesCorrected() {
    // this subroutine assumes that there is no change in the structure of the particles
    // any creation, deletion or reorganization of particles requires a different subroutine

    //    #pragma omp parallel for
    for (int p = 0; p < stdPartNumber; ++p) {
        if (particles[p].active) {
            //getting belonging element index
            const unsigned int clusterIndex = particles[p].clusterIndex;
            particles[p].updateCorrected(elmts[clusterIndex], prototypes);
        }
    }

    if (ghosts.size() != 0) {
        // updating ghost particles
        for (int g = 0; g < ghosts.size(); ++g) {
            //getting original particle index
            const unsigned int originParticleIndex = ghosts[g].ghostIndex;
            const unsigned int ghostParticleIndex = stdPartNumber + g;
            // reconstructing ghost particle
            particles[ghostParticleIndex].ghostUpdate(particles[originParticleIndex], ghosts[g].pbcVector);
        }
    }
}

double DEM::criticalTimeStep() const {
    // determines the critical time step based on the stiffness and mass of the elements
    // we use YADE documentation, see https://yade-dem.org/doc/formulation.html

    // maximum damping coefficient to avoid having too large deltat
    static const double maxDampCoef = 0.9;

    double minRad = elmts[0].radius;
    double minMass = elmts[0].m;
    for (int n = 0; n < elmts.size(); ++n) {
        minRad = std::min(minRad, elmts[n].radius);
        minMass = std::min(minMass, elmts[n].m);
    }

    // double const k=8.0/15.0*sphereMat.youngMod/(1-sphereMat.poisson*sphereMat.poisson)*sqrt(minRad);

    double deltaTCrit = 0.0;
    switch (sphereMat.contactModel) {
        case HERTZIAN:
        {
            // maximum length scale
            const double maxDist = std::max(demSize[0], std::max(demSize[1], demSize[2]));
            // modulus of acceleration
            const double maxAccel = demF.norm();
            // estimate maximum velocity (assuming governed by acceleration field)
            const double maxVel = std::sqrt(2.0 * maxAccel * maxDist);
            // see Landau & Lifshitz, or better Antypov & Elliott
            deltaTCrit = 2.214 * (2.0 * minRad) * std::pow(sphereMat.density / sphereMat.youngMod, 2.0 / 5.0) * std::pow(maxVel, 1.0 / 5.0);
            cout << "Computed duration of collisions: E=" << sphereMat.youngMod << " r=" << minRad << ", then t_c=" << deltaTCrit << endl;
            // deltaTCrit=minRad*sqrt(sphereMat.density/sphereMat.youngMod);
            break;
        }
        case LINEAR:
        {
            // compute dampCoef to use here
            double dampCoefDeltaT = 0.0;
            if (sphereMat.dampCoeff > maxDampCoef) {
                cout << "reducing dampCoeff (alpha_n) from " << sphereMat.dampCoeff << " to " << maxDampCoef << " for the purpose of deltaT computation" << endl;
                dampCoefDeltaT = maxDampCoef;
            } else {
                dampCoefDeltaT = sphereMat.dampCoeff;
            }
            deltaTCrit = M_PI / sqrt(sphereMat.linearStiff / minMass * (1.0 - dampCoefDeltaT * dampCoefDeltaT));
            cout << "Computed duration of collisions: k=" << sphereMat.linearStiff << " alpha_n=" << dampCoefDeltaT << " m=" << minMass << ", then t_c=" << deltaTCrit << endl;
            break;
        }
    }

    return deltaTCrit;
}

// neighbor list functions

void DEM::saveObjectForces() {

    const tVect direction = Xp;

    // save maximum local force
    for (int o = 0; o < objects.size(); ++o) {
        objects[o].updateMax(direction, demTime);
    }

    // compute total force on objects
    tVect totalForce = Zero;
    for (int o = 0; o < objects.size(); ++o) {
        totalForce += objects[o].FHydro + objects[o].FParticle;
    }


    if (objMaxTotalForce.dot(direction) < totalForce.dot(direction)) {
        objMaxTotalForce = totalForce;
        // save forces in this instant
        for (int o = 0; o < objects.size(); ++o) {

            objects[o].saveForces();
        }
    }




}

void DEM::initNeighborParameters() {
    // initializes all parameters useful for neighbor list algorithm

    cout << "Initialize neighbor list\n";

    const double radiusMultiplier = 2.5;

    // maximum radius
    double maxRad = 0.0;
    for (int i = 0; i < elmts.size(); ++i) {
        if (maxRad < elmts[i].radius) {
            maxRad = elmts[i].radius;
        }
    }
    for (int o = 0; o < objects.size(); ++o) {
        if (maxRad < objects[o].r) {
            maxRad = objects[o].r;
        }
    }
    // minimum radius
    double minRad = 1.0e9;
    for (int i = 0; i < elmts.size(); ++i) {
        if (minRad > elmts[i].radius) {
            minRad = elmts[i].radius;
        }
    }
    cout << "Max radius=" << maxRad << endl;
    cout << "Min radius=" << minRad << endl;
    // if there are no particles, just to avoid nonsense numbers, we use an only cell
    // therefore the cell width is the same as the simulation box (which comes from the LB settings)
    if (elmts.size() == 0) {
        cellWidth[0] = demSize[0];
        cellWidth[1] = demSize[1];
        cellWidth[2] = demSize[2];
    }// otherwise the cell size is determined by the size of the particles
    else {
        cellWidth[0] = std::min(maxRad * radiusMultiplier, demSize[0]);
        cellWidth[1] = std::min(maxRad * radiusMultiplier, demSize[1]);
        cellWidth[2] = std::min(maxRad * radiusMultiplier, demSize[2]);
    }

    for (int k = 0; k < 3; ++k) {
        // ncells = numbers of cells for the linked cell algorithm
        nCells[k] = int( floor(demSize[k] / cellWidth[k]));
        // width of the cells (actual)
        cellWidth[k] = demSize[k] / double(nCells[k]);
        // increase by two to give a container for ghost cells, in case of periodicity
        nCells[k] += 2;
    }

    // may need a revision
    nebrRange = std::max(maxRad * radiusMultiplier, 0.5 * std::min(cellWidth[0], std::min(cellWidth[1], cellWidth[2])));
    maxDisp = 100.0 * nebrRange;
    cout << "Neighbor list parameters\n";
    cout << "DEM size " << demSize[0] << " " << demSize[1] << " " << demSize[2] << "\n";
    cout << "Number of Cells " << nCells[0] << " " << nCells[1] << " " << nCells[2] << "\n";
    cout << "Cell width " << cellWidth[0] << " " << cellWidth[1] << " " << cellWidth[2] << "\n";
    cout << "Range " << nebrRange << "\n";
}

void DEM::evalMaxDisp() {

    double maxVel = 0.0;
    for (int a = 0; a < elmts.size(); ++a) {
        unsigned int n = activeElmts[a];
        //if (n==elmts.size()-1) {
        //cout<<"a="<<a<<" ("<<activeElmts.size()<<")"<<endl;
        //cout<<"n="<<n<<" ("<<elmts.size()<<")"<<endl;
        //if (a!=n) {
        //    for (int x=0; x<activeElmts.size(); ++x) {
        //        cout<<"a="<<activeElmts[x]<<" n="<<x<<endl;
        //    }
        //}
        //}
        const double velNorm2 = elmts[n].x1.norm2();
        //cout<<"end"<<endl;
        if (velNorm2 > maxVel) {
            maxVel = velNorm2;
        }
    }
    maxVel = sqrt(maxVel);
    maxDisp += maxVel*deltat;
}

void DEM::evalCellTable() {
    // updates cellTable, a table which associates every particle to the belonging cell

    // delete precedent cell table
    cellTable.clear();
    // create a new cell table of size: number of particles + number of cells
    cellTable.resize(particles.size() + nCells[0] * nCells[1] * nCells[2]);

    // assigns value -1 to all elements of cellTable (the reason becomes clear when looking at the last two lines of this function)
    for (int n = 0; n < cellTable.size(); ++n) {
        // -1 is the default value for unassigned node
        cellTable[n] = -1;
    }

    // cycle through n = cellTable size (number of particles + number of cells)$
    // only active particles take part into neighbor list creation
    for (int a = 0; a < activeParticles.size(); a++) {
        unsigned int p = activeParticles[a];

        // c is a identifier for the belonging cell of a particle. If particle A belongs to cell (2,5,7) in a cell
        // system with 4*6*9 cells, then it will be c = particleSize + 2 + (5*4 + 7)*6 = particleSize + 164 <- bullshit :D
        // basically is a system to store the cell coordinate of a particle in a line vector
        // Here the code: floor(x/cellWidth[0]) + nCells[0]*( floor(y/cellWidth[1]) + (nCells[1]*floor(z/cellWidth[2])) );
        const int c = particles[p].x0.linearizePosition(cellWidth, nCells) + particles.size();
        particles[p].tableCell = c;
        // if the coordinate exceeds the borders of the box, a message is displayed
        if (c > cellTable.size() || c < 0) { // just put control over ghost, here !!!!!!!!!!!!!!!!!!!!!!
            cout << "#neighborList, " << c << " initCellTable: particle " << p << " outside box, ignoring for force calculation." << endl;
            cout << "Position: (" << elmts[p].x0.dot(Xp) << " " << elmts[p].x0.dot(Yp) << " " << elmts[p].x0.dot(Zp) << ")" << endl;
            cout << "Velocity: (" << elmts[p].x1.dot(Xp) << " " << elmts[p].x1.dot(Yp) << " " << elmts[p].x1.dot(Zp) << ")" << endl;
            elmts[p].resetVelocity();
            elmts[p].radius = 0.0;
            elmts[p].x0 = Zero;
            //exit(0);
            //            continue;
        }

        // cellTable is a structure to contain data in an efficient way
        // every element of the array contains the pointer to another particle of the cell, until -1 is found
        // this explains why it needs to have particle.size() + cell.size() elements:
        // every particle points to another particle, and we need cell.size() elements with value -1
        cellTable[p] = cellTable[c];
        cellTable[c] = p;
    }


}

void DEM::destroyElements() {


    if (problemName == HONGKONG || problemName == STVINCENT || problemName == ESERCITAZIONE) {

        double destroyWallPosition = 0.0;
        if (problemName == HONGKONG) {
            destroyWallPosition = 2.0;
        } else if (problemName == STVINCENT) {
            destroyWallPosition = 7.0;
        } else if (problemName == ESERCITAZIONE) {
            destroyWallPosition = 30.0;
        }

        for (int n = 0; n < elmts.size(); ++n) {
            if (elmts[n].active == true && elmts[n].x0.dot(Xp) > destroyWallPosition) {
                elmts[n].active = false;
                elmts[n].resetVelocity();
                elmts[n].FParticle.reset();
                elmts[n].FWall.reset();
                elmts[n].MParticle.reset();
                elmts[n].MWall.reset();
                elmts[n].FSpringP.reset();
                elmts[n].FSpringW.reset();
                elmts[n].MRolling.reset();
                // reset also apparent accellerations
                elmts[n].ACoriolis.reset();
                elmts[n].ACentrifugal.reset();
                elmts[n].maxDtOverlap = 0.0;
                elmts[n].maxOverlap = 0.0;
                elmts[n].slippingCase = 0;
                for (int p = 0; p < elmts[n].components.size(); ++p) {
                    const int index = elmts[n].components[p];
                    particles[index].active = false;
                    // erase springs
                    for (int t = 0; t < 4; ++t) {
                        particles[p].springs[t].clear();
                    }
                    --actvPartNumber;
                }
            }
        }
    } else if (problemName == FILIPPO_SILOS) {

        double destroyWallPosition = 0.3;
        for (int n = 0; n < elmts.size(); ++n) {
            if (elmts[n].active == true && elmts[n].x0.dot(Zp) < destroyWallPosition) {
                elmts[n].active = false;
                elmts[n].resetVelocity();
                elmts[n].FParticle.reset();
                elmts[n].FWall.reset();
                elmts[n].MParticle.reset();
                elmts[n].MWall.reset();
                elmts[n].FSpringP.reset();
                elmts[n].FSpringW.reset();
                elmts[n].MRolling.reset();
                // reset also apparent accellerations
                elmts[n].ACoriolis.reset();
                elmts[n].ACentrifugal.reset();
                elmts[n].maxDtOverlap = 0.0;
                elmts[n].maxOverlap = 0.0;
                elmts[n].slippingCase = 0;
                for (int p = 0; p < elmts[n].components.size(); ++p) {
                    const int index = elmts[n].components[p];
                    particles[index].active = false;
                    // erase springs
                    for (int t = 0; t < 4; ++t) {
                        particles[p].springs[t].clear();
                    }
                    --actvPartNumber;
                }
            }
        }
    }


}

void DEM::evalNeighborTable() {
    //    cout<<"NEW NEIGHBOR TABLE"<<endl;
    //    cout<<"Neighbor table evaluation"<<endl;
    // prototype neighbors

    // cout<<"NEW("<<neighborTable.size()<<")";
    //cout<<"new neighbor table"<<endl;
    const unsigned int neighborCellNumber = 14;

    const static int shiftCell[neighborCellNumber][3] = {
        {0, 0, 0},
        {1, 0, 0},
        {1, 1, 0},
        {0, 1, 0},
        {-1, 1, 0},
        {0, 0, 1},
        {1, 0, 1},
        {1, 1, 1},
        {0, 1, 1},
        {-1, 1, 1},
        {-1, 0, 1},
        {-1, -1, 1},
        {0, -1, 1},
        {1, -1, 1}
    };

    const double nebrRange2 = nebrRange*nebrRange;

    destroyElements();

    // erase ghost particles
    while (particles.size() > stdPartNumber) {
        particles.pop_back();
    }
    //particles.resize(stdPartNumber);
    // if periodicity is active: shift particles, identify ghosts and set the update signal for the LBM
    if (pbcs.size() != 0) {
        // cleaning up
        ghosts.clear();
        //resizing components of elements
        for (int n = 0; n < elmts.size(); ++n) {
            elmts[n].components.resize(elmts[n].size);
        }
        // periodicity shift
        pbcShift();
        // regenerate particles
        updateParticlesCorrected();
        // identify ghosts
        createGhosts();
        // update particle listing to get ready for neighborList
        updateParticlesCorrected();
        // update signal for LBM
        newNeighborList = true;
    }


    activeParticles.clear();
    for (int p = 0; p < particles.size(); ++p) {
        if (particles[p].active) {
            activeParticles.push_back(p);
        }
    }
    //cout<<"tot particles "<<particles.size()<<endl;
    //cout<<"tot active particles "<<activeParticles.size()<<endl;
    activeElmts.clear();
    for (int n = 0; n < elmts.size(); ++n) {
        if (elmts[n].active) {
            activeElmts.push_back(n);
        }
    }
    //cout<<"tot elements "<<elmts.size()<<endl;
    //cout<<"tot active elements "<<activeElmts.size()<<endl;

    // every particle in its cell
    evalCellTable();
    // delete precedent neighbor table
    neighborTable.clear();
    // elements size
    const unsigned int particleSize = particles.size();


    // cycle through all cells ->i0
    // variables for cell coordinate
    unsigned int i0[3];
    for (i0[0] = 0; i0[0] < nCells[0]; ++i0[0]) {
        for (i0[1] = 0; i0[1] < nCells[1]; ++i0[1]) {
            for (i0[2] = 0; i0[2] < nCells[2]; ++i0[2]) {

                // cycle through neighbors vectors -> [s][k]
                for (int s = 0; s < neighborCellNumber; s++) {
                    // variables for cell coordinate
                    unsigned int i1[3];
                    for (int k = 0; k < 3; k++) {
                        // determines the neighbor cells starting from the considered cell and the
                        // prototype for neighbors (shifCell)
                        i1[k] = i0[k] + shiftCell[s][k];
                        // if neighbor cell is in upper border
                        if (i1[k] == nCells[k]) {
                            i1[k] = 0;
                        }
                        // if neighbor cell is in lower border
                        if (i1[k] < 0) {
                            i1[k] = nCells[k] - 1;
                        }
                    }

                    // linearized coordinate of cell0 (first cell of couples)
                    const unsigned int cell0 = (i0[2] * nCells[1] + i0[1]) * nCells[0] + i0[0] + particleSize;
                    // linearized coordinate of cell1 (second cell of couples)
                    const unsigned int cell1 = (i1[2] * nCells[1] + i1[1]) * nCells[0] + i1[0] + particleSize;

                    // this cycles through the elements of a cell, checking all neighbors
                    // the storage system is efficient but not trivial
                    int n0 = cellTable[cell0];
                    while (n0>-1) {
                        int n1 = cellTable[cell1];
                        while (n1>-1) {
                            // we save only the pair n1<n0 (otherwise twice the number of pairs), the pair n0<n1 (which is just the same) is excluded 
                            if (cell0 != cell1 || n1 < n0) {
                                // we exclude particles belonging to the same cluster
                                if (particles[n0].clusterIndex != particles[n1].clusterIndex) {
                                    // we also exclude too far away pairs
                                    const tVect r10 = particles[n1].x0 - particles[n0].x0;
                                    if (r10.norm2() < nebrRange2) {
                                        if (n1 > n0) {
                                            neighborTable.push_back(n1);
                                            neighborTable.push_back(n0);
                                        } else {
                                            neighborTable.push_back(n0);
                                            neighborTable.push_back(n1);
                                        }
                                    }
                                }
                            }
                            // to the next element (cell1)
                            n1 = cellTable[n1];
                        }
                        // to the next element (cell0)
                        n0 = cellTable[n0];
                    }
                }
            }
        }
    }

    if (walls.size() != 0) {
        evalNearWallTable();
    }
    if (cylinders.size() != 0) {
        evalNearCylinderTable();
    }
    if (objects.size() != 0) {
        evalNearObjectTable();
    }
    //cout<<"Size neighbor list:"<<neighborTable.size()<<endl;

}

void DEM::evalNearWallTable() {
    // evaluate the distance of all particles to the walls. Creates a list of all particles
    // inside a certain range (nebrRange)

    nearWallTable.clear();

    for (int n = 0; n < stdPartNumber; ++n) {
        if (particles[n].active) {
            for (int w = 0; w < walls.size(); ++w) {
                if (walls[w].dist(particles[n].x0) < nebrRange) {
                    nearWallTable.push_back(n);
                    nearWallTable.push_back(w);
                }
            }
        }

    }

}

void DEM::evalNearObjectTable() {
    // evaluate the distance of all particles to the walls. Creates a list of all particles
    // inside a certain range (nebrRange)


    nearObjectTable.clear();

    const double nebrRange2 = nebrRange*nebrRange;

    // this loop
    for (int n = 0; n < stdPartNumber; ++n) {
        if (particles[n].active) {
            const tVect posPartHere = particles[n].x0;
            const double radiusPartHere = particles[n].r;
            for (int o = 0; o < objects.size(); ++o) {
                const double radiusObjHere = objects[o].r;
                const tVect x0ij = posPartHere - objects[o].x0;
                if (x0ij.norm2() < nebrRange2) {
                    nearObjectTable.push_back(n);
                    nearObjectTable.push_back(o);
                }
            }
        }
    }

}

void DEM::evalNearCylinderTable() {
    // evaluate the distance of all particles to the cylinders. Creates a list of all particles
    // inside a certain range (nebrRange)

    nearCylinderTable.clear();

    for (int n = 0; n < stdPartNumber; ++n) {
        if (particles[n].active) {
            for (int c = 0; c < cylinders.size(); ++c) {
                if (cylinders[c].dist(particles[n].x0) < nebrRange) {
                    nearCylinderTable.push_back(n);
                    nearCylinderTable.push_back(c);
                }
            }
        }
    }


}

// periodicity functions

void DEM::pbcShift() {
    //shifts elements and particle according to boundary conditions

    // this subroutine works on elements: particles need be generated/updated elsewhere
    for (int b = 0; b < pbcs.size(); ++b) {
        const tVect pbcVector = pbcs[b].v;
        // cycle through elements
        for (int n = 0; n < elmts.size(); ++n) {
            const double leftDist = pbcs[b].pl1.dist(elmts[n].x0);
            const double rightDist = pbcs[b].pl2.dist(elmts[n].x0);
            // first plane of couple (left)
            if (leftDist < 0.0) {
                // cout<<"Shift element number"<<elmts[n].index<<"\n";
                elmts[n].translate(pbcVector);
                // shift particles (not the ghosts)
                //                for (int j=0; j<elmts[n].size; ++j) {
                //                    particles[elmts[n].components[j]].x0+=pbcVector;
                //                }
            }
            // second plane of couple (right)
            if (rightDist < 0.0) {
                //                cout<<"Shift element number"<<elmts[n].index<<"\n";
                elmts[n].translate(-1.0 * pbcVector);
                // shift particles (not the ghosts)
                //                for (int j=0; j<elmts[n].size; ++j) {
                //                    particles[elmts[n].components[j]].x0-=pbcVector;
                //                }
            }
        }
    }
}

void DEM::createGhosts() {

    ghosts.clear();

    // this subroutine implies that elements (and their non-ghost particles) have already been shifted
    for (int b = 0; b < pbcs.size(); ++b) {
        const tVect pbcVector = pbcs[b].v;
        // cycle through elements
        for (int n = 0; n < elmts.size(); ++n) {
            // cycle through standard particles
            for (int j = 0; j < elmts[n].size; ++j) {
                const unsigned int p = elmts[n].components[j];
                // distances from the periodic walls
                const double leftDist = pbcs[b].pl1.dist(particles[p].x0);
                const double rightDist = pbcs[b].pl2.dist(particles[p].x0);
                // first plane of couple (left)
                if (leftDist < nebrRange) {
                    ghost dummyGhost;
                    dummyGhost.ghostIndex = particles[p].particleIndex;
                    dummyGhost.pbcVector = pbcVector;
                    ghosts.push_back(dummyGhost);
                }// second plane of couple (right), we copy only one time
                else if (rightDist < nebrRange) {
                    ghost dummyGhost;
                    dummyGhost.ghostIndex = particles[p].particleIndex;
                    dummyGhost.pbcVector = -1.0 * pbcVector;
                    ghosts.push_back(dummyGhost);
                }
            }
        }
    }

    // now we need to check for particles at corners
    // these are particles have already been shifted by two periodic walls IMPORTANT: (the case of three periodic conditions is missing)
    ghostList cornerGhosts;
    cornerGhosts.clear();
    for (int g1 = 0; g1 < ghosts.size(); ++g1) {
        for (int g2 = g1 + 1; g2 < ghosts.size(); ++g2) {
            if (ghosts[g1].ghostIndex == ghosts[g2].ghostIndex) {
                // we need to create a corner ghost
                ghost dummyGhost;
                dummyGhost.ghostIndex = ghosts[g1].ghostIndex;
                dummyGhost.pbcVector = ghosts[g1].pbcVector + ghosts[g2].pbcVector;
                cornerGhosts.push_back(dummyGhost);
            }
        }
    }
    //    cout<<cornerGhosts.size()<<endl;
    ghosts.insert(ghosts.end(), cornerGhosts.begin(), cornerGhosts.end());

    //resizing components of elements
    for (int t = 0; t < elmts.size(); ++t) {
        elmts[t].components.resize(elmts[t].size);
    }

    // creating particles from ghosts
    // index for ghost particles
    unsigned int ind = stdPartNumber;
    // adding ghost particles to particle list
    for (int g = 0; g < ghosts.size(); ++g) {
        //getting original particle index
        const unsigned int originParticleIndex = ghosts[g].ghostIndex;
        // reconstructing ghost particle
        particle dummyParticle = particles[originParticleIndex];
        //        dummyParticle.x0+=ghosts[g].pbcVector;
        dummyParticle.particleIndex = ind;
        // set ghost
        dummyParticle.isGhost = true;
        // adding ghost particle to components of the mother element
        elmts[dummyParticle.clusterIndex].components.push_back(dummyParticle.particleIndex);
        // updating particle list
        particles.push_back(dummyParticle);
        // updating particle index
        ind++;
    }
}

// force computation functions

void DEM::computeApparentForces() {

    // unit vector of the rotation axes
    const tVect unitRot = demRot / demRot.norm();

    // computes forces (rather, accelerations) that appear if the reference system is rotating
    for (int n = 0; n < elmts.size(); ++n) {
        //local coordinates and velocity

        // Calculate the centrifugal acceleration
        if (solveCentrifugal) {
            // centrifugal acceleration
            elmts[n].ACentrifugal = computeCentrifugal(elmts[n].xp0, demRotCenter, demRot);
        } else {
            elmts[n].ACentrifugal.reset();
        }

        // Calculate Coriolis acceleration
        if (solveCoriolis) {
            // Get Coriolis acceleration: 
            elmts[n].ACoriolis = computeCoriolis(elmts[n].xp1, demRot);
            //elmts[n].ACoriolis.show();
        } else {
            elmts[n].ACoriolis.reset();
        }

    }
}

Elongation* DEM::findSpring(const unsigned int& t, const unsigned int& indexI, particle* partJ) {

    // dummy, inactive spring
    Elongation dummyElongation;
    dummyElongation.reset();
    Elongation* elongation_here_new = &dummyElongation;

    // check if spring already exists
    for (int i = 0; i < partJ->springs[t].size(); ++i) {
        const unsigned int springIndexI = partJ->springs[t][i].indexI;
        // if spring already exist, get it (otherwise it remains empty)
        if (springIndexI == indexI) { ////////////////////
            //cout<<" spring found"<<partj->particleIndex<<endl;
            elongation_here_new = &partJ->springs[t][i];
            elongation_here_new->active = true;
            break;
        }
    }

    // if no spring has been found,
    if (elongation_here_new->active == false) {
        dummyElongation.active = true;
        dummyElongation.indexJ = partJ->particleIndex;
        dummyElongation.indexI = indexI; //////////////////
        partJ->springs[t].push_back(dummyElongation);
        elongation_here_new = &partJ->springs[t].back();
    }

    return elongation_here_new;

}

void DEM::particleParticleContacts() {

    for (unsIntList::iterator ipi = neighborTable.begin(); ipi != neighborTable.end(); ipi = ipi + 2) {
        // couple of contact candidates
        unsIntList::iterator ipj = ipi + 1;
        // pointers to particles
        ASSERT(*ipi>*ipj);

        particle *partI = &particles[*ipi];
        particle *partJ = &particles[*ipj];

        // checking for overlap
        const double ri = partI->r;
        const double rj = partJ->r;
        const double sigij = ri + rj;
        const double sigij2 = sigij*sigij;
        // distance between centers
        const tVect vectorDistance = partJ->x0 - partI->x0;
        const double distance2 = vectorDistance.norm2();

        // check for contact
        if (distance2 < sigij2) {

            // pointer to elongation, initially pointing to an empty spring
            Elongation* elongation_here_new;
            if (staticFrictionSolve) {
                const unsigned int indexI = partI->particleIndex;
                elongation_here_new = findSpring(0, indexI, partJ);
            }

            particleParticleCollision(partI, partJ, vectorDistance, elongation_here_new);

        }
    }
}

void DEM::wallParticleContacts() {

    // to keep conventions, the index i refers to the wall, and j the particle
    //    unsigned int elongIndex = 0;

    for (unsIntList::iterator ip = nearWallTable.begin(); ip != nearWallTable.end(); ip = ip + 2) {
        // couple of contact candidates
        unsIntList::iterator iwall = ip + 1;
        // particle
        particle *partJ = &particles[*ip];
        // wall
        wall *wallI = &walls[*iwall];
        
        // distance from wall (norm)
        const double distance = wallI->dist(partJ->x0);

        bool isInsideLimits = true;

        if (wallI->limited) {
            if (partJ->x0.dot(Xp) < wallI->xMin) {
                isInsideLimits = false;
            } else if (partJ->x0.dot(Xp) > wallI->xMax) {
                isInsideLimits = false;
            } else if (partJ->x0.dot(Yp) < wallI->yMin) {
                isInsideLimits = false;
            } else if (partJ->x0.dot(Yp) > wallI->yMax) {
                isInsideLimits = false;
            } else if (partJ->x0.dot(Zp) < wallI->zMin) {
                isInsideLimits = false;
            } else if (partJ->x0.dot(Zp) > wallI->zMax) {
                isInsideLimits = false;
            }

            // check if we are above the plane (if limited))
            if (isInsideLimits && distance < 0.0) {
                isInsideLimits = false;
            }

        }


        if (isInsideLimits) {
            // radius
            const double rj = partJ->r;

            // distance before contact
            const double overlap = rj - distance;

            if (overlap > 0.0) {

                // pointer to elongation, initially pointing to an empty spring
                Elongation* elongation_here_new;
                if (staticFrictionSolve) {
                    const unsigned int indexI = wallI->index;
                    elongation_here_new = findSpring(1, indexI, partJ);
                }
                wallParticleCollision(wallI, partJ, overlap, elongation_here_new);

            }
        }
    }
}

void DEM::cylinderParticelContacts() {
    //    unsigned int elongIndex = 0;
    // to keep conventions, the index i refers to the cylinder, and j the particle
    // cycling through the cylinders
    for (unsIntList::iterator ip = nearCylinderTable.begin(); ip != nearCylinderTable.end(); ip = ip + 2) {
        // couple of contact candidates
        unsIntList::iterator icyl = ip + 1;
        // particle
        particle *partJ = &particles[*ip];
        // cylinder
        cylinder *cylinderI = &cylinders[*icyl];

        // radius
        const double rj = partJ->r;
        // distance from wall (norm)
        const double distance = cylinderI->dist(partJ->x0);
        // distance before contact
        const double overlap = rj - distance;

        // distance to point 1 of axis
        const tVect p1dist = partJ->x0 - cylinderI->p1;
        // same but projected on the axis
        const tVect p1distax = (p1dist.dot(cylinderI->naxes)) * cylinderI->naxes;
        // distance of point from cylinder axis
        const tVect p1distcylinder = p1distax - p1dist;

        // check for contact
        if (overlap > 0.0) {
            // pointer to elongation, initially pointing to an empty spring
            Elongation* elongation_here_new;
            if (staticFrictionSolve) {
                const unsigned int indexI = cylinderI->index;
                elongation_here_new = findSpring(2, indexI, partJ);
            }
            cylinderParticleCollision(cylinderI, partJ, overlap, elongation_here_new);

        }
    }
}

void DEM::objectParticleContacts() {
    //    unsigned int elongIndex = 0;
    // to keep conventions, the index i refers to the object, and j the particle
    for (unsIntList::iterator ip = nearObjectTable.begin(); ip != nearObjectTable.end(); ip = ip + 2) {
        // couple of contact candidates
        unsIntList::iterator iobj = ip + 1;
        // particle
        particle *partJ = &particles[*ip];
        // object
        object *objectI = &objects[*iobj];

        // radius
        const double rj = partJ->r;
        // distance from object (vector)
        const tVect vectorDistance = partJ->x0 - objectI->x0;
        // distance from object (norm)
        const double distance = vectorDistance.norm();
        // distance before contact
        const double overlap = rj + objectI->r - distance;

        if (overlap > 0.0) {
            // pointer to elongation, initially pointing to an empty spring
            Elongation* elongation_here_new;
            if (staticFrictionSolve) {
                const unsigned int indexI = objectI->index;
                elongation_here_new = findSpring(3, indexI, partJ);
            }
            objectParticleCollision(objectI, partJ, vectorDistance, elongation_here_new);
        }
    }
}

inline void DEM::particleParticleCollision(const particle *partI, const particle *partJ, const tVect& vectorDistance, Elongation* elongation_new) {

    // pointers to elements
    elmt *elmtI = &elmts[partI->clusterIndex];
    elmt *elmtJ = &elmts[partJ->clusterIndex];

    // NORMAL FORCE ///////////////////////////////////////////////////////////////////

    // geometry /////////////////////////////
    const double radI = partI->r;
    const double radJ = partJ->r;
    // distance norm
    const double distance = vectorDistance.norm();
    // overlap
    const double overlap = radI + radJ - distance;
    // relative velocity
    const tVect relVel = partJ->x1 - partI->x1;
    // first local unit vector (normal)
    const tVect en = vectorDistance / distance;
    // relative normal velocity (modulus)
    const double normRelVel = relVel.dot(en);
    // relative normal velocity
    const tVect normalRelVel = en*normRelVel;
    // effective mass
    const double effMass = elmtI->m * elmtJ->m / (elmtI->m + elmtJ->m);
    // effective radius
    const double effRad = radI * radJ / (radI + radJ);
    /*cout<<"overlap "<<overlap<<endl;
    cout<<"effRad "<<effRad<<" effMass "<<effMass<<" normRelVel "<<normRelVel<<" en";
    en.show();
    vectorDistance.show();
    cout<<distance<<" "<<partI->particleIndex<<" "<<partJ->particleIndex;*/

    // force computation /////////////////////////////////
    const double normNormalForce = normalContact(overlap, normRelVel, effRad, effMass);
    ASSERT(normNormalForce >= 0.0);

    //    // Overlap elastic potential energy
    //    energy.elastic+=0.5*fn*xi*xi;

    // force application ///////////////////////////////////
    // vectorial normal force
    const tVect normalForce = en*normNormalForce;

    // cluster geometry
    // vectorized radii
    const tVect vecRadI = radI*en;
    const tVect vecRadj = -radJ*en;
    // vectorized distance contactPoint-center of cluster
    const tVect centerDistI = vecRadI + (partI->radiusVec);
    const tVect centerDistJ = vecRadj + (partJ->radiusVec);

    if (!partI->isGhost) {
        elmtI->FParticle = elmtI->FParticle - normalForce;
        elmtI->solidIntensity += normalForce.abs();
        //  moment generated in non-spherical particles
        if (elmtI->size > 1) {
            elmtI->MParticle = elmtI->MParticle - centerDistI.cross(normalForce);
        }
    }
    if (!partJ->isGhost) {
        elmtJ->FParticle = elmtJ->FParticle + normalForce;
        elmtJ->solidIntensity += normalForce.abs();
        //  moment generated in non-spherical particles
        if (elmtJ->size > 1) {
            elmtJ->MParticle = elmtJ->MParticle + centerDistJ.cross(normalForce);
        }
    }

    // TANGENTIAL FORCE ///////////////////////////////////////////////////////////////////

    // angular velocities (global reference frame)
    const tVect wI = elmtI->wpGlobal; //2.0*quat2vec( elmtI->qp1.multiply( elmtI->qp0.adjoint() ) );
    const tVect wJ = elmtJ->wpGlobal; //2.0*quat2vec( elmtJ->qp1.multiply( elmtJ->qp0.adjoint() ) );
    // relative velocity at contact point
    const tVect relVelContact = relVel - wI.cross(vecRadI) + wJ.cross(vecRadj);
    // tangential component of relative velocity
    const tVect tangRelVelContact = relVelContact - normalRelVel;
    // norm of tangential velocity
    const double normTangRelVelContact = tangRelVelContact.norm();
    // checking if there is any tangential motion
    if (normTangRelVelContact != 0.0) {

        // update spring
        if (staticFrictionSolve) {

            const double elong_old_norm = elongation_new->e.norm();

            double normElong = elongation_new->e.dot(en);
            const tVect normalElong = en*normElong;

            elongation_new->e = elongation_new->e - normalElong;

            if (elongation_new->e.norm() != 0) {
                const double scaling = elong_old_norm / elongation_new->e.norm();
                elongation_new->e = elongation_new->e*scaling;
            } else {
                const double scaling = 0.0;
                elongation_new->e = elongation_new->e*scaling;
            }
        }

        //elongation.e.show();
        tVect tangForce = FRtangentialContact(tangRelVelContact, normNormalForce, overlap, effRad, effMass, elongation_new, sphereMat.frictionCoefPart,
                sphereMat.linearStiff, sphereMat.viscTang);

        // torque updating
        if (!partI->isGhost) {
            elmtI->MParticle = elmtI->MParticle + centerDistI.cross(tangForce);
            elmtI->FSpringP = elmtI->FSpringP + tangForce;
            elmtI->FParticle = elmtI->FParticle + tangForce;
            elmtI->solidIntensity += tangForce.abs();
            if (staticFrictionSolve) {
                elmtI->slippingCase = elongation_new->slippingCase;
            }

        }
        if (!partJ->isGhost) {
            elmtJ->MParticle = elmtJ->MParticle - centerDistJ.cross(tangForce);
            elmtJ->FSpringP = elmtJ->FSpringP - tangForce;
            elmtJ->FParticle = elmtJ->FParticle - tangForce;
            elmtJ->solidIntensity += tangForce.abs();
            if (staticFrictionSolve) {
                elmtJ->slippingCase = elongation_new->slippingCase;
            }
        }

    }


    //ROLLING

    tVect rollingMoment = rollingContact(wI, wJ, effRad, normNormalForce,
            sphereMat.rollingCoefPart);


    if (!partI->isGhost) {
        elmtI->MRolling = elmtI->MRolling - rollingMoment;

    }
    if (!partJ->isGhost) {
        elmtJ->MRolling = elmtJ->MRolling + rollingMoment;
    }

    // save overlap
    if (overlap > elmtI->maxOverlap) {
        elmtI->maxOverlap = overlap;
    }
    if (overlap > elmtJ->maxOverlap) {
        elmtJ->maxOverlap = overlap;
    }

    // updating connectivity
    elmtI->coordination += 1;
    elmtJ->coordination += 1;
}

inline void DEM::wallParticleCollision(wall *wallI, const particle *partJ, const double& overlap, Elongation* elongation_new) {

    // pointers to element
    elmt *elmtJ = &elmts[partJ->clusterIndex];

    // NORMAL FORCE ///////////////////////////////////////////////////////////////////

    // geometry ///////////////
    // particle radius
    const double radJ = partJ->r;
    // first local unit vector (normal)
    const tVect en = wallI->n;
    // speed of the wall at contact point
    const tVect contactPointVelocity = wallI->getSpeed(partJ->x0); // fix this, contact point not defined
    // relative velocity
    const tVect relVel = partJ->x1 - contactPointVelocity;
    // relative normal velocity (modulus)
    const double normRelVel = relVel.dot(en);
    // relative normal velocity
    const tVect normalRelVel = en*normRelVel;

    // force computation /////////////////////////////////
    double normNormalForce = normalContact(overlap, normRelVel, radJ, elmtJ->m); // removed 2.0 *

    //    switch (problemName) {
    //        case TRIAXIAL:
    //        {
    //            if (wallI->index > 5) {
    //                normNormalForce /= 10;
    //            }
    //        }
    //    }


    // Overlap elastic potential energy
    //                    energy.elastic+=fn*xj/2.5;

    // force application ///////////////////////////////////
    // vectorial normal force
    const tVect normalForce = en*normNormalForce;

    // cluster geometry
    // vectorized radius
    const tVect vecRadJ = -radJ*en;
    // vectorized distance contactPoint-center of cluster
    tVect centerDistJ = vecRadJ;
    if (elmtJ->size > 1) {
        // vectorized distance contactPoint-center of cluster
        centerDistJ = centerDistJ + (partJ->x0 - elmtJ->xp0);
    }

    // force updating
    elmtJ->FWall = elmtJ->FWall + normalForce;
    wallI->FParticle = wallI->FParticle - normalForce;
    elmtJ->solidIntensity += normalForce.abs();
    // torque updating
    if (elmtJ->size > 1) {
        elmtJ->MWall = elmtJ->MWall + centerDistJ.cross(normalForce);
    }

    // TANGENTIAL FORCE ///////////////////////////////////////////////////////////////////

    // angular velocities (global reference frame)
    const tVect wJ = elmtJ->wpGlobal; //2.0*quat2vec( elmtJ->qp1.multiply( elmtJ->qp0.adjoint() ) );
    // relative velocity at contact point
    const tVect relVelContact = relVel + wJ.cross(vecRadJ);
    // tangential component of relative velocity
    const tVect tangRelVelContact = relVelContact - normalRelVel;
    // norm of tangential velocity
    const double normTangRelVelContact = tangRelVelContact.norm();


    // checking if there is any tangential motion
    if (normTangRelVelContact != 0.0) {

        // update spring
        if (staticFrictionSolve) {
            const double elong_old_norm = elongation_new->e.norm();

            double normElong = elongation_new->e.dot(en);
            const tVect normalElong = en*normElong;

            elongation_new->e = elongation_new->e - normalElong;

            if (elongation_new->e.norm() != 0) {
                const double scaling = elong_old_norm / elongation_new->e.norm();
                elongation_new->e = elongation_new->e*scaling;
            } else {
                const double scaling = 0.0;
                elongation_new->e = elongation_new->e*scaling;
            }
        }

        tVect tangForce = FRtangentialContact(tangRelVelContact, normNormalForce, overlap, radJ, elmtJ->m, elongation_new, sphereMat.frictionCoefWall,
                sphereMat.linearStiff, sphereMat.viscTang);

        // torque updating
        elmtJ->MWall = elmtJ->MWall - centerDistJ.cross(tangForce);
        // force updating
        elmtJ->FSpringW = elmtJ->FSpringW - tangForce;
        elmtJ->FWall = elmtJ->FWall - tangForce;
        elmtJ->solidIntensity += tangForce.abs();
        wallI->FParticle = wallI->FParticle + tangForce;

        if (staticFrictionSolve) {
            elmtJ->slippingCase = elongation_new->slippingCase;
        }

    }
    //ROLLING

    tVect rollingMoment = rollingContact(Zero, wJ, radJ, normNormalForce,
            2.0 * sphereMat.rollingCoefPart);

    if (!partJ->isGhost) {
        elmtJ->MRolling = elmtJ->MRolling + rollingMoment;
    }

    // save overlap
    if (overlap > elmtJ->maxOverlap) {
        elmtJ->maxOverlap = overlap;
    }

    // updating connectivity
    elmtJ->coordination += 1;

}

inline void DEM::cylinderParticleCollision(cylinder *cylinderI, const particle *partJ, const double& overlap, Elongation* elongation_new) {

    // pointers to element
    elmt *elmtJ = &elmts[partJ->clusterIndex];

    // NORMAL FORCE ///////////////////////////////////////////////////////////////////

    // geometry ///////////////
    // particle radius
    const double radJ = partJ->r;
    // vectorial distance
    const tVect vecDistance = cylinderI->vecDist(partJ->x0);
    // contact point
    const tVect contactPoint = partJ->x0 - vecDistance;
    // first local unit vector (normal)
    const tVect en = vecDistance / (radJ - overlap);
    // speed of the cylinder at contact point
    const tVect contactPointVelocity = cylinderI->getSpeed(contactPoint);
    // relative velocity
    const tVect relVel = partJ->x1 - contactPointVelocity;
    // relative normal velocity (modulus)
    const double normRelVel = relVel.dot(en);
    // relative normal velocity
    const tVect normalRelVel = en*normRelVel;

    // force computation /////////////////////////////////
    const double normNormalForce = normalContact(overlap, normRelVel, radJ, elmtJ->m); // was 2.0 * overlap

    // Overlap elastic potential energy
    //                    energy.elastic+=fn*xj/2.5;

    // force application ///////////////////////////////////
    // vectorial normal force
    const tVect normalForce = en*normNormalForce;

    // cluster geometry
    // vectorized radius
    const tVect vecRadJ = -radJ*en;
    // vectorized distance contactPoint-center of cluster
    const tVect centerDistJ = vecRadJ + (partJ->x0 - elmtJ->xp0);

    // force updating
    elmtJ->FCylinder = elmtJ->FCylinder + normalForce;
    cylinderI->FParticle = cylinderI->FParticle-normalForce;
    elmtJ->solidIntensity += normalForce.abs();
    // wallI->FParticle=wallI->FParticle-fnv;
    // torque updating
    if (elmtJ->size > 1) {
        elmtJ->MCylinder = elmtJ->MCylinder + centerDistJ.cross(normalForce);
    }

    // TANGENTIAL FORCE ///////////////////////////////////////////////////////////////////

    // angular velocities (global reference frame)
    const tVect wJ = elmtJ->wpGlobal; //2.0*quat2vec( elmtJ->qp1.multiply( elmtJ->qp0.adjoint() ) );
    // relative velocity at contact point
    const tVect relVelContact = relVel + wJ.cross(vecRadJ); // couldn't we just use elmtJ.w?
    // tangential component of relative velocity
    const tVect tangRelVelContact = relVelContact - normalRelVel;
    // norm of tangential velocity
    const double normTangRelVelContact = tangRelVelContact.norm();
    // checking if there is any tangential motion
    if (normTangRelVelContact != 0.0) {

        // update spring
        if (staticFrictionSolve) {
            const double elong_old_norm = elongation_new->e.norm();

            double normElong = elongation_new->e.dot(en);
            const tVect normalElong = en*normElong;

            elongation_new->e = elongation_new->e - normalElong;

            if (elongation_new->e.norm() != 0) {
                const double scaling = elong_old_norm / elongation_new->e.norm();
                elongation_new->e = elongation_new->e*scaling;
            } else {
                const double scaling = 0.0;
                elongation_new->e = elongation_new->e*scaling;
            }
        }

        tVect tangForce = FRtangentialContact(tangRelVelContact, normNormalForce, overlap, radJ, elmtJ->m, elongation_new, sphereMat.frictionCoefWall,
                sphereMat.linearStiff, sphereMat.viscTang);

        // torque updating
        elmtJ->MCylinder = elmtJ->MCylinder - centerDistJ.cross(tangForce);
        // force updating
        elmtJ->FCylinder = elmtJ->FCylinder - tangForce;
        elmtJ->solidIntensity += tangForce.abs();
        cylinderI->FParticle = cylinderI->FParticle + tangForce;
        
        if (staticFrictionSolve) {
            elmtJ->slippingCase = elongation_new->slippingCase;
        }
        
    }
    //ROLLING

    tVect rollingMoment = rollingContact(Zero, wJ, radJ, normNormalForce,
            2.0 * sphereMat.rollingCoefPart);



    if (!partJ->isGhost) {
        elmtJ->MRolling = elmtJ->MRolling + rollingMoment;
    }

    // save overlap
    if (overlap > elmtJ->maxOverlap) {
        elmtJ->maxOverlap = overlap;
    }

    // updating connectivity
    elmtJ->coordination += 1;
}

inline void DEM::objectParticleCollision(object *objectI, const particle *partJ, const tVect& vectorDistance, Elongation* elongation_new) {

    // pointers to element
    elmt *elmtJ = &elmts[partJ->clusterIndex];

    // NORMAL FORCE ///////////////////////////////////////////////////////////////////

    // geometry ///////////////
    // particle radius
    const double radJ = partJ->r;
    // distance from object (norm)
    const double distance = vectorDistance.norm();
    // distance before contact
    double overlap = radJ + objectI->r - distance;
    // first local unit vector (normal)
    const tVect en = (1.0 / distance) * vectorDistance;
    // speed of the wall at contact point
    const tVect contactPointVelocity = objectI->x1; // fix this, contact point not defined
    // relative velocity
    const tVect relVel = partJ->x1 - contactPointVelocity;
    // relative normal velocity (modulus)
    const double normRelVel = relVel.dot(en);
    // relative normal velocity
    const tVect normalRelVel = en*normRelVel;

    // force computation /////////////////////////////////
    const double normNormalForce = normalContact(overlap, normRelVel, radJ, elmtJ->m); // was 2.0 * overlap

    // Overlap elastic potential energy
    //                    energy.elastic+=fn*xj/2.5;

    // force application ///////////////////////////////////
    // vectorial normal force
    const tVect normalForce = en*normNormalForce;

    // cluster geometry
    // vectorized radius
    const tVect vecRadJ = -radJ*en;
    // vectorized distance contactPoint-center of cluster
    const tVect centerDistJ = vecRadJ + (partJ->x0 - elmtJ->xp0);

    // force updating
    elmtJ->FWall = elmtJ->FWall + normalForce;
    objectI->FParticle = objectI->FParticle - normalForce;
    elmtJ->solidIntensity += normalForce.abs();

    // torque updating
    if (elmtJ->size > 1) {
        elmtJ->MWall = elmtJ->MWall + centerDistJ.cross(normalForce);
    }

    // TANGENTIAL FORCE ///////////////////////////////////////////////////////////////////

    // angular velocities (global reference frame)
    const tVect wJ = elmtJ->wpGlobal; //2.0*quat2vec( elmtJ->qp1.multiply( elmtJ->qp0.adjoint() ) );
    // relative velocity at contact point
    const tVect relVelContact = relVel + wJ.cross(vecRadJ);
    // tangential component of relative velocity
    const tVect tangRelVelContact = relVelContact - normalRelVel;
    // norm of tangential velocity
    const double normTangRelVelContact = tangRelVelContact.norm();
    // checking if there is any tangential motion
    if (normTangRelVelContact != 0.0) {
        // update spring
        if (staticFrictionSolve) {

            const double elong_old_norm = elongation_new->e.norm();

            double normElong = elongation_new->e.dot(en);
            const tVect normalElong = en*normElong;

            elongation_new->e = elongation_new->e - normalElong;

            if (elongation_new->e.norm() != 0) {
                const double scaling = elong_old_norm / elongation_new->e.norm();
                elongation_new->e = elongation_new->e*scaling;
            } else {
                const double scaling = 0.0;
                elongation_new->e = elongation_new->e*scaling;
            }
        }

        tVect tangForce = FRtangentialContact(tangRelVelContact, normNormalForce, overlap, radJ, elmtJ->m, elongation_new, sphereMat.frictionCoefObj,
                sphereMat.linearStiff, sphereMat.viscTang);

        // torque updating
        elmtJ->MWall = elmtJ->MWall - centerDistJ.cross(tangForce);
        // force updating
        elmtJ->FSpringW = elmtJ->FSpringW - tangForce;
        elmtJ->FWall = elmtJ->FWall - tangForce;
        elmtJ->solidIntensity += tangForce.abs();
        objectI->FParticle = objectI->FParticle + tangForce;

    }
    //ROLLING

    tVect rollingMoment = rollingContact(Zero, wJ, radJ, normNormalForce,
            2.0 * sphereMat.rollingCoefPart);


    if (!partJ->isGhost) {
        elmtJ->MRolling = elmtJ->MRolling + rollingMoment;
    }

    // save overlap
    if (overlap > elmtJ->maxOverlap) {
        elmtJ->maxOverlap = overlap;
    }

    // updating connectivity
    elmtJ->coordination += 1;
}

double DEM::normalContact(const double& overlap, const double& vrelnnorm, const double& effRad, const double& effMass) const {

    // total normal force
    double fn = 0.0;

    switch (sphereMat.contactModel) {
        case HERTZIAN:
        {
            // square root of effective radius
            const double sqrtEffRad = sqrt(effRad);
            // square root of overlap
            const double sqrtOverlap = sqrt(overlap);
            // normal stiffness = 2/3*(Reff^1/2)*Y/(1-nu²)
            const double kn = sphereMat.knConst * sqrtEffRad*sqrtOverlap;
            // damping
            const double gamman = 2.0 * sphereMat.dampCoeff * sqrt(kn * effMass); // 2 was missing
            const double dumpfn = -gamman*vrelnnorm;
            // elastic normal force (Hertzian contact)
            const double elasticfn = kn*overlap;
            // total normal force (limit to repulsive force)
            fn = std::max(elasticfn + dumpfn, 0.0);
            break;
        }
        case LINEAR:
        {
            // elastic normal force (linear contact)
            const double elasticfn = sphereMat.linearStiff*overlap;
            // damping
            const double gamman = 2.0 * sphereMat.dampCoeff * sqrt(sphereMat.linearStiff * effMass);
            const double dumpfn = -gamman*vrelnnorm;
            // total normal force (limit to repulsive force)
            fn = std::max(elasticfn + dumpfn, 0.0);
            break;
        }
    }

    return fn;

}

tVect DEM::FRtangentialContact(const tVect& tangRelVelContact, const double& fn, const double& overlap, const double& effRad,
        const double& effMass, Elongation* elongation_new, const double& friction, const double& tangStiff, const double& viscTang) {
    // tangential force
    tVect fs = tVect(0.0, 0.0, 0.0);
    // tangent stiffness
    double ks = 0.0;
    switch (sphereMat.contactModel) {
        case HERTZIAN:
        {
            static const double power = 1.0 / 3.0;
            // square root of effective radius
            const double sqrtEffRad = sqrt(effRad);
            // square root of overlap
            const double sqrtOverlap = sqrt(overlap);
            // normal stiffness (from normal contac)
            const double kn = sphereMat.knConst * sqrtEffRad*sqrtOverlap;
            // tangent stiffness -> not physically sound, static friction is actually missing
            ks = 2.0 / 7.0 * kn; //sphereMat.ksConst * sqrtEffRad * std::pow(std::abs(fn), power);
            break;
        }
        case LINEAR:
        {
            ks = 2.0 / 7.0 * tangStiff;
            break;
        }
    }

    // maximum force due to static friction
    const double fsMax = friction*fn;
    const double fdMax = 0.9 * friction*fn;

    // viscous force
    const tVect viscousForce = 2.0 * viscTang * sqrt(effMass * ks) * tangRelVelContact;

    // model with static friction
    if (staticFrictionSolve) {

        const tVect F_control = elongation_new->e * ks + viscousForce;

        const double F_control_norm = F_control.norm();

        tVect et = tVect(0.0, 0.0, 0.0);
        if (F_control_norm != 0) {
            et = F_control / F_control_norm;
        }

        if (elongation_new->slipping) {
            if (F_control_norm < fdMax) {
                // was slipping, goes back to static
                fs = F_control;
                elongation_new->e = elongation_new->e + tangRelVelContact*deltat;
                elongation_new->slipping = false;
                elongation_new->slippingCase = 1;
            } else {
                // keeps slipping
                fs = fdMax*et;
                const tVect f = fdMax * et - viscousForce;
                elongation_new->e = f / ks;
                elongation_new->slippingCase = 2;
            }

        } else {

            if (F_control_norm < fsMax) {
                // keeps being static
                fs = F_control;
                elongation_new->e = elongation_new->e + tangRelVelContact*deltat;
                elongation_new->slippingCase = 3;

            } else {
                // was static, begins slipping
                fs = fdMax*et;
                elongation_new->slipping = true;
                const tVect f = fdMax * et - viscousForce;
                elongation_new->e = f / ks;
                elongation_new->slippingCase = 4;
            }

        }


    }// model without static friction
    else {
        const double viscousForceNorm = viscousForce.norm();
        const double viscousForceReal = std::min(viscousForceNorm, fsMax);
        fs = viscousForce / viscousForceNorm*viscousForceReal;
    }



    return fs;
}

double DEM::tangentialContact(const double& vreltNorm, const double& fn, const double& effRad, const double& effMass, const double& friction) const {

    // tangential force
    double fs = 0.0;

    switch (sphereMat.contactModel) {
        case HERTZIAN:
        {
            static const double power = 1.0 / 3.0;
            // square root of effective radius
            const double sqrtEffRad = sqrt(effRad);
            // tangent stiffness -> not physically sound, static friction is actually missing
            const double ks = sphereMat.ksConst * sqrtEffRad * std::pow(std::abs(fn), power);
            // maximum force due to dynamic friction
            const double fsMax = friction*fn;
            // viscous force
            const double gammas = 2.0 * sphereMat.viscTang * sqrt(effMass * ks);
            const double dumpfs = gammas*vreltNorm;
            // old version
            // const double viscousForce=2.0*sphereMat.viscTang*sqrt(effMass*ks)*vreltNorm;
            // capping with dynamic friction
            fs = std::min(dumpfs, fsMax);
            //ASSERT(fsMax>=0.0);
            break;
        }
        case LINEAR:
        {
            // tangent stiffness -> not physically sound, static friction is actually missing
            const double ks = sphereMat.linearStiff;
            // maximum force due to dynamic friction
            const double fsMax = friction*fn;
            // viscous force
            const double gammas = 2.0 * sphereMat.viscTang * sqrt(effMass * ks);
            const double dumpfs = gammas*vreltNorm;
            // old version
            //const double viscousForce=2.0*sphereMat.viscTang*sqrt(effMass*ks)*vreltNorm;
            // capping with dynamic friction
            fs = std::min(dumpfs, fsMax);
            break;
        }
    }

    return fs;

}

tVect DEM::rollingContact(const tVect& wI, const tVect& wJ, const double& effRad, const double& fn, const double& rolling) {

    const tVect wrel = wI - wJ;
    const double wrel_norm = wrel.norm();

    tVect wIJ = tVect(0.0, 0.0, 0.0);
    if (wrel_norm != 0.0) {
        wIJ = wrel / wrel_norm;
    }

    const tVect fr = wIJ * rolling * fn*effRad;

    return fr;
}

// energy functions

void DEM::updateEnergy(double& totalKineticEnergy) {

    //resetting energy
    particleEnergy.reset();

    // kinetic energy
    double tKin(0.0), rKin(0.0), mass(0.0);
    for (int n = 0; n < elmts.size(); ++n) {
        if (elmts[n].active) {
            tKin += 0.5 * elmts[n].m * elmts[n].x1.norm2();
            // adjoint of orientation quaternion
            //const tQuat q0adj=elmts[n].q0.adjoint();
            // rotational velocity (body-fixed reference frame)
            //const tVect w=wLocal;//2.0*quat2vec( q0adj.multiply( elmts[n].q1 ) );
            //        Tvect w=2.0*quat2vec( elmts[i].q1.multiply( elmts[i].q0.adjoint() ) );
            const tVect wSquare = elmts[n].wLocal.compProd(elmts[n].wLocal);
            rKin += elmts[n].I.dot(wSquare);
            mass += elmts[n].m;
        }
    }

    particleEnergy.rotKin = rKin;
    particleEnergy.trKin = tKin;
    particleEnergy.mass = mass;

    // potential energy
    // defining reference plane
    wall zeroWall;
    double g(0.0);
    const double gravityNorm = demF.norm();
    if (gravityNorm != 0.0) {
        zeroWall.n = -1.0 * demF / gravityNorm;
        zeroWall.p = tVect(0.0, 0.0, 0.0);
        if (zeroWall.n.dot(Xp) < 0) {
            zeroWall.p += tVect(demSize[0], 0.0, 0.0);
        }
        if (zeroWall.n.dot(Yp) < 0) {
            zeroWall.p += tVect(0.0, demSize[1], 0.0);
        }
        if (zeroWall.n.dot(Zp) < 0) {
            zeroWall.p += tVect(0.0, 0.0, demSize[2]);
        }
        for (int n = 0; n < elmts.size(); n++) {
            const double heigth = zeroWall.dist(elmts[n].x0);
            g += elmts[n].m * heigth*gravityNorm;
        }
    }
    particleEnergy.grav = g;

    // elastic is missing

    // total
    particleEnergy.updateTotal();
    totalKineticEnergy = totalKineticEnergy + particleEnergy.rotKin + particleEnergy.trKin;

}
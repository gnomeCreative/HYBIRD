#include "DEMParams.h"

#include <string>

#include "getpot.h"
#include "macros.h"
#include "Object2.h"
#include "Element2.h"

void DEMParams::discreteElementGet(GetPot& configFile, GetPot& commandLine, Element2& elmts, Object2& objects) {
    // getting material properties
    PARSE_CLASS_MEMBER(configFile, sphereMat.density, "particleDensity", 0.0);
    ASSERT(sphereMat.density > 0.0);

    std::string contactModelString;
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
    std::string particleFile;
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
    particleFileID>>elmts.count;
    // Allocate memory for elements
    elmts.memoryAlloc<CPU>(elmts.count);
    elmts.componentsIndex[0] = 0;
    for (int n = 0; n < elmts.count; ++n) {
        elmt dummyElmt;

        // import variables
        particleFileID >> elmts.index[n];
        particleFileID >> elmts.size[n];
        // Calculate components index
        elmts.componentsIndex[n + 1] = elmts.componentsIndex[n] + elmts.size[n];
        particleFileID >> elmts.radius[n];
        elmts.radius[n] = elmts.radius[n] * scale;
        // position
        particleFileID>>elmts.x0[n].x;
        particleFileID>>elmts.x0[n].y;
        particleFileID>>elmts.x0[n].z;
        elmts.x0[n] *= scale;
        elmts.x0[n] += translate;
        // translational velocity
        particleFileID>>elmts.x1[n].x;
        particleFileID>>elmts.x1[n].y;
        particleFileID>>elmts.x1[n].z;
        // rotational velocity
        particleFileID>>elmts.w0[n].x;
        particleFileID>>elmts.w0[n].y;
        particleFileID>>elmts.w0[n].z;
        // orientation
        particleFileID>>elmts.q0[n].q0;
        particleFileID>>elmts.q0[n].q1;
        particleFileID>>elmts.q0[n].q2;
        particleFileID>>elmts.q0[n].q3;
        // translational velocity (in quaternion rates))
        particleFileID>>elmts.q1[n].q0;
        particleFileID>>elmts.q1[n].q1;
        particleFileID>>elmts.q1[n].q2;
        particleFileID>>elmts.q1[n].q3;
        elmts.active[n] = true;
    }
    // Allocate memory for componentsData
    elmts.allocComponentsData();
    cout << " done" << endl;

    // objects initial state //////////////////////
    string objectFile;
    PARSE_CLASS_MEMBER(configFile, objectFile, "objectFile", "objects.dat");
    ifstream objectFileID;
    objectFileID.open(objectFile.c_str(), ios::in);
    ASSERT(objectFileID.is_open());
    cout << "Reading " << objectFile.c_str() << "...";
    objectFileID>>objects.count;
    // Allocate memory for elements
    objects.allocObjects<CPU>(objects.count);
    for (int n = 0; n < objects.count; ++n) {
        // import variables
        objectFileID >> objects.index[n];
        // this is used to identify objects belonging to different groups
        objectFileID >> objects.ElID[n]; // must be one
        objectFileID >> objects.r[n];
        objectFileID >> objects.x0[n].x;
        objectFileID >> objects.x0[n].y;
        objectFileID >> objects.x0[n].z;
        objectFileID >> objects.x1[n].x;
        objectFileID >> objects.x1[n].y;
        objectFileID >> objects.x1[n].z;
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
        objects.originalIndex[n] = objects.index[n];
    }
    cout << " done" << endl;

    // numerical viscosity for stability
    PARSE_CLASS_MEMBER(configFile, numVisc, "numVisc", 0.0);
    // set multiplication number (sets the number of DEM steps between two fluid steps)
    PARSE_CLASS_MEMBER(configFile, multiStep, "multiStep", 1);
    // set ratio between time step and estimated duration of contacts (only if multiStep=0)
    PARSE_CLASS_MEMBER(configFile, criticalRatio, "criticalRatio", 0.1);
}

void DEMParams::compositeProperties() {

    vecList prototype1, prototype2, prototype3, prototype4;


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

}
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <csignal>
#include <limits>

#include "getpot.h"

#include "IO.h"
#include "DEM.h"
#include "LB.h"

/*
 *      HYBIRD
 *      DEM-LBM coupling code
 *          
 *      author: Alessandro Leonardi
 *
 */


enum ExitCode {
    UNFINISHED = -1, SUCCESS = 0, TIME_LIMIT_REACHED, SIGNAL_CAUGHT, ERROR
};




ExitCode exit_code = UNFINISHED;
ProblemName problemName = NONE;

//void catchSignal(int sig) {
//    std::cout << "Received signal " << sig << ": " << strsignal(sig) << " " << std::endl;
//    exit_code = (SIGNAL_CAUGHT > exit_code ? SIGNAL_CAUGHT : exit_code);
//}

void print_help() {
    cout << "Usage: " << " [-c configFile] [-d datadir] [-n resultFolderName (opt. time)] [options}\n";
}

void goCycle(IO& io, DEM& dem, LB& lb) {

    // advance one step in time
    io.realTime += lb.unit.Time;
    ++io.currentTimeStep;
    ++lb.time;

    dem.evolveBoundaries();
//    dem.evolveObj();
    //cout<<"1"<<endl;
    if (io.demSolver) {

        dem.discreteElementStep();
    }

    if (io.lbmSolver && dem.demTime >= dem.demInitialRepeat) {

        if (io.demSolver) {

            lb.latticeBoltzmannCouplingStep(dem.newNeighborList, dem.elmts, dem.particles);
        }

        lb.latticeBolzmannStep(dem.elmts, dem.particles, dem.walls, dem.objects);

        // Lattice Boltzmann core steps

        if (lb.freeSurface) {

            lb.latticeBoltzmannFreeSurfaceStep();
        }

    } else if (io.lbmSolver && dem.demTime <= dem.demInitialRepeat) {
        if (io.demSolver) {

            lb.latticeBoltzmannCouplingStep(dem.newNeighborList, dem.elmts, dem.particles);
        }
    }

    io.outputStep(lb, dem);


}

void parseCommandLine(IO& io, GetPot& commandLine) {

    //io.base_directory = command_line.follow("./", "--directory");

    // print help and exit if requested
    if (commandLine.search("-h") || commandLine.search("--help")) {
        print_help();
        exit(0);
    }

    // create the data directory in case it doesn't exist
    io.workDirectory = "data";
    if (commandLine.search("-d")) {
        io.workDirectory = commandLine.next(io.workDirectory);
    }

    // if the name of the simulation is specified create subfolder
    std::string programName = "name";
    if (commandLine.search("-n")) {
        // creating work directory
        programName = commandLine.next(programName);
        // if the specified name is "time" then create a folder named with the initial time
        if (programName == "time") {
            ostringstream convertTime;
            convertTime << std::setfill('0') << std::setw(4) << (io.now->tm_year + 1900)
                    << std::setfill('0') << std::setw(2) << (io.now->tm_mon + 1)
                    << std::setfill('0') << std::setw(2) << io.now->tm_mday
                    << "_"
                    << std::setfill('0') << std::setw(2) << io.now->tm_hour
                    << std::setfill('0') << std::setw(2) << io.now->tm_min
                    << std::setfill('0') << std::setw(2) << io.now->tm_sec;
            string simulationTime = convertTime.str();
            io.workDirectory = io.workDirectory + "/" + simulationTime;
            io.workDirectory = io.workDirectory.substr(0, io.workDirectory.size());
        } else {
            io.workDirectory = io.workDirectory + "/" + programName;
        }

        //        io.workDirectory = command_line.next(io.workDirectory);
    }

    // lbm configuration file
    io.configFileName = "config.cfg";
    if (commandLine.search("-c")) {
        io.configFileName = commandLine.next(io.configFileName.c_str());
        cout << "Using " << io.configFileName << " for LBM and DEM parameters\n";
        // copy the lbm configuration file into the data folder for later reproducibility
        int a;
        a = filesystem::create_directories(io.workDirectory);
        cout << "Work directory created = " << io.workDirectory << ". Result: " << a << "\n";

        std::system(("cp '" + io.configFileName + "' '" + io.workDirectory + "'").c_str());
    }
    // make sure the config files can be read
    std::ifstream configFile(io.configFileName.c_str());
    if (!configFile) {
        cout << "ERROR: Can't open config file \"" << io.configFileName << "\" for reading!\n";
        //        return ERROR;
    }
    configFile.close();
}

void parseConfigFile(IO& io, DEM& dem, LB& lb, GetPot& configFile, GetPot& commandLine) {

    cout << "Parsing input file" << endl;
    // PROBLEM NAME //////////////
    // necessary for hard coded sections of the code
    string problemNameString;
    PARSE_CLASS_MEMBER(configFile, problemNameString, "problemName", "none");
    if (problemNameString == "DRUM") problemName = DRUM;
    else if (problemNameString == "SHEARCELL") problemName = SHEARCELL;
    else if (problemNameString == "NONE") problemName = NONE;
    else if (problemNameString == "AVALANCHE") problemName = AVALANCHE;
    else if (problemNameString == "NET") problemName = NET;
    else if (problemNameString == "ZHOU") problemName = ZHOU;
    else if (problemNameString == "BARRIER") problemName = BARRIER;
    else if (problemNameString == "OPENBARRIER") problemName = OPENBARRIER;
    else if (problemNameString == "HONGKONG") problemName = HONGKONG;
    else if (problemNameString == "ESERCITAZIONE") problemName = ESERCITAZIONE;
    else if (problemNameString == "FILIPPO_SILOS") problemName = FILIPPO_SILOS;
    else if (problemNameString == "STVINCENT") problemName = STVINCENT;
    else if (problemNameString == "STAVA") problemName = STAVA;
    else if (problemNameString == "NIGRO") problemName = NIGRO;
    else if (problemNameString == "DAMBREAK") problemName = DAMBREAK;
    else if (problemNameString == "GRAY_DAMBREAK") problemName = GRAY_DAMBREAK;
    else if (problemNameString == "GRAY_DAMBREAK_2D") problemName = GRAY_DAMBREAK_2D;
    else if (problemNameString == "CAROLINE") problemName = CAROLINE;
    else if (problemNameString == "INCLINEFLOW") problemName = INCLINEFLOW;
    else if (problemNameString == "JOP") problemName = JOP;
    else if (problemNameString == "MANGENEY") problemName = MANGENEY;
    else if (problemNameString == "WILL") problemName = WILL;
    else if (problemNameString == "WILL_SETTLING") problemName = WILL_SETTLING;
    else if (problemNameString == "HK_SMALL") problemName = HK_SMALL;
    else if (problemNameString == "HK_LARGE") problemName = HK_LARGE;
    else if (problemNameString == "KELVIN") problemName = KELVIN;
    else if (problemNameString == "GRAY") problemName = GRAY;
    else if (problemNameString == "HOURGLASS") problemName = HOURGLASS;
    else if (problemNameString == "IERVOLINO") problemName = IERVOLINO;
    else if (problemNameString == "IERVOLINO_2D") problemName = IERVOLINO_2D;
    else if (problemNameString == "IERVOLINO_CYLINDERTEST") problemName = IERVOLINO_CYLINDERTEST;
    else if (problemNameString == "HEAP") problemName = HEAP;
    else if (problemNameString == "TRIAXIAL") problemName = TRIAXIAL;
    else if (problemNameString == "SHEARCELL2022") problemName = SHEARCELL2023;
    else if (problemNameString == "INTRUDER") problemName = INTRUDER;
    else if (problemNameString == "SEGUIN") problemName = SEGUIN;



    // GETTING SIMULATION PARAMETERS  /////////
    // DEM initial iterations
    PARSE_CLASS_MEMBER(configFile, dem.demInitialRepeat, "demInitialRepeat", 0.0);
    ASSERT(dem.demInitialRepeat >= 0);
    // LB iteration without coupling (for initial stability) - > time is frozen here
    PARSE_CLASS_MEMBER(configFile, lb.lbmInitialRepeat, "lbmInitialRepeat", 0);
    ASSERT(lb.lbmInitialRepeat >= 0);
    // maximum time variable value
    PARSE_CLASS_MEMBER(configFile, io.maxTime, "maxTime", 0.0);
    ASSERT(io.maxTime >= 0);
    // threshold energy for adaptive stop
    PARSE_CLASS_MEMBER(configFile, io.energyStopThreshold, "energyStopThreshold", -10.0);
    PARSE_CLASS_MEMBER(configFile, io.energyStop, "energyStop", false);

    // for time integration and output
    PARSE_CLASS_MEMBER(configFile, io.maximumTimeSteps, "maximumTimeSteps", 0);
    ASSERT(io.maximumTimeSteps >= 0);
    PARSE_CLASS_MEMBER(configFile, io.screenExpTime, "screenExpTime", 0.0);
    ASSERT(io.screenExpTime >= 0);
    PARSE_CLASS_MEMBER(configFile, io.fluidExpTime, "fluidExpTime", 0.0);
    ASSERT(io.fluidExpTime >= 0);
    PARSE_CLASS_MEMBER(configFile, io.fluidLagrangianExpTime, "fluidLagrangianExpTime", 0.0);
    ASSERT(io.fluidLagrangianExpTime >= 0);
    PARSE_CLASS_MEMBER(configFile, io.fluid2DExpTime, "fluid2DExpTime", 0.0);
    ASSERT(io.fluid2DExpTime >= 0);
    PARSE_CLASS_MEMBER(configFile, io.partExpTime, "partExpTime", 0.0);
    ASSERT(io.partExpTime >= 0);
    PARSE_CLASS_MEMBER(configFile, io.fluidRecycleExpTime, "fluidRecycleExpTime", 0.0);
    ASSERT(io.fluidRecycleExpTime >= 0);
    PARSE_CLASS_MEMBER(configFile, io.partRecycleExpTime, "partRecycleExpTime", 0.0);
    ASSERT(io.partRecycleExpTime >= 0);
    PARSE_CLASS_MEMBER(configFile, io.objectExpTime, "objectExpTime", 0.0);
    ASSERT(io.objectExpTime >= 0);
    PARSE_CLASS_MEMBER(configFile, io.cylinderExpTime, "cylinderExpTime", 0.0);
    ASSERT(io.cylinderExpTime >= 0);
    // single objects
    unsigned int totSingleObjects = configFile.vector_variable_size("singleObjects");
    if (commandLine.vector_variable_size("-singleObjects") > 0)
        totSingleObjects = commandLine.vector_variable_size("-singleObjects");
    cout << "Single objects (" << totSingleObjects << "): ";
    for (int index = 0; index < totSingleObjects; index++) {
        int singleObjectHere = 0;
        PARSE_CLASS_MEMBER_VEC(configFile, singleObjectHere, "singleObjects", index, 0);
        ASSERT(singleObjectHere >= 0);
        cout << singleObjectHere << " ";
        io.singleObjects.push_back(singleObjectHere);
    }
    cout << endl;
    ASSERT(io.singleObjects.size() == totSingleObjects);


    // flow level sensors
    unsigned int totFlowLevelBegin = configFile.vector_variable_size("flowLevelSensorBegin");
    unsigned int totFlowLevelEnd = configFile.vector_variable_size("flowLevelSensorEnd");
    if (commandLine.vector_variable_size("-flowLevelSensorBegin") > 0)
        totFlowLevelBegin = commandLine.vector_variable_size("-flowLevelSensorBegin");
    if (commandLine.vector_variable_size("-flowLevelSensorEnd") > 0)
        totFlowLevelEnd = commandLine.vector_variable_size("-flowLevelSensorEnd");
    ASSERT(totFlowLevelBegin == totFlowLevelEnd);
    const int totFlowLevel = totFlowLevelBegin;
    cout << "Number of flow level files (" << totFlowLevel << "): ";
    for (int index = 0; index < totFlowLevel; index++) {
        double flowLevelBeginHere = 0.0;
        double flowLevelEndHere = 0.0;
        PARSE_CLASS_MEMBER_VEC(configFile, flowLevelBeginHere, "flowLevelSensorBegin", index, 0.0);
        PARSE_CLASS_MEMBER_VEC(configFile, flowLevelEndHere, "flowLevelSensorEnd", index, 0.0);
        cout << "[" << flowLevelBeginHere << ", " << flowLevelEndHere << "] ";
        ASSERT(flowLevelBeginHere >= 0.0);
        ASSERT(flowLevelEndHere >= 0.0);
        ASSERT(flowLevelEndHere > flowLevelBeginHere);
        io.flowLevelBegin.push_back(flowLevelBeginHere);
        io.flowLevelEnd.push_back(flowLevelEndHere);
    }
    cout << endl;
    ASSERT(io.flowLevelBegin.size() == totFlowLevel);
    ASSERT(io.flowLevelEnd.size() == totFlowLevel);


    // object groups
    unsigned int totObjectGroupBegin = configFile.vector_variable_size("objectGroupBegin");
    unsigned int totObjectGroupEnd = configFile.vector_variable_size("objectGroupEnd");
    if (commandLine.vector_variable_size("-objectGroupBegin") > 0)
        totObjectGroupBegin = commandLine.vector_variable_size("-objectGroupBegin");
    if (commandLine.vector_variable_size("-objectGroupEnd") > 0)
        totObjectGroupEnd = commandLine.vector_variable_size("-objectGroupEnd");
    ASSERT(totObjectGroupBegin == totObjectGroupEnd);
    const unsigned int totObjectGroups = totObjectGroupBegin;
    if (totObjectGroups) {
        cout << "Object groups (" << totObjectGroups << "): ";
        for (int index = 0; index < totObjectGroups; index++) {
            unsigned int objectGroupBeginHere = 0.0;
            unsigned int objectGroupEndHere = 0.0;
            PARSE_CLASS_MEMBER_VEC(configFile, objectGroupBeginHere, "objectGroupBegin", index, 0);
            PARSE_CLASS_MEMBER_VEC(configFile, objectGroupEndHere, "objectGroupEnd", index, 0);
            cout << "[" << objectGroupBeginHere << ", " << objectGroupEndHere << "] ";
            ASSERT(objectGroupBeginHere >= 0);
            ASSERT(objectGroupEndHere >= 0);
            ASSERT(objectGroupEndHere > objectGroupBeginHere);
            io.objectGroupBegin.push_back(objectGroupBeginHere);
            io.objectGroupEnd.push_back(objectGroupEndHere);
        }
        cout << endl;
        ASSERT(io.objectGroupBegin.size() == totObjectGroups);
        ASSERT(io.objectGroupEnd.size() == totObjectGroups);
    }

    // solving particles?
    PARSE_CLASS_MEMBER(configFile, io.demSolver, "demSolver", 0);
    // solving fluid?
    PARSE_CLASS_MEMBER(configFile, io.lbmSolver, "lbSolver", 0);
    // solving rotational reference system (= adding extra apparent accelerations?)
    PARSE_CLASS_MEMBER(configFile, io.coriolisSolver, "coriolisSolver", 0);
    PARSE_CLASS_MEMBER(configFile, io.centrifugalSolver, "centrifugalSolver", 0);
    // are we using tstatic friction? (shutting down improves performances)
    PARSE_CLASS_MEMBER(configFile, dem.staticFrictionSolve, "staticFrictionSolver", 1);
    // solving free surface?
    PARSE_CLASS_MEMBER(configFile, lb.freeSurface, "freeSurfaceSolver", 0);
    // is there a force field? (shutting down improves performances)
    PARSE_CLASS_MEMBER(configFile, lb.forceField, "forceFieldSolver", 0);

    // getting LBM parameters
    lb.latticeBoltzmannGet(configFile, commandLine);
    // getting DEM parameters and particle initial status
    dem.discreteElementGet(configFile, commandLine);

    // problem specific parameters
    switch (problemName) {
        case DRUM:
        {
            // drum rotational speed
            PARSE_CLASS_MEMBER(configFile, dem.drumSpeed, "drumSpeed", 0.0);
            cout << "DRUM SPEED =" << dem.drumSpeed << "\n";
            PARSE_CLASS_MEMBER(configFile, lb.fluidMass, "fluidMass", 0.0);
            cout << "INPUT FLUID MASS=" << lb.fluidMass << "\n";
            break;
        }
        case SHEARCELL:
        {
            PARSE_CLASS_MEMBER(configFile, lb.maxVisc, "maxVisc", 0.0);
            PARSE_CLASS_MEMBER(configFile, lb.maxPlasticVisc, "maxPlasticVisc", 0.0);
            PARSE_CLASS_MEMBER(configFile, lb.maxYieldStress, "maxYieldStress", 0.0);
            PARSE_CLASS_MEMBER(configFile, lb.viscSteps, "viscSteps", 0.0);
            PARSE_CLASS_MEMBER(configFile, lb.maxShearRate, "maxShearRate", 0.0);
            PARSE_CLASS_MEMBER(configFile, lb.shearRateSteps, "shearRateSteps", 0.0);
            break;
        }
        case NET:
        case BARRIER:
        {
            PARSE_CLASS_MEMBER(configFile, lb.avalanchePosit, "avalanchePosit", 0.0);
            break;
        }
        case HONGKONG:
        {
            PARSE_CLASS_MEMBER(configFile, dem.hongkongSmoothWall, "hongkongSmoothWall", 0);
            PARSE_CLASS_MEMBER(configFile, dem.hongkongSlitSize, "hongkongSlitSize", 0.0);
            break;
        }
        case HK_LARGE:
        {
            PARSE_CLASS_MEMBER(configFile, lb.largeFlumeFlowLevel, "largeFlumeFlowLevel", 0.0);
            PARSE_CLASS_MEMBER(configFile, dem.depositArea, "depositArea", false);
            break;
        }
        case HOURGLASS:
        {
            PARSE_CLASS_MEMBER(configFile, dem.hourglassOutletSize, "hourglassOutletSize", 0.0);
            PARSE_CLASS_MEMBER(configFile, dem.hourglassOutletHeight, "hourglassOutletHeight", 0.0);
            lb.hourglassOutletHeight = dem.hourglassOutletHeight;
            break;
        }
        case HEAP:
        {
            PARSE_CLASS_MEMBER(configFile, dem.heapBaseLevel, "heapBaseLevel", 0.0);
            lb.heapBaseLevel = dem.heapBaseLevel;
            break;
        }
        case TRIAXIAL:
        {
            PARSE_CLASS_MEMBER(configFile, dem.triIsopressure, "triIsopressure", 0.0);
            PARSE_CLASS_MEMBER(configFile, dem.triBeginX, "triBeginX", 0.0);
            PARSE_CLASS_MEMBER(configFile, dem.triBeginY, "triBeginY", 0.0);
            PARSE_CLASS_MEMBER(configFile, dem.triBeginZ, "triBeginZ", 0.0);
            PARSE_CLASS_MEMBER(configFile, dem.triDefSpeed, "triDefSpeed", 0.0);
            break;
        }

    }

}

void printUfo(GetPot& command_line, GetPot& configFile) {
    // warn about unused parameters
    std::vector<std::string> ufos = configFile.unidentified_variables();
    if (ufos.size() > 0) {
        cout << INFO << "WARNING: Unrecognized lbm config file parameter(s):";
        for (unsigned int i = 0; i < ufos.size(); ++i) {
            cout << " " << ufos[i];
        }
        cout << endl;
    }
    ufos = command_line.unidentified_arguments();
    if (ufos.size() > 0) {
        cout << INFO << "WARNING: Unrecognized command line argument(s):";
        for (unsigned int i = 0; i < ufos.size(); ++i) {
            cout << " " << ufos[i];
        }
        cout << endl;
    }
}

int main(int argc, char** argv) {

    // Checking number of processes involved
    cout << "Program starts with threads:";
#pragma omp parallel
    {
#pragma omp critical
        cout << " X ";
    }
    cout << "\n";

    // DECLARATION OF VARIABLES - Input-Output ///////////////
    IO io;

    // DECLARATION OF VARIABLES - DEM ///////////////
    DEM dem;

    // DECLARATION OF VARIABLES - LB ///////////////
    LB lb;

    // print some info for restorability
    time_t t = time(0); // get time now
    io.now = localtime(&t);
    cout << "Binary was compiled on " << __DATE__ << " at " << __TIME__ << endl;
    cout << "Invocation on " << (io.now->tm_year + 1900) << "-" << (io.now->tm_mon + 1) << "-" << io.now->tm_mday <<
            ", at " << io.now->tm_hour << ":" << io.now->tm_min << ":" << io.now->tm_sec << endl;
    cout << "Invocation line: " << argv[0];
    for (int i = 1; i < argc; ++i) {
        cout << " " << argv[i];
    }
    cout << endl;

    // parsing command line and configuration file
    GetPot commandLine(argc, argv);
    commandLine.print();
    parseCommandLine(io, commandLine);

    // parsing LBM input file
    GetPot configFile(io.configFileName);
    parseConfigFile(io, dem, lb, configFile, commandLine);

    printUfo(commandLine, configFile);

    //    // bind signal handlers
    //    signal(SIGHUP, SIG_IGN); // for easy running through ssh
    //    signal(SIGQUIT, catchSignal);
    //    signal(SIGTERM, catchSignal); // important for job killing on a cluster

    io.currentTimeStep = 0;

    /* /////////////////////////////////////////////////////////////////////////
    // PROGRAM CORE  ///////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////*/
    cout << "PROBLEM ID: " << problemName << "\n";

    //  GLOBAL INITIALIZATION ///////////////////////////////////////////////
    // initializing input-output
    io.initialize();

    // initializing lattice
    lb.latticeDefinition();
    lb.LBShow();

    // initializing DEM parameters
    const tVect externalForce = lb.lbF * lb.unit.Accel;
    const tVect externalRotation = lb.rotationSpeed * lb.unit.AngVel;
    const tVect externalRotationCenter = lb.rotationCenter * lb.unit.Length;

    dem.discreteElementInit(lb.boundary, lb.lbPhysicalSize, lb.lbBoundaryLocation, externalForce,
            externalRotation, externalRotationCenter, io.coriolisSolver, io.centrifugalSolver, lb.unit.Time);

    if (io.lbmSolver) {
        lb.latticeBolzmannInit(dem.cylinders, dem.walls, dem.particles, dem.objects, io.coriolisSolver, io.centrifugalSolver);
    }

    // setting time
    lb.time = 0;
    dem.demTime = 0.0;
    io.realTime = 0.0;
    dem.demTimeStep = 0;

    // initial output
    io.outputStep(lb, dem);

    // CYCLE /////////////////////////////
    // integrate in time
    while (true) {

        if (io.realTime != 0.0 && io.realTime > io.maxTime) {
            exit_code = SUCCESS;
        }// exit normally if the maximum simulation time has been reached
        else if (io.maximumTimeSteps && io.currentTimeStep >= io.maximumTimeSteps) {
            exit_code = SUCCESS;
        } else if (io.energyExit) {
            exit_code = SUCCESS;
        } else {
            // core of the code, performs time steps
            goCycle(io, dem, lb);

            //            // exit abnormally if a serious problem has occurred
            //            if (io.problem) {
            //                exit_code = ERROR;
            //            }
        }

        if (exit_code > UNFINISHED) {
            break;
        }
    }
    return exit_code;
}

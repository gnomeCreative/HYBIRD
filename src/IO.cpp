#include "IO.h"

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// PUBLIC FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/


void IO::initialize() {

    partDirectory = workDirectory + "/particleData";
    fluidDirectory = workDirectory + "/fluidData";

    bool a;
    //    if (demSolver) {
    a = filesystem::create_directories(partDirectory.c_str());
    cout << "Work directory created = " << partDirectory << ". Result: " << a << "\n";
    //    }
    //    if (lbmSolver) {
    a = filesystem::create_directories(fluidDirectory.c_str());
    cout << "Work directory created = " << fluidDirectory << ". Result: " << a << "\n";
    //    }

    if (singleObjects.size() > 0) {
        // create directory
        singleObjectDirectory = workDirectory + "/singleObjectData";
        int a = filesystem::create_directories(singleObjectDirectory.c_str());
        cout << "Work directory created = " << singleObjectDirectory << ". Result: " << a << "\n";
    }

    // tracking specific elements
    if (singleElements.size() > 0) {
        // create directory
        singleElementDirectory = workDirectory + "/singleElementData";
        int a = filesystem::create_directories(singleElementDirectory.c_str());
        cout << "Work directory created = " << singleElementDirectory << ". Result: " << a << "\n";
    }

    //  initialising energy file
    energyFileName = workDirectory + "/energy.dat";
    //  initializing object force files
    obstacleFileName = workDirectory + "/objectForce.dat";

    switch (problemName) {
        case HONGKONG:
        {
            hongKongFlowFileName = workDirectory + "/hongKongFlow.dat";
            hongKongForceFileName = workDirectory + "/hongKongForce.dat";
            break;
        }
        case HEAP:
        {
            heapFileName = workDirectory + "/heapHeight.dat";
            break;
        }
        case INCLINEFLOW:
        {
            inclineFlowFileName = workDirectory + "/inclineFlow.dat";
            break;
        }
        case TRIAXIAL:
        {
            triaxialFileName = workDirectory + "/triaxial.dat";
            break;
        }
        case HK_SMALL:
        case HK_LARGE:
        case KELVIN:
        {
            frontFileName = workDirectory + "/flowFront.dat";
            break;
        }
        case MANGENEY:
        {
            frontFileName = workDirectory + "/flowFront.dat";
            topFileName = workDirectory + "/flowTop.dat";
            break;
        }
    }


    // build data file names 
    fluidFileFormat = fluidDirectory + "/fluid%010u.vti";
    fluid2DFileFormat = fluidDirectory + "/fluid2D%010u.txt";
    fluidLagrangianFileFormat = fluidDirectory + "/fluidLag%010u.vtu";
    partFileFormat = partDirectory + "/part%010u.vtu";
    fluidRecycleFileFormat = fluidDirectory + "/fluidRecycle%010u.dat";
    partRecycleFileFormat = partDirectory + "/partRecycle%010u.dat";
    objectFileFormat = partDirectory + "/object%010u.vtu";
    cylinderFileFormat = partDirectory + "/cylinder%010u.vtu";

    //  initializing output file
    exportFileName = workDirectory + "/export.dat";


    if (lbmSolver) {

        //  initializing fluid flow rate file
        fluidFlowRateFileName = workDirectory + "/fluidFlowRate.dat";
        fluidFlowRateFile.open(fluidFlowRateFileName.c_str(), ios::app);
        fluidFlowRateFile << "time fluidFlowRateX, fluidFlowRateY, fluidFlowRateZ\n";
        fluidFlowRateFile.close();
        
        //  initializing free surface extent file
        freeSurfaceExtentFileName = workDirectory + "/freeSurfaceExtent.dat";
        freeSurfaceExtentFile.open(freeSurfaceExtentFileName.c_str(), ios::app);
        freeSurfaceExtentFile << "time maxX y(maxX) z(maxX) minX y(minX) z(minX) maxY z(maxY) x(maxY) minY z(minY) x(minY) maxZ x(maxZ) y(maxZ) minZ x(minZ) y(minZ)\n";
        freeSurfaceExtentFile.close();

        //  initializing max speed file
        maxFluidSpeedFileName = workDirectory + "/maxFluidVel.dat";
        maxFluidSpeedFile.open(maxFluidSpeedFileName.c_str(), ios::app);
        maxFluidSpeedFile << "time maxFluidSpeed\n";
        maxFluidSpeedFile.close();

        //  initializing plasticity file
        plasticityFileName = workDirectory + "/plasticity.dat";
        plasticityFile.open(plasticityFileName.c_str(), ios::app);
        plasticityFile << "time percPlastic\n";
        plasticityFile.close();

        //  initializing mass file
        fluidMassFileName = workDirectory + "/fluidMass.dat";
        fluidMassFile.open(fluidMassFileName.c_str(), ios::app);
        fluidMassFile.close();

        initializedPlanarFile = false;

    }

    if (demSolver) {

        //  initializing particle flow rate file
        particleFlowRateFileName = workDirectory + "/particleFlowRate.dat";
        particleFlowRateFile.open(particleFlowRateFileName.c_str(), ios::app);
        particleFlowRateFile << "time particleFlowRateX, particleFlowRateY, particleFlowRateZ\n";
        particleFlowRateFile.close();

        //  initializing force file
        forceFileName = workDirectory + "/force.dat";
        forceFile.open(forceFileName.c_str(), ios::app);
        forceFile << "time collisionForces hydraulicForces\n";
        forceFile.close();

        //  initializing max speed file
        maxParticleSpeedFileName = workDirectory + "/maxParticleVel.dat";
        maxParticleSpeedFile.open(maxParticleSpeedFileName.c_str(), ios::app);
        maxParticleSpeedFile << "time maxParticleSpeed maxParticleSpin\n";
        maxParticleSpeedFile.close();


        // initializing particle center of mass file 
        particleCenterOfMassFileName = workDirectory + "/particleCenterOfMass.dat";
        particleCenterOfMassFile.open(particleCenterOfMassFileName.c_str(), ios::app);
        particleCenterOfMassFile << "time particleCenterOfMass_X particleCenterOfMass_Y particleCenterOfMass_Z\n";
        particleCenterOfMassFile.close();

        // initializing particle center of mass file 
        particleCoordinationFileName = workDirectory + "/particleCoordination.dat";

        // initializing fluid center of mass file 
        fluidCenterOfMassFileName = workDirectory + "/fluidCenterOfMass.dat";
        fluidCenterOfMassFile.open(fluidCenterOfMassFileName.c_str(), ios::app);
        fluidCenterOfMassFile << "time fluidCenterOfMass_X fluidCenterOfMass_Y fluidCenterOfMass_Z\n";
        fluidCenterOfMassFile.close();

        // initializing max overlap files
        // instantaneous
        overlapFileName = workDirectory + "/maxOverlap.dat";
        overlapFile.open(overlapFileName.c_str(), ios::app);
        overlapFile << "time meanMaxOverlap maxOverlap meanRelativeMaxOverlap maxRelativeOverlap\n";
        overlapFile.close();
        // between output steps
        dtOverlapFileName = workDirectory + "/maxDtOverlap.dat";
        dtOverlapFile.open(dtOverlapFileName.c_str(), ios::app);
        dtOverlapFile << "time meanMaxDtOverlap maxDtOverlap meanRelativeMaxDtOverlap maxRelativeDtOverlap\n";
        dtOverlapFile.close();
    }

    // initializing maximum force-on-the-wall file
    maxWallForceFileName = workDirectory + "/maxWallForce.dat";
    maxWallForceFile.open(maxWallForceFileName.c_str(), ios::app);
    maxWallForceFile << "time n_walls x (maxFParticle_X maxFParticle_Y maxFParticle_Z maxFHydro_X maxFHydro_Y maxFHydro_Z)\n";
    maxWallForceFile.close();

    // initializing force-on-the-wall file
    wallForceFileName = workDirectory + "/wallForce.dat";
    wallForceFile.open(wallForceFileName.c_str(), ios::app);
    wallForceFile << "time n_walls x (FParticle_X FParticle_Y FParticle_Z FHydro_X FHydro_Y FHydro_Z)\n";
    wallForceFile.close();
    
    // initialising forces on cylinders file
    cylinderForceFileName = workDirectory + "/cylinderForce.dat";
    cylinderForceFile.open(cylinderForceFileName.c_str(), ios::app);
    cylinderForceFile << "time n_cylinders x (FParticle_X FParticle_Y FParticle_Z FHydro_X FHydro_Y FHydro_Z)\n";
    cylinderForceFile.close();

    lastScreenExp = 0;
    lastFluidExp = 0;
    lastFluidLagrangianExp = 0;
    lastFluid2DExp = 0;
    lastPartExp = 0;
    lastObjectExp = 0;
    lastCylinderExp =0;
    lastFluidRecycleExp = 0;
    lastPartRecycleExp = 0;


}

void IO::outputStep(LB& lb, DEM& dem) {


    // PLOTTING PHASE  ////////////////////////////////////////////////////////////////
    const unsigned int screenExpCounter = (screenExpTime > 0 ? static_cast<unsigned int> (realTime / screenExpTime) + 1 : 0);
    if (screenExpCounter > lastScreenExp) {

        lastScreenExp = screenExpCounter;

        exportFile.open(exportFileName.c_str(), ios::app);

        // current iteration and real time
        cout << currentTimeStep << " time=" << std::scientific << std::setprecision(2) << realTime << " ";
        exportFile << currentTimeStep << " time=" << std::scientific << std::setprecision(2) << realTime << " ";

        if (dem.elmts.size()) {
            cout << "np=" << dem.actvPartNumber << " ";
            cout << "ns=" << dem.totSprings << " ";
            cout << "nc=" << dem.neighborTable.size() << "/" << dem.nearWallTable.size() << "/" << dem.nearObjectTable.size() << "/" << dem.nearCylinderTable.size() << " ";
        }
        if (lbmSolver) {


            const double deltaLB = std::chrono::duration<double, std::micro>(lb.endLBStep - lb.startLBStep).count();
            cout << "t=" << deltaLB << " ";
            if (lb.freeSurface) {
                const double deltaFreeSurface = std::chrono::duration<double, std::micro>(lb.endFreeSurfaceStep - lb.startFreeSurfaceStep).count();
                const double deltaUpdateMass = std::chrono::duration<double, std::micro>(lb.endUpdateMassStep - lb.startUpdateMassStep).count();
                const double deltaUpdateInterface = std::chrono::duration<double, std::micro>(lb.endUpdateInterfaceStep - lb.startUpdateInterfaceStep).count();
                const double deltaFindMutants = std::chrono::duration<double, std::micro>(lb.endFindMutantsStep - lb.startFindMutantsStep).count();
                const double deltaSmoothenInterface_1 = std::chrono::duration<double, std::micro>(lb.endSmoothenInterfaceStep_1 - lb.startSmoothenInterfaceStep_1).count();
                const double deltaSmoothenInterface_2 = std::chrono::duration<double, std::micro>(lb.endSmoothenInterfaceStep_2 - lb.startSmoothenInterfaceStep_2).count();
                const double deltaUpdateMutants = std::chrono::duration<double, std::micro>(lb.endUpdateMutantsStep - lb.startUpdateMutantsStep).count();
                const double deltaRemoveIsolated = std::chrono::duration<double, std::micro>(lb.endRemoveIsolatedStep - lb.startRemoveIsolatedStep).count();
                const double deltaRedistributeMass = std::chrono::duration<double, std::micro>(lb.endRedistributeMassStep - lb.startRedistributeMassStep).count();
                //                cout << "t_fs=" << deltaFreeSurface << " ";
                //                cout << "t_um=" << deltaUpdateMass << " ";
                //                cout << "t_ui=" << deltaUpdateInterface << " ";
                //                cout << "t_fm=" << deltaFindMutants << " ";
                //                cout << "t_si1=" << deltaSmoothenInterface_1 << " ";
                //                cout << "t_si2=" << deltaSmoothenInterface_2 << " ";
                //                cout << "n_fn=" << lb.filledNodes.size()<< " ";
                //                cout << "n_en=" << lb.emptiedNodes.size()<< " ";
                //                cout << "n_nin=" << lb.newInterfaceNodes.size()<< " ";
                //                cout << "t_um=" << deltaUpdateMutants << " ";
                //                cout << "t_ri=" << deltaRemoveIsolated << " ";
                //                cout << "t_rm=" << deltaRedistributeMass << " ";
                //                cout << "n_fs=" << lb.interfaceNodes.size() << " ";
            }
            if (demSolver) {
                const double deltaCoupling = std::chrono::duration<double, std::micro>(lb.endCouplingStep - lb.startCouplingStep).count();
                cout << "t_c=" << deltaCoupling << " ";
            }

            exportMaxSpeedFluid(lb);
            exportFreeSurfaceExtent(lb);
            exportFluidFlowRate(lb);
            exportFluidMass(lb);
            exportFluidCenterOfMass(lb);
            switch (lb.fluidMaterial.rheologyModel) {
                case BINGHAM:
                case FRICTIONAL:
                case VOELLMY:
                {
                    exportPlasticity(lb);
                    break;
                }
            }
            exportMeanViscosity(lb);
        }

        if (dem.elmts.size()) {
            exportParticleFlowRate(dem);
            exportMaxSpeedParticles(dem);
            exportForces(dem);
            exportParticleCenterOfMass(dem);
            exportParticleCoordination(dem);
            exportParticleOverlap(dem);
        }

        if (dem.walls.size() > 0) {
            exportWallForce(dem);
        }
        if (dem.cylinders.size() > 0){
            exportCylinderForces(dem);
        }
        //
        // update energies
        totalKineticEnergy = 0.0;
        energyExit = false;
        if (dem.elmts.size()) {
            dem.updateEnergy(totalKineticEnergy);
        }
        if (lbmSolver) {
            lb.updateEnergy(totalKineticEnergy);
        }
        exportEnergy(dem, lb);
        if (totalKineticEnergy < energyStopThreshold && lb.time > minimumIterations) {
            energyExit = true;
        }
        if (dem.objects.size()) {
            exportForceObstacle(dem.objects);
        }
        switch (problemName) {
            case SHEARCELL:
            {
                exportShearCell(lb, dem);
                break;
            }
            case OPENBARRIER:
            {
                if (dem.objects.size()) {
                    exportForceObstacleElement(dem);
                    exportMomentObstacleElement(dem);
                    exportHPartObstacle(dem);
                    exportKmPartObstacle(dem);
                }
                break;
            }
            case HONGKONG:
            {
                exportHongKongFlow(dem);
                exportHongKongBarrier(dem);
                break;
            }
            case INCLINEFLOW:
            {
                exportInclineFlow(lb);
                break;
            }
            case HEAP:
            {
                exportHeapHeight(lb);
                break;
            }
            case TRIAXIAL:
            {
                exportTriaxial(dem);
                break;
            }
            case MANGENEY:
            {
                exportFront(lb);
                exportTop(lb);
                break;
            }
            case HK_LARGE:
            case HK_SMALL:
            case KELVIN:
            {
                exportFront(lb);
                break;
            }


        }

        if (lbmSolver && fluid2DExpTime > 0) {

            if (!initializedPlanarFile) {
                initialize2DFile(lb);
                initializedPlanarFile = true;
            }

            update2DFile(lb);

        }

        if (singleObjects.size() > 0) {
            exportSingleObjects(dem.objects);
        }

        if (singleElements.size() > 0) {
            exportSingleElements(dem.elmts);
        }

        if (objectGroupBegin.size() > 0) {
            exportGroupForce(dem.objects);
        }

        if (flowLevelBegin.size() > 0) {
            exportFlowLevel(lb);
        }


        // closing file
        cout << endl;
        exportFile << endl;
        exportFile.close();
        cout.flush();


    }

    // FILE CREATION PHASE  ////////////////////////////////////////////////////////////////
    createFiles(lb, dem);

    if (lbmSolver) {
        const unsigned int fluid2DExpCounter = (fluid2DExpTime > 0 ? static_cast<unsigned int> (realTime / fluid2DExpTime) + 1 : 0);
        if (fluid2DExpCounter > lastFluid2DExp) {
            lastFluid2DExp = fluid2DExpCounter;
            char filePathBuffer [1024];
            sprintf(filePathBuffer, fluid2DFileFormat.c_str(), currentTimeStep);
            create2DFile(lb, filePathBuffer);
        }
    }



}

void IO::outputFinal() {
    // drag file closure
    exportFile.close();
}

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// PRIVATE FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/

// paraview files

void IO::createFiles(const LB& lb, const DEM& dem) {
    // write vtk at regular interval defined with the input file

    if (lbmSolver) {

        const unsigned int fluidExpCounter = (fluidExpTime > 0 ? static_cast<unsigned int> (realTime / fluidExpTime) + 1 : 0);
        if (fluidExpCounter > lastFluidExp) {
            lastFluidExp = fluidExpCounter;
            char filePathBuffer [1024];
            sprintf(filePathBuffer, fluidFileFormat.c_str(), currentTimeStep);
            if (fluidFormat == ParaviewFormat::Ascii) {
                exportEulerianParaviewFluid(lb, filePathBuffer);  // Benchmark: 2994s
            } else if (fluidFormat == ParaviewFormat::BinaryLowMem) {
                // exportEulerianParaviewFluid_binary(lb, filePathBuffer);  // Benchmark: 315s
                exportEulerianParaviewFluid_binaryv2(lb, filePathBuffer);  // Benchmark: 211s
            } else {
                exportEulerianParaviewFluid_binaryv3(lb, filePathBuffer);  // Benchmark: 120s (but requires ~200mb extra memory)
            }
        }

        const unsigned int fluidLagrangianExpCounter = (fluidLagrangianExpTime > 0 ? static_cast<unsigned int> (realTime / fluidLagrangianExpTime) + 1 : 0);
        if (fluidLagrangianExpCounter > lastFluidLagrangianExp) {
            lastFluidLagrangianExp = fluidLagrangianExpCounter;
            char filePathBuffer [1024];
            sprintf(filePathBuffer, fluidLagrangianFileFormat.c_str(), currentTimeStep);
            //if (currentTimeStep>120) {
            if (fluidLagrangianFormat == ParaviewFormat::Ascii) {
                exportLagrangianParaviewFluid(lb, filePathBuffer);
            } else {
                exportLagrangianParaviewFluid_binaryv3(lb, filePathBuffer);
            }
            //}
        }

        const unsigned int fluidRecycleExpCounter = (fluidRecycleExpTime > 0 ? static_cast<unsigned int> (realTime / fluidRecycleExpTime) + 1 : 0);
        if (fluidRecycleExpCounter > lastFluidRecycleExp) {
            lastFluidRecycleExp = fluidRecycleExpCounter;
            char filePathBuffer [1024];
            sprintf(filePathBuffer, fluidRecycleFileFormat.c_str(), currentTimeStep);
            // requires the pbcShift, contained int he neighborList function
            exportRecycleFluid(lb, filePathBuffer);
        }

    }

    if (demSolver) {

        const unsigned int partExpCounter = (partExpTime > 0 ? static_cast<unsigned int> (realTime / partExpTime) + 1 : 0);
        if (partExpCounter > lastPartExp) {
            lastPartExp = partExpCounter;
            char filePathBuffer [1024];
            sprintf(filePathBuffer, partFileFormat.c_str(), currentTimeStep);
            if (partExpFormat == ParaviewFormat::Ascii) {
                exportParaviewParticles(dem.elmts, dem.particles, filePathBuffer);
            } else {
                exportParaviewParticles_binaryv3(dem.elmts, dem.particles, filePathBuffer);
            }
        }

        const unsigned int partRecycleExpCounter = (partRecycleExpTime > 0 ? static_cast<unsigned int> (realTime / partRecycleExpTime) + 1 : 0);
        if (partRecycleExpCounter > lastPartRecycleExp) {
            lastPartRecycleExp = partRecycleExpCounter;
            char filePathBuffer [1024];
            sprintf(filePathBuffer, partRecycleFileFormat.c_str(), currentTimeStep);
            // requires the pbcShift, contained int he neighborList function
            exportRecycleParticles(dem.elmts, dem.pbcs, filePathBuffer);
        }
    }
    if (dem.objects.size() > 0) {

        const unsigned int objectExpCounter = (objectExpTime > 0 ? static_cast<unsigned int> (realTime / objectExpTime) + 1 : 0);
        if (objectExpCounter > lastObjectExp) {
            lastObjectExp = objectExpCounter;
            char filePathBuffer [1024];
            sprintf(filePathBuffer, objectFileFormat.c_str(), currentTimeStep);
            exportParaviewObjects(dem.objects, filePathBuffer);
        }
    }
    
    if (dem.cylinders.size() > 0) {
        
        const unsigned int CylinderExpCounter = (cylinderExpTime > 0 ? static_cast<unsigned int> (realTime / cylinderExpTime) + 1 : 0);
        if (CylinderExpCounter > lastCylinderExp) {
            lastCylinderExp = CylinderExpCounter;
            char filePathBuffer [1024];
            sprintf(filePathBuffer, cylinderFileFormat.c_str(), currentTimeStep);
            exportParaviewCylinders(dem.cylinders, filePathBuffer);
        }
        
    }
          
}

void IO::initialize2DFile(const LB& lb) {

    // initialize containers for 2D files
    planarHeight.resize(lb.lbSize[0]);
    planarVel.resize(lb.lbSize[0]);
    planarPointVel.resize(lb.lbSize[0]);
    planarLevel.resize(lb.lbSize[0]);
    maxPlanarHeight.resize(lb.lbSize[0]);
    maxPlanarVel.resize(lb.lbSize[0]);
    maxPlanarPointVel.resize(lb.lbSize[0]);

    for (int i = 0; i < lb.lbSize[0]; ++i) {
        planarHeight[i].resize(lb.lbSize[1]);
        planarVel[i].resize(lb.lbSize[1]);
        planarPointVel[i].resize(lb.lbSize[1]);
        planarLevel[i].resize(lb.lbSize[1]);
        maxPlanarHeight[i].resize(lb.lbSize[1]);
        maxPlanarVel[i].resize(lb.lbSize[1]);
        maxPlanarPointVel[i].resize(lb.lbSize[1]);
        for (int j = 0; j < lb.lbSize[1]; ++j) {
            planarLevel[i][j] = 0.0;
            maxPlanarHeight[i][j] = 0.0;
            maxPlanarVel[i][j] = 0.0;
            maxPlanarPointVel[i][j] = 0.0;
        }
    }

    // compute topographic quantities
    for (nodeList::const_iterator it = lb.wallNodes.begin(); it != lb.wallNodes.end(); ++it) {
        const node* nodeHere = *it;
        const unsigned int index = nodeHere->coord;

        if (nodeHere->isTopography()) {
            const unsigned int iHere = lb.getX(index);
            const unsigned int jHere = lb.getY(index);
            const double kHere = lb.getPositionZ(index);
            planarLevel[iHere][jHere] = kHere;
        }
    }
}

void IO::update2DFile(const LB& lb) {

    // initialize containers for 2D files
    for (int i = 0; i < lb.lbSize[0]; ++i) {
        for (int j = 0; j < lb.lbSize[1]; ++j) {
            planarHeight[i][j] = 0.0;
            planarVel[i][j] = 0.0;
            planarPointVel[i][j] = 0.0;
        }
    }

    // compute 2D quantities
    for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
        const node* nodeHere = *it;
        const unsigned int index = nodeHere->coord;
        const unsigned int iHere = lb.getX(index);
        const unsigned int jHere = lb.getY(index);
        const double velHere = nodeHere->u.norm();
        planarHeight[iHere][jHere] += nodeHere->mass;
        planarVel[iHere][jHere] += nodeHere->mass * velHere;
        planarPointVel[iHere][jHere] = std::max(planarPointVel[iHere][jHere], velHere);
    }


    // average velocity over height
    for (int i = 0; i < lb.lbSize[0]; ++i) {
        for (int j = 0; j < lb.lbSize[1]; ++j) {
            if (planarHeight[i][j] > 0.0) {
                planarVel[i][j] = planarVel[i][j] / planarHeight[i][j];
            }
        }
    }

    // update global maxima
    for (int i = 0; i < lb.lbSize[0]; ++i) {
        for (int j = 0; j < lb.lbSize[1]; ++j) {
            if (maxPlanarHeight[i][j] < planarHeight[i][j]) {
                maxPlanarHeight[i][j] = planarHeight[i][j];
            }
            if (maxPlanarVel[i][j] < planarVel[i][j]) {
                maxPlanarVel[i][j] = planarVel[i][j];
            }
            if (maxPlanarPointVel[i][j] < planarPointVel[i][j]) {
                maxPlanarPointVel[i][j] = planarPointVel[i][j];
            }
        }
    }
}

void IO::create2DFile(const LB& lb, const string& planarFile) {

    // constants for scaling (sorry, hard coded)
    double xScaling = -1.0 * lb.translateTopographyX;
    double yScaling = -1.0 * lb.translateTopographyY;
    double zScaling = -1.0 * lb.translateTopographyZ;
    switch (problemName) {
        case STAVA:
        {
            xScaling = 691859.0;
            yScaling = 5128805.0;
            zScaling = 860.0;
            break;
        }
    }

    ofstream fluid2DFile;
    fluid2DFile.open(planarFile.c_str());

    fluid2DFile << "X, Y, Z, H, U, Hmax, Umax, ULocalMax, ULocalMaxMax" << endl;
    for (int i = 0; i < lb.lbSize[0]; ++i) {
        const double xHere = double(i) * lb.unit.Length + xScaling;
        for (int j = 0; j < lb.lbSize[1]; ++j) {
            const double maxHeightHere = maxPlanarHeight[i][j] * lb.unit.Length;
            //if (maxHeightHere>0.0) {
            const double yHere = double(j) * lb.unit.Length + yScaling;
            const double levelHere = planarLevel[i][j] * lb.unit.Length + zScaling;
            const double heightHere = planarHeight[i][j] * lb.unit.Length;
            const double velHere = planarVel[i][j] * lb.unit.Speed;
            const double maxVelHere = maxPlanarVel[i][j] * lb.unit.Speed;
            const double pointVelHere = planarPointVel[i][j] * lb.unit.Speed;
            const double maxPointVelHere = maxPlanarPointVel[i][j] * lb.unit.Speed;

            fluid2DFile << std::setprecision(3) << std::fixed << xHere << ", "
                    << std::setprecision(3) << std::fixed << yHere << ", "
                    << std::setprecision(3) << std::fixed << levelHere << ", "
                    << std::setprecision(3) << std::fixed << heightHere << ", "
                    << std::setprecision(3) << std::fixed << velHere << ", "
                    << std::setprecision(3) << std::fixed << maxHeightHere << ", "
                    << std::setprecision(3) << std::fixed << maxVelHere << ", "
                    << std::setprecision(3) << std::fixed << pointVelHere << ", "
                    << std::setprecision(3) << std::fixed << maxPointVelHere << endl;
            //}
        }
    }
    fluid2DFile.close();
}

void IO::exportSingleObjects(const objectList& objects) {

    std::vector<string> xyzString;
    xyzString.resize(3);
    xyzString[0] = 'x';
    xyzString[1] = 'y';
    xyzString[2] = 'z';
    for (int i = 0; i < singleObjects.size(); i++) {
        const unsigned int indexHere = singleObjects[i];
        const tVect forceHere = objects[indexHere].FHydro + objects[indexHere].FParticle;
        double f[3] = {forceHere.dot(Xp), forceHere.dot(Yp), forceHere.dot(Zp)};
        for (int xyz = 0; xyz < 3; xyz++) {
            string singleObjectFileName = singleObjectDirectory + "/singleObject_" + std::to_string(indexHere) + "_" + xyzString[xyz] + ".dat";
            ofstream singleObjectFile;
            singleObjectFile.open(singleObjectFileName.c_str(), ios::app);
            singleObjectFile << realTime << " " << f[xyz] << "\n";
            singleObjectFile.close();
        }
    }


}

void IO::exportSingleElements(const elmtList& elmts) {
    // export specific element from the particle.dat file

    for (int i = 0; i < singleElements.size(); i++) {
        // cycle through all the elemnts to be tracked
        const unsigned int indexHere = singleElements[i];
        // get the position of the tracked element
        const double xNow = elmts[indexHere].x0.x;
        const double yNow = elmts[indexHere].x0.y;
        const double zNow = elmts[indexHere].x0.z;
        // get the forces (particle and hydro) of the tracked element
        const tVect forceHydro = elmts[indexHere].FHydro;
        const tVect forcePart = elmts[indexHere].FParticle;
        double FParticleX = forcePart.dot(Xp);
        double FParticleY = forcePart.dot(Yp);
        double FParticleZ = forcePart.dot(Zp);
        double FHydroX = forceHydro.dot(Xp);
        double FHydroY = forceHydro.dot(Yp);
        double FHydroZ = forceHydro.dot(Zp);

        // write the results in a .dat file
        string singleObjectFileName = singleElementDirectory + "/singleElement_" + std::to_string(indexHere) + ".dat";
        ofstream singleElementFile;
        singleElementFile.open(singleObjectFileName.c_str(), ios::app);
        singleElementFile << std::scientific << std::setprecision(6); // Set scientific format and precision

        if (realTime<=1e-10) {
            singleElementFile <<
                "realTime" << " " << "x" << " " << "y" << " " << "z" << " " <<
                "FParticleX" << " " << "FParticleY" << " " <<  "FParticleZ" << " " <<
                "FHydroX" << " " << "FHydroY" << " " << "FHydroZ" << "\n" <<
                realTime << " " <<
                xNow << " " << yNow << " " <<  zNow << " " <<
                FParticleX << " " << FParticleY << " " <<  FParticleZ << " " <<
                FHydroX << " " << FHydroY << " " << FHydroZ << "\n";
        } else {
            singleElementFile <<
                realTime << " " <<
                xNow << " " << yNow << " " <<  zNow << " " <<
                FParticleX << " " << FParticleY << " " <<  FParticleZ << " " <<
                FHydroX << " " << FHydroY << " " << FHydroZ << "\n";
            singleElementFile.close();
        }
    }

}

void IO::exportFlowLevel(const LB& lb) {

    doubleList volumeCount;
    volumeCount.resize(flowLevelBegin.size());
    for (int i = 0; i < flowLevelBegin.size(); i++) {
        volumeCount[i] = 0.0;
    }

    if (lbmSolver) {
        for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
            const node* nodeHere = *it;
            const double xCoordHere = lb.getPositionX(nodeHere->coord) * lb.unit.Length;
            for (int i = 0; i < flowLevelBegin.size(); i++) {
                const double beginHere = flowLevelBegin[i];
                const double endHere = flowLevelEnd[i];
                if (xCoordHere <= endHere && xCoordHere >= beginHere) {
                    //cout<<lb.getPositionX(nodeHere->coord)*lb.unit.Length<<endl;
                    volumeCount[i] += nodeHere->mass;
                }
            }
        }
    }

    for (int i = 0; i < flowLevelBegin.size(); i++) {
        // compute window length in lattice units (to round to exact window measurement in )
        const unsigned int latticeBegin = ceil(flowLevelBegin[i] / lb.unit.Length);
        const unsigned int latticeEnd = floor(flowLevelEnd[i] / lb.unit.Length);
        const unsigned int latticewindowSpan = latticeEnd - latticeBegin + 1;
        ASSERT(latticewindowSpan >= 1);

        const double windowSpan = double(latticewindowSpan) * lb.unit.Length;

        const double windowDepth = double(lb.lbSize[2] - 2) * lb.unit.Length;

        const double flowVolume = volumeCount[i] * lb.unit.Volume;

        const double flowLevelHere = flowVolume / (windowDepth * windowSpan);

        string flowLevelFileName = workDirectory + "/flowLevel_" + std::to_string(i) + ".dat";
        ofstream flowLevelFile;
        flowLevelFile.open(flowLevelFileName.c_str(), ios::app);
        flowLevelFile << realTime << " " << flowLevelHere << "\n";
        flowLevelFile.close();
    }
}

void IO::exportGroupForce(const objectList& objects) {
    // printing total force on object group
    for (int i = 0; i < objectGroupBegin.size(); i++) {
        // compute window length in lattice units (to round to exact window measurement in )
        const unsigned int beginHere = objectGroupBegin[i];
        const unsigned int endHere = objectGroupEnd[i];

        tVect groupObstacleForce = totForceObject(objects, beginHere, endHere);

        string objectGroupForceFileName = workDirectory + "/objectGroupForce_" + std::to_string(i) + ".dat";
        ofstream objectGroupForceFile;
        objectGroupForceFile.open(objectGroupForceFileName.c_str(), ios::app);
        objectGroupForceFile << realTime << " ";
        groupObstacleForce.printLine(objectGroupForceFile);
        objectGroupForceFile.close();
    }

}

void IO::exportParaviewParticles(const elmtList& elmts, const particleList& particles, const string& particleFile) {

    const int one = 1;

    unsigned int activeNumber = 0;
    for (int i = 0; i < particles.size(); ++i) {
        if (particles[i].active) {
            ++activeNumber;
        }
    }

    const int Pnumber = particles.size();
    const tVect nx(1.0, 0.0, 0.0), ny(0.0, 1.0, 0.0), nz(0.0, 0.0, 1.0);
    tVect n1, n2, n3;

    std::cout.precision(10);
    std::cout.fixed;


    // file opening
    ofstream paraviewParticleFile;
    paraviewParticleFile.open(particleFile.c_str());
    // writing on header file
    paraviewParticleFile << "<?xml version=\"1.0\"?>\n";
    paraviewParticleFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
    paraviewParticleFile << " <UnstructuredGrid GhostLevel=\"0\">\n";
    paraviewParticleFile << "  <Piece NumberOfPoints=\"" << activeNumber << "\" NumberOfCells=\"" << activeNumber << "\">\n";
    paraviewParticleFile << "   <PointData>\n";
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"radius\"/>\n";
    for (int i = 0; i < Pnumber; ++i) {
        if (particles[i].active) {
            paraviewParticleFile << particles[i].r << "\n";
        }
    }
    //    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"mass\"/>\n";
    //    for (int i = 0; i < Pnumber; ++i) {
    //        if (particles[i].active) {
    //            paraviewParticleFile << elmts[particles[i].clusterIndex].m << "\n";
    //        }
    //    }
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"particleIndex\"/>\n";
    for (int i = 0; i < Pnumber; ++i) {
        if (particles[i].active) {
            paraviewParticleFile << particles[i].particleIndex << "\n";
        }
    }
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"clusterIndex\"/>\n";
    for (int i = 0; i < Pnumber; ++i) {
        if (particles[i].active) {
            paraviewParticleFile << particles[i].clusterIndex << "\n";
        }
    }
    //    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"tableCell\"/>\n";
    //    for (int i = 0; i < Pnumber; ++i) {
    //        if (particles[i].active) {
    //            paraviewParticleFile << particles[i].tableCell << "\n";
    //        }
    //    }
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"coordinationNumber\"/>\n";
    for (int i = 0; i < Pnumber; ++i) {
        if (particles[i].active) {
            paraviewParticleFile << elmts[particles[i].clusterIndex].coordination << "\n";
        }
    }
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\"/>\n";
    for (int i = 0; i < Pnumber; ++i) {
        if (particles[i].active) {
            particles[i].x1.printFixedLine(paraviewParticleFile);
        }
    }
    //    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"slipCase\" NumberOfComponents=\"1\"/>\n";
    //    for (int i = 0; i < Pnumber; ++i) {
    //        if (particles[i].active) {
    //            paraviewParticleFile << elmts[particles[i].clusterIndex].slippingCase << "\n";
    //        }
    //    }
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"w\" NumberOfComponents=\"3\"/>\n";
    for (int i = 0; i < Pnumber; ++i) {
        if (particles[i].active) {
            elmts[particles[i].clusterIndex].wGlobal.printFixedLine(paraviewParticleFile);
        }
    }
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"FParticle\" NumberOfComponents=\"3\"/>\n";
    for (int i = 0; i < Pnumber; ++i) {
        if (particles[i].active) {
            elmts[particles[i].clusterIndex].FParticle.printFixedLine(paraviewParticleFile);
        }
    }
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"FGrav\" NumberOfComponents=\"3\"/>\n";
    for (int i = 0; i < Pnumber; ++i) {
        if (particles[i].active) {
            elmts[particles[i].clusterIndex].FGrav.printFixedLine(paraviewParticleFile);
        }
    }
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"solidIntensity\" NumberOfComponents=\"3\"/>\n";
    for (int i = 0; i < Pnumber; ++i) {
        if (particles[i].active) {
            elmts[particles[i].clusterIndex].solidIntensity.printFixedLine(paraviewParticleFile);
        }
    }
    //    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"FCoriolis\" NumberOfComponents=\"3\"/>\n";
    //    for (int i = 0; i < Pnumber; ++i) {
    //        if (particles[i].active) {
    //            const tVect FCoriolis = elmts[particles[i].clusterIndex].ACoriolis * elmts[particles[i].clusterIndex].m;
    //            FCoriolis.printFixedLine(paraviewParticleFile);
    //        }
    //    }
    //    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"FCentrifugal\" NumberOfComponents=\"3\"/>\n";
    //    for (int i = 0; i < Pnumber; ++i) {
    //        if (particles[i].active) {
    //            const tVect FCentrifugal = elmts[particles[i].clusterIndex].ACentrifugal * elmts[particles[i].clusterIndex].m;
    //            FCentrifugal.printFixedLine(paraviewParticleFile);
    //        }
    //    }
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"MParticle\" NumberOfComponents=\"3\"/>\n";
    for (int i = 0; i < Pnumber; ++i) {
        if (particles[i].active) {
            elmts[particles[i].clusterIndex].MParticle.printFixedLine(paraviewParticleFile);
        }
    }
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"FWall\" NumberOfComponents=\"3\"/>\n";
    for (int i = 0; i < Pnumber; ++i) {
        if (particles[i].active) {
            elmts[particles[i].clusterIndex].FWall.printFixedLine(paraviewParticleFile);
        }
    }
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"FCylinder\" NumberOfComponents=\"3\"/>\n";
    for (int i = 0; i < Pnumber; ++i) {
        if (particles[i].active) {
            elmts[particles[i].clusterIndex].FCylinder.printFixedLine(paraviewParticleFile);
        }
    }
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"MWall\" NumberOfComponents=\"3\"/>\n";
    for (int i = 0; i < Pnumber; ++i) {
        if (particles[i].active) {
            elmts[particles[i].clusterIndex].MWall.printFixedLine(paraviewParticleFile);
        }
    }
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"MCylinder\" NumberOfComponents=\"3\"/>\n";
    for (int i = 0; i < Pnumber; ++i) {
        if (particles[i].active) {
            elmts[particles[i].clusterIndex].MCylinder.printFixedLine(paraviewParticleFile);
        }
    }
    if (lbmSolver) {
        paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"FHydro\" NumberOfComponents=\"3\"/>\n";
        for (int i = 0; i < Pnumber; ++i) {
            if (particles[i].active) {
                elmts[particles[i].clusterIndex].FHydro.printFixedLine(paraviewParticleFile);
            }
        }
        paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"MHydro\" NumberOfComponents=\"3\"/>\n";
        for (int i = 0; i < Pnumber; ++i) {
            if (particles[i].active) {
                elmts[particles[i].clusterIndex].MHydro.printLine(paraviewParticleFile);
            }
        }
        paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"fluidIntensity\" NumberOfComponents=\"3\"/>\n";
        for (int i = 0; i < Pnumber; ++i) {
            if (particles[i].active) {
                elmts[particles[i].clusterIndex].FHydro.abs().printFixedLine(paraviewParticleFile);
            }
        }
    }
    //    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"n1\" NumberOfComponents=\"3\"/>\n";
    //    for (int i = 0; i < Pnumber; ++i) {
    //        if (particles[i].active) {
    //            n1 = project(nx, elmts[particles[i].clusterIndex].q0);
    //            n1.printLine(paraviewParticleFile);
    //        }
    //    }
    //    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"n2\" NumberOfComponents=\"3\"/>\n";
    //    for (int i = 0; i < Pnumber; ++i) {
    //        if (particles[i].active) {
    //            n2 = project(ny, elmts[particles[i].clusterIndex].q0);
    //            n2.printLine(paraviewParticleFile);
    //        }
    //    }
    //    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"n3\" NumberOfComponents=\"3\"/>\n";
    //    for (int i = 0; i < Pnumber; ++i) {
    //        if (particles[i].active) {
    //            n3 = project(nz, elmts[particles[i].clusterIndex].q0);
    //            n3.printLine(paraviewParticleFile);
    //        }
    //    }
    paraviewParticleFile << "   </PointData>\n";
    paraviewParticleFile << "   <CellData>\n";
    paraviewParticleFile << "   </CellData>\n";
    paraviewParticleFile << "   <Points>\n";
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n";
    for (int i = 0; i < Pnumber; ++i) {
        if (particles[i].active) {
            particles[i].x0.printFixedLine(paraviewParticleFile);
        }
    }
    paraviewParticleFile << "   </Points>\n";
    paraviewParticleFile << "   <Cells>\n";
    paraviewParticleFile << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int i = 1; i < activeNumber + 1; ++i) {
        paraviewParticleFile << i - 1 << "\n";
    }
    paraviewParticleFile << "    </DataArray>\n";
    paraviewParticleFile << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 1; i < activeNumber + 1; ++i) {
        paraviewParticleFile << i << "\n";
    }
    paraviewParticleFile << "    </DataArray>\n";
    paraviewParticleFile << "    <DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < activeNumber; ++i) {
        paraviewParticleFile << one << "\n";
    }
    paraviewParticleFile << "    </DataArray>\n";
    paraviewParticleFile << "   </Cells>\n";
    paraviewParticleFile << "  </Piece>\n";
    //    paraviewHeaderFile<<"  <Piece Source=\""<<pieceFile<<"\"/>\n";
    paraviewParticleFile << " </UnstructuredGrid>\n";
    paraviewParticleFile << "</VTKFile>";
    // header file closing
    paraviewParticleFile.close();
}

void IO::exportParaviewParticles_binaryv3(const elmtList& elmts, const particleList& particles, const string& particleFile) {
    // Build a list of active particles for faster access
    std::vector<unsigned int> active_particles;
    for (int i = 0; i < particles.size(); ++i) {
        if (particles[i].active) {
            active_particles.push_back(i);
        }
    }

    // file opening
    ofstream paraviewParticleFile;
    paraviewParticleFile.open(particleFile.c_str(), std::ios::binary);
    // writing on header file
    static const uint16_t m_endianCheck(0x00ff);
    const bool is_big_endian ( *((const uint8_t*)&m_endianCheck) == 0x0);
    paraviewParticleFile << "<?xml version=\"1.0\"?>\n";
    paraviewParticleFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" << (is_big_endian ? "BigEndian" : "LittleEndian") << "\" header_type=\"UInt32\">\n";
    paraviewParticleFile << " <UnstructuredGrid GhostLevel=\"0\">\n";
    paraviewParticleFile << "  <Piece NumberOfPoints=\"" << active_particles.size() << "\" NumberOfCells=\"" << active_particles.size() << "\">\n";
    paraviewParticleFile << "   <PointData>\n";
    unsigned int offset = 0;
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"radius\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += active_particles.size() * sizeof(double) + sizeof(unsigned int);
    paraviewParticleFile << "    <DataArray type=\"UInt32\" Name=\"particleIndex\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += active_particles.size() * sizeof(unsigned int) + sizeof(unsigned int);
    paraviewParticleFile << "    <DataArray type=\"UInt32\" Name=\"clusterIndex\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += active_particles.size() * sizeof(unsigned int) + sizeof(unsigned int);
    paraviewParticleFile << "    <DataArray type=\"UInt32\" Name=\"coordinationNumber\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += active_particles.size() * sizeof(unsigned int) + sizeof(unsigned int);
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += active_particles.size() * 3 * sizeof(double) + sizeof(unsigned int);
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"w\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += active_particles.size() * 3 * sizeof(double) + sizeof(unsigned int);
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"FParticle\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += active_particles.size() * 3 * sizeof(double) + sizeof(unsigned int);
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"FGrav\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += active_particles.size() * 3 * sizeof(double) + sizeof(unsigned int);
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"solidIntensity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += active_particles.size() * 3 * sizeof(double) + sizeof(unsigned int);
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"MParticle\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += active_particles.size() * 3 * sizeof(double) + sizeof(unsigned int);
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"FWall\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += active_particles.size() * 3 * sizeof(double) + sizeof(unsigned int);
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"MWall\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += active_particles.size() * 3 * sizeof(double) + sizeof(unsigned int);
    if (lbmSolver) {
        paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"FHydro\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\" />\n";
        offset += active_particles.size() * 3 * sizeof(double) + sizeof(unsigned int);
        paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"MHydro\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\" />\n";
        offset += active_particles.size() * 3 * sizeof(double) + sizeof(unsigned int);
        paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"fluidIntensity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\" />\n";
        offset += active_particles.size() * 3 * sizeof(double) + sizeof(unsigned int);
    }
    paraviewParticleFile << "   </PointData>\n";
    paraviewParticleFile << "   <CellData>\n";
    paraviewParticleFile << "   </CellData>\n";
    paraviewParticleFile << "   <Points>\n";
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += active_particles.size() * 3 * sizeof(double) + sizeof(unsigned int);
    paraviewParticleFile << "   </Points>\n";
    paraviewParticleFile << "   <Cells>\n";
    paraviewParticleFile << "    <DataArray type=\"UInt32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += active_particles.size() * sizeof(unsigned int) + sizeof(unsigned int);
    paraviewParticleFile << "    <DataArray type=\"UInt32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += active_particles.size() * sizeof(unsigned int) + sizeof(unsigned int);
    paraviewParticleFile << "    <DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += active_particles.size() * sizeof(unsigned char) + sizeof(unsigned int);
    paraviewParticleFile << "   </Cells>\n";
    paraviewParticleFile << "  </Piece>\n";
    paraviewParticleFile << " </UnstructuredGrid>\n";
    paraviewParticleFile << " <AppendedData encoding=\"raw\">\n  _";
    // Allocate a buffer equal to size of the largest data array
    // Allocate once rather than allocating and freeing per export
    static char *const t_buffer = static_cast<char*>(malloc(active_particles.size() * 3 * sizeof(double)));
    static double *const d_buffer = reinterpret_cast<double*>(t_buffer);
    static tVect *const v_buffer = reinterpret_cast<tVect*>(t_buffer);
    static unsigned int *const u_buffer = reinterpret_cast<unsigned int*>(t_buffer);
    static unsigned char *const uc_buffer = reinterpret_cast<unsigned char*>(t_buffer);
    // radius
    offset = active_particles.size() * sizeof(double);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (int i = 0; i < active_particles.size(); ++i) {
        d_buffer[i] = particles[active_particles[i]].r;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // particleIndex
    offset = active_particles.size() * sizeof(unsigned int);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (int i = 0; i < active_particles.size(); ++i) {
        u_buffer[i] = particles[active_particles[i]].particleIndex;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // clusterIndex
    offset = active_particles.size() * sizeof(unsigned int);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (int i = 0; i < active_particles.size(); ++i) {
        u_buffer[i] = particles[active_particles[i]].clusterIndex;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // coordinationNumber
    offset = active_particles.size() * sizeof(unsigned int);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (int i = 0; i < active_particles.size(); ++i) {
        u_buffer[i] = elmts[particles[active_particles[i]].clusterIndex].coordination;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // v
    offset = active_particles.size() * 3 * sizeof(double);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (int i = 0; i < active_particles.size(); ++i) {
        v_buffer[i] = particles[active_particles[i]].x1;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // w
    offset = active_particles.size() * 3 * sizeof(double);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (int i = 0; i < active_particles.size(); ++i) {
        v_buffer[i] = elmts[particles[active_particles[i]].clusterIndex].wGlobal;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // FParticle
    offset = active_particles.size() * 3 * sizeof(double);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (int i = 0; i < active_particles.size(); ++i) {
        v_buffer[i] = elmts[particles[active_particles[i]].clusterIndex].FParticle;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // FGrav
    offset = active_particles.size() * 3 * sizeof(double);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (int i = 0; i < active_particles.size(); ++i) {
        v_buffer[i] = elmts[particles[active_particles[i]].clusterIndex].FGrav;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // solidIntensity
    offset = active_particles.size() * 3 * sizeof(double);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (int i = 0; i < active_particles.size(); ++i) {
        v_buffer[i] = elmts[particles[active_particles[i]].clusterIndex].solidIntensity;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // MParticle
    offset = active_particles.size() * 3 * sizeof(double);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (int i = 0; i < active_particles.size(); ++i) {
        v_buffer[i] = elmts[particles[active_particles[i]].clusterIndex].MParticle;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // FWall
    offset = active_particles.size() * 3 * sizeof(double);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (int i = 0; i < active_particles.size(); ++i) {
        v_buffer[i] = elmts[particles[active_particles[i]].clusterIndex].FWall;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // MWall
    offset = active_particles.size() * 3 * sizeof(double);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (int i = 0; i < active_particles.size(); ++i) {
        v_buffer[i] = elmts[particles[active_particles[i]].clusterIndex].MWall;
    }
    paraviewParticleFile.write(t_buffer, offset);
    if (lbmSolver) {
        // FHydro
        offset = active_particles.size() * 3 * sizeof(double);
        paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        for (int i = 0; i < active_particles.size(); ++i) {
            v_buffer[i] = elmts[particles[active_particles[i]].clusterIndex].FHydro;
        }
        paraviewParticleFile.write(t_buffer, offset);
        // MHydro
        offset = active_particles.size() * 3 * sizeof(double);
        paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        for (int i = 0; i < active_particles.size(); ++i) {
            v_buffer[i] = elmts[particles[active_particles[i]].clusterIndex].MHydro;
        }
        paraviewParticleFile.write(t_buffer, offset);
        // fluidIntensity
        offset = active_particles.size() * 3 * sizeof(double);
        paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        for (int i = 0; i < active_particles.size(); ++i) {
            v_buffer[i] = elmts[particles[active_particles[i]].clusterIndex].FHydro;
        }
        paraviewParticleFile.write(t_buffer, offset);
    }

    // Points
    offset = active_particles.size() * 3 * sizeof(double);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (int i = 0; i < active_particles.size(); ++i) {
        v_buffer[i] = particles[active_particles[i]].x0;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // connectivity
    offset = active_particles.size() * sizeof(unsigned int);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < active_particles.size(); ++i) {
        u_buffer[i] = i;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // offsets
    offset = active_particles.size() * sizeof(unsigned int);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 1; i < active_particles.size() + 1; ++i) {
        u_buffer[i] = i;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // types
    offset = active_particles.size() * sizeof(unsigned char);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    std::fill(uc_buffer, uc_buffer + active_particles.size(), 1u);
    paraviewParticleFile.write(t_buffer, offset);
    paraviewParticleFile << "</AppendedData>";
    paraviewParticleFile << "</VTKFile>";
    // header file closing
    paraviewParticleFile.close();
}

void IO::exportParaviewObjects(const objectList& objects, const string& objectFile) {

    const int one = 1;
    const int Onumber = objects.size();

    // file opening
    ofstream paraviewObjectFile;
    paraviewObjectFile.open(objectFile.c_str());
    // writing on header file
    paraviewObjectFile << "<?xml version=\"1.0\"?>\n";
    paraviewObjectFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
    paraviewObjectFile << " <UnstructuredGrid GhostLevel=\"0\">\n";
    paraviewObjectFile << "  <Piece NumberOfPoints=\"" << Onumber << "\" NumberOfCells=\"" << Onumber << "\">\n";
    paraviewObjectFile << "   <PointData>\n";
    paraviewObjectFile << "    <DataArray type=\"Float64\" Name=\"radius\"/>\n";
    for (int i = 0; i < Onumber; ++i) {
        paraviewObjectFile << objects[i].r << "\n";
    }
    paraviewObjectFile << "    <DataArray type=\"Float64\" Name=\"objectIndex\"/>\n";
    for (int i = 0; i < Onumber; ++i) {
        paraviewObjectFile << objects[i].index << "\n";
    }
    paraviewObjectFile << "    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\"/>\n";
    for (int i = 0; i < Onumber; ++i) {
        objects[i].x1.printLine(paraviewObjectFile);
    }
    paraviewObjectFile << "    <DataArray type=\"Float64\" Name=\"FParticle\" NumberOfComponents=\"3\"/>\n";
    for (int i = 0; i < Onumber; ++i) {
        objects[i].FParticle.printLine(paraviewObjectFile);
    }
    paraviewObjectFile << "    <DataArray type=\"Float64\" Name=\"FHydro\" NumberOfComponents=\"3\"/>\n";
    for (int i = 0; i < Onumber; ++i) {
        objects[i].FHydro.printLine(paraviewObjectFile);
    }
    paraviewObjectFile << "    <DataArray type=\"Float64\" Name=\"maxFParticle\" NumberOfComponents=\"3\"/>\n";
    for (int i = 0; i < Onumber; ++i) {
        objects[i].maxFParticle.printLine(paraviewObjectFile);
    }
    paraviewObjectFile << "    <DataArray type=\"Float64\" Name=\"savedFParticle\" NumberOfComponents=\"3\"/>\n";
    for (int i = 0; i < Onumber; ++i) {
        objects[i].savedFParticle.printLine(paraviewObjectFile);
    }
    paraviewObjectFile << "    <DataArray type=\"Float64\" Name=\"timeMaxFParticle\"/>\n";
    for (int i = 0; i < Onumber; ++i) {
        paraviewObjectFile << objects[i].timeMaxFParticle << "\n";
    }
    paraviewObjectFile << "   </PointData>\n";
    paraviewObjectFile << "   <CellData>\n";
    paraviewObjectFile << "   </CellData>\n";
    paraviewObjectFile << "   <Points>\n";
    paraviewObjectFile << "    <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n";
    for (int i = 0; i < Onumber; ++i) {
        objects[i].x0.printLine(paraviewObjectFile);
    }
    paraviewObjectFile << "   </Points>\n";
    paraviewObjectFile << "   <Cells>\n";
    paraviewObjectFile << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int i = 1; i < Onumber + 1; ++i) {
        paraviewObjectFile << i - 1 << "\n";
    }
    paraviewObjectFile << "    </DataArray>\n";
    paraviewObjectFile << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 1; i < Onumber + 1; ++i) {
        paraviewObjectFile << i << "\n";
    }
    paraviewObjectFile << "    </DataArray>\n";
    paraviewObjectFile << "    <DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < Onumber; ++i) {
        paraviewObjectFile << one << "\n";
    }
    paraviewObjectFile << "    </DataArray>\n";
    paraviewObjectFile << "   </Cells>\n";
    paraviewObjectFile << "  </Piece>\n";
    //    paraviewHeaderFile<<"  <Piece Source=\""<<pieceFile<<"\"/>\n";
    paraviewObjectFile << " </UnstructuredGrid>\n";
    paraviewObjectFile << "</VTKFile>";
    // header file closing
    paraviewObjectFile.close();
}

void IO::exportParaviewCylinders(const cylinderList& cylinders, const string& cylinderFile) {
    
    const int vtkLineType = 3; // vtk line type
    const int Cnumber = cylinders.size(); // number of cylinders
    
    // file opening
    ofstream paraviewCylinderFile;
    paraviewCylinderFile.open(cylinderFile.c_str());
    
    // write the header of the file
    paraviewCylinderFile << "<?xml version=\"1.0\"?>\n";
    paraviewCylinderFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
    paraviewCylinderFile << " <UnstructuredGrid GhostLevel=\"0\">\n";
    paraviewCylinderFile << "  <Piece NumberOfPoints=\"" << 2 * Cnumber << "\" NumberOfCells=\"" << Cnumber << "\">\n";
    paraviewCylinderFile << "   <PointData>\n";
    
    // Write the radius of the cylinder
    paraviewCylinderFile << "    <DataArray type=\"Float64\" Name=\"radius\" format=\"ascii\">\n";
    for (int i = 0; i < Cnumber; ++i) {
        paraviewCylinderFile << cylinders[i].R << "\n";
        paraviewCylinderFile << cylinders[i].R << "\n"; // Each cylinder has two points with the same radius
    }
    paraviewCylinderFile << "    </DataArray>\n";

    // Write the index
    paraviewCylinderFile << "    <DataArray type=\"Float64\" Name=\"cylinderIndex\" format=\"ascii\">\n";
    for (int i = 0; i < Cnumber; ++i) {
        paraviewCylinderFile << cylinders[i].index << "\n";
        paraviewCylinderFile << cylinders[i].index << "\n"; // Each cylinder has two points with the same index
    }
    paraviewCylinderFile << "    </DataArray>\n";

    // Write velocity
    paraviewCylinderFile << "    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < Cnumber; ++i) {
        cylinders[i].trans.printLine(paraviewCylinderFile);
        cylinders[i].trans.printLine(paraviewCylinderFile); // Each cylinder has two points with the same velocity
    }
    
    // add the vtk data for the forces (needed algorithm to comoute forces)
//        paraviewCylinderFile << "    <DataArray type=\"Float64\" Name=\"FParticle\" NumberOfComponents=\"3\"/>\n";
//    for (int i = 0; i < Onumber; ++i) {
//        cylinders[i].FParticle.printLine(paraviewObjectFile);
//    }
//    paraviewCylinderFile << "    <DataArray type=\"Float64\" Name=\"FHydro\" NumberOfComponents=\"3\"/>\n";
//    for (int i = 0; i < Onumber; ++i) {
//        cylinders[i].FHydro.printLine(paraviewObjectFile);
//    }
    
    // Write the position 
    paraviewCylinderFile << "    </DataArray>\n";

    paraviewCylinderFile << "   </PointData>\n";
    paraviewCylinderFile << "   <CellData>\n";
    paraviewCylinderFile << "   </CellData>\n";
    paraviewCylinderFile << "   <Points>\n";
    paraviewCylinderFile << "    <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < Cnumber; ++i) {
        paraviewCylinderFile << cylinders[i].p1.x << " " << cylinders[i].p1.y << " " << cylinders[i].p1.z << "\n";
        paraviewCylinderFile << cylinders[i].p2.x << " " << cylinders[i].p2.y << " " << cylinders[i].p2.z << "\n";
    }
    
    // Write all the shit that Paraview needs for vtk files (connectivity matrix etc..)
    paraviewCylinderFile << "    </DataArray>\n";
    paraviewCylinderFile << "   </Points>\n";
    paraviewCylinderFile << "   <Cells>\n";
    paraviewCylinderFile << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int i = 0; i < Cnumber; ++i) {
        paraviewCylinderFile << 2 * i << " " << 2 * i + 1 << "\n";
    }
    paraviewCylinderFile << "    </DataArray>\n";
    paraviewCylinderFile << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 1; i <= Cnumber; ++i) {
        paraviewCylinderFile << 2 * i << "\n";
    }
    paraviewCylinderFile << "    </DataArray>\n";
    paraviewCylinderFile << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < Cnumber; ++i) {
        paraviewCylinderFile << vtkLineType << "\n"; // VTK_LINE type
    }
    paraviewCylinderFile << "    </DataArray>\n";
    paraviewCylinderFile << "   </Cells>\n";
    paraviewCylinderFile << "  </Piece>\n";
    paraviewCylinderFile << " </UnstructuredGrid>\n";
    paraviewCylinderFile << "</VTKFile>";
    
    // Close the file
    paraviewCylinderFile.close();
    
}

void IO::exportRecycleFluid(const LB& lb, const string& fluidRecycleFile) {
    // exports a format readable by this code itself, allowing the re-use of computed particle data.
    // needs further work...

    ofstream recycleFluidFile;
    recycleFluidFile.open(fluidRecycleFile.c_str());

    recycleFluidFile << lb.lbSize[0] << " " << lb.lbSize[1] << " " << lb.lbSize[2] << " " << lb.activeNodes.size() << "\n";

    //cout << endl;
    for (nodeMap::const_iterator imap = lb.nodes.begin(); imap != lb.nodes.end(); imap++) {
        const unsigned int it = imap->first;
        const node* nodeHere = &imap->second;
        // to gain precision, scale density and mass by 1 (initDensity)
        if (nodeHere->isActive()) {
            recycleFluidFile << it << " " << nodeHere->type << " " << nodeHere->n - lb.fluidMaterial.initDensity << " " << nodeHere->u.dot(Xp) << " " << nodeHere->u.dot(Yp) << " " << nodeHere->u.dot(Zp) << " " << nodeHere->mass - lb.fluidMaterial.initDensity << " " << nodeHere->visc << " ";
            for (int j = 0; j < lbmDirec; ++j) {
                recycleFluidFile << nodeHere->f[j] << " ";
            }
            recycleFluidFile << endl;
        }
    }
    recycleFluidFile.close();
}

//void IO::exportRecycleFluid(const LB& lb, const string& fluidRecycleFile) {
//    // exports a format readable by this code itself, allowing the re-use of computed particle data.
//    // needs further work...
//
//    ofstream recycleFluidFile;
//    recycleFluidFile.open(fluidRecycleFile.c_str());
//
//    recycleFluidFile << lb.lbSize[0] << " " << lb.lbSize[1] << " " << lb.lbSize[2] << "\n";
//
//    for (int it = 0; it < lb. totPossibleNodes; ++it) {
//        recycleFluidFile << lb.types[it].getType() << " ";
//    }
//    recycleFluidFile << endl;
//    //cout << endl;
//    for (int it = 0; it < lb.totPossibleNodes; ++it) {
//        if (lb.types[it].isActive()) {
//            // to gain precision, scale density and mass by 1 (initDensity)
//            recycleFluidFile << lb.nodes[it]->n - lb.fluidMaterial.initDensity << " " << lb.nodes[it]->u.dot(Xp) << " " << lb.nodes[it]->u.dot(Yp) << " " << lb.nodes[it]->u.dot(Zp) << " " << lb.nodes[it]->mass - lb.fluidMaterial.initDensity << " " << lb.nodes[it]->visc << " ";
//            for (int j = 0; j < lbmDirec; ++j) {
//                recycleFluidFile << lb.nodes[it]->f[j] << " ";
//            }
//            recycleFluidFile << endl;
//        }
//    }
//    recycleFluidFile.close();
//}

void IO::exportRecycleParticles(const elmtList& elmts, const pbcList& pbcs, const string& partRecycleFile) {
    // exports a format readable by this code itself, allowing the re-use of computed particle data.
    // needs further work...


    unsigned int activeNumber = 0;
    for (int i = 0; i < elmts.size(); ++i) {
        if (elmts[i].active) {
            ++activeNumber;
        }
    }

    ofstream recycleParticleFile;
    recycleParticleFile.open(partRecycleFile.c_str());
    recycleParticleFile << std::setprecision(8);
    recycleParticleFile << std::scientific;
    // total number of element

    recycleParticleFile << activeNumber << "\n";

    unsigned int index = 0;
    for (int i = 0; i < elmts.size(); ++i) {
        if (elmts[i].active) {
            // import variables
            recycleParticleFile << index << "\t";
            recycleParticleFile << elmts[i].size << "\t";
            recycleParticleFile << elmts[i].radius << "\t";
            // element center could be out of domain of there are periodic walls, fix this.
            tVect printPosition = elmts[i].x0;
            for (int b = 0; b < pbcs.size(); ++b) {
                const tVect pbcVector = pbcs[b].v;
                const double leftDist = pbcs[b].pl1.dist(elmts[i].x0);
                const double rightDist = pbcs[b].pl2.dist(elmts[i].x0);
                // first plane of couple (left)
                if (leftDist < 0.0) {
                    printPosition = printPosition + pbcVector;
                }
                if (rightDist < 0.0) {
                    printPosition = printPosition - pbcVector;
                }
            }
            printPosition.print(recycleParticleFile);
            elmts[i].x1.print(recycleParticleFile);
            elmts[i].w0.print(recycleParticleFile);
            elmts[i].q0.print(recycleParticleFile);
            elmts[i].q1.printLine(recycleParticleFile);
            ++index;
        }
    }
    recycleParticleFile.close();
}

void IO::exportLagrangianParaviewFluid(const LB& lb, const string& fluidFile) {

    const int one = 1;

    // start printing all the crap required for Paraview
    // header file opening
    ofstream paraviewFluidFile;
    paraviewFluidFile.open(fluidFile.c_str());
    // writing on header file
    paraviewFluidFile << std::scientific << std::setprecision(4);

    paraviewFluidFile << "<?xml version=\"1.0\"?>\n";
    paraviewFluidFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
    paraviewFluidFile << " <UnstructuredGrid GhostLevel=\"0\">\n";
    paraviewFluidFile << "  <Piece NumberOfPoints=\"" << lb.activeNodes.size() << "\" NumberOfCells=\"" << lb.activeNodes.size() << "\">\n";
    paraviewFluidFile << "   <PointData>\n";
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\" format=\"ascii\"/>\n";
    for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
        const node* nodeHere = *it;
        tVect uPhysic = nodeHere->u * lb.unit.Speed;
        uPhysic.printLine(paraviewFluidFile);
    }
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"pressure\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
        const node* nodeHere = *it;
        paraviewFluidFile << 0.3333333 * (nodeHere->n - lb.fluidMaterial.initDensity) * lb.unit.Pressure << "\n";
    }
    //    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"smoothedPressure\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    //    for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
    //        const node* nodeHere = *it;
    //        paraviewFluidFile << nodeHere->smoothedPressure * lb.unit.Pressure << "\n";
    //    }
    //    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"n\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    //    for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
    //        const node* nodeHere = *it;
    //        paraviewFluidFile << nodeHere->n << "\n";
    //    }
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"dynVisc\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
        const node* nodeHere = *it;
        paraviewFluidFile << nodeHere->visc * lb.unit.DynVisc << "\n";
    }
    //    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"applySlip\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    //    for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
    //        const node* nodeHere = *it;
    //        paraviewFluidFile << nodeHere->applySlip << "\n";
    //    }
//    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"viscosity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
//    for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
//        const node* nodeHere = *it;
//        paraviewFluidFile << nodeHere->visc * lb.unit.KinVisc << "\n";
//    }
    if (lb.fluidMaterial.rheologyModel == MUI || lb.fluidMaterial.rheologyModel == FRICTIONAL || lb.fluidMaterial.rheologyModel == VOELLMY) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"friction\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
        for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
            const node* nodeHere = *it;
            paraviewFluidFile << nodeHere->friction << "\n";
        }
    }
    //    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"friction\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    //    for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
    //        const node* nodeHere = *it;
    //        paraviewFluidFile << nodeHere->friction << "\n";
    //    }
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"mass\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
        const node* nodeHere = *it;
        paraviewFluidFile << nodeHere->mass * lb.unit.Density << "\n";
    }
    paraviewFluidFile << "    <DataArray type=\"Int8\" Name=\"type\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"8\"/>\n";
    for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
        const node* nodeHere = *it;
        paraviewFluidFile << static_cast<unsigned>(nodeHere->type) << "\n";
    }
    //    for (int j = 1; j < lbmDirec; ++j) {
    //        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"delta" << j << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    //        for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
    //            const node* nodeHere = *it;
    //            if (nodeHere->curved != 0) {
    //                const tVect deltaVector = nodeHere->curved->delta[j] * lb.unit.Length * v[j];
    //                deltaVector.printLine(paraviewFluidFile);
    //            } else {
    //                Zero.printLine(paraviewFluidFile);
    //
    //            }
    //        }
    //        paraviewFluidFile << "    </DataArray>\n";
    //    }
#ifdef DEBUG 
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"age\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"8\"/>\n";
    for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
        const node* nodeHere = *it;
        paraviewFluidFile << nodeHere->age << "\n";
    }
    //    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"vv0\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    //    for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
    //        const node* nodeHere = *it;
    //        const double deltaZero = (nodeHere->f[0] - coeff[0]);
    //        paraviewFluidFile << deltaZero << "\n";
    //    }
    //    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"surfNormal\" NumberOfComponents=\"3\" format=\"ascii\"/>\n";
    //    for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
    //        const node* nodeHere = *it;
    //        nodeHere->surfaceNormal.printLine(paraviewFluidFile);
    //    }
    //    for (int j = 1; j < lbmDirec2D; ++j) {
    //        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"vv" << j << "\" NumberOfComponents=\"3\" format=\"ascii\"/>\n";
    //        for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
    //            const node* nodeHere = *it;
    //            const tVect deltaVector = (nodeHere->f[two_dim[j]] - coeff[two_dim[j]]) * v[two_dim[j]];
    //            deltaVector.printLine(paraviewFluidFile);
    //        }
    //    }
#endif
    paraviewFluidFile << "   </PointData>\n";
    paraviewFluidFile << "   <CellData>\n";
    paraviewFluidFile << "   </CellData>\n";
    paraviewFluidFile << "   <Points>\n";
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n";
    for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
        const node* nodeHere = *it;
        const unsigned int index = nodeHere->coord;
        const tVect positionHere = lb.getPosition(index) * lb.unit.Length;
        positionHere.printLine(paraviewFluidFile);
    }
    paraviewFluidFile << "   </Points>\n";
    paraviewFluidFile << "   <Cells>\n";
    paraviewFluidFile << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int i = 1; i < lb.activeNodes.size() + 1; ++i) {
        paraviewFluidFile << i - 1 << "\n";
    }
    paraviewFluidFile << "    </DataArray>\n";
    paraviewFluidFile << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 1; i < lb.activeNodes.size() + 1; ++i) {
        paraviewFluidFile << i << "\n";
    }
    paraviewFluidFile << "    </DataArray>\n";
    paraviewFluidFile << "    <DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < lb.activeNodes.size(); ++i) {
        paraviewFluidFile << one << "\n";
    }
    paraviewFluidFile << "    </DataArray>\n";
    paraviewFluidFile << "   </Cells>\n";
    paraviewFluidFile << "  </Piece>\n";
    //    paraviewHeaderFile<<"  <Piece Source=\""<<pieceFile<<"\"/>\n";
    paraviewFluidFile << " </UnstructuredGrid>\n";
    paraviewFluidFile << "</VTKFile>";

    // data file closing
    paraviewFluidFile.close();
}

void IO::exportLagrangianParaviewFluid_binaryv3(const LB& lb, const string& fluidFile) {
    /**
     * This function is a rewrite of exportLagrangianParaviewFluid() that writes to a binary vtkhdf5 format
     * It is intended to provide much faster fluid export performance
     **/

    // start printing all the crap required for Paraview
    // header file opening
    ofstream paraviewFluidFile;
    paraviewFluidFile.open(fluidFile.c_str(), std::ios::binary);
    // paraviewFluidFile << std::scientific << std::setprecision(4);
    // writing on header file

    // Endianness (the order of bits within a byte) depends on the processor hardware
    // but it's probably LittleEndian, IBM processors are the only default BigEndian you're likely to come across
    static const uint16_t m_endianCheck(0x00ff);
    const bool is_big_endian ( *((const uint8_t*)&m_endianCheck) == 0x0);
    paraviewFluidFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" << (is_big_endian ? "BigEndian" : "LittleEndian") << "\"  header_type=\"UInt32\">\n";
    paraviewFluidFile << " <UnstructuredGrid GhostLevel=\"0\">\n";
    paraviewFluidFile << "  <Piece NumberOfPoints=\"" << lb.activeNodes.size() << "\" NumberOfCells=\"" << lb.activeNodes.size() << "\">\n";
    paraviewFluidFile << "   <PointData>\n";
    unsigned int offset = 0;
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\"/>\n";
    offset += lb.activeNodes.size() * 3 * sizeof(double) + sizeof(unsigned int);
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    if (lb.fluidMaterial.rheologyModel != NEWTONIAN || lb.fluidMaterial.turbulenceOn) {
        offset += lb.activeNodes.size() * sizeof(double) + sizeof(unsigned int);
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"dynVisc\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    }
    offset += lb.activeNodes.size() * sizeof(double) + sizeof(unsigned int);
    if (lb.fluidMaterial.rheologyModel == MUI || lb.fluidMaterial.rheologyModel == FRICTIONAL || lb.fluidMaterial.rheologyModel == VOELLMY) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"friction\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
        offset += lb.activeNodes.size() * sizeof(double) + sizeof(unsigned int);
    }
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"mass\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    offset += lb.activeNodes.size() * sizeof(double) + sizeof(unsigned int);
    paraviewFluidFile << "    <DataArray type=\"UInt8\" Name=\"type\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"16\"/>\n";
    offset += lb.activeNodes.size() * sizeof(unsigned char) + sizeof(unsigned int);
#ifdef DEBUG
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"age\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"8\"/>\n";
    offset += lb.activeNodes.size() * sizeof(double) + sizeof(unsigned int);
#endif
    paraviewFluidFile << "   </PointData>\n";
    paraviewFluidFile << "   <CellData>\n";
    paraviewFluidFile << "   </CellData>\n";
    paraviewFluidFile << "   <Points>\n";
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += lb.activeNodes.size() * 3 * sizeof(double) + sizeof(unsigned int);
    paraviewFluidFile << "   </Points>\n";
    paraviewFluidFile << "   <Cells>\n";
    paraviewFluidFile << "    <DataArray type=\"UInt32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += lb.activeNodes.size() * sizeof(unsigned int) + sizeof(unsigned int);
    paraviewFluidFile << "    <DataArray type=\"UInt32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += lb.activeNodes.size() * sizeof(unsigned int) + sizeof(unsigned int);
    paraviewFluidFile << "    <DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += lb.activeNodes.size() * sizeof(unsigned char) + sizeof(unsigned int);
    paraviewFluidFile << "   </Cells>\n";
    paraviewFluidFile << "  </Piece>\n";
    paraviewFluidFile << " </UnstructuredGrid>\n";
    paraviewFluidFile << " <AppendedData encoding=\"raw\">\n  _";
      /**
     * Based on the sparse documentation at https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
     * and alot of testing.
     * Inside <AppendedData> the binary dump must be preceded by an underscore (_)
     * Each DataArray's binary dump must be preceded by it's length.
     * The length should be exported as the integer type specified as header_type in the opening <VTKFile> tag
     * The offset specified in the above <DataArray> tag refers to the offset from the start of the whole binary dump to the start of the length
     */
    // Allocate a buffer equal to size of the largest data array
    // Allocate once rather than allocating and freeing per export
    char *const t_buffer = static_cast<char*>(malloc(lb.activeNodes.size() * 3 * sizeof(double)));
    unsigned char *const uc_buffer = reinterpret_cast<unsigned char*>(t_buffer);
    double *const d_buffer = reinterpret_cast<double*>(t_buffer);
    tVect *const v_buffer = reinterpret_cast<tVect*>(t_buffer);
    unsigned int* const ui_buffer = reinterpret_cast<unsigned int*>(t_buffer);
    // Velocity
    offset = lb.activeNodes.size() * 3 * sizeof(double);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < lb.activeNodes.size(); ++i) {
        v_buffer[i] = lb.activeNodes[i]->u * lb.unit.Speed;
    }
    paraviewFluidFile.write(t_buffer, offset);
    // Pressure
    offset = lb.activeNodes.size() * sizeof(double);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    const double THIRD_PRESSURE = 0.3333333 * lb.unit.Pressure;
    for (unsigned int i = 0; i < lb.activeNodes.size(); ++i) {
        d_buffer[i] = (lb.activeNodes[i]->n - lb.fluidMaterial.initDensity) * THIRD_PRESSURE;
    }
    paraviewFluidFile.write(t_buffer, offset);
    // Dynamic Viscosity
    if (lb.fluidMaterial.rheologyModel != NEWTONIAN || lb.fluidMaterial.turbulenceOn) {
        offset = lb.activeNodes.size() * sizeof(double);
        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        for (unsigned int i = 0; i < lb.activeNodes.size(); ++i) {
            d_buffer[i] = lb.activeNodes[i]->visc * lb.unit.DynVisc;
        }
        paraviewFluidFile.write(t_buffer, offset);
    }
    // Friction
    if (lb.fluidMaterial.rheologyModel == MUI || lb.fluidMaterial.rheologyModel == FRICTIONAL || lb.fluidMaterial.rheologyModel == VOELLMY) {
        offset = lb.activeNodes.size() * sizeof(double);
        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        for (unsigned int i = 0; i < lb.activeNodes.size(); ++i) {
            d_buffer[i] = lb.activeNodes[i]->friction;
        }
        paraviewFluidFile.write(t_buffer, offset);
    }
    // Mass
    offset = lb.activeNodes.size() * sizeof(double);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < lb.activeNodes.size(); ++i) {
        d_buffer[i] = lb.activeNodes[i]->mass * lb.unit.Density;
    }
    paraviewFluidFile.write(t_buffer, offset);
    // Type
    offset = lb.activeNodes.size() * sizeof(unsigned char);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    std::fill(uc_buffer, uc_buffer + lb.activeNodes.size(), (unsigned char)2);
    for (unsigned int i = 0; i < lb.activeNodes.size(); ++i) {
        uc_buffer[i] = lb.activeNodes[i]->type;
    }
    paraviewFluidFile.write(t_buffer, offset);
#ifdef DEBUG
    // Age
    offset = lb.activeNodes.size() * sizeof(double);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    memset(d_buffer, 0, sizeof(double)* lb.activeNodes.size());
    for (unsigned int i = 0; i < lb.activeNodes.size(); ++i) {
        d_buffer[i] = lb.activeNodes[i]->age;
    }
    paraviewFluidFile.write(t_buffer, offset);
#endif
    // Points
    offset = lb.activeNodes.size() * 3 * sizeof(double);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < lb.activeNodes.size(); ++i) {
        v_buffer[i] = lb.getPosition(lb.activeNodes[i]->coord) * lb.unit.Length;
    }
    paraviewFluidFile.write(t_buffer, offset);
    // Connectivity
    offset = lb.activeNodes.size() * sizeof(unsigned int);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < lb.activeNodes.size(); ++i) {
        ui_buffer[i] = i ;
    }
    paraviewFluidFile.write(t_buffer, offset);
    // Offsets
    offset = lb.activeNodes.size() * sizeof(unsigned int);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < lb.activeNodes.size(); ++i) {
        ui_buffer[i] = i + 1;
    }
    paraviewFluidFile.write(t_buffer, offset);
    // Types
    offset = lb.activeNodes.size() * sizeof(unsigned char);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    std::fill(uc_buffer, uc_buffer + lb.activeNodes.size(), (unsigned char)1);
    paraviewFluidFile.write(t_buffer, offset);
    paraviewFluidFile << "</AppendedData>";
    paraviewFluidFile << "</VTKFile>\n";
    // data file closing
    paraviewFluidFile.close();
    free(t_buffer);    
}

void IO::exportEulerianParaviewFluid(const LB& lb, const string& fluidFile) {

    const double zero = 0.0;

    // start printing all the crap required for Paraview
    // header file opening
    ofstream paraviewFluidFile;

    paraviewFluidFile.open(fluidFile.c_str());
//    paraviewFluidFile << std::scientific << std::setprecision(4);
    // writing on header file
    paraviewFluidFile << "<VTKFile type=\"ImageData\" version=\"0.1\">\n";
    paraviewFluidFile << " <ImageData WholeExtent=\"0 " << lb.lbSize[0] - 1 << " 0 " << lb.lbSize[1] - 1 << " 0 " << lb.lbSize[2] - 1 << "\" "
            << "Origin=\"" << -0.5 * lb.unit.Length << " " << -0.5 * lb.unit.Length << " " << -0.5 * lb.unit.Length << "\" "
            << "Spacing=\"" << lb.unit.Length << " " << lb.unit.Length << " " << lb.unit.Length << "\">\n";
    paraviewFluidFile << "  <Piece Extent=\"0 " << lb.lbSize[0] - 1 << " 0 " << lb.lbSize[1] - 1 << " 0 " << lb.lbSize[2] - 1 << "\">\n";
    paraviewFluidFile << "   <PointData>\n";
    paraviewFluidFile << "    <DataArray type=\"Int8\" Name=\"type\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"8\">\n";
    //    paraviewFluidFile << std::fixed << std::setprecision(0);
    for (int i = 0; i < lb.totPossibleNodes; ++i) {
        if (lb.nodes.count(i) != 0) {
            if (!lb.nodes.at(i).isInsideParticle()) {
                paraviewFluidFile << static_cast<unsigned>(lb.nodes.at(i).type) << " ";
            } else {
                paraviewFluidFile << 1 << " ";
            }
        } else {
            paraviewFluidFile << 2 << " ";
        }
    }
    paraviewFluidFile << "    </DataArray>\n";
    paraviewFluidFile << std::scientific << std::setprecision(4);
    //    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"n\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\">\n";
    //    for (int i = 0; i < lb.totPossibleNodes; ++i) {
    //        if (lb.nodes.count(i)!=0) {
    //            paraviewFluidFile << lb.nodes.at(i).n - 1.0 << " ";
    //        } else {
    //            paraviewFluidFile << 0 << " ";
    //        }
    //    }
    //    paraviewFluidFile << "    </DataArray>\n";
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < lb.totPossibleNodes; ++i) {
        if (lb.nodes.count(i) != 0) {
            const tVect uPhysic = lb.nodes.at(i).u * lb.unit.Speed;
            uPhysic.print(paraviewFluidFile);
        } else {
            Zero.print(paraviewFluidFile);
        }
    }
    paraviewFluidFile << "    </DataArray>\n";
    //    for (int j = 1; j < lbmDirec; ++j) {
    //        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"delta" << j << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    //        for (int i = 0; i < lb.totPossibleNodes; ++i) {
    //            if (lb.nodes.count(i) != 0) {
    //                if (lb.nodes.at(i).curved!=0) {
    //                    const tVect deltaVector = lb.nodes.at(i).curved->delta[j] * lb.unit.Length*vDirec[j];
    //                    deltaVector.print(paraviewFluidFile);
    //                } else {
    //                    Zero.print(paraviewFluidFile);
    //                }
    //            } else {
    //                Zero.print(paraviewFluidFile);
    //            }
    //        }
    //        paraviewFluidFile << "    </DataArray>\n";
    //    }
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"pressure\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\">\n";
    for (int i = 0; i < lb.totPossibleNodes; ++i) {
        if (lb.nodes.count(i) != 0) {
            if (lb.nodes.at(i).n != 0.0) {
                paraviewFluidFile << 0.3333333 * (lb.nodes.at(i).n - 1.0) * lb.unit.Pressure << " ";
            } else {
                paraviewFluidFile << lb.nodes.at(i).n << " ";
            }
        } else {
            paraviewFluidFile << zero << " ";
        }
    }
    paraviewFluidFile << "    </DataArray>\n";
    if (lb.fluidMaterial.rheologyModel != NEWTONIAN || lb.fluidMaterial.turbulenceOn) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"dynVisc\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\">\n";
        for (int i = 0; i < lb.totPossibleNodes; ++i) {
            if (lb.nodes.count(i) != 0) {
                paraviewFluidFile << lb.nodes.at(i).visc * lb.unit.DynVisc << " ";
            } else {
                paraviewFluidFile << zero << " ";
            }
        }
        paraviewFluidFile << "    </DataArray>\n";
    }
    if (lb.fluidMaterial.rheologyModel == MUI || lb.fluidMaterial.rheologyModel == FRICTIONAL || lb.fluidMaterial.rheologyModel == VOELLMY) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"friction\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\">\n";
        for (int i = 0; i < lb.totPossibleNodes; ++i) {
            if (lb.nodes.count(i) != 0) {
                paraviewFluidFile << lb.nodes.at(i).friction << " ";
            } else {
                paraviewFluidFile << zero << " ";
            }
        }
        paraviewFluidFile << "    </DataArray>\n";
    }
    if (lb.freeSurface) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"amass\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\">\n";
        for (int i = 0; i < lb.totPossibleNodes; ++i) {
            if (lb.nodes.count(i) != 0) {
                paraviewFluidFile << lb.nodes.at(i).mass * lb.unit.Density << " ";
            } else {
                paraviewFluidFile << zero << " ";
            }
        }
        paraviewFluidFile << "    </DataArray>\n";
    }
    if (lb.lbTopography) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"topSurface\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\">\n";
        for (int i = 0; i < lb.totPossibleNodes; ++i) {
            if (lb.nodes.count(i) != 0) {
                if (lb.nodes.at(i).isTopography()) {
                    paraviewFluidFile << 1.0 << " ";
                } else {
                    paraviewFluidFile << zero << " ";
                }
            } else {
                paraviewFluidFile << zero << " ";
            }
        }
        paraviewFluidFile << "    </DataArray>\n";
    }
    //    if (demSolver) {
    //        paraviewFluidFile << "    <DataArray type=\"Int16\" Name=\"solidIndex\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"8\">\n";
    //        paraviewFluidFile << std::fixed << std::setprecision(0);
    //        for (int i = 0; i < lb.totPossibleNodes; ++i) {
    //            if (lb.nodes.count(i) != 0) {
    //                paraviewFluidFile << lb.nodes.at(i).getSolidIndex() << " ";
    //            } else {
    //                paraviewFluidFile << zero << " ";
    //            }
    //        }
    //        paraviewFluidFile << "    </DataArray>\n";
    //        paraviewFluidFile << std::scientific << std::setprecision(4);
    //    }
    paraviewFluidFile << "   </PointData>\n";
    paraviewFluidFile << "   <CellData>\n";
    paraviewFluidFile << "   </CellData>\n";
    paraviewFluidFile << "  </Piece>\n";
    paraviewFluidFile << " </ImageData>\n";
    paraviewFluidFile << "</VTKFile>\n";
    // data file closing
    paraviewFluidFile.close();
}

void IO::exportEulerianParaviewFluid_binary(const LB& lb, const string& fluidFile) {
    /**
     * This function is a rewrite of exportEulerianParaviewFluid() that writes to a binary vtkhdf5 format
     * It is intended to provide much faster fluid export performance
     **/

    // start printing all the crap required for Paraview
    // header file opening
    ofstream paraviewFluidFile;

    paraviewFluidFile.open(fluidFile.c_str(), std::ios::binary);
//    paraviewFluidFile << std::scientific << std::setprecision(4);
    // writing on header file

    // Endianness (the order of bits within a byte) depends on the processor hardware
    // but it's probably LittleEndian, IBM processors are the only default BigEndian you're likely to come across
    static const uint16_t m_endianCheck(0x00ff);
    const bool is_big_endian ( *((const uint8_t*)&m_endianCheck) == 0x0);
    paraviewFluidFile << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"" << (is_big_endian ? "BigEndian" : "LittleEndian") << "\"  header_type=\"UInt32\">\n";
    paraviewFluidFile << " <ImageData WholeExtent=\"0 " << lb.lbSize[0] - 1 << " 0 " << lb.lbSize[1] - 1 << " 0 " << lb.lbSize[2] - 1 << "\" "
            << "Origin=\"" << -0.5 * lb.unit.Length << " " << -0.5 * lb.unit.Length << " " << -0.5 * lb.unit.Length << "\" "
            << "Spacing=\"" << lb.unit.Length << " " << lb.unit.Length << " " << lb.unit.Length << "\">\n";
    paraviewFluidFile << "  <Piece Extent=\"0 " << lb.lbSize[0] - 1 << " 0 " << lb.lbSize[1] - 1 << " 0 " << lb.lbSize[2] - 1 << "\">\n";
    paraviewFluidFile << "   <PointData>\n";
    unsigned int offset = 0;
    paraviewFluidFile << "    <DataArray type=\"UInt8\" Name=\"type\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"8\"/>\n";
    offset += lb.totPossibleNodes * sizeof(unsigned char) + sizeof(unsigned int);
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\"/>\n";
    offset += lb.totPossibleNodes * 3 * sizeof(double) + sizeof(unsigned int);
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    offset += lb.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
    if (lb.fluidMaterial.rheologyModel != NEWTONIAN || lb.fluidMaterial.turbulenceOn) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"dynVisc\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
        offset += lb.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
    }
    if (lb.fluidMaterial.rheologyModel == MUI || lb.fluidMaterial.rheologyModel == FRICTIONAL || lb.fluidMaterial.rheologyModel == VOELLMY) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"friction\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
        offset += lb.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
    }
    if (lb.freeSurface) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"amass\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
        offset += lb.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
    }
    if (lb.lbTopography) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"topSurface\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
        offset += lb.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
    }
    paraviewFluidFile << "   </PointData>\n";
    paraviewFluidFile << "   <CellData>\n";
    paraviewFluidFile << "   </CellData>\n";
    paraviewFluidFile << "  </Piece>\n";
    paraviewFluidFile << " </ImageData>\n";
    paraviewFluidFile << " <AppendedData encoding=\"raw\">\n  _";
    /**
     * Based on the sparse documentation at https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
     * and alot of testing.
     * Inside <AppendedData> the binary dump must be preceded by an underscore (_)
     * Each DataArray's binary dump must be preceded by it's length.
     * The length should be exported as the integer type specified as header_type in the opening <VTKFile> tag
     * The offset specified in the above <DataArray> tag refers to the offset from the start of the whole binary dump to the start of the length
     */
    static const double zero = 0;
    static const tVect tVect_zero = {0,0,0};
    // May be faster to iterate the (sorted) map and catch missing items
    // Might be possible to do inline compression, this would slow export but unclear by how much
    offset = lb.totPossibleNodes * sizeof(unsigned char);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (int i = 0; i < lb.totPossibleNodes; ++i) {
        const auto &node = lb.nodes.find(i);
        if (node != lb.nodes.end()) {
            if (!node->second.isInsideParticle()) {
                paraviewFluidFile.write(reinterpret_cast<const char*>(&node->second.type), sizeof(node->second.type));
            } else {
                const unsigned char t = 1;
                paraviewFluidFile.write(reinterpret_cast<const char*>(&t), sizeof(unsigned char));
            }
        } else {
            const unsigned char t = 2;
            paraviewFluidFile.write(reinterpret_cast<const char*>(&t), sizeof(unsigned char));
        }
    }
    offset = lb.totPossibleNodes * 3 * sizeof(double);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (int i = 0; i < lb.totPossibleNodes; ++i) {
        const auto &node = lb.nodes.find(i);
        if (node != lb.nodes.end()) {
            static const tVect uPhysic = lb.nodes.at(i).u * lb.unit.Speed;
            paraviewFluidFile.write(reinterpret_cast<const char*>(&uPhysic), sizeof(tVect));
        } else {
            paraviewFluidFile.write(reinterpret_cast<const char*>(&tVect_zero), sizeof(tVect));
        }
    }
    offset = lb.totPossibleNodes * sizeof(double);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (int i = 0; i < lb.totPossibleNodes; ++i) {
        const auto &node = lb.nodes.find(i);
        if (node != lb.nodes.end()) {
            const double t = 0.3333333 * (node->second.n - 1.0) * lb.unit.Pressure;
            paraviewFluidFile.write(reinterpret_cast<const char*>(&t), sizeof(double));
        } else {
            paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
        }
    }
    if (lb.fluidMaterial.rheologyModel != NEWTONIAN || lb.fluidMaterial.turbulenceOn) {
        offset = lb.totPossibleNodes * sizeof(double);
        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        for (int i = 0; i < lb.totPossibleNodes; ++i) {
            const auto &node = lb.nodes.find(i);
            if (node != lb.nodes.end()) {
                const double t = node->second.visc * lb.unit.DynVisc;
                paraviewFluidFile.write(reinterpret_cast<const char*>(&t), sizeof(double));
            } else {
                paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
            }
        }
    }
    if (lb.fluidMaterial.rheologyModel == MUI || lb.fluidMaterial.rheologyModel == FRICTIONAL || lb.fluidMaterial.rheologyModel == VOELLMY) {
        offset = lb.totPossibleNodes * sizeof(double);
        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        for (int i = 0; i < lb.totPossibleNodes; ++i) {
            const auto &node = lb.nodes.find(i);
            if (node != lb.nodes.end()) {
                paraviewFluidFile.write(reinterpret_cast<const char*>(&node->second.friction), sizeof(double));
            } else {
                paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
            }
        }
    }
    if (lb.freeSurface) {
        offset = lb.totPossibleNodes * sizeof(double);
        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        for (int i = 0; i < lb.totPossibleNodes; ++i) {
            const auto &node = lb.nodes.find(i);
            if (node != lb.nodes.end()) {
                const double t = node->second.mass * lb.unit.Density;
                paraviewFluidFile.write(reinterpret_cast<const char*>(&t), sizeof(double));
            } else {
                paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
            }
        }
    }
    if (lb.lbTopography) {
        offset = lb.totPossibleNodes * sizeof(double);
        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        for (int i = 0; i < lb.totPossibleNodes; ++i) {
            const auto &node = lb.nodes.find(i);
            if (node != lb.nodes.end()) {
                if (node->second.isTopography()) {
                    const double one = 1.0;
                    paraviewFluidFile.write(reinterpret_cast<const char*>(&one), sizeof(double));
                } else {
                    paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
                }
            } else {
                paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
            }
        }
    }
    paraviewFluidFile << "</AppendedData>";
    paraviewFluidFile << "</VTKFile>\n";
    // data file closing
    paraviewFluidFile.close();
}

void IO::exportEulerianParaviewFluid_binaryv2(const LB& lb, const string& fluidFile) {
    /**
     * This function is a rewrite of exportEulerianParaviewFluid() that writes to a binary vtkhdf5 format
     * It is intended to provide much faster fluid export performance
     **/

    // start printing all the crap required for Paraview
    // header file opening
    ofstream paraviewFluidFile;

    paraviewFluidFile.open(fluidFile.c_str(), std::ios::binary);
//    paraviewFluidFile << std::scientific << std::setprecision(4);
    // writing on header file

    // Endianness (the order of bits within a byte) depends on the processor hardware
    // but it's probably LittleEndian, IBM processors are the only default BigEndian you're likely to come across
    static const uint16_t m_endianCheck(0x00ff);
    const bool is_big_endian ( *((const uint8_t*)&m_endianCheck) == 0x0);
    paraviewFluidFile << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"" << (is_big_endian ? "BigEndian" : "LittleEndian") << "\"  header_type=\"UInt32\">\n";
    paraviewFluidFile << " <ImageData WholeExtent=\"0 " << lb.lbSize[0] - 1 << " 0 " << lb.lbSize[1] - 1 << " 0 " << lb.lbSize[2] - 1 << "\" "
            << "Origin=\"" << -0.5 * lb.unit.Length << " " << -0.5 * lb.unit.Length << " " << -0.5 * lb.unit.Length << "\" "
            << "Spacing=\"" << lb.unit.Length << " " << lb.unit.Length << " " << lb.unit.Length << "\">\n";
    paraviewFluidFile << "  <Piece Extent=\"0 " << lb.lbSize[0] - 1 << " 0 " << lb.lbSize[1] - 1 << " 0 " << lb.lbSize[2] - 1 << "\">\n";
    paraviewFluidFile << "   <PointData>\n";
    unsigned int offset = 0;
    paraviewFluidFile << "    <DataArray type=\"UInt8\" Name=\"type\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"8\"/>\n";
    offset += lb.totPossibleNodes * sizeof(unsigned char) + sizeof(unsigned int);
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\"/>\n";
    offset += lb.totPossibleNodes * 3 * sizeof(double) + sizeof(unsigned int);
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    offset += lb.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
    if (lb.fluidMaterial.rheologyModel != NEWTONIAN || lb.fluidMaterial.turbulenceOn) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"dynVisc\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
        offset += lb.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
    }
    if (lb.fluidMaterial.rheologyModel == MUI || lb.fluidMaterial.rheologyModel == FRICTIONAL || lb.fluidMaterial.rheologyModel == VOELLMY) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"friction\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
        offset += lb.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
    }
    if (lb.freeSurface) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"amass\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
        offset += lb.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
    }
    if (lb.lbTopography) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"topSurface\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
        offset += lb.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
    }
    paraviewFluidFile << "   </PointData>\n";
    paraviewFluidFile << "   <CellData>\n";
    paraviewFluidFile << "   </CellData>\n";
    paraviewFluidFile << "  </Piece>\n";
    paraviewFluidFile << " </ImageData>\n";
    paraviewFluidFile << " <AppendedData encoding=\"raw\">\n  _";
    /**
     * Based on the sparse documentation at https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
     * and alot of testing.
     * Inside <AppendedData> the binary dump must be preceded by an underscore (_)
     * Each DataArray's binary dump must be preceded by it's length.
     * The length should be exported as the integer type specified as header_type in the opening <VTKFile> tag
     * The offset specified in the above <DataArray> tag refers to the offset from the start of the whole binary dump to the start of the length
     */
    static const double zero = 0;
    static const tVect tVect_zero = {0,0,0};
    // May be faster to iterate the (sorted) map and catch missing items
    // Might be possible to do inline compression, this would slow export but unclear by how much
    offset = lb.totPossibleNodes * sizeof(unsigned char);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    int i = 0;
    for (const auto &[key, node] : lb.nodes) {
        for (;i<key;++i) {
            const unsigned char t = 2;
            paraviewFluidFile.write(reinterpret_cast<const char*>(&t), sizeof(unsigned char));
        }
        if (!node.isInsideParticle()) {
            paraviewFluidFile.write(reinterpret_cast<const char*>(&node.type), sizeof(unsigned char));
        } else {
            const unsigned char t = 1;
            paraviewFluidFile.write(reinterpret_cast<const char*>(&t), sizeof(unsigned char));
        }
        ++i;
    }
    offset = lb.totPossibleNodes * 3 * sizeof(double);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    i = 0;
    for (const auto &[key, node] : lb.nodes) {
        for (;i<key;++i) {
            paraviewFluidFile.write(reinterpret_cast<const char*>(&tVect_zero), sizeof(tVect));
        }
        const tVect uPhysic = lb.nodes.at(i).u * lb.unit.Speed;
        paraviewFluidFile.write(reinterpret_cast<const char*>(&uPhysic), sizeof(tVect));
        ++i;
    }
    offset = lb.totPossibleNodes * sizeof(double);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    i = 0;
    for (const auto &[key, node] : lb.nodes) {
        for (;i<key;++i) {
            paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
        }
        const double t = 0.3333333 * (node.n - 1.0) * lb.unit.Pressure;
        paraviewFluidFile.write(reinterpret_cast<const char*>(&t), sizeof(double));
        ++i;
    }
    if (lb.fluidMaterial.rheologyModel != NEWTONIAN || lb.fluidMaterial.turbulenceOn) {
        offset = lb.totPossibleNodes * sizeof(double);
        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        i = 0;
        for (const auto &[key, node] : lb.nodes) {
            for (;i<key;++i) {
                paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
            }
            const double t = node.visc * lb.unit.DynVisc;
            paraviewFluidFile.write(reinterpret_cast<const char*>(&t), sizeof(double));
            ++i;
        }
    }
    if (lb.fluidMaterial.rheologyModel == MUI || lb.fluidMaterial.rheologyModel == FRICTIONAL || lb.fluidMaterial.rheologyModel == VOELLMY) {
        offset = lb.totPossibleNodes * sizeof(double);
        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        i = 0;
        for (const auto &[key, node] : lb.nodes) {
            for (;i<key;++i) {
                paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
            }
            paraviewFluidFile.write(reinterpret_cast<const char*>(&node.friction), sizeof(double));
            ++i;
        }
    }
    if (lb.freeSurface) {
        offset = lb.totPossibleNodes * sizeof(double);
        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        i = 0;
        for (const auto &[key, node] : lb.nodes) {
            for (;i<key;++i) {
                paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
            }
            const double t = node.mass * lb.unit.Density;
            paraviewFluidFile.write(reinterpret_cast<const char*>(&t), sizeof(double));
            ++i;
        }
    }
    if (lb.lbTopography) {
        offset = lb.totPossibleNodes * sizeof(double);
        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        i = 0;
        for (const auto &[key, node] : lb.nodes) {
            for (;i<key;++i) {
                paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
            }
            if (node.isTopography()) {
                const double one = 1.0;
                paraviewFluidFile.write(reinterpret_cast<const char*>(&one), sizeof(double));
            } else {
                paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
            }
            ++i;
        }
    }
    paraviewFluidFile << "</AppendedData>";
    paraviewFluidFile << "</VTKFile>\n";
    // data file closing
    paraviewFluidFile.close();
}

void IO::exportEulerianParaviewFluid_binaryv3(const LB& lb, const string& fluidFile) {
    /**
     * This function is a rewrite of exportEulerianParaviewFluid() that writes to a binary vtkhdf5 format
     * It is intended to provide much faster fluid export performance
     **/

    // start printing all the crap required for Paraview
    // header file opening
    ofstream paraviewFluidFile;

    paraviewFluidFile.open(fluidFile.c_str(), std::ios::binary);
//    paraviewFluidFile << std::scientific << std::setprecision(4);
    // writing on header file

    // Endianness (the order of bits within a byte) depends on the processor hardware
    // but it's probably LittleEndian, IBM processors are the only default BigEndian you're likely to come across
    static const uint16_t m_endianCheck(0x00ff);
    const bool is_big_endian ( *((const uint8_t*)&m_endianCheck) == 0x0);
    paraviewFluidFile << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"" << (is_big_endian ? "BigEndian" : "LittleEndian") << "\"  header_type=\"UInt32\">\n";
    paraviewFluidFile << " <ImageData WholeExtent=\"0 " << lb.lbSize[0] - 1 << " 0 " << lb.lbSize[1] - 1 << " 0 " << lb.lbSize[2] - 1 << "\" "
            << "Origin=\"" << -0.5 * lb.unit.Length << " " << -0.5 * lb.unit.Length << " " << -0.5 * lb.unit.Length << "\" "
            << "Spacing=\"" << lb.unit.Length << " " << lb.unit.Length << " " << lb.unit.Length << "\">\n";
    paraviewFluidFile << "  <Piece Extent=\"0 " << lb.lbSize[0] - 1 << " 0 " << lb.lbSize[1] - 1 << " 0 " << lb.lbSize[2] - 1 << "\">\n";
    paraviewFluidFile << "   <PointData>\n";
    unsigned int offset = 0;
    paraviewFluidFile << "    <DataArray type=\"UInt8\" Name=\"type\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"8\"/>\n";
    offset += lb.totPossibleNodes * sizeof(unsigned char) + sizeof(unsigned int);
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\"/>\n";
    offset += lb.totPossibleNodes * 3 * sizeof(double) + sizeof(unsigned int);
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    offset += lb.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
    if (lb.fluidMaterial.rheologyModel != NEWTONIAN || lb.fluidMaterial.turbulenceOn) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"dynVisc\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
        offset += lb.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
    }
    if (lb.fluidMaterial.rheologyModel == MUI || lb.fluidMaterial.rheologyModel == FRICTIONAL || lb.fluidMaterial.rheologyModel == VOELLMY) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"friction\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
        offset += lb.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
    }
    if (lb.freeSurface) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"amass\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
        offset += lb.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
    }
    if (lb.lbTopography) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"topSurface\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
        offset += lb.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
    }
    paraviewFluidFile << "   </PointData>\n";
    paraviewFluidFile << "   <CellData>\n";
    paraviewFluidFile << "   </CellData>\n";
    paraviewFluidFile << "  </Piece>\n";
    paraviewFluidFile << " </ImageData>\n";
    paraviewFluidFile << " <AppendedData encoding=\"raw\">\n  _";
    /**
     * Based on the sparse documentation at https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
     * and alot of testing.
     * Inside <AppendedData> the binary dump must be preceded by an underscore (_)
     * Each DataArray's binary dump must be preceded by it's length.
     * The length should be exported as the integer type specified as header_type in the opening <VTKFile> tag
     * The offset specified in the above <DataArray> tag refers to the offset from the start of the whole binary dump to the start of the length
     */
    // Allocate a buffer equal to size of the largest data array
    // Allocate once rather than allocating and freeing per export
    static char *const t_buffer = static_cast<char*>(malloc(lb.totPossibleNodes * 3 * sizeof(double)));
    static unsigned char *const uc_buffer = reinterpret_cast<unsigned char*>(t_buffer);
    static double *const d_buffer = reinterpret_cast<double*>(t_buffer);
    static tVect *const v_buffer = reinterpret_cast<tVect*>(t_buffer);
    // Type
    offset = lb.totPossibleNodes * sizeof(unsigned char);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    std::fill(uc_buffer, uc_buffer + lb.totPossibleNodes, 2);
    for (const auto &[key, node] : lb.nodes) {
        uc_buffer[key] = node.isInsideParticle() ? 1 : node.type;
    }
    paraviewFluidFile.write(t_buffer, offset);
    // Velocity
    offset = lb.totPossibleNodes * 3 * sizeof(double);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    memset(v_buffer, 0, sizeof(tVect) * lb.totPossibleNodes);
    for (const auto &[key, node] : lb.nodes) {
        v_buffer[key] = node.u * lb.unit.Speed;
    }
    paraviewFluidFile.write(t_buffer, offset);
    // Pressure
    offset = lb.totPossibleNodes * sizeof(double);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    memset(d_buffer, 0, sizeof(double) * lb.totPossibleNodes);
    const double THIRD_PRESSURE = 0.3333333 * lb.unit.Pressure;
    for (const auto &[key, node] : lb.nodes) {
        d_buffer[key] = (node.n - 1.0) * THIRD_PRESSURE;
    }
    paraviewFluidFile.write(t_buffer, offset);
    // Dynamic Viscosity
    if (lb.fluidMaterial.rheologyModel != NEWTONIAN || lb.fluidMaterial.turbulenceOn) {
        offset = lb.totPossibleNodes * sizeof(double);
        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        memset(d_buffer, 0, sizeof(double) * lb.totPossibleNodes);
        for (const auto &[key, node] : lb.nodes) {
            d_buffer[key] = node.visc * lb.unit.DynVisc;
        }
        paraviewFluidFile.write(t_buffer, offset);
    }
    // Friction
    if (lb.fluidMaterial.rheologyModel == MUI || lb.fluidMaterial.rheologyModel == FRICTIONAL || lb.fluidMaterial.rheologyModel == VOELLMY) {
        offset = lb.totPossibleNodes * sizeof(double);
        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        memset(d_buffer, 0, sizeof(double) * lb.totPossibleNodes);
        for (const auto &[key, node] : lb.nodes) {
            d_buffer[key] = node.friction;
        }
        paraviewFluidFile.write(t_buffer, offset);
    }
    // AMass
    if (lb.freeSurface) {
        offset = lb.totPossibleNodes * sizeof(double);
        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        memset(d_buffer, 0, sizeof(double) * lb.totPossibleNodes);
        for (const auto &[key, node] : lb.nodes) {
            d_buffer[key] = node.mass * lb.unit.Density;
        }
        paraviewFluidFile.write(t_buffer, offset);
    }
    // Top Surface
    if (lb.lbTopography) {
        offset = lb.totPossibleNodes * sizeof(double);
        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        memset(d_buffer, 0, sizeof(double) * lb.totPossibleNodes);
        for (const auto &[key, node] : lb.nodes) {
            if (node.isTopography()) {
                d_buffer[key] = 1.0;
            }
        }
        paraviewFluidFile.write(t_buffer, offset);
    }
    paraviewFluidFile << "</AppendedData>";
    paraviewFluidFile << "</VTKFile>\n";
    // data file closing
    paraviewFluidFile.close();
}

// print export stuff

void IO::exportMaxSpeedFluid(const LB& lb) {

    static const double soundSpeed = 1.0 / sqrt(3);
    // fluid max velocity
    double maxFluidSpeed = 0.0;
    for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
        const node* nodeHere = *it;
        maxFluidSpeed = std::max(maxFluidSpeed, nodeHere->u.norm2());
    }
    maxFluidSpeed = sqrt(maxFluidSpeed);
    cout << "MaxFSpeed= " << std::scientific << std::setprecision(2) << maxFluidSpeed * lb.unit.Speed << "(Ma=" << std::scientific << std::setprecision(2) << maxFluidSpeed / soundSpeed << ") ";
    exportFile << "MaxFSpeed= " << std::scientific << std::setprecision(2) << maxFluidSpeed * lb.unit.Speed << "(Ma=" << std::scientific << std::setprecision(2) << maxFluidSpeed / soundSpeed << ") ";

    // printing max speed
    maxFluidSpeedFile.open(maxFluidSpeedFileName.c_str(), ios::app);
    maxFluidSpeedFile << realTime << " " << maxFluidSpeed * lb.unit.Speed << "\n";
    maxFluidSpeedFile.close();
}

void IO::exportFreeSurfaceExtent(const LB& lb) {

    // fluid max velocity
    unsigned int maxX = 0;
    unsigned int maxX_Y = 0;
    unsigned int maxX_Z = 0;
    unsigned int minX = UINT_MAX;
    unsigned int minX_Y = 0;
    unsigned int minX_Z = 0;
    //
    unsigned int maxY = 0;
    unsigned int maxY_Z = 0;
    unsigned int maxY_X = 0;
    unsigned int minY = UINT_MAX;
    unsigned int minY_Z = 0;
    unsigned int minY_X = 0;
    //
    unsigned int maxZ = 0;
    unsigned int maxZ_X = 0;
    unsigned int maxZ_Y = 0;
    unsigned int minZ = UINT_MAX;
    unsigned int minZ_X = 0;
    unsigned int minZ_Y = 0;
    for (nodeList::const_iterator it = lb.interfaceNodes.begin(); it != lb.interfaceNodes.end(); ++it) {
        const node* nodeHere = *it;
        const unsigned int index = nodeHere->coord;
        
        const double xHere = lb.getPositionX(index);
        const double yHere = lb.getPositionY(index);
        const double zHere = lb.getPositionZ(index);

        // max
        if(xHere>maxX) {
            maxX=xHere;
            maxX_Y=yHere;
            maxX_Z=zHere;
        }
        if(yHere>maxY) {
            maxY=yHere;
            maxY_Z=zHere;
            maxY_X=xHere;
        }
        if(zHere>maxZ) {
            maxZ=zHere;
            maxZ_X=xHere;
            maxZ_Y=yHere;
        }
        
        //min
        if(xHere<minX) {
            minX=xHere;
            minX_Y=yHere;
            minX_Z=zHere;
        }
        if(yHere<minY) {
            minY=yHere;
            minY_Z=zHere;
            minY_X=xHere;
        }
        if(zHere<minZ) {
            minZ=zHere;
            minZ_X=xHere;
            minZ_Y=yHere;
        }
                
    }

    // printing max speed
    freeSurfaceExtentFile.open(freeSurfaceExtentFileName.c_str(), ios::app);
    freeSurfaceExtentFile << realTime << " " << maxX * lb.unit.Length << " " << maxX_Y * lb.unit.Length << " " << maxX_Z * lb.unit.Length
                                      << " " << minX * lb.unit.Length << " " << minX_Y * lb.unit.Length << " " << minX_Z * lb.unit.Length
                                      << " " << maxY * lb.unit.Length << " " << maxY_Z * lb.unit.Length << " " << maxY_X * lb.unit.Length
                                      << " " << minY * lb.unit.Length << " " << minY_Z * lb.unit.Length << " " << minY_X * lb.unit.Length
                                      << " " << maxZ * lb.unit.Length << " " << maxZ_X * lb.unit.Length << " " << maxZ_Y * lb.unit.Length
                                      << " " << minZ * lb.unit.Length << " " << minZ_X * lb.unit.Length << " " << minZ_Y * lb.unit.Length<< "\n";
    freeSurfaceExtentFile.close();
}

void IO::exportFluidFlowRate(const LB& lb) {
    // fluid flow rate
    tVect flowRate(0.0, 0.0, 0.0);
    for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
        const node* nodeHere = *it;
        if (!nodeHere->isInsideParticle()) {
            flowRate += nodeHere->u * nodeHere->mass;
        }
    }

    const double flowRateX = flowRate.dot(Xp) / double(lb.lbSize[0] - 2);
    const double flowRateY = flowRate.dot(Yp) / double(lb.lbSize[1] - 2);
    const double flowRateZ = flowRate.dot(Zp) / double(lb.lbSize[2] - 2);

    // printing rate
    fluidFlowRateFile.open(fluidFlowRateFileName.c_str(), ios::app);
    fluidFlowRateFile << realTime << " " << flowRateX * lb.unit.FlowRate << " " << flowRateY * lb.unit.FlowRate << " " << flowRateZ * lb.unit.FlowRate << "\n";
    fluidFlowRateFile.close();
}

void IO::exportParticleFlowRate(const DEM& dem) {
    // calculating speed for real time check
    // particle max velocity
    tVect flowRate(0.0, 0.0, 0.0);
    double totMass = 0.0;
    for (int n = 0; n < dem.elmts.size(); ++n) {
        flowRate += dem.elmts[n].m * dem.elmts[n].x1;
        //totMass += dem.elmts[n].m;
    }
    //if (totMass > 0.0) {
    //flowRate *= dem.sphereMat.density / totMass;
    //}


    const double flowRateX = flowRate.dot(Xp) / dem.demSize[0];
    const double flowRateY = flowRate.dot(Yp) / dem.demSize[1];
    const double flowRateZ = flowRate.dot(Zp) / dem.demSize[2];

    // printing rate
    particleFlowRateFile.open(particleFlowRateFileName.c_str(), ios::app);
    particleFlowRateFile << realTime << " " << flowRateX << " " << flowRateY << " " << flowRateZ << "\n";
    particleFlowRateFile.close();


}

void IO::exportMaxSpeedParticles(const DEM& dem) {
    // calculating speed for real time check
    // particle max velocity
    double maxPartSpeed = 0.0;
    for (int n = 0; n < dem.elmts.size(); ++n) {
        maxPartSpeed = std::max(maxPartSpeed, dem.elmts[n].x1.norm2());
    }
    maxPartSpeed = sqrt(maxPartSpeed);
    cout << "MaxPSpeed= " << std::scientific << std::setprecision(2) << maxPartSpeed << " ";
    exportFile << "MaxPSpeed= " << std::scientific << std::setprecision(2) << maxPartSpeed << " ";

    // particle max rotational velocity
    double maxPartSpin = 0.0;
    for (int n = 0; n < dem.elmts.size(); ++n) {
        maxPartSpin = std::max(maxPartSpin, dem.elmts[n].wGlobal.norm2());
    }
    maxPartSpin = sqrt(maxPartSpin);
    cout << "MaxPSpin= " << std::scientific << std::setprecision(2) << maxPartSpin << " ";
    exportFile << "MaxPSpin= " << std::scientific << std::setprecision(2) << maxPartSpin << " ";

    // printing max speed
    maxParticleSpeedFile.open(maxParticleSpeedFileName.c_str(), ios::app);
    maxParticleSpeedFile << realTime << " " << maxPartSpeed << " " << maxPartSpin << "\n";
    maxParticleSpeedFile.close();

}

void IO::exportParticleCenterOfMass(const DEM& dem) {

    // particle center of mass
    const tVect particleCenter = particleCenterOfMass(dem.elmts);

    // printing particle center of mass
    particleCenterOfMassFile.open(particleCenterOfMassFileName.c_str(), ios::app);
    particleCenterOfMassFile << realTime << " " << particleCenter.dot(Xp) << " " << particleCenter.dot(Yp) << " " << particleCenter.dot(Zp) << "\n";
    particleCenterOfMassFile.close();

}

void IO::exportParticleCoordination(const DEM& dem) {

    // particle coordination number and force intensity
    tVect meanSolidIntensity(0.0, 0.0, 0.0);
    tVect meanFluidIntensity(0.0, 0.0, 0.0);
    const unsigned int maxCoordination = 16;
    unsIntList totCoordination;
    totCoordination.resize(maxCoordination);
    for (int c = 0; c < maxCoordination; ++c) {
        totCoordination[c] = 0;
    }
    for (int a = 0; a < dem.activeElmts.size(); ++a) {
        unsigned int p = dem.activeElmts[a];
        const tVect solidIntensityHere = dem.elmts[p].solidIntensity;
        const tVect fluidIntensityHere = dem.elmts[p].FHydro.abs();
        meanSolidIntensity += solidIntensityHere;
        meanFluidIntensity += fluidIntensityHere;
        unsigned int coordinationHere = dem.elmts[p].coordination;
        if (coordinationHere > maxCoordination - 1) {
            coordinationHere = maxCoordination - 1;
        }
        totCoordination[coordinationHere] += 1;
    }
    meanSolidIntensity /= dem.elmts.size();
    meanFluidIntensity /= dem.elmts.size();
    // printing particle center of mass
    particleCoordinationFile.open(particleCoordinationFileName.c_str(), ios::app);
    particleCoordinationFile << realTime << " ";
    particleCoordinationFile << meanSolidIntensity.dot(Xp) << " " << meanSolidIntensity.dot(Yp) << " " << meanSolidIntensity.dot(Zp) << " ";
    particleCoordinationFile << meanFluidIntensity.dot(Xp) << " " << meanFluidIntensity.dot(Yp) << " " << meanFluidIntensity.dot(Zp);
    for (int c = 0; c < maxCoordination; ++c) {
        particleCoordinationFile << " " << totCoordination[c];
    }
    particleCoordinationFile << endl;
    particleCoordinationFile.close();



}

void IO::exportFluidCenterOfMass(const LB& lb) {

    // particle center of mass
    const tVect fluidCenter = fluidCenterOfMass(lb) * lb.unit.Length;

    // printing particle center of mass
    fluidCenterOfMassFile.open(fluidCenterOfMassFileName.c_str(), ios::app);
    fluidCenterOfMassFile << realTime << " " << fluidCenter.dot(Xp) << " " << fluidCenter.dot(Yp) << " " << fluidCenter.dot(Zp) << "\n";
    fluidCenterOfMassFile.close();

}

void IO::exportParticleOverlap(DEM& dem) { // used to be passed as const
    // checking the mean and max overlap to see if system is stiff
    // prints the total mass in the free fluid domain
    double meanOverlap(0.0), meanRelOverlap(0.0);
    double maxOverlap(0.0), maxRelOverlap(0.0);
    double meanDtOverlap(0.0), meanRelDtOverlap(0.0);
    double maxDtOverlap(0.0), maxRelDtOverlap(0.0);
    int totContactParticles(0), totDtContactParticles(0);

    for (int p = 0; p < dem.elmts.size(); ++p) {
        const double radiusHere = dem.elmts[p].radius;
        // istantaneous overlap
        const double overlapHere = dem.elmts[p].maxOverlap;
        if (overlapHere > 0.0) {
            ++totContactParticles;

            const double relOverlapHere = overlapHere / radiusHere;
            if (overlapHere > maxOverlap) {
                maxOverlap = overlapHere;
            }
            if (relOverlapHere > maxRelOverlap) {
                maxRelOverlap = relOverlapHere;
            }
            meanOverlap += overlapHere;
            meanRelOverlap += relOverlapHere;
        }
        // max overlap between output steps
        const double dtOverlapHere = dem.elmts[p].maxDtOverlap;
        if (dtOverlapHere > 0.0) {
            ++totDtContactParticles;
            const double relDtOverlapHere = dtOverlapHere / radiusHere;
            if (dtOverlapHere > maxDtOverlap) {
                maxDtOverlap = dtOverlapHere;
            }
            if (relDtOverlapHere > maxRelDtOverlap) {
                maxRelDtOverlap = relDtOverlapHere;
            }
            meanDtOverlap += dtOverlapHere;
            meanRelDtOverlap += relDtOverlapHere;
        }
        // reset also maximum overlap
        dem.elmts[p].maxDtOverlap = 0.0;
    }

    if (totContactParticles > 0) {
        meanOverlap /= double(totContactParticles);
        meanRelOverlap /= double(totContactParticles);
    }

    if (totDtContactParticles > 0) {
        meanDtOverlap /= double(totDtContactParticles);
        meanRelDtOverlap /= double(totDtContactParticles);
    }


    // printing instantaneous overlap
    overlapFile.open(overlapFileName.c_str(), ios::app);
    overlapFile << realTime << " " << meanOverlap << " " << maxOverlap << " " << meanRelOverlap << " " << maxRelOverlap << "\n";
    overlapFile.close();
    // printing max overlap between output steps
    dtOverlapFile.open(dtOverlapFileName.c_str(), ios::app);
    dtOverlapFile << realTime << " " << meanDtOverlap << " " << maxDtOverlap << " " << meanRelDtOverlap << " " << maxRelDtOverlap << "\n";
    dtOverlapFile.close();

}

void IO::exportWallForce(const DEM& dem) {
    // printing forces acting on all walls
    const int NDIGITS = 6;
    // this specific time step
    wallForceFile.open(wallForceFileName.c_str(), ios::app);
    wallForceFile << std::scientific << std::setprecision(NDIGITS) << realTime;
    wallForceFile << " " << dem.walls.size();
    for (int w = 0; w < dem.walls.size(); ++w) {
        const tVect partF = dem.walls[w].FParticle;
        const tVect hydroF = dem.walls[w].FHydro;
        wallForceFile << " " 
                << std::scientific << std::setprecision(NDIGITS) << partF.dot(Xp) << " " 
                << std::scientific << std::setprecision(NDIGITS) << partF.dot(Yp) << " " 
                << std::scientific << std::setprecision(NDIGITS) << partF.dot(Zp) << " " 
                << std::scientific << std::setprecision(NDIGITS) << hydroF.dot(Xp) << " " 
                << std::scientific << std::setprecision(NDIGITS) << hydroF.dot(Yp) << " " 
                << std::scientific << std::setprecision(NDIGITS) << hydroF.dot(Zp);
    }
    wallForceFile << endl;
    wallForceFile.close();

    // the maximum since last output
    maxWallForceFile.open(maxWallForceFileName.c_str(), ios::app);
    maxWallForceFile << std::scientific << std::setprecision(NDIGITS) << realTime;
    wallForceFile << " " << dem.walls.size();
    for (int w = 0; w < dem.walls.size(); ++w) {
        const tVect partF = dem.walls[w].maxFParticle;
        const tVect hydroF = dem.walls[w].maxFHydro;
        maxWallForceFile << " " 
                << std::scientific << std::setprecision(NDIGITS) << partF.dot(Xp) << " " 
                << std::scientific << std::setprecision(NDIGITS) << partF.dot(Yp) << " " 
                << std::scientific << std::setprecision(NDIGITS) << partF.dot(Zp) << " " 
                << std::scientific << std::setprecision(NDIGITS) << hydroF.dot(Xp) << " " 
                << std::scientific << std::setprecision(NDIGITS) << hydroF.dot(Yp) << " " 
                << std::scientific << std::setprecision(NDIGITS) << hydroF.dot(Zp);
    }
    maxWallForceFile << endl;
    maxWallForceFile.close();

}

void IO::exportCylinderForces(const DEM& dem){
    // printing forces acting on all cylinders
    const int NDIGITS = 6;

    // this specific time step
    cylinderForceFile.open(cylinderForceFileName.c_str(), ios::app);
    cylinderForceFile << std::scientific << std::setprecision(NDIGITS) << realTime;
    cylinderForceFile << " " << dem.cylinders.size();
    for (int c = 0; c < dem.cylinders.size(); ++c) {
        const tVect partF = dem.cylinders[c].FParticle;
        const tVect hydroF = dem.cylinders[c].FHydro;
        cylinderForceFile << " "
                 << std::scientific << std::setprecision(NDIGITS) << partF.dot(Xp) << " "
                 << std::scientific << std::setprecision(NDIGITS) << partF.dot(Yp) << " " 
                 << std::scientific << std::setprecision(NDIGITS) << partF.dot(Zp) << " " 
                 << std::scientific << std::setprecision(NDIGITS) << hydroF.dot(Xp) << " " 
                 << std::scientific << std::setprecision(NDIGITS) << hydroF.dot(Yp) << " " 
                 << std::scientific << std::setprecision(NDIGITS) << hydroF.dot(Zp);
    }
    cylinderForceFile << endl;
    cylinderForceFile.close();
}

void IO::exportFluidMass(const LB& lb) {
    // total fluid mass
    double massTot = totFluidMass(lb);
    cout << "Volume=" << std::scientific << std::setprecision(2) << massTot * lb.unit.Volume << "; Mass = " << std::scientific << std::setprecision(2) << massTot * lb.unit.Mass << " ";
    exportFile << "Volume=" << std::scientific << std::setprecision(2) << massTot * lb.unit.Volume << "; Mass = " << std::scientific << std::setprecision(2) << massTot * lb.unit.Mass << " ";

    // printing fluid mass
    fluidMassFile.open(fluidMassFileName.c_str(), ios::app);
    fluidMassFile << realTime << " " << massTot * lb.unit.Mass << "\n";
    fluidMassFile.close();
}

void IO::exportPlasticity(const LB& lb) {
    // fluid plasticity state
    const double percPlastic = totPlastic(lb);
    cout << "Plastic =" << int(percPlastic) << "% ";
    exportFile << "Plastic =" << int(percPlastic) << "% ";
    // printing plasticity level
    plasticityFile.open(plasticityFileName.c_str(), ios::app);
    plasticityFile << realTime << " " << percPlastic << "\n";
    plasticityFile.close();

}

void IO::exportMeanViscosity(const LB& lb) {
    // fluid plasticity state
    const double meanVisc = meanViscosity(lb);
    cout << "MeanVisc =" << std::scientific << std::setprecision(2) << meanVisc * lb.unit.DynVisc << " ";
    exportFile << "MeanVisc =" << std::scientific << std::setprecision(2) << meanVisc * lb.unit.DynVisc << " ";
}

void IO::exportShearCell(const LB& lb, const DEM& dem) {
    // apparent viscosity from shear cell
    double appVisc = 0.0;
    double externalShear = 0.0;
    double wallStress = 0.0;
    apparentViscosity(lb, dem.walls, externalShear, wallStress, appVisc);
    cout << "App Visc= " << std::scientific << std::setprecision(2) << appVisc << " ";
    exportFile << "App Visc= " << std::scientific << std::setprecision(2) << appVisc << " ";

    tVect xDirec = tVect(1.0, 0.0, 0.0);
    cout << "wallDown = " << std::scientific << std::setprecision(2) << dem.walls[0].FHydro.dot(xDirec) << " wallUp = " << std::scientific << std::setprecision(2) << dem.walls[1].FHydro.dot(xDirec) << " ";
    cout << "wallDown = " << std::scientific << std::setprecision(2) << dem.walls[0].FParticle.dot(xDirec) << " wallUp = " << std::scientific << std::setprecision(2) << dem.walls[1].FParticle.dot(xDirec) << " ";
}

void IO::exportEnergy(const DEM& dem, const LB& lb) {

    if (dem.elmts.size()) {
        cout << "Energy (DEM): ";
        exportFile << "Energy (DEM): ";
        cout << "eKin = " << std::scientific << std::setprecision(2) << dem.particleEnergy.kin << " ";
        exportFile << "eKin = " << std::scientific << std::setprecision(2) << dem.particleEnergy.kin << " ";
        cout << "eGrav = " << std::scientific << std::setprecision(2) << dem.particleEnergy.grav << " ";
        exportFile << "eGrav = " << std::scientific << std::setprecision(2) << dem.particleEnergy.grav << " ";
        cout << "eTot = " << std::scientific << std::setprecision(2) << dem.particleEnergy.total << " ";
        exportFile << "eTot = " << std::scientific << std::setprecision(2) << dem.particleEnergy.total << " ";
    }
    if (lbmSolver) {
        cout << "Energy (LBM): ";
        cout << "eKin = " << std::scientific << std::setprecision(2) << lb.fluidEnergy.kin + lb.fluidImmersedEnergy.kin << " ";
        exportFile << "eKin = " << std::scientific << std::setprecision(2) << lb.fluidEnergy.kin + lb.fluidImmersedEnergy.kin << " ";
        cout << "eGrav = " << std::scientific << std::setprecision(2) << lb.fluidEnergy.grav + lb.fluidImmersedEnergy.grav << " ";
        exportFile << "eGrav = " << std::scientific << std::setprecision(2) << lb.fluidEnergy.grav + lb.fluidImmersedEnergy.grav << " ";
        cout << "eTot = " << std::scientific << std::setprecision(2) << lb.fluidEnergy.total + lb.fluidImmersedEnergy.total << " ";
        exportFile << "eTot = " << std::scientific << std::setprecision(2) << lb.fluidEnergy.total + lb.fluidImmersedEnergy.total << " ";
    }

    ofstream energyFile;
    energyFile.open(energyFileName.c_str(), ios::app);
    // Set energyFile header
    energyFile << std::scientific << std::setprecision(6) << realTime << " ";
    energyFile << std::scientific << std::setprecision(10) << dem.particleEnergy.mass << " " << dem.particleEnergy.trKin << " " << dem.particleEnergy.rotKin << " " << dem.particleEnergy.grav << " ";
    energyFile << std::scientific << std::setprecision(10) << lb.fluidEnergy.mass * lb.unit.Mass << " " << lb.fluidEnergy.trKin * lb.unit.Energy << " " << lb.fluidEnergy.grav * lb.unit.Energy << " ";
    energyFile << std::scientific << std::setprecision(10) << lb.fluidImmersedEnergy.mass * lb.unit.Mass << " " << lb.fluidImmersedEnergy.trKin * lb.unit.Energy << " " << lb.fluidImmersedEnergy.grav * lb.unit.Energy << endl;
    energyFile.close();

}

void IO::exportForces(const DEM& dem) {
    // printing force
    // hydraulic and collision forces
    forceFile.open(forceFileName.c_str(), ios::app);
    const double collisionForce = collisionForceTot(dem.elmts);
    const double hydraulicForce = hydraulicForceTot(dem.elmts);
    forceFile << realTime << " " << collisionForce << " " << hydraulicForce << "\n";
    forceFile.close();
}

void IO::exportForceObstacle(const objectList& objects) {
    // printing total force on obstacle
    ofstream obstacleFile;
    obstacleFile.open(obstacleFileName.c_str(), ios::app);
    tVect totObstacleForce = totForceObject(objects, 0, objects.size() - 1);
    if (realTime <=0) {
        obstacleFile << "t" << "," << "fx" << "," << "fy" << "," << "fz" << endl;
        obstacleFile << realTime << "," <<
        totObstacleForce.x << "," << totObstacleForce.y << "," << totObstacleForce.z << endl;
    } else {
        obstacleFile << realTime << "," <<
        totObstacleForce.x << "," << totObstacleForce.y << "," << totObstacleForce.z << endl;
    }
//    obstacleFile << realTime << " ";
//    totObstacleForce.printLine(obstacleFile);
    obstacleFile.close();
}


// data elaboration

double IO::totPlastic(const LB& lb) const {
    // prints the total mass in the free fluid domain
    unsigned int totPlastic = 0;
    unsigned int totActive = 0;

    for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
        const node* nodeHere = *it;
        ++totActive;
        if (nodeHere->visc > 0.95 * lb.fluidMaterial.lbMaxVisc) {
            ++totPlastic;
        }
    }
    const double perc = 100.0 * double(totPlastic) / double(totActive);
    return perc;
}

double IO::totFluidMass(const LB& lb) const {
    // returns the total mass in the free fluid domain
    double mass = 0.0;

    if (lbmSolver) {
        for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
            const node* nodeHere = *it;

            if (!nodeHere->isInsideParticle()) {
                mass += nodeHere->mass;
            }
        }
    }
    return mass;
}

double IO::meanViscosity(const LB& lb) const {
    // prints the total mass in the free fluid domain
    double meanVisc = 0.0;
    unsigned int counter = 0;

    if (lbmSolver) {
        for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
            const node* nodeHere = *it;

            ++counter;
            meanVisc += nodeHere->visc;
        }
        meanVisc = meanVisc / double(counter);
        return meanVisc;
    }
}

double IO::totParticleMass(const elmtList& elmts) const {
    // prints the total mass in the free fluid domain
    double mass = 0.0;

    for (int n = 0; n < elmts.size(); ++n) {
        // calculate mass
        mass += elmts[n].m;
    }

    return mass;
}

void IO::apparentViscosity(const LB& lb, const wallList& walls, double& externalShear, double& wallStress, double& appVisc) const {
    // computes the apparent viscosity of a shear cell, as the ratio between shear at the wall and imposed shear rate
    // the cell is sheared in direction x and has moving wall in z
    appVisc = 0.0;
    tVect xDirec = tVect(1.0, 0.0, 0.0);

    double cellSize = double(lb.lbSize[2] - 2) * lb.unit.Length;

    externalShear = 0.5 * (walls[1].vel.dot(xDirec) - walls[0].vel.dot(xDirec)) / cellSize;

    double plateSize = double(lb.lbSize[0] - 2) * double(lb.lbSize[1] - 2) * lb.unit.Length * lb.unit.Length;

    wallStress = -0.5 * (walls[1].FHydro.dot(xDirec) + walls[1].FParticle.dot(xDirec) - walls[0].FHydro.dot(xDirec) - walls[0].FParticle.dot(xDirec)) / plateSize;

    appVisc = 0.5 * wallStress / externalShear;
}

tVect IO::fluidCenterOfMass(const LB& lb) const {
    // prints the total mass in the free fluid domain
    tVect center(0.0, 0.0, 0.0);
    double totMass(0.0);

    if (lbmSolver) {
        for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
            const node* nodeHere = *it;
            const unsigned int index = nodeHere->coord;
            if (!nodeHere->isInsideParticle()) {
                center += nodeHere->mass * lb.getPosition(index);
                totMass += nodeHere->mass;
            }
        }
        return center / totMass;
    } else {
        return Zero;
    }
}

tVect IO::particleCenterOfMass(const elmtList& elmts) const {
    // prints the total mass in the free fluid domain
    tVect center(0.0, 0.0, 0.0);
    double totMass(0.0);

    for (int p = 0; p < elmts.size(); ++p) {
        center += elmts[p].m * elmts[p].x0;
        totMass += elmts[p].m;
    }

    return center / totMass;
}

double IO::hydraulicForceTot(const elmtList& elmts) const {
    double totForce = 0.0;

    for (int n = 0; n < elmts.size(); ++n) {
        totForce += elmts[n].FHydro.norm();
    }

    return totForce;
}

double IO::collisionForceTot(const elmtList& elmts) const {
    double totForce = 0.0;

    for (int n = 0; n < elmts.size(); ++n) {
        totForce += elmts[n].FParticle.norm();
        totForce += elmts[n].FWall.norm();
    }

    return totForce;
}

tVect IO::totForceObject(const objectList& objects, const unsigned int& groupBegin, const unsigned int& groupEnd) const {
    // prints the total mass in the free fluid domain
    tVect totForce(0.0, 0.0, 0.0);

    for (int o = 0; o < objects.size(); ++o) {
        if (o >= groupBegin && o <= groupEnd) {
            totForce += objects[o].FParticle + objects[o].FHydro;
        }
    }

    return totForce;
}

// OPEN BARRIER (MADDALENA)

void IO::exportForceObstacleElement(const DEM& dem) {
    unsigned int totObstacles = 0;
    for (int o = 0; o < dem.objects.size(); ++o) {
        if (dem.objects[o].ElID + 1 > totObstacles) {
            totObstacles = dem.objects[o].ElID + 1;
        }
    }
    // printing total force on obstacle
    ofstream obstacleFile;
    obstacleFile.open(obstacleFileName.c_str(), ios::app);
    obstacleFile << realTime << " ";
    for (int i = 0; i < totObstacles; i++) {
        tVect totObstacleForce = totForceObjectElement(dem.objects, i);
        totObstacleForce.print(obstacleFile);
    }
    obstacleFile << "\n";
    obstacleFile.close();
}

void IO::exportMomentObstacleElement(const DEM& dem) {
    unsigned int totObstacles = 0;
    for (int o = 0; o < dem.objects.size(); ++o) {
        if (dem.objects[o].ElID + 1 > totObstacles) {
            totObstacles = dem.objects[o].ElID + 1;
        }
    }
    // printing total force on obstacle
    ofstream obstacleMomentFile;
    obstacleMomentFile.open(obstacleMomentFileName.c_str(), ios::app);
    obstacleMomentFile << realTime << " ";
    for (int i = 0; i < totObstacles; i++) {
        tVect totObstacleMoment = totMomentObjectElement(dem.walls, dem.objects, i);
        totObstacleMoment.print(obstacleMomentFile);
    }
    obstacleMomentFile << "\n";
    obstacleMomentFile.close();
}

void IO::exportHPartObstacle(const DEM& dem) {

    ofstream obstacleHeightFile;
    obstacleHeightFile.open(obstacleHeightFileName.c_str(), ios::app);
    obstacleHeightFile << realTime << " ";

    double radMax = 0.0;
    for (int n = 0; n < dem.elmts.size(); ++n) {
        if (radMax < dem.elmts[n].radius) {
            radMax = dem.elmts[n].radius;
        }
    }
    double xObMax = 0.0;
    for (int o = 0; o < dem.objects.size(); ++o) {
        if (xObMax < dem.objects[o].x0.dot(Xp)) {
            xObMax = dem.objects[o].x0.dot(Xp);
        }
    }


    double xObMin = xObMax;
    for (int o = 0; o < dem.objects.size(); ++o) {
        if (xObMin > dem.objects[o].x0.dot(Xp)) {
            xObMin = dem.objects[o].x0.dot(Xp);
        }
    }
    xObMin = xObMin - 4.0 * radMax;

    double zObMax = 0.0;
    for (int o = 0; o < dem.objects.size(); ++o) {
        if (zObMax < dem.objects[o].x0.dot(Zp)) {
            zObMax = dem.objects[o].x0.dot(Zp);
        }
    }

    unsigned int totObstacles = 0;
    for (int o = 0; o < dem.objects.size(); ++o) {
        if (dem.objects[o].ElID + 1 > totObstacles) {
            totObstacles = dem.objects[o].ElID + 1;
        }
    }

    for (int i = 0; i < totObstacles - 1; i++) {

        double yObMin = 0.0;
        double yObMax = dem.walls[3].p.dot(Yp);
        for (int o = 0; o < dem.objects.size(); ++o) {
            if (dem.objects[o].ElID == i) {
                if (yObMin < dem.objects[o].x0.dot(Yp))
                    yObMin = dem.objects[o].x0.dot(Yp);
            } else if (dem.objects[o].ElID == i + 1) {
                if (yObMax > dem.objects[o].x0.dot(Yp))
                    yObMax = dem.objects[o].x0.dot(Yp);
            }

        }

        unsigned int numPart = 0;
        double maxPartHeight = 0.0;
        doubleList heightPart;
        heightPart.clear();
        //unsigned int indexvis=0.0;
        for (int n = 0; n < dem.elmts.size(); ++n) {
            if ((dem.elmts[n].x0.dot(Xp) < xObMax) && (dem.elmts[n].x0.dot(Xp) > xObMin)) {
                if ((dem.elmts[n].x0.dot(Yp) < yObMax) && (dem.elmts[n].x0.dot(Yp) > yObMin)) {
                    numPart += 1;
                    heightPart.push_back(dem.elmts[n].x0.dot(Zp));
                    if (maxPartHeight < dem.elmts[n].x0.dot(Zp)) {
                        maxPartHeight = dem.elmts[n].x0.dot(Zp);
                    }

                }

            }

        }
        obstacleHeightFile << maxPartHeight << " ";
        std::sort(heightPart.rbegin(), heightPart.rend());

        double meanThree = 0.0;
        double meanTen = 0.0;
        double meanHeight = 0.0;

        if (heightPart.size() >= 3) {
            for (int n = 0; n < 3; ++n) {
                meanThree += heightPart[n];
            }
        }

        meanThree = meanThree / 3.0;
        if (heightPart.size() >= 10) {
            for (int n = 0; n < 10; ++n) {
                meanTen += heightPart[n];
            }
        }
        meanTen = meanTen / 10.0;

        for (int n = 0; n < heightPart.size(); ++n) {
            meanHeight += heightPart[n];
        }
        meanHeight = meanHeight / numPart;

        obstacleHeightFile << meanThree << " ";
        obstacleHeightFile << meanTen << " ";
        obstacleHeightFile << meanHeight << " ";
        obstacleHeightFile << numPart << " ";
    }

    obstacleHeightFile << "\n";
    obstacleHeightFile.close();
}

void IO::exportKmPartObstacle(const DEM& dem) {
    ofstream obstacleKmFile;
    obstacleKmFile.open(obstacleKmFileName.c_str(), ios::app);
    obstacleKmFile << realTime << " ";

    double radMax = 0.0;
    for (int n = 0; n < dem.elmts.size(); ++n) {
        if (radMax < dem.elmts[n].radius) {
            radMax = dem.elmts[n].radius;
        }
    }

    double xObMax = 0.0;
    for (int o = 0; o < dem.objects.size(); ++o) {
        if (xObMax < dem.objects[o].x0.dot(Xp)) {
            xObMax = dem.objects[o].x0.dot(Xp);
        }
    }


    double xObMin = xObMax;
    for (int o = 0; o < dem.objects.size(); ++o) {
        if (xObMin > dem.objects[o].x0.dot(Xp)) {
            xObMin = dem.objects[o].x0.dot(Xp);
        }
    }
    xObMin = xObMin - 4.0 * radMax;



    double zObMax = 0.0;
    for (int o = 0; o < dem.objects.size(); ++o) {
        if (zObMax < dem.objects[o].x0.dot(Zp)) {
            zObMax = dem.objects[o].x0.dot(Zp);
        }
    }


    unsigned int totObstacles = 0;
    for (int o = 0; o < dem.objects.size(); ++o) {
        if (dem.objects[o].ElID + 1 > totObstacles) {
            totObstacles = dem.objects[o].ElID + 1;
        }
    }

    for (int i = 0; i < totObstacles - 1; i++) {

        double yObMin = 0.0;
        double yObMax = dem.walls[3].p.dot(Yp);
        for (int o = 0; o < dem.objects.size(); ++o) {
            if (dem.objects[o].ElID == i) {
                if (yObMin < dem.objects[o].x0.dot(Yp))
                    yObMin = dem.objects[o].x0.dot(Yp);
            } else if (dem.objects[o].ElID == i + 1) {
                if (yObMax > dem.objects[o].x0.dot(Yp))
                    yObMax = dem.objects[o].x0.dot(Yp);
            }

        }


        double subdZ = zObMax / 10.0;
        doubleList enKinZ;
        enKinZ.clear();
        double zInt = 0.0;
        double zInf = 0.0;
        double enKinM = 0.0;
        double numPartInt;
        unsIntList numPartList;
        numPartList.clear();


        for (int m = 0; m < 10; ++m) {
            zInf = zInt;
            zInt += subdZ;
            for (int n = 0; n < dem.elmts.size(); ++n) {
                if ((dem.elmts[n].x0.dot(Xp) < xObMax) && (dem.elmts[n].x0.dot(Xp) > xObMin)) {
                    if ((dem.elmts[n].x0.dot(Yp) < yObMax) && (dem.elmts[n].x0.dot(Yp) > yObMin)) {
                        if ((dem.elmts[n].x0.dot(Zp) < zInt) && (dem.elmts[n].x0.dot(Zp) >= zInf)) {
                            numPartInt += 1;
                            const tVect wSquare = dem.elmts[n].wLocal.compProd(dem.elmts[n].wLocal);
                            enKinM += (0.5 * dem.elmts[n].m * dem.elmts[n].x1.norm2() + 0.5 * dem.elmts[n].I.dot(wSquare));
                            //cout<<"enKinM"<<enKinM<<" "<<0.5*dem.elmts[n].m * dem.elmts[n].x1.norm2()<<" "<<0.5*dem.elmts[n].I.dot(wSquare);

                        }
                    }
                }
            }
            enKinM = enKinM / numPartInt;
            enKinZ.push_back(enKinM);
            numPartList.push_back(numPartInt);
            numPartInt = 0;
            enKinM = 0.0;

        }
        for (int n = 0; n < enKinZ.size(); ++n) {
            obstacleKmFile << enKinZ[n] << " ";
        }
        for (int n = 0; n < numPartList.size(); ++n) {
            obstacleKmFile << numPartList[n] << " ";
        }
    }
    obstacleKmFile << "\n";
    obstacleKmFile.close();
}

tVect IO::totForceObjectElement(const objectList& objects, const unsigned int& obstacleIndex) const {
    // prints the total mass in the free fluid domain
    tVect totForce(0.0, 0.0, 0.0);

    for (int o = 0; o < objects.size(); ++o) {
        if (objects[o].ElID == obstacleIndex) {
            totForce += objects[o].FParticle;
        }
    }

    return totForce;
}

tVect IO::totMomentObjectElement(const wallList& walls, const objectList& objects, const unsigned int& obstacleIndex) const {
    // prints the total mass in the free fluid domain
    tVect totMoment(0.0, 0.0, 0.0);

    for (int o = 0; o < objects.size(); ++o) {
        if (objects[o].ElID == obstacleIndex) {
            double distance = walls[4].dist(objects[o].x0);
            totMoment += objects[o].FParticle*distance;
        }
    }

    return totMoment;
}

void IO::exportHongKongFlow(DEM& dem) {
    // creates a file that contains the variables necessary for flow characterization
    // the barrier is at x=1.4 and is 0.01 thick
    static const double wallPosition = 1.4;
    static const double barrierBuffer = dem.meanObjRadius;

    // we create grain-coarsed variables in two intervals: upstream and downstream
    static const double windowSize = 0.1;
    // flow characterization (mean velocity and flow level) -> both close to walls and global
    static const double wallBuffer = 2.5 * dem.meanPartRadius;

    // window limits
    static const double upBegin = wallPosition - barrierBuffer - windowSize;
    static const double upEnd = wallPosition - barrierBuffer;
    static const double downBegin = wallPosition + barrierBuffer;
    static const double downEnd = wallPosition + barrierBuffer + windowSize;
    // lateral limits
    static const double leftBegin = dem.walls[2].p.dot(Yp);
    static const double leftEnd = dem.walls[2].p.dot(Yp) + wallBuffer;
    static const double rightBegin = dem.walls[3].p.dot(Yp) - wallBuffer;
    static const double rightEnd = dem.walls[3].p.dot(Yp);
    //cout<<"(("<<leftBegin<<" "<<leftEnd<<", "<<rightBegin<<" "<<rightEnd<<"))";

    // flow velocity (mean velocity of all particles in the interval)
    double uUp(0.0), uUpLeft(0.0), uUpRight(0.0);
    double uDown(0.0), uDownLeft(0.0), uDownRight(0.0);

    // mass passed at control point
    double flowedMassUp(0.0), flowedMassDown(0.0);
    // empty arays to collect heights
    doubleList insideUp, insideUpLeft, insideUpRight;
    doubleList insideDown, insideDownLeft, insideDownRight;
    insideUp.clear();
    insideUpLeft.clear();
    insideUpRight.clear();
    insideDown.clear();
    insideDownLeft.clear();
    insideDownRight.clear();

    // collect data
    for (int i = 0; i < dem.elmts.size(); ++i) {
        const double elmtX = dem.elmts[i].x0.dot(Xp);
        const double elmtY = dem.elmts[i].x0.dot(Yp);
        const double elmtZ = dem.elmts[i].x0.dot(Zp) + dem.elmts[i].radius;
        const double elmtU = dem.elmts[i].x1.dot(Xp);
        const double elmtMass = dem.elmts[i].m;
        // check if we are inside observation window (upstream)
        if (elmtX > upBegin && elmtX < upEnd) {
            uUp += elmtU;
            insideUp.push_back(elmtZ);
            if (elmtY > leftBegin && elmtY < leftEnd) {
                uUpLeft += elmtU;
                insideUpLeft.push_back(elmtZ);
            }
            if (elmtY > rightBegin && elmtY < rightEnd) {
                uUpRight += elmtU;
                insideUpRight.push_back(elmtZ);
            }
        } else if (elmtX > downBegin && elmtX < downEnd) {
            uDown += elmtU;
            insideDown.push_back(elmtZ);
            if (elmtY > leftBegin && elmtY < leftEnd) {
                uDownLeft += elmtU;
                insideDownLeft.push_back(elmtZ);
            }
            if (elmtY > rightBegin && elmtY < rightEnd) {
                uDownRight += elmtU;
                insideDownRight.push_back(elmtZ);
            }
        }

        if (elmtX > upBegin) {
            flowedMassUp += elmtMass;
        }
        if (elmtX > downBegin) {
            flowedMassDown += elmtMass;
        }

    }

    // compute mean velocities and flow height (mean position of the top 10% particles)
    double hUp(0.0);
    if (insideUp.size() > 0) {
        hUp = getSurface(insideUp);
        uUp /= double(insideUp.size());
    }

    double hUpLeft(0.0);
    if (insideUpLeft.size() > 0) {
        hUpLeft = getSurface(insideUpLeft);
        uUpLeft /= double(insideUpLeft.size());
    }

    double hUpRight(0.0);
    if (insideUpRight.size() > 0) {
        hUpRight = getSurface(insideUpRight);
        uUpRight /= double(insideUpRight.size());
    }

    double hDown(0.0);
    if (insideDown.size() > 0) {
        hDown = getSurface(insideDown);
        uDown /= double(insideDown.size());
    }

    double hDownLeft(0.0);
    if (insideDownLeft.size() > 0) {
        hDownLeft = getSurface(insideDownLeft);
        uDownLeft /= double(insideDownLeft.size());
    }

    double hDownRight(0.0);
    if (insideDownRight.size() > 0) {
        hDownRight = getSurface(insideDownRight);
        uDownRight /= double(insideDownRight.size());
    }

    cout << "[[" << insideUp.size() << " " << insideUpLeft.size() << " " << insideUpRight.size() << " " << insideDown.size() << " " << insideDownLeft.size() << " " << insideDownRight.size() << "]]";

    ofstream hongKongFlowFile;
    hongKongFlowFile.open(hongKongFlowFileName.c_str(), ios::app);
    //hongKongFlowFile << "time percPlastic\n";
    hongKongFlowFile << realTime << " " << uUp << " " << uUpLeft << " " << uUpRight << " "
            << uDown << " " << uDownLeft << " " << uDownRight << " "
            << hUp << " " << hUpLeft << " " << hUpRight << " "
            << hDown << " " << hDownLeft << " " << hDownRight << " "
            << flowedMassUp << " " << flowedMassDown << endl;
    hongKongFlowFile.close();


}

void IO::exportFront(const LB& lb) {

    const double flowFront = lb.maxHeight(0) * lb.unit.Length;

    ofstream frontFile;
    frontFile.open(frontFileName.c_str(), ios::app);
    //hongKongFlowFile << "time percPlastic\n";
    frontFile << realTime << " " << flowFront << endl;
    frontFile.close();
}

void IO::exportTop(const LB& lb) {

    const double flowTop = lb.maxHeight(1) * lb.unit.Length;

    ofstream topFile;
    topFile.open(topFileName.c_str(), ios::app);
    //hongKongFlowFile << "time percPlastic\n";
    topFile << realTime << " " << flowTop << endl;
    topFile.close();
}

void IO::exportHeapHeight(const LB& lb) {

    const double heapHeight = lb.maxHeight(2) * lb.unit.Length;

    ofstream heapFile;
    heapFile.open(heapFileName.c_str(), ios::app);
    //hongKongFlowFile << "time percPlastic\n";
    heapFile << realTime << " " << heapHeight << endl;
    heapFile.close();
}

void IO::exportInclineFlow(const LB& lb) {

    double heapHeight = lb.maxHeight(2) * lb.unit.Length;

    const unsigned int sizeY = lb.lbSize[1];
    doubleList inclineFlowVelocity, inclineFlowPressure;
    doubleList inclineFlowViscosity, inclineFlowHeight, inclineFlowMass;
    inclineFlowHeight.resize(sizeY);
    inclineFlowVelocity.resize(sizeY);
    inclineFlowPressure.resize(sizeY);
    inclineFlowViscosity.resize(sizeY);
    inclineFlowMass.resize(sizeY);

    for (int iy = 0; iy < sizeY; iy++) {
        inclineFlowVelocity[iy] = 0.0;
        inclineFlowPressure[iy] = 0.0;
        inclineFlowViscosity[iy] = 0.0;
        inclineFlowHeight[iy] = 0.0;
        inclineFlowMass[iy] = 0.0;
    }


    for (int it = 0; it < lb.totPossibleNodes; ++it) {
        const unsigned int xHere = lb.getX(it);
        const unsigned int zHere = lb.getZ(it);
        if (xHere == 1 && zHere == 1) {
            const unsigned int yHere = lb.getY(it);
            if (lb.nodes.count(it) != 0) {
                if (lb.nodes.at(it).isActive()) {
                    inclineFlowVelocity[yHere] = lb.nodes.at(it).u.dot(Xp) * lb.unit.Speed;
                    inclineFlowPressure[yHere] = 0.3333333 * (lb.nodes.at(it).n - lb.fluidMaterial.initDensity) * lb.unit.Pressure;
                    inclineFlowViscosity[yHere] = lb.nodes.at(it).visc * lb.unit.KinVisc;
                    inclineFlowMass[yHere] = lb.nodes.at(it).mass;
                }
            }
            inclineFlowHeight[yHere] = (double(yHere) - 0.5) * lb.unit.Length;
        }
    }


    ofstream inclineFlowFile;
    inclineFlowFile.open(inclineFlowFileName.c_str(), ios::app);


    //hongKongFlowFile << "time percPlastic\n";
    inclineFlowFile << realTime;
    for (int iy = 0; iy < sizeY; iy++) {
        inclineFlowFile << " " << inclineFlowHeight[iy];
    }
    for (int iy = 0; iy < sizeY; iy++) {
        inclineFlowFile << " " << inclineFlowVelocity[iy];
    }
    for (int iy = 0; iy < sizeY; iy++) {
        inclineFlowFile << " " << inclineFlowPressure[iy];
    }
    for (int iy = 0; iy < sizeY; iy++) {
        inclineFlowFile << " " << inclineFlowViscosity[iy];
    }
    for (int iy = 0; iy < sizeY; iy++) {
        inclineFlowFile << " " << inclineFlowMass[iy];
    }
    inclineFlowFile << endl;
    inclineFlowFile.close();

}

double IO::getSurface(doubleList& surfaceParticles) {
    static const double percentile = 0.9;

    double surfaceLevel = 0.0;
    const unsigned int surface = std::max(1, int(double(surfaceParticles.size()) * percentile));
    sort(surfaceParticles.begin(), surfaceParticles.end());
    for (int i = surface - 1; i < surfaceParticles.size(); ++i) {
        surfaceLevel += surfaceParticles[i];
    }
    surfaceLevel /= double(surfaceParticles.size() - surface + 1);
    return surfaceLevel;
}

void IO::exportHongKongBarrier(DEM& dem) {

    static const double barrierHeigth = 0.5;
    static const double midChannel = 0.5 * (dem.walls[2].p.dot(Yp) + dem.walls[3].p.dot(Yp));

    tVect forcesPartLeft(0.0, 0.0, 0.0), forcesPartRight(0.0, 0.0, 0.0);
    tVect forcesHydroLeft(0.0, 0.0, 0.0), forcesHydroRight(0.0, 0.0, 0.0);

    // collect data
    for (int o = 0; o < dem.objects.size(); ++o) {
        const double objY = dem.objects[o].x0.dot(Yp);
        const tVect objPartForce = dem.objects[o].FParticle;
        const tVect objHydroForce = dem.objects[o].FHydro;
        if (objY < midChannel) {
            forcesPartLeft += objPartForce;
            forcesHydroLeft += objHydroForce;
        } else {
            forcesPartRight += objPartForce;
            forcesHydroRight += objHydroForce;
        }

    }

    // adding wall forces
    if (dem.hongkongSmoothWall) {
        forcesPartLeft += dem.walls[6].FParticle;
        forcesHydroLeft += dem.walls[6].FHydro;
        forcesPartRight += dem.walls[7].FParticle;
        forcesHydroRight += dem.walls[7].FHydro;
    }

    ofstream hongKongForceFile;
    hongKongForceFile.open(hongKongForceFileName.c_str(), ios::app);
    //hongKongFlowFile << "time percPlastic\n";
    hongKongForceFile << realTime << " ";
    forcesPartLeft.print(hongKongForceFile);
    forcesHydroLeft.print(hongKongForceFile);
    forcesPartRight.print(hongKongForceFile);
    forcesHydroRight.print(hongKongForceFile);
    hongKongForceFile << endl;
    hongKongForceFile.close();

    //        static const double sensorSize=0.06;
    //    static const double barrierHeigth=0.5;
    //    static const unsigned int intervals=std::ceil(barrierHeigth/sensorSize);
    //    static const double midChannel=0.5*(dem.walls[2].p.dot(Yp)+dem.walls[3].p.dot(Yp));
    //    
    //    vecList forcesPartLeft,forcesPartRight;
    //    vecList forcesHydroLeft,forcesHydroRight;
    //    forcesPartLeft.resize(intervals);
    //    forcesPartRight.resize(intervals);
    //    forcesHydroLeft.resize(intervals);
    //    forcesHydroRight.resize(intervals);
    //    for (int i=0; i<intervals; ++i) {
    //        forcesPartLeft[i]=Zero;
    //        forcesPartRight[i]=Zero;
    //        forcesHydroLeft[i]=Zero;
    //        forcesHydroRight[i]=Zero;
    //    }
    //    
    //    // collect data
    //    for (int o = 0; o < dem.objects.size(); ++o) {
    //        const double objY=dem.objects[o].x0.dot(Yp);
    //        const double objZ=dem.objects[o].x0.dot(Zp);
    //        const unsigned int objSector=std::floor(objZ/sensorSize);
    //        const tVect objPartForce=dem.objects[o].FParticle;
    //        const tVect objHydroForce=dem.objects[o].FHydro;
    //        if (objY<midChannel) {
    //            forcesPartLeft[objSector]+=objPartForce;
    //            forcesHydroLeft[objSector]+=objHydroForce;
    //        }
    //        else {
    //            forcesPartRight[objSector]+=objPartForce;
    //            forcesHydroRight[objSector]+=objHydroForce;
    //        }
    //
    //    }
    //    
    //    
    //    hongKongForceFile.open(hongKongForceFileName.c_str(), ios::app);
    //    //hongKongFlowFile << "time percPlastic\n";
    //    hongKongForceFile << realTime << " " << intervals<< " ";
    //    for (int i=0; i<intervals; ++i) {
    //        forcesPartLeft[i].print(hongKongForceFile);
    //    }
    //    for (int i=0; i<intervals; ++i) {
    //        forcesHydroLeft[i].print(hongKongForceFile);
    //    }
    //    for (int i=0; i<intervals; ++i) {
    //        forcesPartRight[i].print(hongKongForceFile);
    //    }
    //    for (int i=0; i<intervals; ++i) {
    //        forcesHydroRight[i].print(hongKongForceFile);
    //    }
    //    hongKongForceFile<<endl;
    //    hongKongForceFile.close();



}

void IO::exportTriaxial(DEM& dem) {

    ofstream triaxialFile;
    triaxialFile.open(triaxialFileName.c_str(), ios::app);
    //hongKongFlowFile << "time percPlastic\n";
    triaxialFile << realTime << " ";
    triaxialFile << dem.walls[6].p.dot(Xp) << " ";
    triaxialFile << dem.walls[7].p.dot(Yp) << " ";
    triaxialFile << dem.walls[8].p.dot(Zp) << " ";
    triaxialFile << dem.pressureX << " ";
    triaxialFile << dem.pressureY << " ";
    triaxialFile << dem.pressureZ << " ";
    triaxialFile << endl;
    triaxialFile.close();


}

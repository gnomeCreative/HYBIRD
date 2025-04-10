
#include "IO2.h"
#include "gpu/LB2.h"

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// PUBLIC FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/

void IO2::outputStep(LB2& lb, DEM& dem) {


    //// PLOTTING PHASE  ////////////////////////////////////////////////////////////////
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


    //        const double deltaLB = std::chrono::duration<double, std::micro>(lb.endLBStep - lb.startLBStep).count();
    //        cout << "t=" << deltaLB << " ";
    //        if (PARAMS.freeSurface) {
    //            const double deltaFreeSurface = std::chrono::duration<double, std::micro>(lb.endFreeSurfaceStep - lb.startFreeSurfaceStep).count();
    //            const double deltaUpdateMass = std::chrono::duration<double, std::micro>(lb.endUpdateMassStep - lb.startUpdateMassStep).count();
    //            const double deltaUpdateInterface = std::chrono::duration<double, std::micro>(lb.endUpdateInterfaceStep - lb.startUpdateInterfaceStep).count();
    //            const double deltaFindMutants = std::chrono::duration<double, std::micro>(lb.endFindMutantsStep - lb.startFindMutantsStep).count();
    //            const double deltaSmoothenInterface_1 = std::chrono::duration<double, std::micro>(lb.endSmoothenInterfaceStep_1 - lb.startSmoothenInterfaceStep_1).count();
    //            const double deltaSmoothenInterface_2 = std::chrono::duration<double, std::micro>(lb.endSmoothenInterfaceStep_2 - lb.startSmoothenInterfaceStep_2).count();
    //            const double deltaUpdateMutants = std::chrono::duration<double, std::micro>(lb.endUpdateMutantsStep - lb.startUpdateMutantsStep).count();
    //            const double deltaRemoveIsolated = std::chrono::duration<double, std::micro>(lb.endRemoveIsolatedStep - lb.startRemoveIsolatedStep).count();
    //            const double deltaRedistributeMass = std::chrono::duration<double, std::micro>(lb.endRedistributeMassStep - lb.startRedistributeMassStep).count();
    //            //                cout << "t_fs=" << deltaFreeSurface << " ";
    //            //                cout << "t_um=" << deltaUpdateMass << " ";
    //            //                cout << "t_ui=" << deltaUpdateInterface << " ";
    //            //                cout << "t_fm=" << deltaFindMutants << " ";
    //            //                cout << "t_si1=" << deltaSmoothenInterface_1 << " ";
    //            //                cout << "t_si2=" << deltaSmoothenInterface_2 << " ";
    //            //                cout << "n_fn=" << lb.filledNodes.size()<< " ";
    //            //                cout << "n_en=" << lb.emptiedNodes.size()<< " ";
    //            //                cout << "n_nin=" << lb.newInterfaceNodes.size()<< " ";
    //            //                cout << "t_um=" << deltaUpdateMutants << " ";
    //            //                cout << "t_ri=" << deltaRemoveIsolated << " ";
    //            //                cout << "t_rm=" << deltaRedistributeMass << " ";
    //            //                cout << "n_fs=" << lb.interfaceNodes.size() << " ";
            }
    //        if (demSolver) {
    //            const double deltaCoupling = std::chrono::duration<double, std::micro>(lb.endCouplingStep - lb.startCouplingStep).count();
    //            cout << "t_c=" << deltaCoupling << " ";
    //        }

            //exportMaxSpeedFluid(lb);
    //        exportFreeSurfaceExtent(lb);
    //        exportFluidFlowRate(lb);
    //        exportFluidMass(lb);
    //        exportFluidCenterOfMass(lb);
            //switch (PARAMS.fluidMaterial.rheologyModel) {
            //    case BINGHAM:
            //    case FRICTIONAL:
            //    case VOELLMY:
            //    {
            //        exportPlasticity(lb);
            //        break;
            //    }
            //}
    //        exportMeanViscosity(lb);
    //    }

    //    if (dem.elmts.size()) {
    //        exportParticleFlowRate(dem);
    //        exportMaxSpeedParticles(dem);
    //        exportForces(dem);
    //        exportParticleCenterOfMass(dem);
    //        exportParticleCoordination(dem);
    //        exportParticleOverlap(dem);
    //    }

    //    if (dem.walls.size() > 0) {
    //        exportWallForce(dem);
    //    }
    //    //
    //    // update energies
    //    totalKineticEnergy = 0.0;
    //    energyExit = false;
    //    if (dem.elmts.size()) {
    //        dem.updateEnergy(totalKineticEnergy);
    //    }
    //    if (lbmSolver) {
    //        lb.updateEnergy(totalKineticEnergy);
    //    }
    //    exportEnergy(dem, lb);
    //    if (totalKineticEnergy < energyStopThreshold && PARAMS.time > minimumIterations) {
    //        energyExit = true;
    //    }
    //    if (dem.objects.size()) {
    //        exportForceObstacle(dem.objects);
    //    }
    //    switch (problemName) {
    //        case SHEARCELL:
    //        {
    //            exportShearCell(lb, dem);
    //            break;
    //        }
    //        case OPENBARRIER:
    //        {
    //            if (dem.objects.size()) {
    //                exportForceObstacleElement(dem);
    //                exportMomentObstacleElement(dem);
    //                exportHPartObstacle(dem);
    //                exportKmPartObstacle(dem);
    //            }
    //            break;
    //        }
    //        case HONGKONG:
    //        {
    //            exportHongKongFlow(dem);
    //            exportHongKongBarrier(dem);
    //            break;
    //        }
    //        case INCLINEFLOW:
    //        {
    //            exportInclineFlow(lb);
    //            break;
    //        }
    //        case HEAP:
    //        {
    //            exportHeapHeight(lb);
    //            break;
    //        }
    //        case TRIAXIAL:
    //        {
    //            exportTriaxial(dem);
    //            break;
    //        }
    //        case MANGENEY:
    //        {
    //            exportFront(lb);
    //            exportTop(lb);
    //            break;
    //        }
    //        case HK_LARGE:
    //        case HK_SMALL:
    //        case KELVIN:
    //        {
    //            exportFront(lb);
    //            break;
    //        }


    //    }

    //    if (lbmSolver && fluid2DExpTime > 0) {

    //        if (!initializedPlanarFile) {
    //            initialize2DFile(lb);
    //            initializedPlanarFile = true;
    //        }

    //        update2DFile(lb);

    //    }

    //    if (singleObjects.size() > 0) {
    //        exportSingleObjects(dem.objects);
    //    }

    //    if (objectGroupBegin.size() > 0) {
    //        exportGroupForce(dem.objects);
    //    }

    //    if (flowLevelBegin.size() > 0) {
    //        exportFlowLevel(lb);
    //    }


    //    // closing file
        cout << endl;
    //    exportFile << endl;
    //    exportFile.close();
    //    cout.flush();


    }

    // FILE CREATION PHASE  ////////////////////////////////////////////////////////////////
    createFiles(lb, dem);

    //if (lbmSolver) {
    //    const unsigned int fluid2DExpCounter = (fluid2DExpTime > 0 ? static_cast<unsigned int> (realTime / fluid2DExpTime) + 1 : 0);
    //    if (fluid2DExpCounter > lastFluid2DExp) {
    //        lastFluid2DExp = fluid2DExpCounter;
    //        char filePathBuffer [1024];
    //        sprintf(filePathBuffer, fluid2DFileFormat.c_str(), currentTimeStep);
    //        create2DFile(lb, filePathBuffer);
    //    }
    //}
}

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// PRIVATE FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/

// paraview files

void IO2::createFiles(LB2& lb, const DEM& dem) {
    // write vtk at regular interval defined with the input file

    if (lbmSolver) {

        const unsigned int fluidExpCounter = (fluidExpTime > 0 ? static_cast<unsigned int> (realTime / fluidExpTime) + 1 : 0);
        if (fluidExpCounter > lastFluidExp) {
            lastFluidExp = fluidExpCounter;
            char filePathBuffer [1024];
            sprintf(filePathBuffer, fluidFileFormat.c_str(), currentTimeStep);
            //if (fluidLagrangianFormat == ParaviewFormat::Ascii) {
            //    exportEulerianParaviewFluid(lb, filePathBuffer);  // Benchmark: 2994s
            //} else if (fluidLagrangianFormat == ParaviewFormat::BinaryLowMem) {
            //    // exportEulerianParaviewFluid_binary(lb, filePathBuffer);  // Benchmark: 315s
            //    exportEulerianParaviewFluid_binaryv2(lb, filePathBuffer);  // Benchmark: 211s
            //} else {
                exportEulerianParaviewFluid_binaryv3(lb, filePathBuffer);  // Benchmark: 120s (but requires ~200mb extra memory)
            //}
        }

        const unsigned int fluidLagrangianExpCounter = (fluidLagrangianExpTime > 0 ? static_cast<unsigned int> (realTime / fluidLagrangianExpTime) + 1 : 0);
        if (fluidLagrangianExpCounter > lastFluidLagrangianExp) {
            lastFluidLagrangianExp = fluidLagrangianExpCounter;
            char filePathBuffer [1024];
            sprintf(filePathBuffer, fluidLagrangianFileFormat.c_str(), currentTimeStep);
            //if (currentTimeStep>120) {
                // if (fluidLagrangianFormat == ParaviewFormat::Ascii) {
                    //exportLagrangianParaviewFluid(lb, filePathBuffer);
                // } else {
                    exportLagrangianParaviewFluid_binaryv3(lb, filePathBuffer);
                // }
            //}
        }

        //const unsigned int fluidRecycleExpCounter = (fluidRecycleExpTime > 0 ? static_cast<unsigned int> (realTime / fluidRecycleExpTime) + 1 : 0);
        //if (fluidRecycleExpCounter > lastFluidRecycleExp) {
        //    lastFluidRecycleExp = fluidRecycleExpCounter;
        //    char filePathBuffer [1024];
        //    sprintf(filePathBuffer, fluidRecycleFileFormat.c_str(), currentTimeStep);
        //    // requires the pbcShift, contained in the neighborList function
        //    exportRecycleFluid(lb, filePathBuffer);
        //}

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

            //const unsigned int partRecycleExpCounter = (partRecycleExpTime > 0 ? static_cast<unsigned int> (realTime / partRecycleExpTime) + 1 : 0);
            //if (partRecycleExpCounter > lastPartRecycleExp) {
            //    lastPartRecycleExp = partRecycleExpCounter;
            //    char filePathBuffer [1024];
            //    sprintf(filePathBuffer, partRecycleFileFormat.c_str(), currentTimeStep);
            //    // requires the pbcShift, contained in the neighbourList function
            //    exportRecycleParticles(dem.elmts, dem.pbcs, filePathBuffer);
            //}
    }
    //if (dem.objects.size() > 0) {

    //    const unsigned int objectExpCounter = (objectExpTime > 0 ? static_cast<unsigned int> (realTime / objectExpTime) + 1 : 0);
    //    if (objectExpCounter > lastObjectExp) {
    //        lastObjectExp = objectExpCounter;
    //        char filePathBuffer [1024];
    //        sprintf(filePathBuffer, objectFileFormat.c_str(), currentTimeStep);
    //        exportParaviewObjects(dem.objects, filePathBuffer);
    //    }
    //}
    //
    //if (dem.cylinders.size() > 0) {
    //    
    //    const unsigned int CylinderExpCounter = (cylinderExpTime > 0 ? static_cast<unsigned int> (realTime / cylinderExpTime) + 1 : 0);
    //    if (CylinderExpCounter > lastCylinderExp) {
    //        lastCylinderExp = CylinderExpCounter;
    //        char filePathBuffer [1024];
    //        sprintf(filePathBuffer, cylinderFileFormat.c_str(), currentTimeStep);
    //        exportParaviewCylinders(dem.cylinders, filePathBuffer);
    //    }
    //    
    //}
    
}
/*
void IO2::initialize2DFile(const LB2& lb) {

    // initialize containers for 2D files
    planarHeight.resize(PARAMS.lbSize[0]);
    planarVel.resize(PARAMS.lbSize[0]);
    planarPointVel.resize(PARAMS.lbSize[0]);
    planarLevel.resize(PARAMS.lbSize[0]);
    maxPlanarHeight.resize(PARAMS.lbSize[0]);
    maxPlanarVel.resize(PARAMS.lbSize[0]);
    maxPlanarPointVel.resize(PARAMS.lbSize[0]);

    for (int i = 0; i < PARAMS.lbSize[0]; ++i) {
        planarHeight[i].resize(PARAMS.lbSize[1]);
        planarVel[i].resize(PARAMS.lbSize[1]);
        planarPointVel[i].resize(PARAMS.lbSize[1]);
        planarLevel[i].resize(PARAMS.lbSize[1]);
        maxPlanarHeight[i].resize(PARAMS.lbSize[1]);
        maxPlanarVel[i].resize(PARAMS.lbSize[1]);
        maxPlanarPointVel[i].resize(PARAMS.lbSize[1]);
        for (int j = 0; j < PARAMS.lbSize[1]; ++j) {
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

void IO2::update2DFile(const LB2& lb) {

    // initialize containers for 2D files
    for (int i = 0; i < PARAMS.lbSize[0]; ++i) {
        for (int j = 0; j < PARAMS.lbSize[1]; ++j) {
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
    for (int i = 0; i < PARAMS.lbSize[0]; ++i) {
        for (int j = 0; j < PARAMS.lbSize[1]; ++j) {
            if (planarHeight[i][j] > 0.0) {
                planarVel[i][j] = planarVel[i][j] / planarHeight[i][j];
            }
        }
    }

    // update global maxima
    for (int i = 0; i < PARAMS.lbSize[0]; ++i) {
        for (int j = 0; j < PARAMS.lbSize[1]; ++j) {
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

void IO2::create2DFile(const LB2& lb, const string& planarFile) {

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
    for (int i = 0; i < PARAMS.lbSize[0]; ++i) {
        const double xHere = double(i) * PARAMS.unit.Length + xScaling;
        for (int j = 0; j < PARAMS.lbSize[1]; ++j) {
            const double maxHeightHere = maxPlanarHeight[i][j] * PARAMS.unit.Length;
            //if (maxHeightHere>0.0) {
            const double yHere = double(j) * PARAMS.unit.Length + yScaling;
            const double levelHere = planarLevel[i][j] * PARAMS.unit.Length + zScaling;
            const double heightHere = planarHeight[i][j] * PARAMS.unit.Length;
            const double velHere = planarVel[i][j] * PARAMS.unit.Speed;
            const double maxVelHere = maxPlanarVel[i][j] * PARAMS.unit.Speed;
            const double pointVelHere = planarPointVel[i][j] * PARAMS.unit.Speed;
            const double maxPointVelHere = maxPlanarPointVel[i][j] * PARAMS.unit.Speed;

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

void IO2::exportFlowLevel(const LB2& lb) {

    doubleList volumeCount;
    volumeCount.resize(flowLevelBegin.size());
    for (int i = 0; i < flowLevelBegin.size(); i++) {
        volumeCount[i] = 0.0;
    }

    if (lbmSolver) {
        for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
            const node* nodeHere = *it;
            const double xCoordHere = lb.getPositionX(nodeHere->coord) * PARAMS.unit.Length;
            for (int i = 0; i < flowLevelBegin.size(); i++) {
                const double beginHere = flowLevelBegin[i];
                const double endHere = flowLevelEnd[i];
                if (xCoordHere <= endHere && xCoordHere >= beginHere) {
                    //cout<<lb.getPositionX(nodeHere->coord)*PARAMS.unit.Length<<endl;
                    volumeCount[i] += nodeHere->mass;
                }
            }
        }
    }

    for (int i = 0; i < flowLevelBegin.size(); i++) {
        // compute window length in lattice units (to round to exact window measurement in )
        const unsigned int latticeBegin = ceil(flowLevelBegin[i] / PARAMS.unit.Length);
        const unsigned int latticeEnd = floor(flowLevelEnd[i] / PARAMS.unit.Length);
        const unsigned int latticewindowSpan = latticeEnd - latticeBegin + 1;
        ASSERT(latticewindowSpan >= 1);

        const double windowSpan = double(latticewindowSpan) * PARAMS.unit.Length;

        const double windowDepth = double(PARAMS.lbSize[2] - 2) * PARAMS.unit.Length;

        const double flowVolume = volumeCount[i] * PARAMS.unit.Volume;

        const double flowLevelHere = flowVolume / (windowDepth * windowSpan);

        string flowLevelFileName = workDirectory + "/flowLevel_" + std::to_string(i) + ".dat";
        ofstream flowLevelFile;
        flowLevelFile.open(flowLevelFileName.c_str(), ios::app);
        flowLevelFile << realTime << " " << flowLevelHere << "\n";
        flowLevelFile.close();
    }
}
*/
void IO2::exportParaviewParticles(const elmtList& elmts, const particleList& particles, const string& particleFile) {

    const int one = 1;

    unsigned int activeNumber = 0;
    for (int i = 0; i < particles.size(); ++i) {
        if (particles[i].active) {
            ++activeNumber;
        }
    }

    const int Pnumber = (int)particles.size();
    // const tVect nx(1.0, 0.0, 0.0), ny(0.0, 1.0, 0.0), nz(0.0, 0.0, 1.0);
    // tVect n1, n2, n3;

    std::cout.precision(10);
    // std::cout << std::fixed;


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
    paraviewParticleFile << "    <DataArray type=\"Float64\" Name=\"MWall\" NumberOfComponents=\"3\"/>\n";
    for (int i = 0; i < Pnumber; ++i) {
        if (particles[i].active) {
            elmts[particles[i].clusterIndex].MWall.printFixedLine(paraviewParticleFile);
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
    for (unsigned int i = 1; i < activeNumber + 1; ++i) {
        paraviewParticleFile << i - 1 << "\n";
    }
    paraviewParticleFile << "    </DataArray>\n";
    paraviewParticleFile << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (unsigned int i = 1; i < activeNumber + 1; ++i) {
        paraviewParticleFile << i << "\n";
    }
    paraviewParticleFile << "    </DataArray>\n";
    paraviewParticleFile << "    <DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n";
    for (unsigned int i = 0; i < activeNumber; ++i) {
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

void IO2::exportParaviewParticles_binaryv3(const elmtList& elmts, const particleList& particles, const string& particleFile) {
    // Build a list of active particles for faster access
    std::vector<unsigned int> active_particles;
    for (unsigned int i = 0; i < particles.size(); ++i) {
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
    size_t offset = 0;
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
    for (unsigned int i = 0; i < active_particles.size(); ++i) {
        d_buffer[i] = particles[active_particles[i]].r;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // particleIndex
    offset = active_particles.size() * sizeof(unsigned int);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < active_particles.size(); ++i) {
        u_buffer[i] = particles[active_particles[i]].particleIndex;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // clusterIndex
    offset = active_particles.size() * sizeof(unsigned int);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < active_particles.size(); ++i) {
        u_buffer[i] = particles[active_particles[i]].clusterIndex;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // coordinationNumber
    offset = active_particles.size() * sizeof(unsigned int);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < active_particles.size(); ++i) {
        u_buffer[i] = elmts[particles[active_particles[i]].clusterIndex].coordination;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // v
    offset = active_particles.size() * 3 * sizeof(double);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < active_particles.size(); ++i) {
        v_buffer[i] = particles[active_particles[i]].x1;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // w
    offset = active_particles.size() * 3 * sizeof(double);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < active_particles.size(); ++i) {
        v_buffer[i] = elmts[particles[active_particles[i]].clusterIndex].wGlobal;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // FParticle
    offset = active_particles.size() * 3 * sizeof(double);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < active_particles.size(); ++i) {
        v_buffer[i] = elmts[particles[active_particles[i]].clusterIndex].FParticle;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // FGrav
    offset = active_particles.size() * 3 * sizeof(double);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < active_particles.size(); ++i) {
        v_buffer[i] = elmts[particles[active_particles[i]].clusterIndex].FGrav;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // solidIntensity
    offset = active_particles.size() * 3 * sizeof(double);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < active_particles.size(); ++i) {
        v_buffer[i] = elmts[particles[active_particles[i]].clusterIndex].solidIntensity;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // MParticle
    offset = active_particles.size() * 3 * sizeof(double);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < active_particles.size(); ++i) {
        v_buffer[i] = elmts[particles[active_particles[i]].clusterIndex].MParticle;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // FWall
    offset = active_particles.size() * 3 * sizeof(double);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < active_particles.size(); ++i) {
        v_buffer[i] = elmts[particles[active_particles[i]].clusterIndex].FWall;
    }
    paraviewParticleFile.write(t_buffer, offset);
    // MWall
    offset = active_particles.size() * 3 * sizeof(double);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < active_particles.size(); ++i) {
        v_buffer[i] = elmts[particles[active_particles[i]].clusterIndex].MWall;
    }
    paraviewParticleFile.write(t_buffer, offset);
    if (lbmSolver) {
        // FHydro
        offset = active_particles.size() * 3 * sizeof(double);
        paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        for (unsigned int i = 0; i < active_particles.size(); ++i) {
            v_buffer[i] = elmts[particles[active_particles[i]].clusterIndex].FHydro;
        }
        paraviewParticleFile.write(t_buffer, offset);
        // MHydro
        offset = active_particles.size() * 3 * sizeof(double);
        paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        for (unsigned int i = 0; i < active_particles.size(); ++i) {
            v_buffer[i] = elmts[particles[active_particles[i]].clusterIndex].MHydro;
        }
        paraviewParticleFile.write(t_buffer, offset);
        // fluidIntensity
        offset = active_particles.size() * 3 * sizeof(double);
        paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        for (unsigned int i = 0; i < active_particles.size(); ++i) {
            v_buffer[i] = elmts[particles[active_particles[i]].clusterIndex].FHydro;
        }
        paraviewParticleFile.write(t_buffer, offset);
    }

    // Points
    offset = active_particles.size() * 3 * sizeof(double);
    paraviewParticleFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < active_particles.size(); ++i) {
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
    std::fill(uc_buffer, uc_buffer + active_particles.size(), (unsigned char)1);
    paraviewParticleFile.write(t_buffer, offset);
    paraviewParticleFile << "</AppendedData>";
    paraviewParticleFile << "</VTKFile>";
    // header file closing
    paraviewParticleFile.close();
}
/*
void IO2::exportRecycleFluid(const LB2& lb, const string& fluidRecycleFile) {
    // exports a format readable by this code itself, allowing the re-use of computed particle data.
    // needs further work...

    ofstream recycleFluidFile;
    recycleFluidFile.open(fluidRecycleFile.c_str());

    recycleFluidFile << PARAMS.lbSize[0] << " " << PARAMS.lbSize[1] << " " << PARAMS.lbSize[2] << " " << lb.activeNodes.size() << "\n";

    //cout << endl;
    for (nodeMap::const_iterator imap = lb.nodes.begin(); imap != lb.nodes.end(); imap++) {
        const unsigned int it = imap->first;
        const node* nodeHere = &imap->second;
        // to gain precision, scale density and mass by 1 (initDensity)
        if (nodeHere->isActive()) {
            recycleFluidFile << it << " " << nodeHere->type << " " << nodeHere->n - PARAMS.fluidMaterial.initDensity << " " << nodeHere->u.dot(Xp) << " " << nodeHere->u.dot(Yp) << " " << nodeHere->u.dot(Zp) << " " << nodeHere->mass - PARAMS.fluidMaterial.initDensity << " " << nodeHere->visc << " ";
            for (int j = 0; j < lbmDirec; ++j) {
                recycleFluidFile << nodeHere->f[j] << " ";
            }
            recycleFluidFile << endl;
        }
    }
    recycleFluidFile.close();
}
*/
void IO2::exportLagrangianParaviewFluid(LB2& lb, const string& fluidFile) {

    const Node2 nodes = lb.getNodes();
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
    paraviewFluidFile << "  <Piece NumberOfPoints=\"" << nodes.activeCount << "\" NumberOfCells=\"" << nodes.activeCount << "\">\n";
    paraviewFluidFile << "   <PointData>\n";
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\" format=\"ascii\"/>\n";
    for (unsigned int i = 0; i < nodes.activeCount; ++i) {
        tVect uPhysic = nodes.u[nodes.activeI[i]] * PARAMS.unit.Speed;
        uPhysic.printLine(paraviewFluidFile);
    }
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"pressure\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    for (unsigned int i = 0; i < nodes.activeCount; ++i) {
        paraviewFluidFile << 0.3333333 * (nodes.n[nodes.activeI[i]] - PARAMS.fluidMaterial.initDensity) * PARAMS.unit.Pressure << "\n";
    }
    //    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"smoothedPressure\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    //    for (int i = 0; i < nodes.activeCount; ++i) {
    //        paraviewFluidFile << nodes.smoothedPressure[nodes.activeI[i]] * PARAMS.unit.Pressure << "\n";
    //    }
    //    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"n\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    //    for (int i = 0; i < nodes.activeCount; ++i) {
    //        paraviewFluidFile << nodes.n[nodes.activeI[i]] << "\n";
    //    }
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"dynVisc\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    for (unsigned int i = 0; i < nodes.activeCount; ++i) {
        paraviewFluidFile << nodes.visc[nodes.activeI[i]] * PARAMS.unit.DynVisc << "\n";
    }
    //    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"applySlip\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    //    for (int i = 0; i < nodes.activeCount; ++i) {
    //        paraviewFluidFile << nodes.applySlip[nodes.activeI[i]] << "\n";
    //    }
//    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"viscosity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
//    for (int i = 0; i < nodes.activeCount; ++i) {
//        paraviewFluidFile << nodes.visc[nodes.activeI[i]] * PARAMS.unit.KinVisc << "\n";
//    }
    if (PARAMS.fluidMaterial.rheologyModel == MUI || PARAMS.fluidMaterial.rheologyModel == FRICTIONAL || PARAMS.fluidMaterial.rheologyModel == VOELLMY) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"friction\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
        for (unsigned int i = 0; i < nodes.activeCount; ++i) {
            paraviewFluidFile << nodes.friction[nodes.activeI[i]] << "\n";
        }
    }
    //    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"friction\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    //    for (int i = 0; i < nodes.activeCount; ++i) {
    //        paraviewFluidFile << nodes.friction[nodes.activeI[i]] << "\n";
    //    }
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"mass\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    for (unsigned int i = 0; i < nodes.activeCount; ++i) {
        paraviewFluidFile << nodes.mass[nodes.activeI[i]] * PARAMS.unit.Density << "\n";
    }
    paraviewFluidFile << "    <DataArray type=\"Int8\" Name=\"type\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"8\"/>\n";
    for (unsigned int i = 0; i < nodes.activeCount; ++i) {
        paraviewFluidFile << static_cast<unsigned>(nodes.type[nodes.activeI[i]]) << "\n";
    }
    //    for (unsigned int j = 1; j < lbmDirec; ++j) {
    //        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"delta" << j << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    //        for (unsigned int i = 0; i < nodes.activeCount; ++i) {
    //            if (nodes.curved[nodes.activeI[i]] != std::numeric_limits<unsigned int>::max())) {
    //                const tVect deltaVector = nodes.curves[nodes.curved[nodes.activeI[i]]].delta[j] * PARAMS.unit.Length * v[j];
    //                deltaVector.printLine(paraviewFluidFile);
    //            } else {
    //                Zero.printLine(paraviewFluidFile);    //
    //            }
    //        }
    //        paraviewFluidFile << "    </DataArray>\n";
    //    }
#ifdef DEBUG 
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"age\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"8\"/>\n";
    for (unsigned int i = 0; i < nodes.activeCount; ++i) {
        paraviewFluidFile << nodes.age[nodes.activeI[i]] << "\n";
    }
    //    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"vv0\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    //    for (unsigned int i = 0; i < nodes.activeCount; ++i) {
    //        const double deltaZero = (nodes.f[nodes.activeI[i]*lbmDirec + 0] - coeff[0]);
    //        paraviewFluidFile << deltaZero << "\n";
    //    }
    //    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"surfNormal\" NumberOfComponents=\"3\" format=\"ascii\"/>\n";
    //    for (unsigned int i = 0; i < nodes.activeCount; ++i) {
    //        nodeHere->surfaceNormal.printLine(paraviewFluidFile);
    //    }
    //    for (unsigned int j = 1; j < lbmDirec2D; ++j) {
    //        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"vv" << j << "\" NumberOfComponents=\"3\" format=\"ascii\"/>\n";
    //        for (unsigned int i = 0; i < nodes.activeCount; ++i) {
    //            const tVect deltaVector = (nodes.f[nodes.activeI[i]*lbmDirec + two_dim[j]] - coeff[two_dim[j]]) * v[two_dim[j]];
    //            deltaVector.printLine(paraviewFluidFile);
    //        }
    //    }
#endif
    paraviewFluidFile << "   </PointData>\n";
    paraviewFluidFile << "   <CellData>\n";
    paraviewFluidFile << "   </CellData>\n";
    paraviewFluidFile << "   <Points>\n";
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n";
    for (unsigned int i = 0; i < nodes.activeCount; ++i) {
        const tVect positionHere = nodes.getPosition(nodes.activeI[i]) * PARAMS.unit.Length;
        positionHere.printLine(paraviewFluidFile);
    }
    paraviewFluidFile << "   </Points>\n";
    paraviewFluidFile << "   <Cells>\n";
    paraviewFluidFile << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (unsigned int i = 1; i < nodes.activeCount + 1; ++i) {
        paraviewFluidFile << i - 1 << "\n";
    }
    paraviewFluidFile << "    </DataArray>\n";
    paraviewFluidFile << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (unsigned int i = 1; i < nodes.activeCount + 1; ++i) {
        paraviewFluidFile << i << "\n";
    }
    paraviewFluidFile << "    </DataArray>\n";
    paraviewFluidFile << "    <DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n";
    for (unsigned int i = 0; i < nodes.activeCount; ++i) {
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

void IO2::exportLagrangianParaviewFluid_binaryv3(LB2& lb, const string& fluidFile) {
    /**
     * This function is a rewrite of exportLagrangianParaviewFluid() that writes to a binary vtkhdf5 format
     * It is intended to provide much faster fluid export performance
     **/
    const Node2 nodes = lb.getNodes();

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
    paraviewFluidFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" << (is_big_endian ? "BigEndian" : "LittleEndian") << "\"  header_type=\"UInt32\">\n";
    paraviewFluidFile << " <UnstructuredGrid GhostLevel=\"0\">\n";
    paraviewFluidFile << "  <Piece NumberOfPoints=\"" << nodes.activeCount << "\" NumberOfCells=\"" << nodes.activeCount << "\">\n";
    paraviewFluidFile << "   <PointData>\n";
    unsigned int offset = 0;
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\"/>\n";
    offset += nodes.activeCount * 3 * sizeof(double) + sizeof(unsigned int);
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    if (PARAMS.fluidMaterial.rheologyModel != NEWTONIAN || PARAMS.fluidMaterial.turbulenceOn) {
        offset += nodes.activeCount * sizeof(double) + sizeof(unsigned int);
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"dynVisc\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    }
    offset += nodes.activeCount * sizeof(double) + sizeof(unsigned int);
    if (PARAMS.fluidMaterial.rheologyModel == MUI || PARAMS.fluidMaterial.rheologyModel == FRICTIONAL || PARAMS.fluidMaterial.rheologyModel == VOELLMY) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"friction\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
        offset += nodes.activeCount * sizeof(double) + sizeof(unsigned int);
    }
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"mass\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    offset += nodes.activeCount * sizeof(double) + sizeof(unsigned int);
    paraviewFluidFile << "    <DataArray type=\"UInt8\" Name=\"type\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"16\"/>\n";
    offset += nodes.activeCount * sizeof(unsigned char) + sizeof(unsigned int);
#ifdef DEBUG
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"age\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"8\"/>\n";
    offset += nodes.activeCount * sizeof(double) + sizeof(unsigned int);
#endif
    paraviewFluidFile << "   </PointData>\n";
    paraviewFluidFile << "   <CellData>\n";
    paraviewFluidFile << "   </CellData>\n";
    paraviewFluidFile << "   <Points>\n";
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += nodes.activeCount * 3 * sizeof(double) + sizeof(unsigned int);
    paraviewFluidFile << "   </Points>\n";
    paraviewFluidFile << "   <Cells>\n";
    paraviewFluidFile << "    <DataArray type=\"UInt32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += nodes.activeCount * sizeof(unsigned int) + sizeof(unsigned int);
    paraviewFluidFile << "    <DataArray type=\"UInt32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += nodes.activeCount * sizeof(unsigned int) + sizeof(unsigned int);
    paraviewFluidFile << "    <DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += nodes.activeCount * sizeof(unsigned char) + sizeof(unsigned int);
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
    char *const t_buffer = static_cast<char*>(malloc(nodes.activeCount * 3 * sizeof(double)));
    unsigned char *const uc_buffer = reinterpret_cast<unsigned char*>(t_buffer);
    double *const d_buffer = reinterpret_cast<double*>(t_buffer);
    tVect *const v_buffer = reinterpret_cast<tVect*>(t_buffer);
    unsigned int* const ui_buffer = reinterpret_cast<unsigned int*>(t_buffer);
    // Velocity
    offset = nodes.activeCount * 3 * sizeof(double);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < nodes.activeCount; ++i) {
        v_buffer[i] = nodes.u[nodes.activeI[i]] * PARAMS.unit.Speed;
    }
    paraviewFluidFile.write(t_buffer, offset);
    // Pressure
    offset = nodes.activeCount * sizeof(double);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    const double THIRD_PRESSURE = 0.3333333 * PARAMS.unit.Pressure;
    for (unsigned int i = 0; i < nodes.activeCount; ++i) {
        d_buffer[i] = (nodes.n[nodes.activeI[i]] - PARAMS.fluidMaterial.initDensity) * THIRD_PRESSURE;
    }
    paraviewFluidFile.write(t_buffer, offset);
    // Dynamic Viscosity
    if (PARAMS.fluidMaterial.rheologyModel != NEWTONIAN || PARAMS.fluidMaterial.turbulenceOn) {
        offset = nodes.activeCount * sizeof(double);
        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        for (unsigned int i = 0; i < nodes.activeCount; ++i) {
            d_buffer[i] = nodes.visc[nodes.activeI[i]] * PARAMS.unit.DynVisc;
        }
        paraviewFluidFile.write(t_buffer, offset);
    }
    // Friction
    if (PARAMS.fluidMaterial.rheologyModel == MUI || PARAMS.fluidMaterial.rheologyModel == FRICTIONAL || PARAMS.fluidMaterial.rheologyModel == VOELLMY) {
        offset = nodes.activeCount * sizeof(double);
        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        for (unsigned int i = 0; i < nodes.activeCount; ++i) {
            d_buffer[i] = nodes.friction[nodes.activeI[i]];
        }
        paraviewFluidFile.write(t_buffer, offset);
    }
    // Mass
    offset = nodes.activeCount * sizeof(double);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < nodes.activeCount; ++i) {
        d_buffer[i] = nodes.mass[nodes.activeI[i]] * PARAMS.unit.Density;
    }
    paraviewFluidFile.write(t_buffer, offset);
    // Type
    offset = nodes.activeCount * sizeof(unsigned char);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    std::fill(uc_buffer, uc_buffer + nodes.activeCount, (unsigned char)2);
    for (unsigned int i = 0; i < nodes.activeCount; ++i) {
        uc_buffer[i] = nodes.type[nodes.activeI[i]];
    }
    paraviewFluidFile.write(t_buffer, offset);
#ifdef DEBUG
    // Age
    offset = nodes.activeCount * sizeof(double);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    memset(d_buffer, 0, sizeof(double)* nodes.activeCount);
    for (unsigned int i = 0; i < nodes.activeCount; ++i) {
        d_buffer[i] = nodes.age[nodes.activeI[i]];
    }
    paraviewFluidFile.write(t_buffer, offset);
#endif
    // Points
    offset = nodes.activeCount * 3 * sizeof(double);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < nodes.activeCount; ++i) {
        v_buffer[i] = nodes.getPosition(nodes.activeI[i]) * PARAMS.unit.Length;
    }
    paraviewFluidFile.write(t_buffer, offset);
    // Connectivity
    offset = nodes.activeCount * sizeof(unsigned int);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < nodes.activeCount; ++i) {
        ui_buffer[i] = i ;
    }
    paraviewFluidFile.write(t_buffer, offset);
    // Offsets
    offset = nodes.activeCount * sizeof(unsigned int);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    for (unsigned int i = 0; i < nodes.activeCount; ++i) {
        ui_buffer[i] = i + 1;
    }
    paraviewFluidFile.write(t_buffer, offset);
    // Types
    offset = nodes.activeCount * sizeof(unsigned char);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    std::fill(uc_buffer, uc_buffer + nodes.activeCount, (unsigned char)1);
    paraviewFluidFile.write(t_buffer, offset);
    paraviewFluidFile << "</AppendedData>";
    paraviewFluidFile << "</VTKFile>\n";
    // data file closing
    paraviewFluidFile.close();
    free(t_buffer);
}

//void IO2::exportEulerianParaviewFluid(const LB2& lb, const string& fluidFile) {
//
//    const double zero = 0.0;
//
//    // start printing all the crap required for Paraview
//    // header file opening
//    ofstream paraviewFluidFile;
//
//    paraviewFluidFile.open(fluidFile.c_str());
////    paraviewFluidFile << std::scientific << std::setprecision(4);
//    // writing on header file
//    paraviewFluidFile << "<VTKFile type=\"ImageData\" version=\"0.1\">\n";
//    paraviewFluidFile << " <ImageData WholeExtent=\"0 " << PARAMS.lbSize[0] - 1 << " 0 " << PARAMS.lbSize[1] - 1 << " 0 " << PARAMS.lbSize[2] - 1 << "\" "
//            << "Origin=\"" << -0.5 * PARAMS.unit.Length << " " << -0.5 * PARAMS.unit.Length << " " << -0.5 * PARAMS.unit.Length << "\" "
//            << "Spacing=\"" << PARAMS.unit.Length << " " << PARAMS.unit.Length << " " << PARAMS.unit.Length << "\">\n";
//    paraviewFluidFile << "  <Piece Extent=\"0 " << PARAMS.lbSize[0] - 1 << " 0 " << PARAMS.lbSize[1] - 1 << " 0 " << PARAMS.lbSize[2] - 1 << "\">\n";
//    paraviewFluidFile << "   <PointData>\n";
//    paraviewFluidFile << "    <DataArray type=\"Int8\" Name=\"type\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"8\">\n";
//    //    paraviewFluidFile << std::fixed << std::setprecision(0);
//    for (int i = 0; i < PARAMS.totPossibleNodes; ++i) {
//        if (lb.nodes.count(i) != 0) {
//            if (!lb.nodes.at(i).isInsideParticle()) {
//                paraviewFluidFile << static_cast<unsigned>(lb.nodes.at(i).type) << " ";
//            } else {
//                paraviewFluidFile << 1 << " ";
//            }
//        } else {
//            paraviewFluidFile << 2 << " ";
//        }
//    }
//    paraviewFluidFile << "    </DataArray>\n";
//    paraviewFluidFile << std::scientific << std::setprecision(4);
//    //    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"n\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\">\n";
//    //    for (int i = 0; i < PARAMS.totPossibleNodes; ++i) {
//    //        if (lb.nodes.count(i)!=0) {
//    //            paraviewFluidFile << lb.nodes.at(i).n - 1.0 << " ";
//    //        } else {
//    //            paraviewFluidFile << 0 << " ";
//    //        }
//    //    }
//    //    paraviewFluidFile << "    </DataArray>\n";
//    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\" format=\"ascii\">\n";
//    for (int i = 0; i < PARAMS.totPossibleNodes; ++i) {
//        if (lb.nodes.count(i) != 0) {
//            const tVect uPhysic = lb.nodes.at(i).u * PARAMS.unit.Speed;
//            uPhysic.print(paraviewFluidFile);
//        } else {
//            Zero.print(paraviewFluidFile);
//        }
//    }
//    paraviewFluidFile << "    </DataArray>\n";
//    //    for (int j = 1; j < lbmDirec; ++j) {
//    //        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"delta" << j << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
//    //        for (int i = 0; i < PARAMS.totPossibleNodes; ++i) {
//    //            if (lb.nodes.count(i) != 0) {
//    //                if (lb.nodes.at(i).curved!=0) {
//    //                    const tVect deltaVector = lb.nodes.at(i).curved->delta[j] * PARAMS.unit.Length*vDirec[j];
//    //                    deltaVector.print(paraviewFluidFile);
//    //                } else {
//    //                    Zero.print(paraviewFluidFile);
//    //                }
//    //            } else {
//    //                Zero.print(paraviewFluidFile);
//    //            }
//    //        }
//    //        paraviewFluidFile << "    </DataArray>\n";
//    //    }
//    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"pressure\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\">\n";
//    for (int i = 0; i < PARAMS.totPossibleNodes; ++i) {
//        if (lb.nodes.count(i) != 0) {
//            if (lb.nodes.at(i).n != 0.0) {
//                paraviewFluidFile << 0.3333333 * (lb.nodes.at(i).n - 1.0) * PARAMS.unit.Pressure << " ";
//            } else {
//                paraviewFluidFile << lb.nodes.at(i).n << " ";
//            }
//        } else {
//            paraviewFluidFile << zero << " ";
//        }
//    }
//    paraviewFluidFile << "    </DataArray>\n";
//    if (PARAMS.fluidMaterial.rheologyModel != NEWTONIAN || PARAMS.fluidMaterial.turbulenceOn) {
//        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"dynVisc\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\">\n";
//        for (int i = 0; i < PARAMS.totPossibleNodes; ++i) {
//            if (lb.nodes.count(i) != 0) {
//                paraviewFluidFile << lb.nodes.at(i).visc * PARAMS.unit.DynVisc << " ";
//            } else {
//                paraviewFluidFile << zero << " ";
//            }
//        }
//        paraviewFluidFile << "    </DataArray>\n";
//    }
//    if (PARAMS.fluidMaterial.rheologyModel == MUI || PARAMS.fluidMaterial.rheologyModel == FRICTIONAL || PARAMS.fluidMaterial.rheologyModel == VOELLMY) {
//        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"friction\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\">\n";
//        for (int i = 0; i < PARAMS.totPossibleNodes; ++i) {
//            if (lb.nodes.count(i) != 0) {
//                paraviewFluidFile << lb.nodes.at(i).friction << " ";
//            } else {
//                paraviewFluidFile << zero << " ";
//            }
//        }
//        paraviewFluidFile << "    </DataArray>\n";
//    }
//    if (lb.freeSurface) {
//        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"amass\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\">\n";
//        for (int i = 0; i < PARAMS.totPossibleNodes; ++i) {
//            if (lb.nodes.count(i) != 0) {
//                paraviewFluidFile << lb.nodes.at(i).mass * PARAMS.unit.Density << " ";
//            } else {
//                paraviewFluidFile << zero << " ";
//            }
//        }
//        paraviewFluidFile << "    </DataArray>\n";
//    }
//    if (lb.lbTopography) {
//        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"topSurface\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\">\n";
//        for (int i = 0; i < PARAMS.totPossibleNodes; ++i) {
//            if (lb.nodes.count(i) != 0) {
//                if (lb.nodes.at(i).isTopography()) {
//                    paraviewFluidFile << 1.0 << " ";
//                } else {
//                    paraviewFluidFile << zero << " ";
//                }
//            } else {
//                paraviewFluidFile << zero << " ";
//            }
//        }
//        paraviewFluidFile << "    </DataArray>\n";
//    }
//    //    if (demSolver) {
//    //        paraviewFluidFile << "    <DataArray type=\"Int16\" Name=\"solidIndex\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"8\">\n";
//    //        paraviewFluidFile << std::fixed << std::setprecision(0);
//    //        for (int i = 0; i < PARAMS.totPossibleNodes; ++i) {
//    //            if (lb.nodes.count(i) != 0) {
//    //                paraviewFluidFile << lb.nodes.at(i).getSolidIndex() << " ";
//    //            } else {
//    //                paraviewFluidFile << zero << " ";
//    //            }
//    //        }
//    //        paraviewFluidFile << "    </DataArray>\n";
//    //        paraviewFluidFile << std::scientific << std::setprecision(4);
//    //    }
//    paraviewFluidFile << "   </PointData>\n";
//    paraviewFluidFile << "   <CellData>\n";
//    paraviewFluidFile << "   </CellData>\n";
//    paraviewFluidFile << "  </Piece>\n";
//    paraviewFluidFile << " </ImageData>\n";
//    paraviewFluidFile << "</VTKFile>\n";
//    // data file closing
//    paraviewFluidFile.close();
//}
//
//void IO2::exportEulerianParaviewFluid_binary(const LB2& lb, const string& fluidFile) {
//    /**
//     * This function is a rewrite of exportEulerianParaviewFluid() that writes to a binary vtkhdf5 format
//     * It is intended to provide much faster fluid export performance
//     **/
//
//    // start printing all the crap required for Paraview
//    // header file opening
//    ofstream paraviewFluidFile;
//
//    paraviewFluidFile.open(fluidFile.c_str(), std::ios::binary);
////    paraviewFluidFile << std::scientific << std::setprecision(4);
//    // writing on header file
//
//    // Endianness (the order of bits within a byte) depends on the processor hardware
//    // but it's probably LittleEndian, IBM processors are the only default BigEndian you're likely to come across
//    static const uint16_t m_endianCheck(0x00ff);
//    const bool is_big_endian ( *((const uint8_t*)&m_endianCheck) == 0x0);
//    paraviewFluidFile << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"" << (is_big_endian ? "BigEndian" : "LittleEndian") << "\"  header_type=\"UInt32\">\n";
//    paraviewFluidFile << " <ImageData WholeExtent=\"0 " << PARAMS.lbSize[0] - 1 << " 0 " << PARAMS.lbSize[1] - 1 << " 0 " << PARAMS.lbSize[2] - 1 << "\" "
//            << "Origin=\"" << -0.5 * PARAMS.unit.Length << " " << -0.5 * PARAMS.unit.Length << " " << -0.5 * PARAMS.unit.Length << "\" "
//            << "Spacing=\"" << PARAMS.unit.Length << " " << PARAMS.unit.Length << " " << PARAMS.unit.Length << "\">\n";
//    paraviewFluidFile << "  <Piece Extent=\"0 " << PARAMS.lbSize[0] - 1 << " 0 " << PARAMS.lbSize[1] - 1 << " 0 " << PARAMS.lbSize[2] - 1 << "\">\n";
//    paraviewFluidFile << "   <PointData>\n";
//    unsigned int offset = 0;
//    paraviewFluidFile << "    <DataArray type=\"UInt8\" Name=\"type\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"8\"/>\n";
//    offset += PARAMS.totPossibleNodes * sizeof(unsigned char) + sizeof(unsigned int);
//    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\"/>\n";
//    offset += PARAMS.totPossibleNodes * 3 * sizeof(double) + sizeof(unsigned int);
//    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
//    offset += PARAMS.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
//    if (PARAMS.fluidMaterial.rheologyModel != NEWTONIAN || PARAMS.fluidMaterial.turbulenceOn) {
//        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"dynVisc\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
//        offset += PARAMS.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
//    }
//    if (PARAMS.fluidMaterial.rheologyModel == MUI || PARAMS.fluidMaterial.rheologyModel == FRICTIONAL || PARAMS.fluidMaterial.rheologyModel == VOELLMY) {
//        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"friction\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
//        offset += PARAMS.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
//    }
//    if (lb.freeSurface) {
//        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"amass\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
//        offset += PARAMS.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
//    }
//    if (lb.lbTopography) {
//        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"topSurface\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
//        offset += PARAMS.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
//    }
//    paraviewFluidFile << "   </PointData>\n";
//    paraviewFluidFile << "   <CellData>\n";
//    paraviewFluidFile << "   </CellData>\n";
//    paraviewFluidFile << "  </Piece>\n";
//    paraviewFluidFile << " </ImageData>\n";
//    paraviewFluidFile << " <AppendedData encoding=\"raw\">\n  _";
//    /**
//     * Based on the sparse documentation at https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
//     * and alot of testing.
//     * Inside <AppendedData> the binary dump must be preceded by an underscore (_)
//     * Each DataArray's binary dump must be preceded by it's length.
//     * The length should be exported as the integer type specified as header_type in the opening <VTKFile> tag
//     * The offset specified in the above <DataArray> tag refers to the offset from the start of the whole binary dump to the start of the length
//     */
//    static const double zero = 0;
//    static const tVect tVect_zero = {0,0,0};
//    // May be faster to iterate the (sorted) map and catch missing items
//    // Might be possible to do inline compression, this would slow export but unclear by how much
//    offset = PARAMS.totPossibleNodes * sizeof(unsigned char);
//    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
//    for (int i = 0; i < PARAMS.totPossibleNodes; ++i) {
//        const auto &node = lb.nodes.find(i);
//        if (node != lb.nodes.end()) {
//            if (node->second.isInsideParticle()) {
//                paraviewFluidFile.write(reinterpret_cast<const char*>(&node->second.type), sizeof(node->second.type));
//            } else {
//                const unsigned char t = 1;
//                paraviewFluidFile.write(reinterpret_cast<const char*>(&t), sizeof(unsigned char));
//            }
//        } else {
//            const unsigned char t = 2;
//            paraviewFluidFile.write(reinterpret_cast<const char*>(&t), sizeof(unsigned char));
//        }
//    }
//    offset = PARAMS.totPossibleNodes * 3 * sizeof(double);
//    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
//    for (int i = 0; i < PARAMS.totPossibleNodes; ++i) {
//        const auto &node = lb.nodes.find(i);
//        if (node != lb.nodes.end()) {
//            static const tVect uPhysic = lb.nodes.at(i).u * PARAMS.unit.Speed;
//            paraviewFluidFile.write(reinterpret_cast<const char*>(&uPhysic), sizeof(tVect));
//        } else {
//            paraviewFluidFile.write(reinterpret_cast<const char*>(&tVect_zero), sizeof(tVect));
//        }
//    }
//    offset = PARAMS.totPossibleNodes * sizeof(double);
//    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
//    for (int i = 0; i < PARAMS.totPossibleNodes; ++i) {
//        const auto &node = lb.nodes.find(i);
//        if (node != lb.nodes.end()) {
//            const double t = 0.3333333 * (node->second.n - 1.0) * PARAMS.unit.Pressure;
//            paraviewFluidFile.write(reinterpret_cast<const char*>(&t), sizeof(double));
//        } else {
//            paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
//        }
//    }
//    if (PARAMS.fluidMaterial.rheologyModel != NEWTONIAN || PARAMS.fluidMaterial.turbulenceOn) {
//        offset = PARAMS.totPossibleNodes * sizeof(double);
//        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
//        for (int i = 0; i < PARAMS.totPossibleNodes; ++i) {
//            const auto &node = lb.nodes.find(i);
//            if (node != lb.nodes.end()) {
//                const double t = node->second.visc * PARAMS.unit.DynVisc;
//                paraviewFluidFile.write(reinterpret_cast<const char*>(&t), sizeof(double));
//            } else {
//                paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
//            }
//        }
//    }
//    if (PARAMS.fluidMaterial.rheologyModel == MUI || PARAMS.fluidMaterial.rheologyModel == FRICTIONAL || PARAMS.fluidMaterial.rheologyModel == VOELLMY) {
//        offset = PARAMS.totPossibleNodes * sizeof(double);
//        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
//        for (int i = 0; i < PARAMS.totPossibleNodes; ++i) {
//            const auto &node = lb.nodes.find(i);
//            if (node != lb.nodes.end()) {
//                paraviewFluidFile.write(reinterpret_cast<const char*>(&node->second.friction), sizeof(double));
//            } else {
//                paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
//            }
//        }
//    }
//    if (lb.freeSurface) {
//        offset = PARAMS.totPossibleNodes * sizeof(double);
//        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
//        for (int i = 0; i < PARAMS.totPossibleNodes; ++i) {
//            const auto &node = lb.nodes.find(i);
//            if (node != lb.nodes.end()) {
//                const double t = node->second.mass * PARAMS.unit.Density;
//                paraviewFluidFile.write(reinterpret_cast<const char*>(&t), sizeof(double));
//            } else {
//                paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
//            }
//        }
//    }
//    if (lb.lbTopography) {
//        offset = PARAMS.totPossibleNodes * sizeof(double);
//        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
//        for (int i = 0; i < PARAMS.totPossibleNodes; ++i) {
//            const auto &node = lb.nodes.find(i);
//            if (node != lb.nodes.end()) {
//                if (node->second.isTopography()) {
//                    const double one = 1.0;
//                    paraviewFluidFile.write(reinterpret_cast<const char*>(&one), sizeof(double));
//                } else {
//                    paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
//                }
//            } else {
//                paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
//            }
//        }
//    }
//    paraviewFluidFile << "</AppendedData>";
//    paraviewFluidFile << "</VTKFile>\n";
//    // data file closing
//    paraviewFluidFile.close();
//}
//
//void IO2::exportEulerianParaviewFluid_binaryv2(const LB2& lb, const string& fluidFile) {
//    /**
//     * This function is a rewrite of exportEulerianParaviewFluid() that writes to a binary vtkhdf5 format
//     * It is intended to provide much faster fluid export performance
//     **/
//
//    // start printing all the crap required for Paraview
//    // header file opening
//    ofstream paraviewFluidFile;
//
//    paraviewFluidFile.open(fluidFile.c_str(), std::ios::binary);
////    paraviewFluidFile << std::scientific << std::setprecision(4);
//    // writing on header file
//
//    // Endianness (the order of bits within a byte) depends on the processor hardware
//    // but it's probably LittleEndian, IBM processors are the only default BigEndian you're likely to come across
//    static const uint16_t m_endianCheck(0x00ff);
//    const bool is_big_endian ( *((const uint8_t*)&m_endianCheck) == 0x0);
//    paraviewFluidFile << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"" << (is_big_endian ? "BigEndian" : "LittleEndian") << "\"  header_type=\"UInt32\">\n";
//    paraviewFluidFile << " <ImageData WholeExtent=\"0 " << PARAMS.lbSize[0] - 1 << " 0 " << PARAMS.lbSize[1] - 1 << " 0 " << PARAMS.lbSize[2] - 1 << "\" "
//            << "Origin=\"" << -0.5 * PARAMS.unit.Length << " " << -0.5 * PARAMS.unit.Length << " " << -0.5 * PARAMS.unit.Length << "\" "
//            << "Spacing=\"" << PARAMS.unit.Length << " " << PARAMS.unit.Length << " " << PARAMS.unit.Length << "\">\n";
//    paraviewFluidFile << "  <Piece Extent=\"0 " << PARAMS.lbSize[0] - 1 << " 0 " << PARAMS.lbSize[1] - 1 << " 0 " << PARAMS.lbSize[2] - 1 << "\">\n";
//    paraviewFluidFile << "   <PointData>\n";
//    unsigned int offset = 0;
//    paraviewFluidFile << "    <DataArray type=\"UInt8\" Name=\"type\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"8\"/>\n";
//    offset += PARAMS.totPossibleNodes * sizeof(unsigned char) + sizeof(unsigned int);
//    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\"/>\n";
//    offset += PARAMS.totPossibleNodes * 3 * sizeof(double) + sizeof(unsigned int);
//    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
//    offset += PARAMS.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
//    if (PARAMS.fluidMaterial.rheologyModel != NEWTONIAN || PARAMS.fluidMaterial.turbulenceOn) {
//        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"dynVisc\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
//        offset += PARAMS.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
//    }
//    if (PARAMS.fluidMaterial.rheologyModel == MUI || PARAMS.fluidMaterial.rheologyModel == FRICTIONAL || PARAMS.fluidMaterial.rheologyModel == VOELLMY) {
//        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"friction\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
//        offset += PARAMS.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
//    }
//    if (lb.freeSurface) {
//        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"amass\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
//        offset += PARAMS.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
//    }
//    if (lb.lbTopography) {
//        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"topSurface\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
//        offset += PARAMS.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
//    }
//    paraviewFluidFile << "   </PointData>\n";
//    paraviewFluidFile << "   <CellData>\n";
//    paraviewFluidFile << "   </CellData>\n";
//    paraviewFluidFile << "  </Piece>\n";
//    paraviewFluidFile << " </ImageData>\n";
//    paraviewFluidFile << " <AppendedData encoding=\"raw\">\n  _";
//    /**
//     * Based on the sparse documentation at https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
//     * and alot of testing.
//     * Inside <AppendedData> the binary dump must be preceded by an underscore (_)
//     * Each DataArray's binary dump must be preceded by it's length.
//     * The length should be exported as the integer type specified as header_type in the opening <VTKFile> tag
//     * The offset specified in the above <DataArray> tag refers to the offset from the start of the whole binary dump to the start of the length
//     */
//    static const double zero = 0;
//    static const tVect tVect_zero = {0,0,0};
//    // May be faster to iterate the (sorted) map and catch missing items
//    // Might be possible to do inline compression, this would slow export but unclear by how much
//    offset = PARAMS.totPossibleNodes * sizeof(unsigned char);
//    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
//    int i = 0;
//    for (const auto &[key, node] : lb.nodes) {
//        for (;i<key;++i) {
//            const unsigned char t = 2;
//            paraviewFluidFile.write(reinterpret_cast<const char*>(&t), sizeof(unsigned char));
//        }
//        if (node.isInsideParticle()) {
//            paraviewFluidFile.write(reinterpret_cast<const char*>(&node.type), sizeof(unsigned char));
//        } else {
//            const unsigned char t = 1;
//            paraviewFluidFile.write(reinterpret_cast<const char*>(&t), sizeof(unsigned char));
//        }
//        ++i;
//    }
//    offset = PARAMS.totPossibleNodes * 3 * sizeof(double);
//    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
//    i = 0;
//    for (const auto &[key, node] : lb.nodes) {
//        for (;i<key;++i) {
//            paraviewFluidFile.write(reinterpret_cast<const char*>(&tVect_zero), sizeof(tVect));
//        }
//        const tVect uPhysic = lb.nodes.at(i).u * PARAMS.unit.Speed;
//        paraviewFluidFile.write(reinterpret_cast<const char*>(&uPhysic), sizeof(tVect));
//        ++i;
//    }
//    offset = PARAMS.totPossibleNodes * sizeof(double);
//    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
//    i = 0;
//    for (const auto &[key, node] : lb.nodes) {
//        for (;i<key;++i) {
//            paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
//        }
//        const double t = 0.3333333 * (node.n - 1.0) * PARAMS.unit.Pressure;
//        paraviewFluidFile.write(reinterpret_cast<const char*>(&t), sizeof(double));
//        ++i;
//    }
//    if (PARAMS.fluidMaterial.rheologyModel != NEWTONIAN || PARAMS.fluidMaterial.turbulenceOn) {
//        offset = PARAMS.totPossibleNodes * sizeof(double);
//        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
//        i = 0;
//        for (const auto &[key, node] : lb.nodes) {
//            for (;i<key;++i) {
//                paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
//            }
//            const double t = node.visc * PARAMS.unit.DynVisc;
//            paraviewFluidFile.write(reinterpret_cast<const char*>(&t), sizeof(double));
//            ++i;
//        }
//    }
//    if (PARAMS.fluidMaterial.rheologyModel == MUI || PARAMS.fluidMaterial.rheologyModel == FRICTIONAL || PARAMS.fluidMaterial.rheologyModel == VOELLMY) {
//        offset = PARAMS.totPossibleNodes * sizeof(double);
//        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
//        i = 0;
//        for (const auto &[key, node] : lb.nodes) {
//            for (;i<key;++i) {
//                paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
//            }
//            paraviewFluidFile.write(reinterpret_cast<const char*>(&node.friction), sizeof(double));
//            ++i;
//        }
//    }
//    if (lb.freeSurface) {
//        offset = PARAMS.totPossibleNodes * sizeof(double);
//        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
//        i = 0;
//        for (const auto &[key, node] : lb.nodes) {
//            for (;i<key;++i) {
//                paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
//            }
//            const double t = node.mass * PARAMS.unit.Density;
//            paraviewFluidFile.write(reinterpret_cast<const char*>(&t), sizeof(double));
//            ++i;
//        }
//    }
//    if (lb.lbTopography) {
//        offset = PARAMS.totPossibleNodes * sizeof(double);
//        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
//        i = 0;
//        for (const auto &[key, node] : lb.nodes) {
//            for (;i<key;++i) {
//                paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
//            }
//            if (node.isTopography()) {
//                const double one = 1.0;
//                paraviewFluidFile.write(reinterpret_cast<const char*>(&one), sizeof(double));
//            } else {
//                paraviewFluidFile.write(reinterpret_cast<const char*>(&zero), sizeof(double));
//            }
//            ++i;
//        }
//    }
//    paraviewFluidFile << "</AppendedData>";
//    paraviewFluidFile << "</VTKFile>\n";
//    // data file closing
//    paraviewFluidFile.close();
//}

void IO2::exportEulerianParaviewFluid_binaryv3(LB2& lb, const string& fluidFile) {
    /**
     * This function is a rewrite of exportEulerianParaviewFluid() that writes to a binary vtkhdf5 format
     * It is intended to provide much faster fluid export performance
     **/
    const Node2 nodes = lb.getNodes();

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
    paraviewFluidFile << " <ImageData WholeExtent=\"0 " << PARAMS.lbSize[0] - 1 << " 0 " << PARAMS.lbSize[1] - 1 << " 0 " << PARAMS.lbSize[2] - 1 << "\" "
            << "Origin=\"" << -0.5 * PARAMS.unit.Length << " " << -0.5 * PARAMS.unit.Length << " " << -0.5 * PARAMS.unit.Length << "\" "
            << "Spacing=\"" << PARAMS.unit.Length << " " << PARAMS.unit.Length << " " << PARAMS.unit.Length << "\">\n";
    paraviewFluidFile << "  <Piece Extent=\"0 " << PARAMS.lbSize[0] - 1 << " 0 " << PARAMS.lbSize[1] - 1 << " 0 " << PARAMS.lbSize[2] - 1 << "\">\n";
    paraviewFluidFile << "   <PointData>\n";
    unsigned int offset = 0;
    paraviewFluidFile << "    <DataArray type=\"UInt8\" Name=\"type\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"8\"/>\n";
    offset += PARAMS.totPossibleNodes * sizeof(unsigned char) + sizeof(unsigned int);
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\"/>\n";
    offset += PARAMS.totPossibleNodes * 3 * sizeof(double) + sizeof(unsigned int);
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    offset += PARAMS.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
    if (PARAMS.fluidMaterial.rheologyModel != NEWTONIAN || PARAMS.fluidMaterial.turbulenceOn) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"dynVisc\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
        offset += PARAMS.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
    }
    if (PARAMS.fluidMaterial.rheologyModel == MUI || PARAMS.fluidMaterial.rheologyModel == FRICTIONAL || PARAMS.fluidMaterial.rheologyModel == VOELLMY) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"friction\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
        offset += PARAMS.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
    }
    if (PARAMS.freeSurface) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"amass\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
        offset += PARAMS.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
    }
    if (PARAMS.lbTopography) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"topSurface\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
        offset += PARAMS.totPossibleNodes * sizeof(double) + sizeof(unsigned int);
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
    static char *const t_buffer = static_cast<char*>(malloc(PARAMS.totPossibleNodes * 3 * sizeof(double)));
    static unsigned char *const uc_buffer = reinterpret_cast<unsigned char*>(t_buffer);
    static double *const d_buffer = reinterpret_cast<double*>(t_buffer);
    static tVect *const v_buffer = reinterpret_cast<tVect*>(t_buffer);
    // Type
    offset = PARAMS.totPossibleNodes * sizeof(unsigned char);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    std::fill(uc_buffer, uc_buffer + PARAMS.totPossibleNodes, (unsigned char)2);
    for (unsigned int i = 0; i < nodes.count; ++i) {
        uc_buffer[i] = nodes.isInsideParticle(i) ? nodes.type[i] : 1;
    }
    paraviewFluidFile.write(t_buffer, offset);
    // Velocity
    offset = PARAMS.totPossibleNodes * 3 * sizeof(double);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    memset(v_buffer, 0, sizeof(tVect) * PARAMS.totPossibleNodes);
    for (unsigned int i = 0; i < nodes.count; ++i) {
        v_buffer[i] = nodes.u[i] * PARAMS.unit.Speed;
    }
    paraviewFluidFile.write(t_buffer, offset);
    // Pressure
    offset = PARAMS.totPossibleNodes * sizeof(double);
    paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
    memset(d_buffer, 0, sizeof(double) * PARAMS.totPossibleNodes);
    const double THIRD_PRESSURE = 0.3333333 * PARAMS.unit.Pressure;
    for (unsigned int i = 0; i < nodes.count; ++i) {
        d_buffer[i] = (nodes.n[i] - 1.0) * THIRD_PRESSURE;
    }
    paraviewFluidFile.write(t_buffer, offset);
    // Dynamic Viscosity
    if (PARAMS.fluidMaterial.rheologyModel != NEWTONIAN || PARAMS.fluidMaterial.turbulenceOn) {
        offset = PARAMS.totPossibleNodes * sizeof(double);
        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        memset(d_buffer, 0, sizeof(double) * PARAMS.totPossibleNodes);
        for (unsigned int i = 0; i < nodes.count; ++i) {
            d_buffer[i] = nodes.visc[i] * PARAMS.unit.DynVisc;
        }
        paraviewFluidFile.write(t_buffer, offset);
    }
    // Friction
    if (PARAMS.fluidMaterial.rheologyModel == MUI || PARAMS.fluidMaterial.rheologyModel == FRICTIONAL || PARAMS.fluidMaterial.rheologyModel == VOELLMY) {
        offset = PARAMS.totPossibleNodes * sizeof(double);
        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        memset(d_buffer, 0, sizeof(double) * PARAMS.totPossibleNodes);
        for (unsigned int i = 0; i < nodes.count; ++i) {
            d_buffer[i] = nodes.friction[i];
        }
        paraviewFluidFile.write(t_buffer, offset);
    }
    // AMass
    if (PARAMS.freeSurface) {
        offset = PARAMS.totPossibleNodes * sizeof(double);
        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        memset(d_buffer, 0, sizeof(double) * PARAMS.totPossibleNodes);
        for (unsigned int i = 0; i < nodes.count; ++i) {
            d_buffer[i] = nodes.mass[i] * PARAMS.unit.Density;
        }
        paraviewFluidFile.write(t_buffer, offset);
    }
    // Top Surface
    if (PARAMS.lbTopography) {
        offset = PARAMS.totPossibleNodes * sizeof(double);
        paraviewFluidFile.write(reinterpret_cast<const char*>(&offset), sizeof(unsigned int));
        memset(d_buffer, 0, sizeof(double) * PARAMS.totPossibleNodes);
        for (unsigned int i = 0; i < nodes.count; ++i) {
            if (nodes.type[i] == TOPO) {
                d_buffer[i] = 1.0;
            }
        }
        paraviewFluidFile.write(t_buffer, offset);
    }
    paraviewFluidFile << "</AppendedData>";
    paraviewFluidFile << "</VTKFile>\n";
    // data file closing
    paraviewFluidFile.close();
}

//// print export stuff
//
//void IO2::exportMaxSpeedFluid(const LB2& lb) {
//
//    static const double soundSpeed = 1.0 / sqrt(3);
//    // fluid max velocity
//    double maxFluidSpeed = 0.0;
//    for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
//        const node* nodeHere = *it;
//        maxFluidSpeed = std::max(maxFluidSpeed, nodeHere->u.norm2());
//    }
//    maxFluidSpeed = sqrt(maxFluidSpeed);
//    cout << "MaxFSpeed= " << std::scientific << std::setprecision(2) << maxFluidSpeed * PARAMS.unit.Speed << "(Ma=" << std::scientific << std::setprecision(2) << maxFluidSpeed / soundSpeed << ") ";
//    exportFile << "MaxFSpeed= " << std::scientific << std::setprecision(2) << maxFluidSpeed * PARAMS.unit.Speed << "(Ma=" << std::scientific << std::setprecision(2) << maxFluidSpeed / soundSpeed << ") ";
//
//    // printing max speed
//    maxFluidSpeedFile.open(maxFluidSpeedFileName.c_str(), ios::app);
//    maxFluidSpeedFile << realTime << " " << maxFluidSpeed * PARAMS.unit.Speed << "\n";
//    maxFluidSpeedFile.close();
//}
//
//void IO2::exportFreeSurfaceExtent(const LB2& lb) {
//
//    // fluid max velocity
//    unsigned int maxX = 0;
//    unsigned int maxX_Y = 0;
//    unsigned int maxX_Z = 0;
//    unsigned int minX = UINT_MAX;
//    unsigned int minX_Y = 0;
//    unsigned int minX_Z = 0;
//    //
//    unsigned int maxY = 0;
//    unsigned int maxY_Z = 0;
//    unsigned int maxY_X = 0;
//    unsigned int minY = UINT_MAX;
//    unsigned int minY_Z = 0;
//    unsigned int minY_X = 0;
//    //
//    unsigned int maxZ = 0;
//    unsigned int maxZ_X = 0;
//    unsigned int maxZ_Y = 0;
//    unsigned int minZ = UINT_MAX;
//    unsigned int minZ_X = 0;
//    unsigned int minZ_Y = 0;
//    for (nodeList::const_iterator it = lb.interfaceNodes.begin(); it != lb.interfaceNodes.end(); ++it) {
//        const node* nodeHere = *it;
//        const unsigned int index = nodeHere->coord;
//        
//        const double xHere = lb.getPositionX(index);
//        const double yHere = lb.getPositionY(index);
//        const double zHere = lb.getPositionZ(index);
//
//        // max
//        if(xHere>maxX) {
//            maxX=xHere;
//            maxX_Y=yHere;
//            maxX_Z=zHere;
//        }
//        if(yHere>maxY) {
//            maxY=yHere;
//            maxY_Z=zHere;
//            maxY_X=xHere;
//        }
//        if(zHere>maxZ) {
//            maxZ=zHere;
//            maxZ_X=xHere;
//            maxZ_Y=yHere;
//        }
//        
//        //min
//        if(xHere<minX) {
//            minX=xHere;
//            minX_Y=yHere;
//            minX_Z=zHere;
//        }
//        if(yHere<minY) {
//            minY=yHere;
//            minY_Z=zHere;
//            minY_X=xHere;
//        }
//        if(zHere<minZ) {
//            minZ=zHere;
//            minZ_X=xHere;
//            minZ_Y=yHere;
//        }
//                
//    }
//
//    // printing max speed
//    freeSurfaceExtentFile.open(freeSurfaceExtentFileName.c_str(), ios::app);
//    freeSurfaceExtentFile << realTime << " " << maxX * PARAMS.unit.Length << " " << maxX_Y * PARAMS.unit.Length << " " << maxX_Z * PARAMS.unit.Length
//                                      << " " << minX * PARAMS.unit.Length << " " << minX_Y * PARAMS.unit.Length << " " << minX_Z * PARAMS.unit.Length
//                                      << " " << maxY * PARAMS.unit.Length << " " << maxY_Z * PARAMS.unit.Length << " " << maxY_X * PARAMS.unit.Length
//                                      << " " << minY * PARAMS.unit.Length << " " << minY_Z * PARAMS.unit.Length << " " << minY_X * PARAMS.unit.Length
//                                      << " " << maxZ * PARAMS.unit.Length << " " << maxZ_X * PARAMS.unit.Length << " " << maxZ_Y * PARAMS.unit.Length
//                                      << " " << minZ * PARAMS.unit.Length << " " << minZ_X * PARAMS.unit.Length << " " << minZ_Y * PARAMS.unit.Length<< "\n";
//    freeSurfaceExtentFile.close();
//}
//
//void IO2::exportFluidFlowRate(const LB2& lb) {
//    // fluid flow rate
//    tVect flowRate(0.0, 0.0, 0.0);
//    for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
//        const node* nodeHere = *it;
//        if (!nodeHere->isInsideParticle()) {
//            flowRate += nodeHere->u * nodeHere->mass;
//        }
//    }
//
//    const double flowRateX = flowRate.dot(Xp) / double(PARAMS.lbSize[0] - 2);
//    const double flowRateY = flowRate.dot(Yp) / double(PARAMS.lbSize[1] - 2);
//    const double flowRateZ = flowRate.dot(Zp) / double(PARAMS.lbSize[2] - 2);
//
//    // printing rate
//    fluidFlowRateFile.open(fluidFlowRateFileName.c_str(), ios::app);
//    fluidFlowRateFile << realTime << " " << flowRateX * PARAMS.unit.FlowRate << " " << flowRateY * PARAMS.unit.FlowRate << " " << flowRateZ * PARAMS.unit.FlowRate << "\n";
//    fluidFlowRateFile.close();
//}
//
//void IO2::exportFluidCenterOfMass(const LB2& lb) {
//
//    // particle center of mass
//    const tVect fluidCenter = fluidCenterOfMass(lb) * PARAMS.unit.Length;
//
//    // printing particle center of mass
//    fluidCenterOfMassFile.open(fluidCenterOfMassFileName.c_str(), ios::app);
//    fluidCenterOfMassFile << realTime << " " << fluidCenter.dot(Xp) << " " << fluidCenter.dot(Yp) << " " << fluidCenter.dot(Zp) << "\n";
//    fluidCenterOfMassFile.close();
//
//}
//
//void IO2::exportFluidMass(const LB2& lb) {
//    // total fluid mass
//    double massTot = totFluidMass(lb);
//    cout << "Volume=" << std::scientific << std::setprecision(2) << massTot * PARAMS.unit.Volume << "; Mass = " << std::scientific << std::setprecision(2) << massTot * PARAMS.unit.Mass << " ";
//    exportFile << "Volume=" << std::scientific << std::setprecision(2) << massTot * PARAMS.unit.Volume << "; Mass = " << std::scientific << std::setprecision(2) << massTot * PARAMS.unit.Mass << " ";
//
//    // printing fluid mass
//    fluidMassFile.open(fluidMassFileName.c_str(), ios::app);
//    fluidMassFile << realTime << " " << massTot * PARAMS.unit.Mass << "\n";
//    fluidMassFile.close();
//}
//
//void IO2::exportPlasticity(const LB2& lb) {
    //// fluid plasticity state
    //const double percPlastic = totPlastic(lb);
    //cout << "Plastic =" << int(percPlastic) << "% ";
    //exportFile << "Plastic =" << int(percPlastic) << "% ";
    //// printing plasticity level
    //plasticityFile.open(plasticityFileName.c_str(), ios::app);
    //plasticityFile << realTime << " " << percPlastic << "\n";
    //plasticityFile.close();

//}
//
//void IO2::exportMeanViscosity(const LB2& lb) {
//    // fluid plasticity state
//    const double meanVisc = meanViscosity(lb);
//    cout << "MeanVisc =" << std::scientific << std::setprecision(2) << meanVisc * PARAMS.unit.DynVisc << " ";
//    exportFile << "MeanVisc =" << std::scientific << std::setprecision(2) << meanVisc * PARAMS.unit.DynVisc << " ";
//}
//
//void IO2::exportShearCell(const LB2& lb, const DEM& dem) {
//    // apparent viscosity from shear cell
//    double appVisc = 0.0;
//    double externalShear = 0.0;
//    double wallStress = 0.0;
//    apparentViscosity(lb, dem.walls, externalShear, wallStress, appVisc);
//    cout << "App Visc= " << std::scientific << std::setprecision(2) << appVisc << " ";
//    exportFile << "App Visc= " << std::scientific << std::setprecision(2) << appVisc << " ";
//
//    tVect xDirec = tVect(1.0, 0.0, 0.0);
//    cout << "wallDown = " << std::scientific << std::setprecision(2) << dem.walls[0].FHydro.dot(xDirec) << " wallUp = " << std::scientific << std::setprecision(2) << dem.walls[1].FHydro.dot(xDirec) << " ";
//    cout << "wallDown = " << std::scientific << std::setprecision(2) << dem.walls[0].FParticle.dot(xDirec) << " wallUp = " << std::scientific << std::setprecision(2) << dem.walls[1].FParticle.dot(xDirec) << " ";
//}
//
//void IO2::exportEnergy(const DEM& dem, const LB2& lb) {
//
//    if (dem.elmts.size()) {
//        cout << "Energy (DEM): ";
//        exportFile << "Energy (DEM): ";
//        cout << "eKin = " << std::scientific << std::setprecision(2) << dem.particleEnergy.kin << " ";
//        exportFile << "eKin = " << std::scientific << std::setprecision(2) << dem.particleEnergy.kin << " ";
//        cout << "eGrav = " << std::scientific << std::setprecision(2) << dem.particleEnergy.grav << " ";
//        exportFile << "eGrav = " << std::scientific << std::setprecision(2) << dem.particleEnergy.grav << " ";
//        cout << "eTot = " << std::scientific << std::setprecision(2) << dem.particleEnergy.total << " ";
//        exportFile << "eTot = " << std::scientific << std::setprecision(2) << dem.particleEnergy.total << " ";
//    }
//    if (lbmSolver) {
//        cout << "Energy (LBM): ";
//        cout << "eKin = " << std::scientific << std::setprecision(2) << lb.fluidEnergy.kin + lb.fluidImmersedEnergy.kin << " ";
//        exportFile << "eKin = " << std::scientific << std::setprecision(2) << lb.fluidEnergy.kin + lb.fluidImmersedEnergy.kin << " ";
//        cout << "eGrav = " << std::scientific << std::setprecision(2) << lb.fluidEnergy.grav + lb.fluidImmersedEnergy.grav << " ";
//        exportFile << "eGrav = " << std::scientific << std::setprecision(2) << lb.fluidEnergy.grav + lb.fluidImmersedEnergy.grav << " ";
//        cout << "eTot = " << std::scientific << std::setprecision(2) << lb.fluidEnergy.total + lb.fluidImmersedEnergy.total << " ";
//        exportFile << "eTot = " << std::scientific << std::setprecision(2) << lb.fluidEnergy.total + lb.fluidImmersedEnergy.total << " ";
//    }
//
//    ofstream energyFile;
//    energyFile.open(energyFileName.c_str(), ios::app);
//    // Set energyFile header
//    energyFile << std::scientific << std::setprecision(6) << realTime << " ";
//    energyFile << std::scientific << std::setprecision(10) << dem.particleEnergy.mass << " " << dem.particleEnergy.trKin << " " << dem.particleEnergy.rotKin << " " << dem.particleEnergy.grav << " ";
//    energyFile << std::scientific << std::setprecision(10) << lb.fluidEnergy.mass * PARAMS.unit.Mass << " " << lb.fluidEnergy.trKin * PARAMS.unit.Energy << " " << lb.fluidEnergy.grav * PARAMS.unit.Energy << " ";
//    energyFile << std::scientific << std::setprecision(10) << lb.fluidImmersedEnergy.mass * PARAMS.unit.Mass << " " << lb.fluidImmersedEnergy.trKin * PARAMS.unit.Energy << " " << lb.fluidImmersedEnergy.grav * PARAMS.unit.Energy << endl;
//    energyFile.close();
//
//}
//
//// data elaboration
//
//double IO2::totPlastic(const LB2& lb) const {
//    // prints the total mass in the free fluid domain
//    unsigned int totPlastic = 0;
//    unsigned int totActive = 0;
//
//    for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
//        const node* nodeHere = *it;
//        ++totActive;
//        if (nodeHere->visc > 0.95 * PARAMS.fluidMaterial.lbMaxVisc) {
//            ++totPlastic;
//        }
//    }
//    const double perc = 100.0 * double(totPlastic) / double(totActive);
//    return perc;
//}
//
//double IO2::totFluidMass(const LB2& lb) const {
//    // returns the total mass in the free fluid domain
//    double mass = 0.0;
//
//    if (lbmSolver) {
//        for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
//            const node* nodeHere = *it;
//
//            if (!nodeHere->isInsideParticle()) {
//                mass += nodeHere->mass;
//            }
//        }
//    }
//    return mass;
//}
//
//double IO2::meanViscosity(const LB2& lb) const {
//    // prints the total mass in the free fluid domain
//    double meanVisc = 0.0;
//    unsigned int counter = 0;
//
//    if (lbmSolver) {
//        for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
//            const node* nodeHere = *it;
//
//            ++counter;
//            meanVisc += nodeHere->visc;
//        }
//        meanVisc = meanVisc / double(counter);
//        return meanVisc;
//    }
//    return 0;
//}
//
//void IO2::apparentViscosity(const LB2& lb, const wallList& walls, double& externalShear, double& wallStress, double& appVisc) const {
//    // computes the apparent viscosity of a shear cell, as the ratio between shear at the wall and imposed shear rate
//    // the cell is sheared in direction x and has moving wall in z
//    appVisc = 0.0;
//    tVect xDirec = tVect(1.0, 0.0, 0.0);
//
//    double cellSize = double(PARAMS.lbSize[2] - 2) * PARAMS.unit.Length;
//
//    externalShear = 0.5 * (walls[1].vel.dot(xDirec) - walls[0].vel.dot(xDirec)) / cellSize;
//
//    double plateSize = double(PARAMS.lbSize[0] - 2) * double(PARAMS.lbSize[1] - 2) * PARAMS.unit.Length * PARAMS.unit.Length;
//
//    wallStress = -0.5 * (walls[1].FHydro.dot(xDirec) + walls[1].FParticle.dot(xDirec) - walls[0].FHydro.dot(xDirec) - walls[0].FParticle.dot(xDirec)) / plateSize;
//
//    appVisc = 0.5 * wallStress / externalShear;
//}
//
//tVect IO2::fluidCenterOfMass(const LB2& lb) const {
//    // prints the total mass in the free fluid domain
//    tVect center(0.0, 0.0, 0.0);
//    double totMass(0.0);
//
//    if (lbmSolver) {
//        for (nodeList::const_iterator it = lb.activeNodes.begin(); it != lb.activeNodes.end(); ++it) {
//            const node* nodeHere = *it;
//            const unsigned int index = nodeHere->coord;
//            if (!nodeHere->isInsideParticle()) {
//                center += nodeHere->mass * lb.getPosition(index);
//                totMass += nodeHere->mass;
//            }
//        }
//        return center / totMass;
//    } else {
//        return Zero;
//    }
//}

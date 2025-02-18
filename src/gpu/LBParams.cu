#include "LBParams.h"

#include "getpot.h"
#include "macros.h"

#ifdef USE_CUDA
__constant__ LBParams d_PARAMS;
#endif
LBParams h_PARAMS;

void LBParams::latticeDefinition() {
    // LATTICE PARAMETERS  ////////////////
    //size-dependent; the others are in the lattice.h

    // domain movement variables
    shift[2] = lbSize[0] * lbSize[1];
    shift[1] = lbSize[0];
    shift[0] = 1;
    //
    domain[2] = lbSize[2] * lbSize[1] * lbSize[0] - 2 * shift[2];
    domain[1] = lbSize[1] * lbSize[0] - 2 * shift[1];
    domain[0] = lbSize[0] - 2 * shift[0];

    // standard neighbors shifting
    // order is O,x,y,z,xy,yz,zx.
    ne[0] = 0;
    //
    ne[1] = (int)shift[0];
    ne[2] = -(int)shift[0];
    //
    ne[3] = (int)shift[1];
    ne[4] = -(int)shift[1];
    //
    ne[5] = (int)shift[2];
    ne[6] = -(int)shift[2];
    //
    ne[7] = (int)shift[0] + shift[1];
    ne[8] = -(int)shift[0] - shift[1];
    ne[9] = -(int)shift[0] + shift[1];
    ne[10] = (int)shift[0] - shift[1];
    //
    ne[11] = (int)shift[1] + shift[2];
    ne[12] = -(int)shift[1] - shift[2];
    ne[13] = -(int)shift[1] + shift[2];
    ne[14] = (int)shift[1] - shift[2];
    //
    ne[15] = (int)shift[2] + shift[0];
    ne[16] = -(int)shift[2] - shift[0];
    ne[17] = -(int)shift[2] + shift[0];
    ne[18] = (int)shift[2] - shift[0];

}

void LBParams::latticeBoltzmannGet(GetPot& configFile, GetPot& commandLine,  LBInitParams &initP) {
    // conversion units //////////////////////////////////////////////////////////////////////////////
    // fluid and granular matter are solved in different measure units
    // for a reference, check Feng, Han, Owen, 2007
    // measure units for DEM solver are in the international system

    cout << "Getting LBM info from config file" << endl;

    // restart
    PARSE_CLASS_MEMBER(configFile, lbRestart, "restartFluid", false);
    if (lbRestart) {
        PARSE_CLASS_MEMBER(configFile, initP.lbRestartFile, "fluidRestartFile", "");
    }
    // imposed volume
    PARSE_CLASS_MEMBER(configFile, imposeFluidVolume, "imposeFluidVolume", false);
    if (imposeFluidVolume) {
        PARSE_CLASS_MEMBER(configFile, imposedFluidVolume, "imposedFluidVolume", 0.0);
    }
    // increase volume
    PARSE_CLASS_MEMBER(configFile, increaseVolume, "increaseVolume", false);
    if (increaseVolume) {
        ASSERT(imposedFluidVolume == 0.0);
        PARSE_CLASS_MEMBER(configFile, deltaVolume, "deltaVolume", 0.0);
        PARSE_CLASS_MEMBER(configFile, deltaTime, "deltaTime", 0.0);
    }
    if (imposeFluidVolume) {
        cout << "Fixed volume:  " << imposedFluidVolume << endl;
    }
    else if (increaseVolume) {
        cout << "Volume increasing by : " << deltaVolume << " in " << deltaTime << " time" << endl;
    }
    // two-relaxation time (TRT) solver
    PARSE_CLASS_MEMBER(configFile, TRTsolver, "TRTsolver", false);
    if (TRTsolver) {
        PARSE_CLASS_MEMBER(configFile, magicNumber, "magicNumber", 0.25);
    }
    // topography
    PARSE_CLASS_MEMBER(configFile, lbTopography, "applyTopography", false);
    if (lbTopography) {
        PARSE_CLASS_MEMBER(configFile, lbTopographySurface, "fluidFromTopography", false);
        PARSE_CLASS_MEMBER(configFile, initP.lbTopographyFile, "topographyFile", "");
        PARSE_CLASS_MEMBER(configFile, translateTopographyX, "translateTopographyX", 0.0);
        PARSE_CLASS_MEMBER(configFile, translateTopographyY, "translateTopographyY", 0.0);
        PARSE_CLASS_MEMBER(configFile, translateTopographyZ, "translateTopographyZ", 0.0);
    }
    if (lbTopography && !lbRestart) {
        cout << "Using topography file for boundary and initial surface:  " << initP.lbTopographyFile << endl;
    }
    else if (lbTopography && lbRestart) {
        cout << "Using topography file for boundary: " << initP.lbTopographyFile << ", and restart file for fluid:  " << initP.lbRestartFile << endl;
    }
    else if (lbRestart) {
        cout << "Computation from LB restart file " << initP.lbRestartFile << endl;
    }
    else {
        cout << "Computation from scratch - for interface geometry check LB::initializeInterface() " << endl;
    }

    // primary
    PARSE_CLASS_MEMBER(configFile, unit.Length, "latticeSpacing", 0.0);
    ASSERT(unit.Length > 0);
    PARSE_CLASS_MEMBER(configFile, unit.Time, "fluidTimeStep", 1.0);
    //ASSERT(unit.Time > 0);
    PARSE_CLASS_MEMBER(configFile, unit.Density, "fluidDensity", 1.0);
    ASSERT(unit.Density > 0);

    // fixing length: domain size
    double lbSizeX, lbSizeY, lbSizeZ;
    PARSE_CLASS_MEMBER(configFile, lbSizeX, "domainSizeX", 0.0);
    ASSERT(lbSizeX > 0.0);
    PARSE_CLASS_MEMBER(configFile, lbSizeY, "domainSizeY", 0.0);
    ASSERT(lbSizeY > 0.0);
    PARSE_CLASS_MEMBER(configFile, lbSizeZ, "domainSizeZ", 0.0);
    ASSERT(lbSizeZ > 0.0);
    // scaling and adapting length unit to get exact domain discretization
    lbSize[0] = int(floor(lbSizeX / unit.Length + 0.001)) + 2;
    lbSize[1] = int(floor(lbSizeY / unit.Length + 0.001)) + 2;
    lbSize[2] = int(floor(lbSizeZ / unit.Length + 0.001)) + 2;



    // free surface size (case NONE)
    double fluidMinX, fluidMaxX, fluidMinY, fluidMaxY, fluidMinZ, fluidMaxZ;
    PARSE_CLASS_MEMBER(configFile, fluidMinX, "fluidMinX", 0.0);
    PARSE_CLASS_MEMBER(configFile, fluidMaxX, "fluidMaxX", lbSizeX);
    PARSE_CLASS_MEMBER(configFile, fluidMinY, "fluidMinY", 0.0);
    PARSE_CLASS_MEMBER(configFile, fluidMaxY, "fluidMaxY", lbSizeY);
    PARSE_CLASS_MEMBER(configFile, fluidMinZ, "fluidMinZ", 0.0);
    PARSE_CLASS_MEMBER(configFile, fluidMaxZ, "fluidMaxZ", lbSizeZ);
    freeSurfaceBorders[0] = int(fluidMinX / unit.Length);
    freeSurfaceBorders[1] = int(fluidMaxX / unit.Length);
    freeSurfaceBorders[2] = int(fluidMinY / unit.Length);
    freeSurfaceBorders[3] = int(fluidMaxY / unit.Length);
    freeSurfaceBorders[4] = int(fluidMinZ / unit.Length);
    freeSurfaceBorders[5] = int(fluidMaxZ / unit.Length);

    // domain in physical units
    lbPhysicalSize[0] = double(lbSize[0] - 2) * unit.Length;
    lbPhysicalSize[1] = double(lbSize[1] - 2) * unit.Length;
    lbPhysicalSize[2] = double(lbSize[2] - 2) * unit.Length;

    // domain inside boundaries
    lbInnerPhysicalSize[0] = double(lbSize[0]) * unit.Length;
    lbInnerPhysicalSize[1] = double(lbSize[1]) * unit.Length;
    lbInnerPhysicalSize[2] = double(lbSize[2]) * unit.Length;

    // location of the domain boundary
    lbBoundaryLocation[0] = tVect(0.0, 0.0, 0.0);
    lbBoundaryLocation[1] = tVect(double(lbSize[0] - 2) * unit.Length, 0.0, 0.0);
    lbBoundaryLocation[2] = tVect(0.0, 0.0, 0.0);
    lbBoundaryLocation[3] = tVect(0.0, double(lbSize[1] - 2) * unit.Length, 0.0);
    lbBoundaryLocation[4] = tVect(0.0, 0.0, 0.0);
    lbBoundaryLocation[5] = tVect(0.0, 0.0, double(lbSize[2] - 2) * unit.Length);

    // material //////////////////////////////////////////////////
    string rheologyModelString;
    PARSE_CLASS_MEMBER(configFile, rheologyModelString, "rheologyModel", "none");
    if (rheologyModelString == "NEWTONIAN") fluidMaterial.rheologyModel = NEWTONIAN;
    else if (rheologyModelString == "BINGHAM") fluidMaterial.rheologyModel = BINGHAM;
    else if (rheologyModelString == "FRICTIONAL") fluidMaterial.rheologyModel = FRICTIONAL;
    else if (rheologyModelString == "VOELLMY") fluidMaterial.rheologyModel = VOELLMY;
    else if (rheologyModelString == "BAGNOLD") fluidMaterial.rheologyModel = BAGNOLD;
    else if (rheologyModelString == "MUI") fluidMaterial.rheologyModel = MUI;

    // getting viscosity
    PARSE_CLASS_MEMBER(configFile, fluidMaterial.initDynVisc, "initVisc", 1.0);
    ASSERT(fluidMaterial.initDynVisc > 0.0);

    switch (fluidMaterial.rheologyModel) {
    case NEWTONIAN:
    {
        break;
    }
    case BINGHAM:
    {
        PARSE_CLASS_MEMBER(configFile, fluidMaterial.plasticVisc, "plasticVisc", 0.0);
        ASSERT(fluidMaterial.plasticVisc >= 0.0);
        PARSE_CLASS_MEMBER(configFile, fluidMaterial.yieldStress, "yieldStress", 0.0);
        ASSERT(fluidMaterial.yieldStress >= 0.0);
        break;
    }
    case FRICTIONAL:
    {
        PARSE_CLASS_MEMBER(configFile, fluidMaterial.frictionCoefFluid, "frictionCoefFluid", 0.0);
        ASSERT(fluidMaterial.frictionCoefFluid >= 0.0 && fluidMaterial.frictionCoefFluid <= 1.0);
        PARSE_CLASS_MEMBER(configFile, fluidMaterial.basalFrictionCoefFluid, "basalFrictionCoefFluid", 0.0);
        if (fluidMaterial.basalFrictionCoefFluid == 0.0) {
            fluidMaterial.basalFrictionCoefFluid = fluidMaterial.frictionCoefFluid;
        }
        else {
            ASSERT(fluidMaterial.basalFrictionCoefFluid >= 0.0 && fluidMaterial.basalFrictionCoefFluid <= 1.0);
        }
        ASSERT(fluidMaterial.frictionCoefFluid >= 0.0 && fluidMaterial.frictionCoefFluid <= 1.0);
        break;
    }
    case VOELLMY:
    {
        PARSE_CLASS_MEMBER(configFile, fluidMaterial.frictionCoefFluid, "frictionCoefFluid", 0.0);
        PARSE_CLASS_MEMBER(configFile, fluidMaterial.basalFrictionCoefFluid, "basalFrictionCoefFluid", 0.0);
        if (fluidMaterial.basalFrictionCoefFluid == 0.0) {
            fluidMaterial.basalFrictionCoefFluid = fluidMaterial.frictionCoefFluid;
        }
        else {
            ASSERT(fluidMaterial.basalFrictionCoefFluid >= 0.0 && fluidMaterial.basalFrictionCoefFluid <= 1.0);
        }
        ASSERT(fluidMaterial.frictionCoefFluid >= 0.0 && fluidMaterial.frictionCoefFluid <= 1.0);
        PARSE_CLASS_MEMBER(configFile, fluidMaterial.particleDiameter, "particleDiameterFluid", 0.0);
        ASSERT(fluidMaterial.particleDiameter > 0.0);
        PARSE_CLASS_MEMBER(configFile, fluidMaterial.particleDensity, "particleDensityFluid", 0.0);
        ASSERT(fluidMaterial.particleDensity > 0.0);
        fluidMaterial.rhod2 = fluidMaterial.particleDensity * fluidMaterial.particleDiameter * fluidMaterial.particleDiameter;
        break;
    }
    case BAGNOLD:
    {
        PARSE_CLASS_MEMBER(configFile, fluidMaterial.particleDiameter, "particleDiameterFluid", 0.0);
        ASSERT(fluidMaterial.particleDiameter > 0.0);
        PARSE_CLASS_MEMBER(configFile, fluidMaterial.particleDensity, "particleDensityFluid", 0.0);
        ASSERT(fluidMaterial.particleDensity > 0.0);
        fluidMaterial.rhod2 = fluidMaterial.particleDensity * fluidMaterial.particleDiameter * fluidMaterial.particleDiameter;
        break;
    }
    case MUI:
    {
        PARSE_CLASS_MEMBER(configFile, fluidMaterial.frictionCoefFluid, "frictionCoefFluid", 0.0);
        PARSE_CLASS_MEMBER(configFile, fluidMaterial.basalFrictionCoefFluid, "basalFrictionCoefFluid", 0.0);
        if (fluidMaterial.basalFrictionCoefFluid == 0.0) {
            fluidMaterial.basalFrictionCoefFluid = fluidMaterial.frictionCoefFluid;
        }
        else {
            ASSERT(fluidMaterial.basalFrictionCoefFluid >= 0.0 && fluidMaterial.basalFrictionCoefFluid <= 1.0);
        }
        ASSERT(fluidMaterial.frictionCoefFluid >= 0.0 && fluidMaterial.frictionCoefFluid <= 1.0);
        // getting other parameters for voellmy model
        PARSE_CLASS_MEMBER(configFile, fluidMaterial.deltaFriction, "deltaFriction", 0.0);
        PARSE_CLASS_MEMBER(configFile, fluidMaterial.baseInertial, "baseInertial", 0.0);
        PARSE_CLASS_MEMBER(configFile, fluidMaterial.particleDiameter, "particleDiameterFluid", 0.0);
        ASSERT(fluidMaterial.particleDiameter > 0.0);
        PARSE_CLASS_MEMBER(configFile, fluidMaterial.particleDensity, "particleDensityFluid", 0.0);
        ASSERT(fluidMaterial.particleDensity > 0.0);
        break;
    }
    }
    // are we using turbulence modeling? (shutting down improves performances)
    PARSE_CLASS_MEMBER(configFile, fluidMaterial.turbulenceOn, "turbulenceSolver", 0);
    PARSE_CLASS_MEMBER(configFile, fluidMaterial.turbConst, "turbConst", 0.0);
    ASSERT(fluidMaterial.turbConst >= 0.0);

    // if given, assign an earth pressure coefficient for static load on structures
    PARSE_CLASS_MEMBER(configFile, fluidMaterial.earthPressureCoeff, "earthPressureCoeff", 1.0);
    ASSERT(fluidMaterial.earthPressureCoeff >= 0.0);

    // minimum and maximum relaxation time
    double minTau, maxTau;
    PARSE_CLASS_MEMBER(configFile, minTau, "minTau", 0.5001);
    PARSE_CLASS_MEMBER(configFile, maxTau, "maxTau", 1.0);
    fluidMaterial.lbMinVisc = (minTau - 0.5) / 3.0;
    fluidMaterial.lbMaxVisc = (maxTau - 0.5) / 3.0;

    // if time step is chosen as optimal, compute it
    // see thesis, §4.1.1
    if (unit.Time == 0.0) {
        const double tauOpt = 1.0;
        switch (fluidMaterial.rheologyModel) {
        case BINGHAM:
        {
            unit.Time = unit.Density * unit.Length * unit.Length / fluidMaterial.plasticVisc * (1.0 / 3.0) * (tauOpt - 0.5);
            break;
        }
        default:
        {
            unit.Time = unit.Density * unit.Length * unit.Length / fluidMaterial.initDynVisc * (1.0 / 3.0) * (tauOpt - 0.5);
            break;
        }
        cout << "Time step automatically chosen: " << unit.Time << " s" << endl;
        }
    }


    // compute all non-primary conversion units
    unit.setComposite();

    double initVelocityX, initVelocityY, initVelocityZ;
    PARSE_CLASS_MEMBER(configFile, initVelocityX, "fluidInitVelocityX", 0.0);
    PARSE_CLASS_MEMBER(configFile, initVelocityY, "fluidInitVelocityY", 0.0);
    PARSE_CLASS_MEMBER(configFile, initVelocityZ, "fluidInitVelocityZ", 0.0);
    initVelocity = tVect(initVelocityX, initVelocityY, initVelocityZ);
    initVelocity /= unit.Speed;

    // external force field (as acceleration)
    double lbFX, lbFY, lbFZ;
    PARSE_CLASS_MEMBER(configFile, lbFX, "forceX", 0.0);
    PARSE_CLASS_MEMBER(configFile, lbFY, "forceY", 0.0);
    PARSE_CLASS_MEMBER(configFile, lbFZ, "forceZ", 0.0);
    lbF = tVect(lbFX, lbFY, lbFZ);
    lbF /= unit.Accel;

    // rotation of the local coordinate system (rad/s)
    double rotationX, rotationY, rotationZ;
    PARSE_CLASS_MEMBER(configFile, rotationX, "rotationX", 0.0);
    PARSE_CLASS_MEMBER(configFile, rotationY, "rotationY", 0.0);
    PARSE_CLASS_MEMBER(configFile, rotationZ, "rotationZ", 0.0);
    rotationSpeed = tVect(rotationX, rotationY, rotationZ);
    rotationSpeed /= unit.AngVel;

    // position of the rotation center of the local coordinate system
    // default is (-1.0,-1.0,-1.0) to avoid computing a zero distance from the center
    double rotationCenterX, rotationCenterY, rotationCenterZ;
    PARSE_CLASS_MEMBER(configFile, rotationCenterX, "rotationCenterX", -1.0);
    PARSE_CLASS_MEMBER(configFile, rotationCenterY, "rotationCenterY", -1.0);
    PARSE_CLASS_MEMBER(configFile, rotationCenterZ, "rotationCenterZ", -1.0);
    rotationCenter = tVect(rotationCenterX, rotationCenterY, rotationCenterZ);
    rotationCenter /= unit.Length;
    
    string boundaryString;
    PARSE_CLASS_MEMBER(configFile, boundaryString, "boundary0", "");
    boundary[0] = boundary2Type(boundaryString);
    PARSE_CLASS_MEMBER(configFile, boundaryString, "boundary1", "");
    boundary[1] = boundary2Type(boundaryString);
    PARSE_CLASS_MEMBER(configFile, boundaryString, "boundary2", "");
    boundary[2] = boundary2Type(boundaryString);
    PARSE_CLASS_MEMBER(configFile, boundaryString, "boundary3", "");
    boundary[3] = boundary2Type(boundaryString);
    PARSE_CLASS_MEMBER(configFile, boundaryString, "boundary4", "");
    boundary[4] = boundary2Type(boundaryString);
    PARSE_CLASS_MEMBER(configFile, boundaryString, "boundary5", "");
    boundary[5] = boundary2Type(boundaryString);


    PARSE_CLASS_MEMBER(configFile, slipCoefficient, "slipCoefficient", 0.0);
    ASSERT(slipCoefficient >= 0.0);

    PARSE_CLASS_MEMBER(configFile, hydrodynamicRadius, "hydrodynamicRadius", 1.0);
    ASSERT(hydrodynamicRadius > 0.0);
    ASSERT(hydrodynamicRadius <= 1.0);
    if (hydrodynamicRadius < 1.0) {
        cout << "Particle radius reduced for hydrodynamic interaction computations. Reduction factor: " << hydrodynamicRadius << endl;
    }


    // scaling rheological parameters and computing derived quantities
    fluidMaterial.initDynVisc /= unit.DynVisc;
    fluidMaterial.plasticVisc /= unit.DynVisc;
    fluidMaterial.yieldStress /= unit.Stress;
    fluidMaterial.particleDiameter /= unit.Length;
    fluidMaterial.rhod2 /= unit.Length * unit.Length * unit.Density;
    cout << "scaling " << fluidMaterial.particleDensity << " by " << unit.Density << " obtaining " << fluidMaterial.particleDensity / unit.Density << endl;
    fluidMaterial.particleDensity /= unit.Density;
    fluidMaterial.minimumPressure = 0.0;//1.0 * lbF.norm();
    //    PARSE_CLASS_MEMBER(lbmCfgFile, initDensity, "initDensity",1.0);
    //    initDensity/=unit.Density;
    // to avoid errors we set this to be just 1
    fluidMaterial.initDensity = 1.0;
    ASSERT(fluidMaterial.initDensity > 0.0);
}
void LBParams::LBShow() const {
    cout << "LATTICE CHARACTERISTICS (lattice units):" << endl;
    cout << "directions=" << lbmDirec << ";" << endl;
    cout << "time unit = " << unit.Time << "; length unit = " << unit.Length << "; density unit = " << unit.Density << ";" << endl;
    cout << "dynVisc unit = " << unit.DynVisc << "; kinVisc unit = " << unit.KinVisc << "; speed unit = " << unit.Speed << ";" << endl;
    cout << "xdim =" << lbSize[0] << "; ydim= " << lbSize[1] << "; zdim= " << lbSize[2] << ";" << endl;
    cout << "Init tau = " << (0.5 + 3.0 * fluidMaterial.initDynVisc) << ";init visc =" << fluidMaterial.initDynVisc << "; init density = " << fluidMaterial.initDensity << ";" << endl;
    cout << "Min tau = " << (0.5 + 3.0 * fluidMaterial.lbMinVisc) << "; Max tau = " << (0.5 + 3.0 * fluidMaterial.lbMaxVisc) << ";" << endl;
    switch (fluidMaterial.rheologyModel) {
        case BINGHAM:
        {
            cout << "Bingham rheology:" << endl;
            cout << "Plastic viscosity = " << fluidMaterial.plasticVisc << "(tau = " << (0.5 + 3.0 * fluidMaterial.plasticVisc) << " ); yield stress=" << fluidMaterial.yieldStress << ";" << endl;
            break;
        }
        case FRICTIONAL:
        {
            cout << "Simple frictional rheology:" << endl;
            cout << "Fluid friction coefficient = " << fluidMaterial.frictionCoefFluid << ";" << endl;
            cout << "Basal friction coefficient = " << fluidMaterial.basalFrictionCoefFluid << ";" << endl;
            break;
        }
        case VOELLMY:
        {
            cout << "Voellmy rheology:" << endl;
            cout << "Continuum particle diameter (NOT DEM) = " << fluidMaterial.particleDiameter << endl;
            cout << "Fluid friction coefficient = " << fluidMaterial.frictionCoefFluid << ";" << endl;
            cout << "Basal friction coefficient = " << fluidMaterial.basalFrictionCoefFluid << ";" << endl;
            break;
        }
        case BAGNOLD:
        {
            cout << "Bagnold rheology:" << endl;
            cout << "Continuum particle diameter (NOT DEM) = " << fluidMaterial.particleDiameter << endl;
            cout << "Continuum particle density (NOT DEM) = " << fluidMaterial.particleDensity << endl;
            break;
        }
        case MUI:
        {
            cout << "Rheology mu(I):" << endl;
            cout << "Fluid friction coefficient = " << fluidMaterial.frictionCoefFluid << "; friction interval = " << fluidMaterial.deltaFriction << endl;
            cout << "Basal friction coefficient = " << fluidMaterial.basalFrictionCoefFluid << ";" << endl;
            cout << "Basic inertial number = " << fluidMaterial.baseInertial << "; continuum particle diameter (NOT DEM) = " << fluidMaterial.particleDiameter << endl;
            cout << "Continuum particle density (NOT DEM) = " << fluidMaterial.particleDensity << endl;
            break;
        }
    }
    if (fluidMaterial.earthPressureCoeff!=1.0) {
        cout << "Earth pressure coefficient = " << fluidMaterial.earthPressureCoeff << endl;
}
    cout << "Initial Velocity=";
    initVelocity.show();
    cout << ";" << endl;
    cout << "F=";
    lbF.show();
    cout << endl;
    cout << "Rotation speed =";
    rotationSpeed.show();
    cout << ";" << endl;
    cout << "Rotation center =";
    rotationCenter.show();
    cout << ";" << endl;
    cout << "X(-)=" << boundary[0] << "; X(+)=" << boundary[1] << "; Y(-)=" << boundary[2] << "; Y(+)=" << boundary[3] << "; Z(-)=" << boundary[4] << "; Z(+)=" << boundary[5] << ";" << endl;
    cout << "PHISICAL CHARACTERISTICS (physical units)" << endl;
    cout << "xdim=" << (lbSize[0] - 2) * unit.Length << "; ydim= " << (lbSize[1] - 2) * unit.Length << "; zdim= " << (lbSize[2] - 2) * unit.Length << ";" << endl;
    cout << "u_sound=" << unit.Speed/ sqrt(3)<<endl;
    cout << "init kin viscosity=" << fluidMaterial.initDynVisc * unit.KinVisc << ", init dyn viscosity=" << fluidMaterial.initDynVisc * unit.DynVisc << "; init density=" << fluidMaterial.initDensity * unit.Density << ";" << endl;
    switch (fluidMaterial.rheologyModel) {
        case BINGHAM:
        {
            cout << "Bingham rheology:" << endl;
            cout << "Plastic kin viscosity=" << fluidMaterial.plasticVisc * unit.KinVisc << ", plastic dyn viscosity=" << fluidMaterial.plasticVisc * unit.DynVisc << "; yield stress=" << fluidMaterial.yieldStress * unit.Stress << ";" << endl;
            break;
        }
        case FRICTIONAL:
        {
            cout << "Simple frictional rheology:" << endl;
            cout << "Fluid friction coefficient = " << fluidMaterial.frictionCoefFluid << ";" << endl;
            cout << "Basal friction coefficient = " << fluidMaterial.basalFrictionCoefFluid << ";" << endl;
            break;
        }
        case VOELLMY:
        {
            cout << "Voellmy rheology:" << endl;
            cout << "Continuum particle diameter (NOT DEM) = " << fluidMaterial.particleDiameter * unit.Length << endl;
            cout << "Fluid friction coefficient = " << fluidMaterial.frictionCoefFluid << ";" << endl;
            cout << "Basal friction coefficient = " << fluidMaterial.basalFrictionCoefFluid << ";" << endl;
            break;
        }
        case BAGNOLD:
        {
            cout << "Bagnold rheology:" << endl;
            cout << "Continuum particle diameter (NOT DEM) = " << fluidMaterial.particleDiameter * unit.Length << endl;
            cout << "Continuum particle density (NOT DEM) = " << fluidMaterial.particleDensity * unit.Density << endl;
            break;
        }
        case MUI:
        {
            cout << "Rheology mu(I):" << endl;
            cout << "Fluid friction coefficient = " << fluidMaterial.frictionCoefFluid << "; friction interval = " << fluidMaterial.deltaFriction << endl;
            cout << "Basal friction coefficient = " << fluidMaterial.basalFrictionCoefFluid << ";" << endl;
            cout << "Basic inertial number = " << fluidMaterial.baseInertial << "; continuum particle diameter (NOT DEM) = " << fluidMaterial.particleDiameter * unit.Length << endl;
            cout << "Continuum particle density (NOT DEM) = " << fluidMaterial.particleDensity * unit.Density << endl;
            break;
        }
    }
    cout << "Min dyn viscosity = " << fluidMaterial.lbMinVisc * unit.DynVisc << "; Max dyn viscosity = " << fluidMaterial.lbMaxVisc * unit.DynVisc << ";" << endl;
    cout << "Initial Velocity=";
    const tVect initVelocityPhysic = initVelocity * unit.Speed;
    initVelocityPhysic.show();
    cout << ";" << endl;
    cout << "F=";
    const tVect FPhysic = lbF * unit.Accel;
    FPhysic.show();
    cout << ";" << endl;
    cout << "Rotation speed =";
    const tVect RotPhysic = rotationSpeed * unit.AngVel;
    RotPhysic.show();
    cout << ";" << endl;
    cout << "Rotation center =";
    const tVect RotCenterPhysic = rotationCenter * unit.Length;
    RotCenterPhysic.show();
    cout << ";" << endl;
}

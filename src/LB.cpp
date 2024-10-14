 
#include "LB.h"

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
 // PUBLIC FUNCTIONS
 //////////////////////////////////////////////////////////////////////////////////////////////////////*/

void LB::LBShow() const {
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

void LB::latticeDefinition() {
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
    ne[1] = shift[0];
    ne[2] = -shift[0];
    //
    ne[3] = shift[1];
    ne[4] = -shift[1];
    //
    ne[5] = shift[2];
    ne[6] = -shift[2];
    //
    ne[7] = shift[0] + shift[1];
    ne[8] = -shift[0] - shift[1];
    ne[9] = -shift[0] + shift[1];
    ne[10] = shift[0] - shift[1];
    //
    ne[11] = shift[1] + shift[2];
    ne[12] = -shift[1] - shift[2];
    ne[13] = -shift[1] + shift[2];
    ne[14] = shift[1] - shift[2];
    //
    ne[15] = shift[2] + shift[0];
    ne[16] = -shift[2] - shift[0];
    ne[17] = -shift[2] + shift[0];
    ne[18] = shift[2] - shift[0];

}

void LB::latticeBoltzmannGet(GetPot& configFile, GetPot& commandLine) {


    // conversion units //////////////////////////////////////////////////////////////////////////////
    // fluid and granular matter are solved in different measure units
    // for a reference, check Feng, Han, Owen, 2007
    // measure units for DEM solver are in the international system

    cout << "Getting LBM info from config file" << endl;

    // restart
    PARSE_CLASS_MEMBER(configFile, lbRestart, "restartFluid", false);
    if (lbRestart) {
        PARSE_CLASS_MEMBER(configFile, lbRestartFile, "fluidRestartFile", "");
    }
    // imposed volume
    PARSE_CLASS_MEMBER(configFile, imposeFluidVolume, "imposeFluidVolume", false);
    if (imposeFluidVolume) {
        PARSE_CLASS_MEMBER(configFile, imposedFluidVolume, "imposedFluidVolume", 0.0);
    }
    // increase volume
    PARSE_CLASS_MEMBER(configFile, increaseVolume, "increaseVolume", false);
    if (increaseVolume) {
        ASSERT(imposedFluidVolume == false);
        PARSE_CLASS_MEMBER(configFile, deltaVolume, "deltaVolume", 0.0);
        PARSE_CLASS_MEMBER(configFile, deltaTime, "deltaTime", 0.0);
    }
    if (imposeFluidVolume) {
        cout << "Fixed volume:  " << imposedFluidVolume << endl;
    } else if (increaseVolume) {
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
        PARSE_CLASS_MEMBER(configFile, lbTopographyFile, "topographyFile", "");
        PARSE_CLASS_MEMBER(configFile, translateTopographyX, "translateTopographyX", 0.0);
        PARSE_CLASS_MEMBER(configFile, translateTopographyY, "translateTopographyY", 0.0);
        PARSE_CLASS_MEMBER(configFile, translateTopographyZ, "translateTopographyZ", 0.0);
    }
    if (lbTopography && !lbRestart) {
        cout << "Using topography file for boundary and initial surface:  " << lbTopographyFile << endl;
    } else if (lbTopography && lbRestart) {
        cout << "Using topography file for boundary: " << lbTopographyFile << ", and restart file for fluid:  " << lbRestartFile << endl;
    } else if (lbRestart) {
        cout << "Computation from LB restart file " << lbRestartFile << endl;
    } else {
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
    lbSize.resize(lbmDim);
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
    lbPhysicalSize.resize(3);
    lbPhysicalSize[0] = double(lbSize[0] - 2) * unit.Length;
    lbPhysicalSize[1] = double(lbSize[1] - 2) * unit.Length;
    lbPhysicalSize[2] = double(lbSize[2] - 2) * unit.Length;

    // domain inside boundaries
    lbInnerPhysicalSize.resize(3);
    lbInnerPhysicalSize[0] = double(lbSize[0]) * unit.Length;
    lbInnerPhysicalSize[1] = double(lbSize[1]) * unit.Length;
    lbInnerPhysicalSize[2] = double(lbSize[2]) * unit.Length;

    // location of the domain boundary
    lbBoundaryLocation.resize(6);
    //    lbBoundaryLocation[0] = tVect(0.5 * unit.Length, 0.0, 0.0);
    //    lbBoundaryLocation[1] = tVect(double(lbSize[0] - 1.5) * unit.Length, 0.0, 0.0);
    //    lbBoundaryLocation[2] = tVect(0.0, 0.5 * unit.Length, 0.0);
    //    lbBoundaryLocation[3] = tVect(0.0, double(lbSize[1] - 1.5) * unit.Length, 0.0);
    //    lbBoundaryLocation[4] = tVect(0.0, 0.0, 0.5 * unit.Length);
    //    lbBoundaryLocation[5] = tVect(0.0, 0.0, double(lbSize[2] - 1.5) * unit.Length);
    lbBoundaryLocation[0] = tVect(0.0, 0.0, 0.0);
    lbBoundaryLocation[1] = tVect(double(lbSize[0] - 2) * unit.Length, 0.0, 0.0);
    lbBoundaryLocation[2] = tVect(0.0, 0.0, 0.0);
    lbBoundaryLocation[3] = tVect(0.0, double(lbSize[1] - 2) * unit.Length, 0.0);
    lbBoundaryLocation[4] = tVect(0.0, 0.0, 0.0);
    lbBoundaryLocation[5] = tVect(0.0, 0.0, double(lbSize[2] - 2) * unit.Length);
    
//    lbBoundaryLocation[0] = tVect(0.5 * unit.Length, 0.0, 0.0);
//    lbBoundaryLocation[1] = tVect(double(lbSize[0] - 2) * unit.Length, 0.0, 0.0);
//    lbBoundaryLocation[2] = tVect(0.0, 0.5 * unit.Length, 0.0);
//    lbBoundaryLocation[3] = tVect(0.0, double(lbSize[1] - 2) * unit.Length, 0.0);
//    lbBoundaryLocation[4] = tVect(0.0, 0.0, 0.5 * unit.Length);
//    lbBoundaryLocation[5] = tVect(0.0, 0.0, double(lbSize[2] - 2) * unit.Length);


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
            } else {
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
            } else {
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
            } else {
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
                unit.Time = unit.Density * unit.Length * unit.Length / fluidMaterial.plasticVisc * (1.0 / 3.0)*(tauOpt - 0.5);
                break;
            }
            default:
            {
                unit.Time = unit.Density * unit.Length * unit.Length / fluidMaterial.initDynVisc * (1.0 / 3.0)*(tauOpt - 0.5);
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


    boundary.resize(6);

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
    if (hydrodynamicRadius<1.0) {
        cout<< "Particle radius reduced for hydrodynamic interaction computations. Reduction factor: "<<hydrodynamicRadius<<endl;
    }


    // scaling rheological parameters and computing derived quantities
    fluidMaterial.initDynVisc /= unit.DynVisc;
    fluidMaterial.plasticVisc /= unit.DynVisc;
    fluidMaterial.yieldStress /= unit.Stress;
    fluidMaterial.particleDiameter /= unit.Length;
    fluidMaterial.rhod2 /= unit.Length * unit.Length * unit.Density;
    cout<<"scaling "<<fluidMaterial.particleDensity<<" by "<<unit.Density<<" obtaining "<<fluidMaterial.particleDensity/unit.Density<<endl;
    fluidMaterial.particleDensity /= unit.Density;
    fluidMaterial.minimumPressure = 0.0;//1.0 * lbF.norm();
    //    PARSE_CLASS_MEMBER(lbmCfgFile, initDensity, "initDensity",1.0);
    //    initDensity/=unit.Density;
    // to avoid errors we set this to be just 1
    fluidMaterial.initDensity = 1.0;
    ASSERT(fluidMaterial.initDensity > 0.0);

}

void LB::latticeBolzmannInit(cylinderList& cylinders, wallList& walls, particleList& particles, objectList& objects, const bool externalSolveCoriolis, const bool externalSolveCentrifugal) {
    //  Lattice Boltzmann initialization steps
    
    // switchers for apparent accelerations
    solveCoriolis = externalSolveCoriolis;
    solveCentrifugal = externalSolveCentrifugal;

    // first comes the initialization of the data structures
    initializeNodes();

    // application of lattice boundaries
    initializeLatticeBoundaries();

    // then the initial node type must be identified for every node (if not specified, it is already Fluid)
    initializeTypes(walls, cylinders, objects);

    ifstream fluidFileID;
    if (lbRestart) {
        // open fluid restart file
        fluidFileID.open(lbRestartFile.c_str(), ios::in);
        ASSERT(fluidFileID.is_open());
        // check if the restart file size is ok
        unsigned int restartX, restartY, restartZ, restartNodes;
        fluidFileID>>restartX;
        fluidFileID>>restartY;
        fluidFileID>>restartZ;
        fluidFileID>>restartNodes;
        ASSERT(restartX == lbSize[0]);
        ASSERT(restartY == lbSize[1]);
        ASSERT(restartZ == lbSize[2]);
        // read active nodes from file and generate
        restartInterface(fluidFileID, restartNodes);
    } else {
        // initialize interface
        initializeInterface(particles.size());
        // initialize variables for active nodes
        initializeVariables();
    }

    // initialize variables for wall nodes
    initializeWalls(walls, cylinders, objects);

    // initializing curved properties
    initializeCurved(cylinders);

    // initialize node properties and lists
    initializeLists();

    //create active list and checks for errors in the list
    cleanLists();

    // application of particle initial position
    initializeParticleBoundaries(particles);

    // in case mass needs to be kept constant, compute it here
    totalMass = 0.0;
    if (imposeFluidVolume) {
        // volume and mass is the same in lattice units
        totalMass = imposedFluidVolume / unit.Volume;
    } else {
        switch (problemName) {
            case DRUM:
            {
                totalMass = fluidMass / unit.Mass;
                break;
            }
            case STAVA:
            {
                totalMass = 200000.0 / unit.Volume;
                break;
            }
            default:
            {
                for (nodeList::iterator it = activeNodes.begin(); it != activeNodes.end(); ++it) {
                    const node* nodeHere = *it;
                    if (!nodeHere->isInsideParticle()) {
                        totalMass += nodeHere->mass;
                    }
                }
                break;
            }
        }
    }
    if (increaseVolume) {
        deltaVolume /= unit.Volume;
        deltaTime /= unit.Time;
    }
    cout << "Done with initialization" << endl;
}

void LB::latticeBolzmannStep(elmtList& elmts, particleList& particles, wallList& walls, objectList& objects) {

    // Lattice Boltzmann core steps

    // measure time for performance check (begin)
    startLBStep = std::chrono::steady_clock::now();

    //create active list and checks for errors in the list
    cleanLists();

    // initializing the elements forces (lattice units)
#pragma omp parallel for
    for (int n = 0; n < elmts.size(); ++n) {
        //initializing this time step hydrodynamic force
        elmts[n].FHydro = tVect(0.0, 0.0, 0.0);
        elmts[n].MHydro = tVect(0.0, 0.0, 0.0);
        // initializing the fluid mass for buoyancy
        elmts[n].fluidVolume = 0.0;
    }

    if (!forceField) {
        lbF.reset();
    }

#pragma omp parallel for
    for (int it = 0; it < activeNodes.size(); ++it) {
        node* nodeHere = activeNodes[it];

        // reconstruction of macroscopic variables from microscopic distribution
        // this step is necessary to proceed to the collision step
        nodeHere->reconstruct();

        // compute interaction forces
        if (elmts.size()) {
            computeHydroForces(nodeHere, elmts, particles);
        }

        //collision operator
        collision(nodeHere);
    }

    // shifting elements forces and torques to physical units
#pragma omp parallel for
    for (int n = 0; n < elmts.size(); ++n) {

        // adding buoyancy
        const tVect buoyancy = (-1.0) * elmts[n].fluidVolume*lbF;

        //        elmts[n].FHydro+=buoyancy;
        elmts[n].FHydro *= unit.Force;
        elmts[n].MHydro *= unit.Torque;
        elmts[n].fluidVolume *= unit.Volume;
    }

    //streaming operator
    streaming(walls, objects);
    
    // measure time for performance check (end)
    endLBStep = std::chrono::steady_clock::now();

}

void LB::latticeBoltzmannFreeSurfaceStep() {

    // free-surface  management functions

    // measure time for performance check (begin)
    startFreeSurfaceStep = std::chrono::steady_clock::now();

    // in case mass needs to be kept constant, call enforcing function here
    if (imposeFluidVolume) {
        enforceMassConservation();
    } else if (increaseVolume) {
        if (double(time) < deltaTime)
            redistributeMass(deltaVolume / deltaTime);
    } else {
        switch (problemName) {
            case DRUM:
            {
                enforceMassConservation();
                break;
            }
            case STAVA:
            {
                enforceMassConservation();
                break;
            }
        }
    }
    // mass and free surface update

    updateMass();

    updateInterface();

    cleanLists();

    // measure time for performance check (end)
    endFreeSurfaceStep = std::chrono::steady_clock::now();

}
//

void LB::latticeBoltzmannCouplingStep(bool& newNeighborList, elmtList& elmts, particleList& particles) {
    // identifies which nodes need to have an update due to particle movement
    // the complexity arises from trying to keep the scaling as close to linear as possible
    // maybe the best idea is to do this in two step:
    // 1) the first is to check for new active nodes and initialise them
    // 2) the second is checking for new solid nodes.
    // this automatically check also for the hideous case of particle to particle double transition

    // measure time for performance check (begin)
    startCouplingStep = std::chrono::steady_clock::now();

    // first we check if a new neighbour table has been defined. In that case, the indexing needs to be reinitialised

    if (newNeighborList) {
        cout<<endl<<"New neighbor list"<<endl;
        initializeParticleBoundaries(particles);
        newNeighborList = false;
    }
    

    // declaring list for mutant nodes
    nodeList newPopUpNodes, newSolidNodes;
    // emptying list of mutant nodes
    newPopUpNodes.clear();
    newSolidNodes.clear();
    // double massSurplus=0.0;

    // SOLID TO ACTIVE CHECK
    findNewActive(newPopUpNodes, elmts, particles);

    //    solidToActive(newPopUpNodes, elmts, massSurplus);

    // ACTIVE TO SOLID CHECK
    findNewSolid(newSolidNodes, elmts, particles);
    
    if (freeSurface && time % 1 == 0) {
        checkNewInterfaceParticles(elmts,particles);
    }

    //activeToSolid(newSolidNodes, elmts, massSurplus);

    // redistributing extra mass due to popping of nodes
    // redistributeMass(massSurplus);
    // measure time for performance check (begin)
    endCouplingStep = std::chrono::steady_clock::now();

}


/*//////////////////////////////////////////////////////////////////////////////////////////////////////
 // PRIVATE FUNCTIONS
 //////////////////////////////////////////////////////////////////////////////////////////////////////*/

// node creation and deletion

void LB::generateNode(const unsigned int& coord, const types& typeHere) {

    //initialize node
    node dummyNode;
    // set coordinate
    dummyNode.coord = coord;

    // set type
    dummyNode.setType(typeHere);
    dummyNode.setOutsideParticle();
    dummyNode.age=0.0;

    //add it to the map;
    nodes.insert(std::pair<unsigned int, node>(coord, dummyNode));

    // get its pointer
    node* nodeHere = &nodes[coord];

    // find neighbor indices
    unsigned int neighborCoord[lbmDirec];
    findNeighbors(neighborCoord, nodeHere);

    nodeHere->basal = false;
    
    // set centrifugal acceleration
    nodeHere->centrifugalForce = computeCentrifugal(getPosition(nodeHere->coord), rotationCenter, rotationSpeed);

    // assign neighbor nodes
    for (int j = 1; j < lbmDirec; ++j) {
        // linearized coordinate of neighbor nodes
        const unsigned int link = neighborCoord[j];
        // check if node at that location exists
        if (nodes.count(link) != 0) {
            // if exists, get pointer to that node
            node* linkNode = &nodes[link];
            // assign neighbor for local node
            nodeHere->d[j] = linkNode;
            // if neighbor node is also active, link it to local node
            if (nodeHere->isActive()) {
                linkNode->d[opp[j]] = nodeHere;
                // if the neighbor is a curved wall, set parameters accordingly
                if (linkNode->isTopography()) {

                    if (nodeHere->curved == 0) {
                        nodeHere->curved = new curve;
                    }
                    // set curved
                    const tVect nodePosHere = unit.Length * getPosition(coord);
                    // xf - xw
                    const double topographyDistance = 1.0 * lbTop.directionalDistance(nodePosHere, vDirec[j]) / unit.Length;
                    // wall normal
                    nodeHere->curved->wallNormal = lbTop.surfaceNormal(nodePosHere);
                    //cout << topographyDistance << endl;
                    const double deltaHere = topographyDistance / vNorm[j];
                    nodeHere->curved->delta[j] = std::min(0.99, std::max(0.01, deltaHere));
                    nodeHere->curved->computeCoefficients();
                }
                if (linkNode->isWall()) {
                    nodeHere->basal = true;
                }
            }
        } else {
            nodeHere->d[j] = nullptr;
        }
    }
}

void LB::eraseNode(const unsigned int& coord) {

    // find neighbor indices

    node* nodeHere = &nodes[coord];
    unsigned int neighborCoord[lbmDirec];
    findNeighbors(neighborCoord, nodeHere);

    //remove from map;
    delete nodeHere->curved;
    //delete nodeHere->d;
    nodes.erase(nodeHere->coord);

    // assign neighbor nodes
    for (int j = 1; j < lbmDirec; ++j) {
        // linearized coordinate of neighbor nodes
        const unsigned int link = neighborCoord[j];
        // check if node at that location exists
        if (nodes.count(link) != 0) {
            // if exists, get pointer to that node
            node* linkNode = &nodes[link];
            // if neighbor node is active, remove link to local node
            if (linkNode->isActive()) {
                linkNode->d[opp[j]] = nullptr;
            }
        }
    }


}



// initialization functions

void LB::initializeTypes(wallList& walls, cylinderList& cylinders, objectList& objects) {

    // application of solid walls
    initializeWallBoundaries(walls);
    // application of solid cylinders
    initializeCylinderBoundaries(cylinders);
    // application of objects
    initializeObjectBoundaries(objects);
    // initializing topography if one is present
    initializeTopography();

}

void LB::initializeNodes() {

    cout << "Initializing nodes containers and types" << endl;

    // total number of nodes
    totPossibleNodes = lbSize[0] * lbSize[1] * lbSize[2];

    // map with nodes   
    nodes.clear();

}

void LB::initializeLatticeBoundaries() {
    // assign boundary characteristic to nodes (see class)
    // if not differently defined, type is 0 (fluid)

    // BOUNDARY CONDITIONS ///////////////////////////
    // solid boundary wins over all in corners, where more than 1 bc is defined
    cout << "Initializing boundaries" << endl;

    unsigned int indexHere = 0;

    // XY
    for (int x = 0; x < lbSize[0]; ++x) {
        for (int y = 0; y < lbSize[1]; ++y) {
            // bottom
            indexHere = getIndex(x, y, 0);
            if (nodes.count(indexHere) == 0) {
                generateNode(indexHere, boundary[4]);
            }
            // top
            indexHere = getIndex(x, y, lbSize[2] - 1);
            if (nodes.count(indexHere) == 0) {
                generateNode(indexHere, boundary[5]);
            }
        }
    }

    // YZ
    for (int y = 0; y < lbSize[1]; ++y) {
        for (int z = 0; z < lbSize[2]; ++z) {
            // bottom
            indexHere = getIndex(0, y, z);
            if (nodes.count(indexHere) == 0) {
                generateNode(indexHere, boundary[0]);
            }
            // top
            indexHere = getIndex(lbSize[0] - 1, y, z);
            if (nodes.count(indexHere) == 0) {
                generateNode(indexHere, boundary[1]);
            }
        }
    }

    // ZX
    for (int z = 0; z < lbSize[2]; ++z) {
        for (int x = 0; x < lbSize[0]; ++x) {
            // bottom
            indexHere = getIndex(x, 0, z);
            if (nodes.count(indexHere) == 0) {
                generateNode(indexHere, boundary[2]);
            }
            // top
            indexHere = getIndex(x, lbSize[1] - 1, z);
            if (nodes.count(indexHere) == 0) {
                generateNode(indexHere, boundary[3]);
            }
        }
    }
}

void LB::findNeighbors(unsigned int neighborCoord[], const node* nodeHere) {
    // assign boundary characteristic to nodes (see class)
    // if not differently defined, type is 0 (fluid)

    const unsigned int it = nodeHere->coord;

    for (int j = 1; j < lbmDirec; ++j) {
        neighborCoord[j] = it + ne[j];
    }

    // BOUNDARY CONDITIONS ///////////////////////////
    // nodes on the boundary have no neighbors
    if (nodeHere->isWall()) {
        const unsigned int xHere = getX(it);
        if (xHere == 0) {
            for (int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Xp) < 0.0) {
                    neighborCoord[j] = it;
                }
            }
        }
        if (xHere == lbSize[0] - 1) {
            for (int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Xp) > 0.0) {
                    neighborCoord[j] = it;
                }
            }
        }
        const unsigned int yHere = getY(it);
        if (yHere == 0) {
            for (int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Yp) < 0.0) {
                    neighborCoord[j] = it;
                }
            }
        }
        if (yHere == lbSize[1] - 1) {
            for (int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Yp) > 0.0) {
                    neighborCoord[j] = it;
                }
            }
        }
        const unsigned int zHere = getZ(it);
        if (zHere == 0) {
            for (int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Zp) < 0.0) {
                    neighborCoord[j] = it;
                }
            }
        }
        if (zHere == lbSize[2] - 1) {
            for (int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Zp) > 0.0) {
                    neighborCoord[j] = it;
                }
            }
        }
    }// PERIODICITY ////////////////////////////////////////////////
        // assigning periodicity conditions (this needs to be done after applying boundary conditions)
        // runs through free cells and identifies neighboring cells. If neighbor cell is
        // a special cell (periodic) then the proper neighboring condition is applied
        // calculates the effect of periodicity
    else if (nodeHere->isActive()) {
        // neighboring and periodicity vector for boundary update
        unsigned int pbc[lbmDirec];
        for (int j = 0; j < lbmDirec; ++j) {
            pbc[j] = 0;
        }
        if (nodes.count(neighborCoord[1]) != 0) {
            if (nodes[neighborCoord[1]].isPeriodic()) {
                for (int j = 1; j < lbmDirec; ++j) {
                    if (v[j].dot(Xp) > 0.0) {
                        pbc[j] -= domain[0];
                    }
                }
            }
        }
        if (nodes.count(neighborCoord[2]) != 0) {
            if (nodes[neighborCoord[2]].isPeriodic()) {
                for (int j = 1; j < lbmDirec; ++j) {
                    if (v[j].dot(Xp) < 0.0) {
                        pbc[j] += domain[0];
                    }
                }
            }
        }
        if (nodes.count(neighborCoord[3]) != 0) {
            if (nodes[neighborCoord[3]].isPeriodic()) {
                for (int j = 1; j < lbmDirec; ++j) {
                    if (v[j].dot(Yp) > 0.0) {
                        pbc[j] -= domain[1];
                    }
                }
            }
        }
        if (nodes.count(neighborCoord[4]) != 0) {
            if (nodes[neighborCoord[4]].isPeriodic()) {
                for (int j = 1; j < lbmDirec; ++j) {
                    if (v[j].dot(Yp) < 0.0) {
                        pbc[j] += domain[1];
                    }
                }
            }
        }
        if (nodes.count(neighborCoord[5]) != 0) {
            if (nodes[neighborCoord[5]].isPeriodic()) {
                for (int j = 1; j < lbmDirec; ++j) {
                    if (v[j].dot(Zp) > 0.0) {
                        pbc[j] -= domain[2];
                    }
                }
            }
        }
        if (nodes.count(neighborCoord[6]) != 0) {
            if (nodes[neighborCoord[6]].isPeriodic()) {
                for (int j = 1; j < lbmDirec; ++j) {
                    if (v[j].dot(Zp) < 0.0) {
                        pbc[j] += domain[2];
                    }
                }
            }
        }

        // apply periodicity
        for (int j = 1; j < lbmDirec; ++j) {
            neighborCoord[j] += pbc[j];
            //cout << "pbc " << j << " " << neighborCoord[j] << endl;
        }

    }
    //
    //    for (int j = 1; j < lbmDirec; ++j) {
    //        cout << "it " << it << " " << getX(it) << " " << getY(it) << " " << getZ(it) << " " << nodeHere->getType() << " " << j << " " << neighborCoord[j] << endl;
    //    }

}

void LB::initializeWallBoundaries(wallList& walls) {

    const double wallThickness = 2.0 * unit.Length;
    // SOLID WALLS ////////////////////////
    cout << "Initializing solid walls" << endl;
    for (int iw = 0; iw < walls.size(); ++iw) {
        const tVect convertedWallp = walls[iw].p / unit.Length;
        const tVect normHere = walls[iw].n;
        const unsigned short int indexHere = walls[iw].index;
        const bool slipHere = walls[iw].slip;
        const bool movingHere = walls[iw].moving;
#pragma omp parallel for
        for (int it = 0; it < totPossibleNodes; ++it) {
            // check if the node is solid
            bool nodeInsideWall = false;
            // all walls have max thickness 2 nodes
            const double wallDistance = getPosition(it).distance2Plane(convertedWallp, normHere);
            if (wallDistance>-2.0 && wallDistance < 0.0) {
                nodeInsideWall = true;
                //check for borders in limted walls
                if (walls[iw].limited) {
                    const double xHere = getPositionX(it) * unit.Length;
                    const double yHere = getPositionY(it) * unit.Length;
                    const double zHere = getPositionZ(it) * unit.Length;
                    // check if beyond limits
                    if (xHere < walls[iw].xMin || xHere > walls[iw].xMax ||
                            yHere < walls[iw].yMin || yHere > walls[iw].yMax ||
                            zHere < walls[iw].zMin || zHere > walls[iw].zMax) {
                        nodeInsideWall = false;

                    }
                }
            }
            if (nodeInsideWall) {
                // generate node (tentatively as static wall)
                generateNode(it, STAT_WALL);
                // setting solidIndex
                nodes[it].setSolidIndex(indexHere);
                // setting type: 5-6=slip, 7-8=no-slip
                if (slipHere) {
                    // setting type for slip: 5=static, 6=moving
                    if (movingHere) {
                        nodes[it].setSlipDynWall();
                    } else {
                        nodes[it].setSlipStatWall();
                    }
                } else {
                    // setting type for no-slip: 7=static, 8=moving
                    if (movingHere) {
                        nodes[it].setDynWall();
                    } else {
                        nodes[it].setStatWall();
                    }
                }

            }
        }
    }
}

void LB::initializeObjectBoundaries(objectList& objects) {

    // SOLID WALLS ////////////////////////
    cout << "Initializing objects" << endl;
    for (int io = 0; io < objects.size(); ++io) {
        const tVect convertedPosition = objects[io].x0 / unit.Length;
        const double convertedRadius = objects[io].r / unit.Length;
        const unsigned int indexHere = objects[io].index;

#pragma omp parallel for
        for (int it = 0; it < totPossibleNodes; ++it) {
            const tVect nodePosition = getPosition(it);
            if (nodePosition.insideSphere(convertedPosition, convertedRadius)) {// setting solidIndex
                // generate node (tentatively as static wall)
                generateNode(it, OBJ);
                // setting solidIndex
                nodes[it].setSolidIndex(indexHere);
            }
        }
    }
}

void LB::initializeCylinderBoundaries(cylinderList& cylinders) {

    // SOLID CYLINDERS ////////////////////////
    cout << "Initializing solid cylinders" << endl;
    for (int ic = 0; ic < cylinders.size(); ++ic) {

        const tVect convertedCylinderp1 = cylinders[ic].p1 / unit.Length;
        const tVect naxesHere = cylinders[ic].naxes;
        const double convertedRadius = cylinders[ic].R / unit.Length;
        const unsigned int indexHere = cylinders[ic].index;
        const bool slipHere = cylinders[ic].slip;
        const bool movingHere = cylinders[ic].moving;
#pragma omp parallel for
        for (int it = 0; it < totPossibleNodes; ++it) {
            bool nodeInsideCylinder = false;
            // creating solid cells
            bool isOutside=getPosition(it).insideCylinder(convertedCylinderp1, naxesHere, convertedRadius, convertedRadius + 3.0);
            bool isInside=getPosition(it).insideCylinder(convertedCylinderp1, naxesHere, max(convertedRadius - 3.0, 0.0), convertedRadius);
            if ((cylinders[ic].type == FULL && isInside) ||
                (cylinders[ic].type == EMPTY && isOutside)) {
                nodeInsideCylinder = true;
                //check for borders in limted walls
                if (cylinders[ic].limited) {
                    const double xHere = getPositionX(it) * unit.Length;
                    const double yHere = getPositionY(it) * unit.Length;
                    const double zHere = getPositionZ(it) * unit.Length;
                    // check if beyond limits
                    if (xHere < cylinders[ic].xMin || xHere > cylinders[ic].xMax ||
                            yHere < cylinders[ic].yMin || yHere > cylinders[ic].yMax ||
                            zHere < cylinders[ic].zMin || zHere > cylinders[ic].zMax) {
                        nodeInsideCylinder = false;

                    }
                }
                if (nodeInsideCylinder) {
                    // tentatively static
                    generateNode(it, STAT_WALL);
                    // setting solidIndex
                    nodes[it].setSolidIndex(indexHere);
                    // setting type: 5-6=slip, 7-8=no-slip
                    if (slipHere) {
                        // setting type for slip: 5=static, 6=moving
                        if (movingHere) {
                            nodes[it].setSlipDynWall();
                        } else {
                            nodes[it].setSlipStatWall();
                        }
                    } else {
                        // setting type for no-slip: 7=static, 8=moving
                        if (movingHere) {
                            nodes[it].setDynWall();
                        } else {
                            nodes[it].setStatWall();
                        }
                    }
                }
            }
        }
    }
}

void LB::initializeTopography() {

    const double surfaceThickness = 1.75 * unit.Length;

    // TOPOGRAPHY ////////////////////////
    if (lbTopography) {
        cout << "Setting topography" << endl;
        lbTop.readFromFile(lbTopographyFile, translateTopographyX, translateTopographyY, translateTopographyZ);
        lbTop.show();
        // check if topography grid contains the fluid domain
        ASSERT(lbTop.coordX[0] < unit.Length);
        ASSERT(lbTop.coordY[0] < unit.Length);

        cout << "lbTop.coordX[lbTop.sizeX - 1]=" << lbTop.coordX[lbTop.sizeX - 1] << endl;
        cout << "lbSize[0]) * unit.Length=" << lbSize[0] * unit.Length << endl;
        ASSERT(lbTop.coordX[lbTop.sizeX - 1] > lbSize[0] * unit.Length);
        cout << "lbTop.coordY[lbTop.sizeY - 1]=" << lbTop.coordY[lbTop.sizeY - 1] << endl;
        cout << "lbSize[1]) * unit.Length=" << lbSize[1] * unit.Length << endl;
        ASSERT(lbTop.coordY[lbTop.sizeY - 1] > lbSize[1] * unit.Length);


        for (int ix = 1; ix < lbSize[0] - 1; ++ix) {
            for (int iy = 1; iy < lbSize[1] - 1; ++iy) {
#pragma omp parallel for
                for (int iz = 1; iz < lbSize[2] - 1; ++iz) {

                    //for (int it = 0; it < totPossibleNodes; ++it) {
                    // control is done in real coordinates
                    //            const unsigned int nodeX = getX(it);
                    //            if (nodeX != 0 && nodeX != lbSize[0]) {
                    //                const unsigned int nodeY = getY(it);
                    //                if (nodeY != 0 && nodeY != lbSize[1]) {
                    //                    const unsigned int nodeZ = getZ(it);
                    //                    if (nodeZ != 0 && nodeZ != lbSize[2]) {

                    const tVect nodePosition = tVect(ix, iy, iz) * unit.Length;
                    const double distanceFromTopography = lbTop.distance(nodePosition);

                    //cout<<distanceFromTopography<<endl;
                    if (distanceFromTopography < 0.0 && distanceFromTopography>-1.0 * surfaceThickness) {// setting solidIndex
                        const unsigned int it = ix + iy * lbSize[0] + iz * lbSize[0] * lbSize[1];
#pragma omp critical
                        {
                            generateNode(it, STAT_WALL);
                        }
                        // setting type: 5-6=slip, 7-8=no-slip
                        nodes[it].setTopography();
                        // initializing info for curved (in this case inclined) boundary
                    }
                }
            }
        }
    }

}

void LB::initializeCurved(cylinderList& cylinders) {
    cout << "Initializing curved boundaries" << endl;
    for (nodeList::iterator iw = wallNodes.begin(); iw != wallNodes.end(); ++iw) {
        node* wallNodeHere = *iw;
        if (wallNodeHere->isCylinderWall()) {
            wallNodeHere->curved = new curve;
            const unsigned int index = wallNodeHere->coord;
            const tVect nodePos = unit.Length * getPosition(index);
            for (int j = 1; j < lbmDirec; ++j) {
                wallNodeHere->curved->delta[j] = 1.0 - cylinders[0].segmentIntercept(nodePos, unit.Length * v[j]);
                wallNodeHere->curved->computeCoefficients();
            }
        }
        //        else if (wallNodeHere->isTopography()) {
        //            wallNodeHere->curved = new curve;
        //            unsigned int neighborCoord[lbmDirec];
        //            findNeighbors(wallNodeHere->coord, neighborCoord, wallNodeHere);
        //            for (int j = 1; j < lbmDirec; ++j) {
        //                const unsigned int neighborIndex = neighborCoord[j];
        //                bool neighIsWall=(nodes.find(neighborIndex)!=nodes.end()
        //                if (!(nodes.find(neighborIndex)!=nodes.end()) {
        //                    const tVect neighborNodePos = unit.Length * getPosition(neighborIndex);
        //                    const double topographyDistance = lbTop.directionalDistance(neighborNodePos, vDirec[opp[j]]) / unit.Length;
        //                    //cout << topographyDistance << endl;
        //                    const double deltaHere = (vNorm[opp[j]] - topographyDistance) / vNorm[opp[j]];
        //                    cout<<topographyDistance<<endl;
        //                    wallNodeHere->curved->delta[j] = std::min(0.9, std::max(0.01, deltaHere));
        //                    wallNodeHere->curved->computeCoefficients();
        //                }
        //                else {
        //                    wallNodeHere->curved->delta[j] = 0.0;
        //                }
        //            }
        //        }
    }
}

void LB::setTopographySurface() {

#pragma omp parallel for
    for (int it = 0; it < totPossibleNodes; ++it) {
        if (nodes.find(it) == nodes.end()) {
            // control is done in real coordinates
            const tVect nodePosition = getPosition(it) * unit.Length;
            const double surfaceIsoparameterHere = lbTop.surfaceIsoparameter(nodePosition);
            //cout<<distanceFromTopography<<endl;
            if (surfaceIsoparameterHere > 0.0 && surfaceIsoparameterHere <= 1.0) {// setting solidIndex
#pragma omp critical
                {
                    generateNode(it, LIQUID);
                }
            }
        }
    }
}

void LB::initializeInterface(double totParticles) {
    // creates an interface electing interface cells from active cells
    // defining a surface

    cout << "Initializing interface" << endl;

    if (lbTopographySurface) {
        cout << "Initializing interface as set in topography file" << endl;
        setTopographySurface();
    } else {
        switch (problemName) {
            case HK_LARGE:
            {
                cout << "Initializing interface as in Hong Kong large scale flume - initial neight controlled  via config file (largeFlumeFlowLevel)" << endl;
                unsigned int yMax = int((largeFlumeFlowLevel / unit.Length));
                unsigned int xMax = int((10.0 / unit.Length));
                // release plane
                wall dummyWall_surface;
                dummyWall_surface.p = tVect(10.0, largeFlumeFlowLevel, 0.0) / unit.Length;
                dummyWall_surface.n = tVect(-0.342020, 0.939693, 0.0);
                wall dummyWall_bottom;
                dummyWall_bottom.p = tVect(10.0, 0.0, 0.0) / unit.Length;
                dummyWall_bottom.n = tVect(0.17365, 0.98481, 0.0);

                for (int it = 0; it < totPossibleNodes; ++it) {
                    if (nodes.count(it) == 0) {
                        const tVect nodePosition = getPosition(it);
                        // creating fluid cells
                        if (nodePosition.x < xMax && nodePosition.y < yMax) {
                            const double surface_distance = dummyWall_surface.dist(nodePosition);
                            const double bottom_distance = dummyWall_bottom.dist(nodePosition);
                            if (surface_distance < 0.0 && bottom_distance > 0.0) {
                                generateNode(it, LIQUID);
                            }
                        }
                    }
                }

                break;
            }
            case KELVIN:
            {
                
                const double totVolume=0.011;
                const double width=0.2;
                const double base=0.45;
                const double triangleHeight=base*tan(26*M_PI/180);
                const double triangleVolume=0.5*triangleHeight*base*width;
                const double rectangleVolume=totVolume-triangleVolume;
                const double rectangleHeight=rectangleVolume/base/width;
                const double totHeight=rectangleHeight+triangleHeight;
                
                cout << "Initializing interface as in Hong Kong small scale flume. Water level = "<<triangleHeight<<" + " <<rectangleHeight<<" = "<< totHeight << endl;
                
                unsigned int yMax = int((totHeight / unit.Length));
                unsigned int xMax = int((base / unit.Length));
                // release plane
                wall dummyWall_surface;
                dummyWall_surface.p = tVect(base, totHeight, 0.0) / unit.Length;
                dummyWall_surface.n = tVect(-0.43837, 0.898794, 0.0);

                for (int it = 0; it < totPossibleNodes; ++it) {
                    if (nodes.count(it) == 0) {
                        const tVect nodePosition = getPosition(it);
                        // creating fluid cells
                        if (nodePosition.x < xMax && nodePosition.y < yMax) {
                            const double surface_distance = dummyWall_surface.dist(nodePosition);
                            if (surface_distance < 0.0) {
                                generateNode(it, LIQUID);
                            }
                        }
                    }
                }

                break;
            }
            case NET:
            case BARRIER:
            {
                cout << "Initializing interface as in flexible barrier case study" << endl;
                // BOX BEFORE THE NET
                unsigned int xMin = int((avalanchePosit / unit.Length));
                unsigned int xMax = lbSize[0];
                unsigned int yMin = 0.0;
                unsigned int yMax = lbSize[1];
                unsigned int zMin = 0.0;
                unsigned int zMax = int((4.0 / unit.Length));
                for (int it = 0; it < totPossibleNodes; ++it) {
                    if (nodes.count(it) == 0) {
                        // creating fluid cells
                        if ((getPositionX(it) > xMin) && (getPositionX(it) < xMax) && (getPositionY(it) > yMin) && (getPositionY(it) < yMax) && (getPositionZ(it) > zMin) && (getPositionZ(it) < zMax)) {
                            generateNode(it, LIQUID);
                        }
                    }
                }
                // PLANE
                tVect p, n;
                p = tVect(int((avalanchePosit / unit.Length)), 0.0, 0.0);
                n = tVect(-1.0, 0.0, 1.0);
                for (int it = 0; it < totPossibleNodes; ++it) {
                    if (nodes.count(it) == 0) {
                        // creating fluid cells
                        const tVect dist = getPosition(it) - p;
                        if (dist.dot(n) < 0.0) {
                            generateNode(it, LIQUID);
                        }
                    }
                }
                break;
            }
            case DRUM:
            {
                cout << "Initializing interface as in drum case study" << endl;
                // curved shape
                cylinder dummyCylinder;
                dummyCylinder.p1 = tVect(0.45, 0.0, -1.0 + 0.0 + 0.04); //tVect(0.7+0.35,0.0,-1.0+0.0);
                dummyCylinder.p2 = tVect(0.45, 1.0, -1.0 + 0.0 + 0.04); //tVect(0.7+0.35,0.0,-1.0+0.0);
                const double fluidRadius = 1.01 + 0.00035298 * sqrt(7891604 - 5713 * (1389.9 - fluidMass));
                // const double fluidRadius=1.01+0.00035298*sqrt(7891604-5666*(1389.9-0.0045*totParticles-fluidMass));
                cout << "Fluid radius = " << fluidRadius << "\n";
                dummyCylinder.R = fluidRadius;
                dummyCylinder.initAxes();
                for (int it = 0; it < totPossibleNodes; ++it) {
                    if (nodes.count(it) == 0) {
                        // creating fluid cells
                        if (getPosition(it).insideCylinder(dummyCylinder.p1 / unit.Length, dummyCylinder.naxes, 0.0, dummyCylinder.R / unit.Length)) {
                            generateNode(it, LIQUID);
                        }
                    }
                }
                break;
            }
            case OPENBARRIER:
            {
                cout << "Initializing interface as in openbarrier case study" << endl;
                unsigned int xMin = 0; //int((0.0375 / unit.Length)); /; // int((avalanchePosit/unit.Length));
                unsigned int xMax = int((4.175 / unit.Length)); //lbSize[0];
                unsigned int zMin = 0;
                unsigned int zMax = int((1.095 / unit.Length));
                for (int it = 0; it < totPossibleNodes; ++it) {
                    if (nodes.count(it) == 0) {
                        // creating fluid cells
                        if ((getPositionX(it) > xMin) && (getPositionX(it) < xMax) && (getPositionZ(it) > zMin) && (getPositionZ(it) < zMax)) {
                            generateNode(it, LIQUID);
                        }
                    }
                }

                break;
            }
            case HONGKONG:
            {
                cout << "Initializing interface as in Hong Kong experiments" << endl;
                unsigned int xMin = 0; //int((0.0375 / unit.Length)); /; // int((avalanchePosit/unit.Length));
                unsigned int xMax = int((0.5 / unit.Length)); //lbSize[0];
                unsigned int zMin = 0;
                unsigned int zMax = int((0.5 / unit.Length));
                for (int it = 0; it < totPossibleNodes; ++it) {
                    if (nodes.count(it) == 0) {
                        // creating fluid cells
                        if ((getPositionX(it) > xMin) && (getPositionX(it) < xMax) && (getPositionZ(it) > zMin) && (getPositionZ(it) < zMax)) {
                            generateNode(it, LIQUID);
                        }
                    }
                }

                break;
            }
            case ESERCITAZIONE:
            {
                cout << "Initializing interface as in Hong Kong experiments" << endl;
                unsigned int xMin = 0; //int((0.0375 / unit.Length)); /; // int((avalanchePosit/unit.Length));
                unsigned int xMax = int((10.0 / unit.Length)); //lbSize[0];
                unsigned int zMin = 0;
                unsigned int zMax = int((9.0 / unit.Length));
                for (int it = 0; it < totPossibleNodes; ++it) {
                    if (nodes.count(it) == 0) {
                        // creating fluid cells
                        if ((getPositionX(it) > xMin) && (getPositionX(it) < xMax) && (getPositionZ(it) > zMin) && (getPositionZ(it) < zMax)) {
                            generateNode(it, LIQUID);
                        }
                    }
                }

                break;
            }
            case STAVA:
            {
                cout << "Initializing interface as in Stava valley study case" << endl;
                unsigned int xMin = int((0.0 / unit.Length));
                unsigned int xMax = int((120.0 / unit.Length));
                unsigned int yMin = int((3670.0 / unit.Length));
                unsigned int yMax = int((3870.0 / unit.Length));
                unsigned int zMin = int((400.0 / unit.Length));
                unsigned int zMax = int((510.0 / unit.Length));
                for (int it = 0; it < totPossibleNodes; ++it) {
                    if (nodes.count(it) == 0) {
                        // creating fluid cells
                        if ((getPositionX(it) > xMin) && (getPositionX(it) < xMax) && (getPositionY(it) > yMin) && (getPositionY(it) < yMax) && (getPositionZ(it) > zMin) && (getPositionZ(it) < zMax)) {
                            generateNode(it, LIQUID);
                        }
                    }
                }

                break;
            }
            case NIGRO:
            {
                cout << "Initializing interface as in Stava valley study case" << endl;
                unsigned int xMin = int((0.25 / unit.Length));
                unsigned int xMax = int((0.75 / unit.Length));
                unsigned int yMin = int((0.25 / unit.Length));
                unsigned int yMax = int((0.75 / unit.Length));
                unsigned int zMin = int((0.0 / unit.Length));
                unsigned int zMax = int((0.9 / unit.Length));
                for (int it = 0; it < totPossibleNodes; ++it) {
                    if (nodes.count(it) == 0) {
                        // creating fluid cells
                        if ((getPositionX(it) > xMin) && (getPositionX(it) < xMax) && (getPositionY(it) > yMin) && (getPositionY(it) < yMax) && (getPositionZ(it) > zMin) && (getPositionZ(it) < zMax)) {
                            generateNode(it, LIQUID);
                        }
                    }
                }

                break;
            }
            case CAROLINE:
            {
                cout << "Initializing interface as in Caroline's experiments" << endl;
                unsigned int zMax = int((3.65e-3 / unit.Length));
                for (int it = 0; it < totPossibleNodes; ++it) {
                    if (nodes.count(it) == 0) {
                        // creating fluid cells
                        if (getPositionZ(it) < zMax) {
                            generateNode(it, LIQUID);
                        }
                    }
                }

                break;
            }
            case DAMBREAK:
            {
                cout << "Initializing interface as in a dam break experiment" << endl;
                unsigned int xMax = int((0.04 / unit.Length));
                unsigned int yMax = int((0.06 / unit.Length));
                for (int it = 0; it < totPossibleNodes; ++it) {
                    if (nodes.count(it) == 0) {
                        // creating fluid cells
                        if (getPositionX(it) < xMax && getPositionY(it) < yMax) {
                            generateNode(it, LIQUID);
                        }
                    }
                }

                break;
            }
            case IERVOLINO:
            case IERVOLINO_2D:
            {
                cout << "Initializing interface as in Iervolino's geometry" << endl;
                const double reservoirX = 0.8;
                const double wallSizeX = 0.025;
                const double outletSizeY = 0.3;
                const double obstacleSizeX = 0.155;
                const double obstacleSizeY = 0.3;
                const double erodibleSizeX = 0.51;
                const double erodibleHeight = 0.02;
                const double edgeRadius = 0.005;
                const double reservoirLevel = 0.1;
                unsigned int xReservoirMax = int(ceil(reservoirX / unit.Length));
                unsigned int zReservoirMin = int(floor(erodibleHeight / unit.Length));
                unsigned int zReservoirMax = int(ceil((erodibleHeight + reservoirLevel) / unit.Length));
                unsigned int xBedMin = int(floor((reservoirX + wallSizeX) / unit.Length));
                unsigned int xBedMax = int(ceil((reservoirX + wallSizeX + erodibleSizeX) / unit.Length));
                unsigned int zBedMax = int(ceil(erodibleHeight / unit.Length));
                for (int it = 0; it < totPossibleNodes; ++it) {
                    if (nodes.count(it) == 0) {
                        unsigned int xHere=getPositionX(it);
                        unsigned int zHere=getPositionZ(it);
                        // creating fluid cells
                        if ((xHere <= xReservoirMax && zHere <= zReservoirMax && zHere >= zReservoirMin)  ||
                        (zHere <= zBedMax && xHere <= xBedMax && xHere >= xBedMin)){
                            generateNode(it, LIQUID);
                        }
                    }
                }

                break;
            }
            case GRAY_DAMBREAK:
            {
                cout << "Initializing interface as in a simplified gray experiment" << endl;
                const double xCenter = (0.5) / unit.Length;
                const double yCenter = (0.65) / unit.Length;
                const double zCenter = (0.0) / unit.Length;
                const tVect center(xCenter, yCenter, zCenter);
                const double radius = 0.35 / unit.Length;
                const double radius2 = radius*radius;
#pragma omp parallel for
                for (int it = 0; it < totPossibleNodes; ++it) {
                    if (nodes.count(it) == 0) {
                        // creating fluid cells
                        const tVect nodePosition = getPosition(it);
                        const tVect centerDistance = center - nodePosition;
                        const double isoParameter = centerDistance.norm2() - radius2;
                        if (isoParameter < 0.0) {
#pragma omp critical
                            {
                                generateNode(it, LIQUID);
                            }
                        }
                    }
                }
                break;
            }
            case GRAY_DAMBREAK_2D:
            {
                cout << "Initializing interface as in a simplified gray experiment" << endl;
                const double xCenter = (0.5) / unit.Length;
                const double yCenter = (0.0) / unit.Length;
                const double zCenter = (0.0) / unit.Length;
                const tVect center(xCenter, yCenter, zCenter);
                const double radius = 0.35 / unit.Length;
                const double radius2 = radius*radius;
#pragma omp parallel for
                for (int it = 0; it < totPossibleNodes; ++it) {
                    if (nodes.count(it) == 0) {
                        // creating fluid cells
                        const tVect nodePosition = getPosition(it);
                        const tVect centerDistance = center - nodePosition;
                        const double isoParameter = centerDistance.norm2() - radius2;
                        if (isoParameter < 0.0) {
#pragma omp critical
                            {
                                generateNode(it, LIQUID);
                            }
                        }
                    }
                }
                break;
            }
            case INCLINEFLOW:
            {
                cout << "Initializing interface as in Lagrée, Staron, Popinet" << endl;
                unsigned int yMax = int((1.0 / unit.Length));
                for (int it = 0; it < totPossibleNodes; ++it) {
                    if (nodes.count(it) == 0) {
                        // creating gas cells
                        if (getPositionY(it) < yMax) {
                            generateNode(it, LIQUID);
                        }
                    }
                }

                break;
            }
            case JOP:
            {
                cout << "Initializing interface as in Jop, Forterre, Pouliquen" << endl;
                unsigned int zMax = int((0.03 / unit.Length)); //0.0265 0.0848  0.075 
                for (int it = 0; it < totPossibleNodes; ++it) {
                    if (nodes.count(it) == 0) {
                        // creating fluid cells
                        if (getPositionZ(it) < zMax) {
                            generateNode(it, LIQUID);
                        }
                    }
                }
                break;
            }
            case WILL:
            {
                cout << "Initializing interface as in Will's Experiments (UniNottingham)" << endl;
                const unsigned int zMax = int((0.05 / unit.Length)); 
                const double xCenter = (0.0 / unit.Length);  
                const double yCenter = (0.0 / unit.Length);
                const double radius2 = (0.054*0.054 / unit.Length/ unit.Length);
                for (int it = 0; it < totPossibleNodes; ++it) {
                if (nodes.count(it) == 0) {
                        const double xHere=double(getPositionX(it));
                        const double yHere=double(getPositionY(it));
                        const unsigned int zHere=getPositionZ(it);
                        const double distance2=(xHere-xCenter)*(xHere-xCenter)+(yHere-yCenter)*(yHere-yCenter);
                        // creating fluid cells
                        if ((zHere <= zMax) && (distance2 <= radius2)){
                            generateNode(it, LIQUID);
                        } 
                }
                }
                break;
            }
            case GRAY:
            {
                cout << "Initializing interface as in Gray, Wieland, Hutter" << endl;
                const double xCenter = (-0.07887 + translateTopographyX) / unit.Length;
                const double yCenter = (-0.00 + translateTopographyY) / unit.Length;
                const double zCenter = (1.2422 + translateTopographyZ) / unit.Length;
                const tVect center(xCenter, yCenter, zCenter);
                const double radius = 0.3427 / unit.Length;
                const double radius2 = radius*radius;
#pragma omp parallel for
                for (int it = 0; it < totPossibleNodes; ++it) {
                    if (nodes.count(it) == 0) {
                        // creating fluid cells
                        const tVect nodePosition = getPosition(it);
                        const tVect centerDistance = center - nodePosition;
                        const double isoParameter = centerDistance.norm2() - radius2;
                        if (isoParameter < 0.0) {
                            const double surfaceIsoparameter = lbTop.surfaceIsoparameter(nodePosition * unit.Length);
                            if (surfaceIsoparameter < 0.0) {
#pragma omp critical
                                {
                                    generateNode(it, LIQUID);
                                }
                            }
                        }
                    }
                }

                break;
            }
            case HOURGLASS:
            {
                cout << "Initializing interface as in hourglass case study" << endl;
                unsigned int zMin = int((hourglassOutletHeight / unit.Length));
                unsigned int zMax = lbSize[2] - 3.0;
                for (int it = 0; it < totPossibleNodes; ++it) {
                    if (nodes.count(it) == 0) {
                        // creating fluid cells
                        if (getPositionZ(it) > zMin && getPositionZ(it) < zMax) {
                            generateNode(it, LIQUID);
                        }
                    }
                }

                break;
            }
            case HEAP:
            {
                cout << "Initializing interface as in heap formation (continuum) case study" << endl;
                unsigned int zMin = int((heapBaseLevel / unit.Length));
                unsigned int zMax = lbSize[2] - 3.0;
                for (int it = 0; it < totPossibleNodes; ++it) {
                    if (nodes.count(it) == 0) {
                        // creating gas cells
                        if (getPositionZ(it) > zMin && getPositionZ(it) < zMax) {
                            generateNode(it, LIQUID);
                        }
                    }
                }

                break;
            }
            case NONE:
            default:
            {
                cout << "Initializing interface using box defined in config file:" << endl;
                cout << "X=(" << double(freeSurfaceBorders[0]) * unit.Length << ", " << double(freeSurfaceBorders[1]) * unit.Length << ")" << endl;
                cout << "Y=(" << double(freeSurfaceBorders[2]) * unit.Length << ", " << double(freeSurfaceBorders[3]) * unit.Length << ")" << endl;
                cout << "Z=(" << double(freeSurfaceBorders[4]) * unit.Length << ", " << double(freeSurfaceBorders[5]) * unit.Length << ")" << endl;
                for (int it = 0; it < totPossibleNodes; ++it) {
                    if (nodes.count(it) == 0) {
                        // creating fluid cells
                        if ((getPositionX(it) > freeSurfaceBorders[0]) &&
                                (getPositionX(it) < freeSurfaceBorders[1]) &&
                                (getPositionY(it) > freeSurfaceBorders[2]) &&
                                (getPositionY(it) < freeSurfaceBorders[3]) &&
                                (getPositionZ(it) > freeSurfaceBorders[4]) &&
                                (getPositionZ(it) < freeSurfaceBorders[5])) {
                            generateNode(it, LIQUID);
                        }
                    }
                }

                break;
                //                cout << "Fluid everywhere" << endl;
                //                for (int it = 0; it < totPossibleNodes; ++it) {
                //                    if (nodes.count(it) == 0) {
                //                        generateNode(it, LIQUID);
                //                    }
                //                }
                //
                //                break;
            }
        }
    }
}

void LB::restartInterface(ifstream& fluidFileID, unsigned int& restartNodes) {
    // creates an interface electing interface cells from active cells
    // defining a surface

    cout << "Restarting interface" << endl;

    // creating list and initialize macroscopic variables for all nodes except walls
    double itHere = 0;
    unsigned int typeHere = 0;
    double nHere = 0.0;
    double uXHere = 0.0;
    double uYHere = 0.0;
    double uZHere = 0.0;
    double massHere = 0.0;
    double viscHere = 0.0;
    double fHere[lbmDirec];
    for (int j = 0; j < lbmDirec; ++j) {
        fHere[j] = 0.0;
    }

    // checking for boundary between gas and fluid and assigning interface properties
    for (int i = 0; i < restartNodes; ++i) {
        fluidFileID>>itHere;
        fluidFileID>>typeHere;
        fluidFileID>>nHere;
        fluidFileID>>uXHere;
        fluidFileID>>uYHere;
        fluidFileID>>uZHere;
        const tVect uHere(uXHere, uYHere, uZHere);
        fluidFileID>>massHere;
        fluidFileID>>viscHere;
        for (int j = 0; j < lbmDirec; ++j) {
            fluidFileID >> fHere[j];
        }

        // creating node
        generateNode(itHere, types(typeHere));

        // setting macroscopic variables
        if (nodes[itHere].isFluid()) {
            nodes[itHere].restart(nHere + fluidMaterial.initDensity, uHere, massHere + fluidMaterial.initDensity, viscHere, fHere);
        } else if (nodes[itHere].isInterface()) {
            // density=initDensity; velocity=initVelocity, mass=0.5*initDensity; viscosity=initVisc; force=lbF
            nodes[itHere].restart(nHere + fluidMaterial.initDensity, uHere, massHere + fluidMaterial.initDensity, viscHere, fHere);
        }
    }
}

void LB::initializeParticleBoundaries(particleList& particles) {

    // INITIAL PARTICLE POSITION ////////////////////////
    cout << "Initializing particle nodes ...";

    // set all outside
    //#pragma omp parallel for
    for (nodeMap::iterator imap = nodes.begin(); imap != nodes.end(); imap++) {
        node* nodeHere = &imap->second;
        unsigned int indexHere = imap->first;
        // reset old list
        nodeHere->setOutsideParticle();
    }

    for (int n = 0; n < particles.size(); ++n) {
        const tVect convertedPosition = particles[n].x0 / unit.Length;
        const double convertedRadius = hydrodynamicRadius * particles[n].r / unit.Length;
        const unsigned int particleIndexHere = particles[n].particleIndex;
#pragma omp parallel for
        for (int it = 0; it < activeNodes.size(); ++it) {
            node* nodeHere = activeNodes[it];
            unsigned int nodeIndexHere = nodeHere->coord;
            // checking if node is inside a particle
            if (getPosition(nodeIndexHere).insideSphere(convertedPosition, convertedRadius)) { //-0.5?
                nodeHere->setInsideParticle();
                nodeHere->setSolidIndex(particleIndexHere);
            }
        }
    }

    //    particleNodes.clear();
    //    // creating list and initialize macroscopic variables for all nodes except walls
    //    for (nodeMap::iterator imap = nodes.begin(); imap != nodes.end(); imap++) {
    //        unsigned int it = imap->first;
    //        node* nodeHere = &imap->second;
    //
    //        // PARTICLE NODES ////
    //        if (imap->second.isInsideParticle()) {
    //            // adding the particle node to the list
    //            particleNodes.push_back(nodeHere);
    //        }
    //    }
}

void LB::checkNewInterfaceParticles(elmtList& elmts, particleList& particles) {

    // INITIAL PARTICLE POSITION ////////////////////////
    //cout << "Checking for new particle nodes at interface ...";

    for (int m = 0; m < elmts.size(); ++m) {
        if (elmts[m].FHydro.norm2() == 0.0) {
            for (int n = 0; n < elmts[m].components.size(); ++n) {

                const tVect convertedPosition = particles[elmts[m].components[n]].x0 / unit.Length;
                const double convertedRadius = hydrodynamicRadius * particles[elmts[m].components[n]].r / unit.Length;
                const unsigned int particleIndexHere = particles[elmts[m].components[n]].particleIndex;
#pragma omp parallel for
                for (int it = 0; it < interfaceNodes.size(); ++it) {
                    node* nodeHere = interfaceNodes[it];
                    if (!nodeHere->isInsideParticle()) {
                        unsigned int nodeIndexHere = nodeHere->coord;
                        // checking if node is inside a particle
                        if (getPosition(nodeIndexHere).insideSphere(convertedPosition, convertedRadius)) { //-0.5?
                            nodeHere->setInsideParticle();
                            nodeHere->setSolidIndex(particleIndexHere);
                            //cout<<"new particle node"<<nodeHere->coord<<" for "<<particleIndexHere<<endl;
                        }
                    }
                }
            }
        }
    }

}

void LB::initializeLists() {

    cout << "Resetting lists ...";

    // note that interface is not defined here. All fluid, interface and gas cells are 0 at the moment
    fluidNodes.clear();
    interfaceNodes.clear();

    // creating list and initialize macroscopic variables for all nodes except walls
    for (nodeMap::iterator imap = nodes.begin(); imap != nodes.end(); imap++) {
        unsigned int it = imap->first;
        node* nodeHere = &imap->second;
        // FLUID NODES ////
        if (imap->second.isFluid()) {
            // adding the free node indexes to the list
            fluidNodes.push_back(nodeHere);
        }// INTERFACE NODES ////
        else if (imap->second.isInterface()) {
            // interface cells list update
            interfaceNodes.push_back(nodeHere);
        }

    }
    cout << " done" << endl;
}

void LB::initializeVariables() {

    cout << "Initializing variables" << endl;
    // note that interface is not defined here. All fluid, interface and gas cells are uninitialized at the moment
    // calculate maximum height of the fluid

    // find "taller" and "deepest" points
    double minProjection = std::numeric_limits<double>::max();
    double maxProjection = -std::numeric_limits<double>::max();
        
    if (!solveCentrifugal) {
        for (nodeMap::iterator imap = nodes.begin(); imap != nodes.end(); imap++) {
            if (imap->second.isActive()) {
                const unsigned int it = imap->first;
                const tVect position = getPosition(it);
                const double projection = position.dot(lbF);
                minProjection = std::min(minProjection, projection);
                maxProjection = std::max(maxProjection, projection);
            }
        }
        cout << "minProjection = " << minProjection << endl;
    } else {
        for (nodeMap::iterator imap = nodes.begin(); imap != nodes.end(); imap++) {
            if (imap->second.isActive()) {
                const unsigned int it = imap->first;
                const tVect position = getPosition(it);
                const double projection = position.dot(nodes[it].centrifugalForce);
                minProjection = std::min(minProjection, projection);
                maxProjection = std::max(maxProjection, projection);
            }
        }
        cout << "minProjection = " << minProjection << endl;
    }
    
    //    double maxX(0.0), maxY(0.0), maxZ(0.0);
    //    for (nodeMap::iterator imap = nodes.begin(); imap != nodes.end(); imap++) {
    //        if (imap->second.isActive()) {
    //            unsigned int it = imap->first;
    //            maxX = std::max(maxX, getPositionX(it));
    //            maxY = std::max(maxY, getPositionY(it));
    //            maxZ = std::max(maxZ, getPositionZ(it));
    //        }
    //    }
    //    const tVect maxP(maxX, maxY, maxZ);

    // checking for boundary between gas and fluid and assigning interface properties
    // at this point fluid cells contain actual fluid cells and potential interface cells, so we create the node anyway
    double massFluid = 0.0;
    double massInterface = 0.0;
    for (nodeMap::iterator imap = nodes.begin(); imap != nodes.end(); imap++) {
        node* nodeHere = &imap->second;
        unsigned int indexHere = imap->first;
        if (nodeHere->isFluid()) {
            // check if it is interface
            for (int j = 1; j < lbmDirec; ++j) {
                node* linkNode = nodeHere->d[j];
                if (linkNode == 0) {
                    nodeHere->setInterface();
                    break;
                }
            }
        }
        // now assign macroscopic quantities accordingly
        // FLUID NODES ////
        if (nodeHere->isFluid()) {
            massFluid += 1.0;
            // setting macroscopic variables
            // density is calculated using hydrostatic profile
            //tVect deltah = getPosition(indexHere) - maxP;
            // density=initDensity; velocity=initVelocity, mass=initDensity; viscosity=initVisc; force=lbF
            const tVect position = getPosition(indexHere);
            // nodeHere->initialize(fluidMaterial.initDensity + 3.0 * fluidMaterial.initDensity * (deltah.dot(lbF+nodeHere->centrifugalForce) + 0.5 * lbF.norm()), initVelocity, fluidMaterial.initDensity, fluidMaterial.initDynVisc, lbF, 1.0, rotationSpeed);
            if (!solveCentrifugal) {
                const double projection = position.dot(lbF);
                nodeHere->initialize(fluidMaterial.initDensity + 3.0 * fluidMaterial.initDensity * (projection-minProjection), initVelocity, fluidMaterial.initDensity, fluidMaterial.initDynVisc, lbF, 1.0, Zero);
            } else {
                const double projection = position.dot(nodeHere->centrifugalForce);
                nodeHere->initialize(fluidMaterial.initDensity + 3.0 * fluidMaterial.initDensity * (projection-minProjection), initVelocity, fluidMaterial.initDensity, fluidMaterial.initDynVisc, lbF, 1.0, rotationSpeed);
            }
        }// INTERFACE NODES ////
        else if (nodeHere->isInterface()) {
            massInterface += 0.5;
            // setting macroscopic variables
            // density=initDensity; velocity=initVelocity, mass=0.5*initDensity; viscosity=initVisc; force=lbF
            nodeHere->initialize(fluidMaterial.initDensity, initVelocity, 0.5 * fluidMaterial.initDensity, fluidMaterial.initDynVisc, lbF, 1.0, rotationSpeed);
        }

    }
    cout << "Approximate volume = " << massFluid * unit.Volume << " (fluid body), " << massInterface * unit.Volume << " (interface), " << (massFluid + massInterface) * unit.Volume << " (tot), " << endl;


}

void LB::initializeWalls(wallList& walls, cylinderList& cylinders, objectList& objects) {
    cout << "Initializing wall nodes" << endl;
    const double zero = 0.0;

    wallNodes.clear();

    // initializing wall nodes
    // note that, in the hypothesis that these walls are not evolving, only nodes at the interface need creation
    for (nodeMap::iterator imap = nodes.begin(); imap != nodes.end(); imap++) {
        node* nodeHere = &imap->second;
        unsigned int indexHere = imap->first;
        if (nodeHere->isWall()) {
            // initialize node
            // STATIC WALL NODES ////
            if (nodeHere->isStatWall() || nodeHere->isSlipStatWall() || nodeHere->isObjectWall() || nodeHere->isTopography()) {
                // reset velocity and mass (useful for plotting)
                // density=0.0; velocity=(0.0,0.0,0.0), mass=0.0; viscosity=0.0; force=(0.0,0.0,0.0)
                nodeHere->initialize(fluidMaterial.initDensity, Zero, zero, zero, Zero, 1.0, Zero);
            }// DYNAMIC WALL NODES ////
            else if (nodeHere->isDynWall() || nodeHere->isSlipDynWall() || nodeHere->isCylinderWall()) {
                // need to define velocity. It could be part of a cylinder or wall, we check both
                tVect solidVelocity;
                const tVect nodePosition = getPosition(indexHere);
                unsigned int solidIndex = nodeHere->getSolidIndex();
                // wall
                if (nodePosition.insidePlane(walls[solidIndex].p / unit.Length, walls[solidIndex].n)) {
                    solidVelocity = walls[solidIndex].getSpeed(nodePosition * unit.Length) / unit.Speed;
                }// cylinder
                else if (!nodePosition.insideCylinder(cylinders[solidIndex].p1 / unit.Length, cylinders[solidIndex].naxes, 0.0, cylinders[solidIndex].R / unit.Length)) {
                    solidVelocity = cylinders[solidIndex].getSpeed(nodePosition * unit.Length) / unit.Speed;
                }// objects
                else if (nodePosition.insideSphere(objects[solidIndex].x0 / unit.Length, objects[solidIndex].r / unit.Length)) {
                    solidVelocity = objects[solidIndex].x1 / unit.Speed;
                }
                // reset velocity and mass (useful for plotting)
                // density=0.0; velocity=solidVelocity, mass=0.0; viscosity=0.0; force=(0.0,0.0,0.0)
                nodeHere->initialize(fluidMaterial.initDensity, solidVelocity, zero, zero, Zero, 1.0, rotationSpeed);
            }
            // add node to list
            wallNodes.push_back(nodeHere);
        }
    }
}


// integration functions

void LB::cleanLists() {
#ifdef _DEBUG
    // sorts the active-node list and removes duplicates
    std::sort(fluidNodes.begin(), fluidNodes.end());
    std::sort(interfaceNodes.begin(), interfaceNodes.end());

    for (int ind = fluidNodes.size() - 1; ind > 0; --ind) {
        node* i = fluidNodes[ind];
        node* j = fluidNodes[ind - 1];
        if (i == j) {
            cout << "duplicate-fluid!" << endl;
            fluidNodes.erase(fluidNodes.begin() + ind);
        }
    }
    for (int ind = interfaceNodes.size() - 1; ind > 0; --ind) {
        node* i = interfaceNodes[ind];
        node* j = interfaceNodes[ind - 1];
        if (i == j) {
            cout << "duplicate-interface!" << endl;
            interfaceNodes.erase(interfaceNodes.begin() + ind);
        }
    }
#endif

    // list with active nodes i.e. nodes where collision and streaming are solved
    // solid nodes, particle nodes and gas nodes are excluded
    activeNodes.clear();
    activeNodes.reserve(fluidNodes.size() + interfaceNodes.size());
    activeNodes.insert(activeNodes.end(), fluidNodes.begin(), fluidNodes.end());
    activeNodes.insert(activeNodes.end(), interfaceNodes.begin(), interfaceNodes.end());

}


void LB::collision(node* nodeHere) {
    if (TRTsolver == false) {
        // equilibrium distributions
        double feq[lbmDirec];
        // force field
        tVect force=lbF;
        if (solveCentrifugal) {
            force+=nodeHere->centrifugalForce;
        }
        if (solveCentrifugal) {
            force+=computeCoriolis(nodeHere->u,rotationSpeed);
        }
        // shift velocity field to F/2
        nodeHere->shiftVelocity(force);
        nodeHere->computeEquilibrium(feq);
        //        if (types[index].isInsideParticle()) {
        //            nodes[index]->resetViscosity(initDynVisc);
        //        }
        if (fluidMaterial.rheologyModel != NEWTONIAN || fluidMaterial.turbulenceOn) {
            // compute shear rate tensor, find invariant calculate viscosity (Bingham)
            nodeHere->computeApparentViscosity(feq, fluidMaterial);
        }
        //if (nodeHere->isInterface()) nodeHere->visc=fluidMaterial.lbMaxVisc;
        //else if (fluidMaterial.rheologyModel == NEWTONIAN)  nodeHere->visc=fluidMaterial.initDynVisc;
        // compute new distributions
        nodeHere->solveCollision(feq);
        // add force term to new distributions
        //nodeHere->addForce(lbF);
        nodeHere->addForce(feq,force);
    } else {
        // equilibrium distributions
        double feqp[lbmDirec];
        double feqm[lbmDirec];
        // force field
        tVect force=lbF;
        if (solveCentrifugal) {
            force+=nodeHere->centrifugalForce;
        }
        if (solveCentrifugal) {
            force+=computeCoriolis(nodeHere->u,rotationSpeed);
        }
        // shift velocity field to F/2
        nodeHere->shiftVelocity(force);
        nodeHere->computeEquilibriumTRT(feqp, feqm);
        //        if (types[index].isInsideParticle()) {
        //            nodes[index]->resetViscosity(initDynVisc);
        //        }
        if (fluidMaterial.rheologyModel != NEWTONIAN || fluidMaterial.turbulenceOn) {
            // compute shear rate tensor, find invariant calculate viscosity (Bingham)
            nodeHere->computeApparentViscosity(feqp, fluidMaterial);
        }
        else nodeHere->visc = fluidMaterial.initDynVisc;
        // compute new distributions
        nodeHere->solveCollisionTRT(feqp, feqm, magicNumber);
        // add force term to new distributions
        nodeHere->addForceTRT(force);
    }

}

double LB::curvedWallReconstruction(const unsigned int& j, const node* nodeHere, const tVect& wallSpeed) const {


    // bounce-back as in Guo, Zhen, Shi (2002)
    // first thing we check if delta is large enough to use u_f, or if we have to use u_ff
    const double deltaHere = nodeHere->curved->delta[j];
    const tVect velHere = nodeHere->u;
    const double tauHere = 0.5 + 3.0 * nodeHere->visc;
    // get value for chi and ubf
    tVect ubf(0.0, 0.0, 0.0);
    double chi = 0.0;
    if (deltaHere >= 0.5) {
        // get fluid reference velocity
        ubf = (2.0 * deltaHere - 3.0) / (2.0 * deltaHere) * velHere + 3.0 / (2.0 * deltaHere) * wallSpeed;
        chi = (2.0 * deltaHere - 1.0) / (tauHere + 0.5);
    } else {
        const node* linkBack = nodeHere->d[opp[j]];
        bool linkBackNodeExists = false;
        if (linkBack != 0) {
            if (linkBack->isActive()) {
                linkBackNodeExists = true;
            }
        }
        if (linkBackNodeExists) {
            // uw2 can be taken exactly
            //ubf=(deltaHere-1.0)/deltaHere*velHere+wallSpeed/deltaHere;
            ubf = linkBack->u; // nodeHere->curved->m3[j] *
        } else {
            // for uw2 we assume uff=uf, since there are not enough fluid nodes.
            ubf = velHere;

        }
        chi = (2.0 * deltaHere - 1.0) / (tauHere - 2.0);
    }

    // calculate missing distribution
    const double usq = velHere.norm2();
    const double vu = velHere.dot(v[j]);
    // coefficients for equilibrium
    static const double C1 = 3.0;
    static const double C2 = 4.5;
    static const double C3 = -1.5;
    // equilibrium
    const double fStar = nodeHere->n * coeff[j]*(1 + C1 * v[j].dot(ubf) + C2 * vu * vu + C3 * usq);

    return (1.0 - chi) * nodeHere->fs[j] + chi * fStar;
    // adding the extra mass to the surplus
    //extraMass += nodes[it]->mass * (chi * nodes[it]->fs[j] - chi * fStar + BBi);
}



//double LB::curvedWallReconstruction(const unsigned int& j, const node* nodeHere) const {
//
//    // coefficients for equilibrium
//    static const double C1 = 3.0;
//    static const double C2 = 4.5;
//    static const double C3 = 1.5;
//    
//    // bounce-back as in Guo, Zhen, Shi (2002)
//    // first thing we check if delta is large enough to use u_f, or if we have to use u_ff
//    const double deltaHere = nodeHere->curved->delta[j];
//    nodeHere->curved->delta[j] = 0.5;
//    // compute uw, the extrapolated fluid velocity at the wall, as uw1 (see paper)
//    tVect uw = nodeHere->curved->m1[j] * nodeHere->u;
//    // tentatively calculate equilibrium of uw1
//    double usq = uw.norm2();
//    double vu = uw.dot(v[opp[j]]);
//    double feq_uw = nodeHere->n * coeff[j]*(1 + C1 * vu + C2 * vu * vu - C3 * usq);
//    // and the non-equilibrium distribution
//    double fneq_uw = nodeHere->fs[j] - feq_uw;
////    if (deltaHere < 0.75) {
////        // let's see if we can use uw2, an extra node is needed
////        const node* linkBack = nodeHere->d[opp[j]];
////        bool linkBackNodeExists = false;
////        if (linkBack != 0) {
////            if (linkBack->isActive()) {
////                linkBackNodeExists = true;
////            }
////        }
////        tVect uw2(0.0, 0.0, 0.0);
////        if (linkBackNodeExists) {
////            // uw2 can be taken exactly
////            uw2 = nodeHere->curved->m3[j] * linkBack->u;
////        } else {
////            // for uw2 we assume uff=uf, since there are not enough fluid nodes.
////            uw2 = nodeHere->curved->m3[j] * nodeHere->u;
////        }
////        // any case, update uw as uw=delta*uw1(1-delta)*uw2
////        uw = deltaHere * uw + (1.0 - deltaHere) * uw2;
////        // recompute equilibrium
////        usq = uw.norm2();
////        vu = uw.dot(v[j]);
////        feq_uw = deltaHere * feq_uw + (1.0 - deltaHere)*(nodeHere->n * coeff[j]*(1 + C1 * vu + C2 * vu * vu - C3 * usq));
////        fneq_uw = deltaHere * fneq_uw + (1.0 - deltaHere)*(nodeHere->fs[j] - feq_uw);
////    }
//    // tau needed
//    const double tauHere = 1.0; //0.5 + 3.0 * nodeHere->visc;
//    // new distribution
//    return feq_uw + (1.0 - 1.0 / tauHere) * fneq_uw;
//    
//}

void LB::streaming(wallList& walls, objectList& objects) {
    // STREAMING STEP

    // coefficient for free-surface
    static const double C2x2 = 9.0;
    static const double C3x2 = 3.0;
    // coefficient for slip conditions
    static const double S1 = slipCoefficient;
    static const double S2 = (1.0 - slipCoefficient);
    // creating list for collision function
    double staticPres[lbmDirec];
    for (int j = 0; j < lbmDirec; j++) {
        staticPres[j] = fluidMaterial.initDensity * coeff[j];
    }
    // coefficient for bounce-back
    static const double BBCoeff = 2.0 * 3.0;
    // extra mass due to bounce-back and moving walls
    double extraMass = 0.0;

    // initializing wall forces
    for (int iw = 0; iw < walls.size(); ++iw) {
        walls[iw].FHydro.reset();
    }
    // initializing object forces
    for (int io = 0; io < objects.size(); ++io) {
        objects[io].FHydro.reset();
    }

    //  Saving in support variables f->fs
#pragma omp parallel for
    for (int it = 0; it < activeNodes.size(); ++it) {
        activeNodes[it]->store();
    }

    //    for (int in = 0; in < wallNodes.size(); ++in) {
    //        wallNodes[in]->surfaceCorner = false;;
    //    }
    //    for (int in = 0; in < interfaceNodes.size(); ++in) {
    //        node* nodeHere = interfaceNodes[in];
    //        // cycling through neighbours
    //        for (int j = 1; j < 7; ++j) {
    //            // getting neighbour index
    //            node* linkNode = nodeHere->d[j];
    //            if (linkNode!=nullptr) {
    //                if (linkNode->isWall()) {
    //                    linkNode->surfaceCorner = true;
    //                }
    //            }
    //        }
    //    }

    //  Streaming
#pragma omp parallel for ordered reduction(+:extraMass)
    // cycling through active nodes
    for (int in = 0; in < activeNodes.size(); ++in) {

        // pointer to active node
        node* nodeHere = activeNodes[in];
        // cycling through neighbours
        for (int j = 1; j < lbmDirec; ++j) {
            // getting neighbour index
            const node* linkNode = nodeHere->d[j];
            // if neighbour is normal fluid cell what follows is true

            if (linkNode == nullptr) { // is gas
                // additional variables for equilibrium f computation
                const double usq = nodeHere->u.norm2();
                const double vuj = nodeHere->u.dot(v[j]);
                //streaming with constant pressure interface
                nodeHere->f[opp[j]] = -nodeHere->fs[j] + coeff[j] * fluidMaterial.initDensity * (2.0 + C2x2 * (vuj * vuj) - C3x2 * usq);

            } else {

                switch (linkNode->type) {
                    case LIQUID:
                    {
#ifdef DEBUG
                        // TEST USING AGE //////////////////////////////////////
                        const double usq = nodeHere->u.norm2();
                        const double vuj = nodeHere->u.dot(v[j]);
                        nodeHere->f[opp[j]] = linkNode->age * linkNode->fs[opp[j]]+
                                (1.0 - linkNode->age)*(-nodeHere->fs[j] + coeff[j] * fluidMaterial.initDensity * (2.0 + C2x2 * (vuj * vuj) - C3x2 * usq));
                        nodeHere->f[opp[j]] = linkNode->fs[opp[j]];
#else
                        nodeHere->f[opp[j]] = linkNode->fs[opp[j]];
#endif
                        break;
                    }
                    case INTERFACE:
                    {
#ifdef DEBUG 
                        // TEST FOR RECONSTRUCTED F DEPENDING ON ACTUAL SUFACE POSITION //////////////////////////////////////
//                        if (nodeHere->isInterface()) {
//                            // nodeHere->f[opp[j]] = linkNode->fs[opp[j]];
//                            tVect distance = unit.Length * (getPosition(linkNode->coord) - getPosition(nodeHere->coord));
//                            distance /= distance.norm();
//                            const double semiSpace = nodeHere->surfaceNormal.dot(distance);
//                            if (semiSpace > 0.0) {
//                                // additional variables for equilibrium f computation
//                                const double usq = nodeHere->u.norm2();
//                                const double vuj = nodeHere->u.dot(v[j]);
//                                //streaming with constant pressure interface
//                                nodeHere->f[opp[j]] = (1.0 - semiSpace) * linkNode->fs[opp[j]] + semiSpace * (-nodeHere->fs[j] + coeff[j] * fluidMaterial.initDensity * (2.0 + C2x2 * (vuj * vuj) - C3x2 * usq));
//                            } else {
//                                nodeHere->f[opp[j]] = linkNode->fs[opp[j]];
//                            }
//                        } else {
//                            nodeHere->f[opp[j]] = linkNode->fs[opp[j]];
//                            const double usq = nodeHere->u.norm2();
//                            const double vuj = nodeHere->u.dot(v[j]);
//                            const double liquidFraction = (linkNode->mass / linkNode->n);
//                            nodeHere->f[opp[j]] = liquidFraction * linkNode->fs[opp[j]]+
//                                    (1.0 - liquidFraction)*(-nodeHere->fs[j] + coeff[j] * fluidMaterial.initDensity * (2.0 + C2x2 * (vuj * vuj) - C3x2 * usq));
//                            nodeHere->f[opp[j]] = linkNode->fs[opp[j]];
//                        }
                        
                        // TEST USING AGE //////////////////////////////////////
                        const double usq = nodeHere->u.norm2();
                        const double vuj = nodeHere->u.dot(v[j]);
                        nodeHere->f[opp[j]] = linkNode->age * linkNode->fs[opp[j]]+
                                (1.0 - linkNode->age)*(-nodeHere->fs[j] + coeff[j] * fluidMaterial.initDensity * (2.0 + C2x2 * (vuj * vuj) - C3x2 * usq));

                        

#else

                        nodeHere->f[opp[j]] = linkNode->fs[opp[j]];

                        //const double usq = nodeHere->u.norm2();
                        //const double vuj = nodeHere->u.dot(v[j]);
                        //const double liquidFraction= (linkNode->mass/linkNode->n);
                        //nodeHere->f[opp[j]] = liquidFraction*linkNode->fs[opp[j]]+ (1.0-liquidFraction)*(-nodeHere->fs[j] + coeff[j] * fluidMaterial.initDensity * (2.0 + C2x2 * (vuj * vuj) - C3x2 * usq));

#endif
                        break;
                        //} else {
                        // simple bounce-back
                        //    nodes[it]->f[opp[j]] = nodes[it]->fs[j];
                        //}
                    } // for walls there is simple bounce-back
                    case STAT_WALL:
                    {
#ifndef DEBUG 
                        if (nodeHere->isInterface()) {
                            // additional variables for equilibrium f computation
                            const double usq = nodeHere->u.norm2();
                            const double vuj = nodeHere->u.dot(v[j]);
                            //streaming with constant pressure interface
                            //nodeHere->f[opp[j]] = nodeHere->fs[j];
                            nodeHere->f[opp[j]] = -nodeHere->fs[j] + coeff[j] * fluidMaterial.initDensity * (2.0 + C2x2 * (vuj * vuj) - C3x2 * usq);

                        } else {

                            // getting the index of the wall to compute force in the right object
                            const int solidIndex = linkNode->getSolidIndex();

                            // static pressure is subtracted in order to correctly compute buoyancy for floating objects
                            const tVect BBforce = nodeHere->bounceBackForce(j, staticPres, 0.0, fluidMaterial.earthPressureCoeff, fluidMaterial.lbMaxVisc);
                            // updating force and torque on the object (lattice units). This point is critical since many nodes update the force on the same object (lattice units)
#pragma omp ordered
                            {
                                walls[solidIndex].FHydro += BBforce;
                            }
                            nodeHere->f[opp[j]] = nodeHere->fs[j];
                        }
#else
                        // getting the index of the wall to compute force in the right object
                        const int solidIndex = linkNode->getSolidIndex();

                        // static pressure is subtracted in order to correctly compute buoyancy for floating objects
                        const tVect BBforce = nodeHere->bounceBackForce(j, staticPres, 0.0, fluidMaterial.earthPressureCoeff, fluidMaterial.lbMaxVisc);
                        // updating force and torque on the object (lattice units). This point is critical since many nodes update the force on the same object (lattice units)
#pragma omp ordered
                        {
                            walls[solidIndex].FHydro += BBforce;
                        }
                        nodeHere->f[opp[j]] = nodeHere->fs[j];
#endif      
                        break;
                    }
                        //                    case TOPO:
                        //                    {
                        //                        // simple bounce-back
                        //                        nodeHere->f[opp[j]] += nodeHere->fs[j];
                        //                        break;
                        //                    }
                        //                    case TOPO:
                        //                    {
                        //                        // simple bounce-back
                        ////                        if (nodeHere->applySlip==false || j < 7) {
                        ////                            nodeHere->f[opp[j]] = nodeHere->fs[j];
                        ////                        }
                        ////                        else {
                        //                            bool active1 = false;
                        //                            bool active2 = false;
                        //                            const node* nodeCheck1 = nodeHere->d[slip1Check[j]];
                        //                            const node* nodeCheck2 = nodeHere->d[slip2Check[j]];
                        //                            // check for the environment
                        //                            if (nodeCheck1 != 0) {
                        //                                if (nodeCheck1->isActive()) {
                        //                                    active1 = true;
                        //                                }
                        //                            }
                        //                            if (nodeCheck2 != 0) {
                        //                                if (nodeCheck2->isActive()) {
                        //                                    active2 = true;
                        //                                }
                        //                            }
                        //                            // given the environment, perform the right operation
                        //                            if (active1 && !active2) {
                        //                                // first
                        //                                nodeHere->f[opp[j]] = S1 * nodeCheck1->fs[slip1[j]] + S2 * nodeHere->fs[j];
                        //
                        //                            } else if (!active1 && active2) {
                        //                                // second
                        //                                nodeHere->f[opp[j]] = S1 * nodeCheck2->fs[slip2[j]] + S2 * nodeHere->fs[j];
                        //                            } else {
                        //                                // standard BB
                        //                                nodeHere->f[opp[j]] = nodeHere->fs[j];
                        //                            }
                        ////                        }
                        //                        break;
                        //                    }
                        // for curved walls there is the rule of Mei-Luo-Shyy
                    case TOPO:
                    {
                        nodeHere->f[opp[j]] = nodeHere->fs[j];

                        //                        // TEST: DOES NOT WORK FOR NIGRO
                        //                        const tVect wallNormalHere = nodeHere->curved->wallNormal;
                        //                        tVect slipVelocity(0.0, 0.0, 0.0);
                        //                        if (wallNormalHere.dot(Zp) < 0.9) { //(nodeHere->applySlip) {
                        //                            slipVelocity = 0.98 * (nodeHere->u - (nodeHere->u.dot(wallNormalHere)) * wallNormalHere);
                        //                        }
                        //
                        //                        const double newDistribution = curvedWallReconstruction(j, nodeHere, slipVelocity);
                        //
                        //                        //const double BBi nodeHere->n*6.0*coeff[j]*v[opp[j]].dot(wallSpeed)
                        //                        // variation in Bounce-Back due to moving object
                        //                        //double velComp=(slipVelocity).dot(v[opp[j]]);
                        //                        //if (j>6){
                        //                        const double BBi = BBCoeff * nodeHere->n * coeff[j] * slipVelocity.dot(v[opp[j]]);
                        //                        //cout<<velComp<<" "<<newDistribution<<" "<<BBi<<endl;
                        //                        //}
                        //
                        //                        nodeHere->f[opp[j]] = newDistribution + BBi;
                        //                        //nodeHere->f[opp[j]] += nodeHere->fs[j];
                        //                        //                        for (int j1 = 0; j1 < lbmDirec; ++j1) {
                        //                        //                            nodeHere->f[j1]+=-coeff[j] *BBi;
                        //                        //                        }
                        break;


                        // adding the extra mass to the surplus
                        //extraMass += nodeHere->mass * (chi * nodeHere->fs[j] - chi * fStar); //+ BBi);
                    }
                        //                    case TOPO: // OLD!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        //                    {
                        //                        // velocity of the wall
                        //                        const tVect vel = nodes[link]->u;
                        //                        // variation in Bounce-Back due to moving object
                        //                        const double BBi = BBCoeff * nodes[it]->n * coeff[j] * vel.dot(v[j]); // mass!!!!!
                        //                        // get value for chi coefficient
                        //                        const double chi = nodes[link]->curved->getChi(opp[j], 0.5 + 3.0 * nodes[it]->visc);
                        //                        // get fluid reference velocity
                        //                        tVect ubf(0.0, 0.0, 0.0);
                        //                        if (nodes[link]->curved->delta[opp[j]] >= 0.5) {
                        //                            ubf = nodes[link]->curved->m1[opp[j]] * nodes[it]->u + nodes[link]->curved->m2[opp[j]] * vel;
                        //                        } else {
                        //                            const unsigned int linkBack = nodes[it]->d[opp[j]];
                        //                            if (types[linkBack].isFluid() && nodes[linkBack] != 0) {
                        //                                ubf = nodes[linkBack]->u;
                        //                            } else {
                        //                                ubf = vel;
                        //                            }
                        //                        }
                        //                        // calculate missing distribution
                        //                        const double usq = nodes[it]->u.norm2();
                        //                        const double vu = nodes[it]->u.dot(v[j]);
                        //                        const double fStar = nodes[it]->n * coeff[j]*(1 + C1 * v[j].dot(ubf) + C2 * vu * vu + C3 * usq);
                        //
                        //                        nodes[it]->f[opp[j]] = (1.0 - chi) * nodes[it]->fs[j] + chi * fStar - BBi;
                        //                        // adding the extra mass to the surplus
                        //                        extraMass += nodes[it]->mass * (chi * nodes[it]->fs[j] - chi * fStar + BBi);
                        //                    }
                    case OUTLET:
                    {
                        //const unsigned int backLink = nodes[it]->d[opp[j]];
                        // simple bounce-back
                        //nodes[it]->f[opp[j]] = nodes[it]->fs[opp[j]];
                        nodeHere->f[opp[j]] = min(nodeHere->fs[opp[j]], nodeHere->fs[j]);
                        break;
                    }
                        // for moving walls there is simple bounce-back with velocity correction
                    case DYN_WALL:
                    {
                        // getting the index of the wall to compute force in the right object
                        const int solidIndex = linkNode->getSolidIndex();
                        // velocity of the wall
                        const tVect vel = linkNode->u;
                        // variation in Bounce-Back due to moving object
                        const double BBi = BBCoeff * nodeHere->n * coeff[j] * vel.dot(v[j]); // mass!!!!!

                        // static pressure is subtracted in order to correctly compute buoyancy for floating objects
                        const tVect BBforce = nodeHere->bounceBackForce(j, staticPres, BBi, fluidMaterial.earthPressureCoeff, fluidMaterial.lbMaxVisc);
                        // updating force and torque on the object (lattice units). This point is critical since many nodes update the force on the same object (lattice units)
#pragma omp ordered
                        {
                            walls[solidIndex].FHydro += BBforce;
                        }
                        nodeHere->f[opp[j]] = nodeHere->fs[j] - BBi;
                        // adding the extra mass to the surplus
                        extraMass = BBi * nodeHere->mass;
                        break;
                    }// for walls there is simple bounce-back
                    case OBJ:
                    {
                        // getting the index of the wall to compute force in the right object
                        const int solidIndex = linkNode->getSolidIndex();
                        //cout<<"solidIndex "<<solidIndex<<endl;
                        // static pressure is subtracted in order to correctly compute buoyancy for floating objects
                        const tVect BBforce = nodeHere->bounceBackForce(j, staticPres, 0.0, fluidMaterial.earthPressureCoeff, fluidMaterial.lbMaxVisc);
                        // updating force and torque on the object (lattice units). This point is critical since many nodes update the force on the same object (lattice units)
#pragma omp ordered
                        {
                            objects[solidIndex].FHydro += BBforce;
                        }
                        nodeHere->f[opp[j]] = nodeHere->fs[j];
                        break;
                    }
                    case SLIP_STAT_WALL:
                    {
                        if (j > 6) {
                            bool active1 = false;
                            bool active2 = false;
                            const node* nodeCheck1 = nodeHere->d[slip1Check[j]];
                            const node* nodeCheck2 = nodeHere->d[slip2Check[j]];
                            // check for the environment
                            if (nodeCheck1 != 0) {
                                if (nodeCheck1->isActive()) {
                                    active1 = true;
                                }
                            }
                            if (nodeCheck2 != 0) {
                                if (nodeCheck2->isActive()) {
                                    active2 = true;
                                }
                            }
                            // given the environment, perform the right operation
                            if (active1 && !active2) {
                                //cout<<"S1 "<<S1<<" S2 "<<S2<<endl;;
                                // first
                                nodeHere->f[opp[j]] = S1 * nodeCheck1->fs[slip1[j]] + S2 * nodeHere->fs[j];

                            } else if (!active1 && active2) {
                                // second
                                nodeHere->f[opp[j]] = S1 * nodeCheck2->fs[slip2[j]] + S2 * nodeHere->fs[j];
                            } else {
                                // standard BB
                                nodeHere->f[opp[j]] = nodeHere->fs[j];
                            }
                        } else {
                            // standard BB
                            nodeHere->f[opp[j]] = nodeHere->fs[j];
                        }
                        break;
                    }
                    case SLIP_DYN_WALL:
                    {

                        // velocity of the wall
                        const tVect vel = linkNode->u;
                        // variation in Bounce-Back due to moving object
                        const double BBi = BBCoeff * nodeHere->n * coeff[j] * vel.dot(v[j]);
                        if (j > 6) {
                            bool active1 = false;
                            bool active2 = false;
                            const node* nodeCheck1 = nodeHere->d[slip1Check[j]];
                            const node* nodeCheck2 = nodeHere->d[slip2Check[j]];
                            // check for the environment
                            if (nodeCheck1->isActive()) {
                                active1 = true;
                            }
                            if (nodeCheck2->isActive()) {
                                active2 = true;
                            }
                            // given the environment, perform the right operation
                            if (active1 && !active2) {
                                // first
                                nodeHere->f[opp[j]] = S1 * nodeCheck1->fs[slip1[j]] + S2 * (nodeHere->fs[j] - BBi);
                                // adding the extra mass to the surplus
                                extraMass += S2 * nodeHere->mass*BBi;
                            } else if (!active1 && active2) {
                                // second
                                nodeHere->f[opp[j]] = S1 * nodeCheck2->fs[slip2[j]] + S2 * (nodeHere->fs[j] - BBi);
                                // adding the extra mass to the surplus
                                extraMass += S2 * nodeHere->mass*BBi;
                            } else {
                                // standard BB
                                nodeHere->f[opp[j]] = nodeHere->fs[j] - BBi;
                                // adding the extra mass to the surplus
                                extraMass += nodeHere->mass*BBi;
                            }
                        } else {
                            // standard BB
                            nodeHere->f[opp[j]] = nodeHere->fs[j] - BBi;
                            // adding the extra mass to the surplus
                            extraMass += nodeHere->mass*BBi;
                        }
                        break;
                    }
                    default:
                    {
#pragma omp ordered
                        {
                            cout << nodeHere->coord << " " << getX(nodeHere->coord) << " " << getY(nodeHere->coord) << " " << getZ(nodeHere->coord) << " " << typeString(nodeHere->type) << " TYPE ERROR:" << j << endl;
                            for (int j = 1; j < lbmDirec; ++j) {
                                cout << "before error: j=" << j << " link=" << nodeHere->d[j]->coord << endl;
                            }
                            cout << linkNode->coord << " " << getX(linkNode->coord) << " " << getY(linkNode->coord) << " " << getZ(linkNode->coord) << " " << typeString(linkNode->type) << " TYPE ERROR" << endl;

                            exit(0);
                        }
                        break;

                    }
                }
            }
        }
    }


    for (int in = 0; in < activeNodes.size(); ++in) {
        // pointer to active node
        node* nodeHere = activeNodes[in];
        for (int j = 1; j < lbmDirec; ++j) {
            if (nodeHere->f[j] == 0) {
                cout << "Error!" << endl;
            }
        }
    }
    // redistributing extra mass due to bounce back to interface cells
    //redistributeMass(extraMass);

    for (int w = 0; w < walls.size(); ++w) {
        walls[w].FHydro *= unit.Force;
    }
    for (int o = 0; o < objects.size(); ++o) {
        objects[o].FHydro *= unit.Force;
    }
}

//// free surface functions

void LB::updateMass() {
    // refer to the article of Miller or of Svec et al. for this

    // measure time for performance check (begin)
    startUpdateMassStep = std::chrono::steady_clock::now();

    // mass for interface nodes is regulated by the evolution equation
#pragma omp parallel for
    for (int it = 0; it < interfaceNodes.size(); ++it) {
        node* nodeHere = interfaceNodes[it];
        nodeHere->newMass = nodeHere->mass;
        // additional mass streaming to/from interface
        double deltaMass = 0.0;
        //cycling through neighbors
        for (int j = 1; j < lbmDirec; ++j) {
            // getting neighbor index
            const node* linkNode = nodeHere->d[j];
            // average liquid fraction
            if (linkNode == 0) {
                // do nothing
            } else if (linkNode->isInterface()) {
                // average liquid fraction
                const double averageMass = 0.5 * (linkNode->mass / linkNode->n + nodeHere->mass / nodeHere->n);
                deltaMass += averageMass * nodeHere->massStream(j);
                //                deltaMass+=averageMass*(nodeHere->f[opp[j]]-nodeHere->fs[j]);
            } else if (linkNode->isFluid()) {
                const double averageMass = 1.0;
                deltaMass += averageMass * nodeHere->massStream(j);
                //                dummyNode.f[opp[j]]=nodes[link].fs[opp[j]];
                //                deltaMass+=averageMass*(nodeHere->f[opp[j]]-nodeHere->fs[j]);
                //} else if (types[link].isGas()) {
                // averageMass = 0.0;
                //                averageMass=1.0;//*nodeHere->mass;
                //                deltaMass+=averageMass*(nodeHere->f[opp[j]]-nodeHere->fs[j]);
            } else if (linkNode->isDynWall()) {
                const double averageMass = 1.0 * nodeHere->mass; //1.0*nodeHere->mass;//0.5*nodeHere->u.norm()/(nodes[link].u.norm()+nodeHere->u.norm());
                deltaMass += averageMass * nodeHere->massStream(j);
                //              double BBi=2.0*3.0*dummyNode.n*coeff[j]*vel.dot(v[j]);
                //              dummyNode.f[opp[j]]=dummyNode.fs[j]-BBi;
                //                deltaMass+=averageMass*(nodeHere->f[opp[j]]-nodeHere->fs[j]);
                //                cout<<"delta mass ="<<averageMass<<","<<deltaMass<<"\n";
            } else if (linkNode->isCylinderWall()) {
                const double averageMass = 1.0 * nodeHere->mass; //1.0*nodeHere->mass;//0.5*nodeHere->u.norm()/(nodes[link].u.norm()+nodeHere->u.norm());
                deltaMass += averageMass * nodeHere->massStream(j);
                //              double BBi=2.0*3.0*dummyNode.n*coeff[j]*vel.dot(v[j]);
                //              dummyNode.f[opp[j]]=dummyNode.fs[j]-BBi;
                //                deltaMass+=averageMass*(nodeHere->f[opp[j]]-nodeHere->fs[j]);
                //                cout<<"delta mass ="<<averageMass<<","<<deltaMass<<"\n";
            } else if (linkNode->isSlipDynWall()) {
                if (j > 6) {
                    bool active1 = false;
                    bool active2 = false;
                    const node* nodeCheck1 = nodeHere->d[slip1Check[j]];
                    const node* nodeCheck2 = nodeHere->d[slip2Check[j]];
                    // check for the environment
                    if (nodeCheck1 != 0) {
                        if (nodeCheck1->isActive()) {
                            active1 = true;
                        }
                    }
                    if (nodeCheck2 != 0) {
                        if (nodeCheck2->isActive()) {
                            active2 = true;
                        }
                    }
                    // given the environment, perform the right operation
                    double averageMass = 0.0;
                    if (active1 && !active2) {
                        // adding the extra mass to the surplus
                        averageMass += 1.0 * (1.0 - slipCoefficient) * nodeHere->mass;
                    } else if (!active1 && active2) {
                        // adding the extra mass to the surplus
                        averageMass += 1.0 * (1.0 - slipCoefficient) * nodeHere->mass;
                    } else {
                        // adding the extra mass to the surplus
                        averageMass += 1.0 * nodeHere->mass;
                    }
                    deltaMass += averageMass * nodeHere->massStream(j);
                } else {
                    // adding the extra mass to the surplus
                    const double averageMass = 1.0 * nodeHere->mass;
                    deltaMass += averageMass * nodeHere->massStream(j);
                }
            }
        }
        nodeHere->newMass += deltaMass;

        interfaceNodes[it]->mass = interfaceNodes[it]->newMass;
        interfaceNodes[it]->age=min(interfaceNodes[it]->age+ageRatio,1.0);
    }
    // mass for fluid nodes is equal to density
#pragma omp parallel for
    for (int it = 0; it < fluidNodes.size(); ++it) {
        fluidNodes[it]->mass = fluidNodes[it]->n;
        fluidNodes[it]->age=min(fluidNodes[it]->age+ageRatio,1.0);
    }

    // measure time for performance check (begin)
    endUpdateMassStep = std::chrono::steady_clock::now();

}

void LB::updateInterface() {
    // updates the location of the interface, creating and deleting interface nodes

    // measure time for performance check (begin)
    startUpdateInterfaceStep = std::chrono::steady_clock::now();

    // variable for storage of mass surplus
    double massSurplus = 0.0;

    // lists for "mutant" nodes
    filledNodes.clear();
    emptiedNodes.clear();
    newInterfaceNodes.clear();

    // filling lists of mutant nodes and changing their type
    findInterfaceMutants();

    // fixing the interface (always one interface between fluid and gas)
    smoothenInterface(massSurplus);

    // updating characteristics of mutant nodes
    updateMutants(massSurplus);

    // remove isolated interface cells (both surrounded by gas and by fluid)
    removeIsolated(massSurplus);

    // distributing surplus to interface cells
    redistributeMass(massSurplus);

#ifdef DEBUG
    // compute surface normal vectors
    computeSurfaceNormal();
#endif

    // measure time for performance check (begin)
    endUpdateInterfaceStep = std::chrono::steady_clock::now();

}

void LB::findInterfaceMutants() {
    // checks for nodes going through the following transitions
    // interface -> fluid (filled node)
    // interface -> gas (emptied node)

    // measure time for performance check (begin)
    startFindMutantsStep = std::chrono::steady_clock::now();

    filledNodes.clear();
    emptiedNodes.clear();

    // CHECKING FOR MUTANT NODES
    for (int ind = interfaceNodes.size() - 1; ind >= 0; --ind) {
        //const unsigned int i = interfaceNodes[ind];
        node* nodeHere = interfaceNodes[ind];
        // CHECKING FOR NEW FLUID NODES from filling
        if (nodeHere->mass > nodeHere->n) {
            // updating mutant list
            filledNodes.push_back(nodeHere);
            // updating type
            nodeHere->setFluid();
            // updating lists (interface -> fluid)
            fluidNodes.push_back(nodeHere);
            interfaceNodes.erase(interfaceNodes.begin() + ind);
        }// CHECKING FOR NEW GAS NODES from emptying
        else if (nodeHere->mass < 0.0) {
            // updating mutant list
            emptiedNodes.push_back(nodeHere);
            // updating lists (interface ->noList)
            interfaceNodes.erase(interfaceNodes.begin() + ind);
        }
    }

    // measure time for performance check (end)
    endFindMutantsStep = std::chrono::steady_clock::now();

}

void LB::smoothenInterface(double& massSurplus) {


    // measure time for performance check (begin)
    startSmoothenInterfaceStep_1 = std::chrono::steady_clock::now();

    static const double marginalMass = 1.0e-2;
    newInterfaceNodes.clear();
    // CHECKING FOR NEW INTERFACE NODES from neighboring a new fluid node
    for (int it = 0; it < filledNodes.size(); ++it) {
        const node* nodeHere = filledNodes[it];
        // neighor indices
        unsigned int neighborCoord[lbmDirec];
        findNeighbors(neighborCoord, nodeHere);
        // cycling through neighbors
        for (int j = 1; j < lbmDirec; ++j) {
            // neighbor index
            const unsigned int linkIndex = neighborCoord[j];
            // checking if node is gas (so to be tranformed into interface)
            if (nodes.count(linkIndex) == 0) {
                // create new interface node
                generateNode(linkIndex, INTERFACE);
                // point at it to start initializing it
                node* linkNode = &nodes[linkIndex];
                // add it to interface node list
                newInterfaceNodes.push_back(linkNode);
                // node is becoming active and needs to be initialized
                //double newNodeN(0.0);
                //tVect newNodeU(0.0,0.0,0.0);
                //double newNodeVisc(0.0);
                //                double newNodeDistributions[lbmDirec]; 
                //                newNodeDistributions[0]=coeff[0];
                //                for (int j =1; j < lbmDirec; ++j) {
                //                    newNodeDistributions[j]=0.0;
                //                    const node* neighbourNode=linkNode->d[j];
                //                    if (neighbourNode!=nullptr) {
                //                        if (neighbourNode->isActive() && neighbourNode!=linkNode) {
                //                            newNodeDistributions[opp[j]]=neighbourNode->f[j];
                //                        }
                //                        else {
                //                            newNodeDistributions[opp[j]]=coeff[opp[j]];
                //                        }
                //                    } else {
                //                            newNodeDistributions[opp[j]]=coeff[opp[j]];
                //                        }
                //                }

                //averageFromNeighbors(linkNode,newNodeN,newNodeU,newNodeVisc,newNodeDistributions);

                double massSurplusHere = -marginalMass * fluidMaterial.initDensity;

                // same density and velocity; 1% of the mass
                //                if (fluidMaterial.rheologyModel==NEWTONIAN) {
                //linkNode->initialize(fluidMaterial.initDensity, nodeHere->u, -1.0 * massSurplusHere, fluidMaterial.initDynVisc, nodeHere->hydroForce + lbF);
                linkNode->copy(nodeHere);
                linkNode->mass = -massSurplusHere;
                //                }
                //                else {
                //                    linkNode->initialize(fluidMaterial.initDensity, nodeHere->u, -1.0 * massSurplusHere, fluidMaterial.lbMaxVisc, nodeHere->hydroForce + lbF);
                //                }
                //                linkNode->restart(nodeHere->n, nodeHere->u, -1.0*massSurplusHere, nodeHere->visc, newNodeDistributions);

                //                for (int j = 1; j < lbmDirec; ++j) {
                //                    const node* neighbourNode = linkNode->d[j];
                //                    if (neighbourNode != nullptr) {
                //                        if (neighbourNode->isActive() && neighbourNode != linkNode) {
                //                            linkNode->f[opp[j]] = neighbourNode->f[opp[j]];
                //                        }
                //                    }
                //                }


                //nodes[linkIndex].initialize(nodeHere->n, nodeHere->u, -1.0*massSurplusHere, nodeHere->visc, nodeHere->hydroForce + lbF);

                //nodes[linkIndex].restart(newNodeN, newNodeU, -1.0*massSurplusHere, newNodeVisc, newNodeDistributions);
                //nodes[linkIndex].initialize(newNodeN, newNodeU, -1.0*massSurplusHere, newNodeVisc, Zero);
                //nodes[link].restart(nodeHere->n, nodeHere->u, -1.0*massSurplusHere, nodeHere->visc, nodeHere->f);
                // nodes[link].restart(nodeHere->n, nodeHere->u, -1.0*massSurplusHere, nodeHere->visc, nodeHere->f);
                //*nodes[link]=*nodes[index];
                //nodes[link]->mass=0.0001 * fluidMaterial.initDensity;
                // the 1% of the mass is taken form the surplus
                nodes[linkIndex].scatterMass(massSurplusHere);
                //massSurplus += massSurplusHere;
            }
        }
    }
    // measure time for performance check (end)
    endSmoothenInterfaceStep_1 = std::chrono::steady_clock::now();

    // measure time for performance check (begin)
    startSmoothenInterfaceStep_2 = std::chrono::steady_clock::now();

    // CHECKING FOR NEW INTERFACE NODES from neighboring a new gas node
    // tested unordered_set, was slower
    std::set<unsigned int> nodesToErase;
    for (int it = 0; it < emptiedNodes.size(); ++it) {
        // empied node
        node* nodeHere = emptiedNodes[it];
        // neighbors of emptied node
        unsigned int neighborCoord[lbmDirec];
        findNeighbors(neighborCoord, nodeHere);
        // cycling through neighbors
        for (int j = 1; j < lbmDirec; ++j) {
            // neighbor node
            const unsigned int link = neighborCoord[j];
            if (nodes.count(link) != 0) {
                node* linkNode = nodeHere->d[j];
                if (linkNode->isFluid()) {
                    linkNode->setInterface();
                    newInterfaceNodes.push_back(linkNode);
                    double massSurplusHere = marginalMass * linkNode->n;
                    // characteristics are inherited by previous fluid cell. Only mass must be updated to 99% of initial mass
                    linkNode->mass = linkNode->n - massSurplusHere;
                    // the remaining 1% of the mass is added to the surplus
                    linkNode->scatterMass(massSurplusHere);
                    //massSurplus += massSurplusHere;
                    // Store the node we want to erase later
                    nodesToErase.insert(link);

                }
            }
        }
    }
    // Process and remove all fluid nodes that have been converted
    for (int ind = fluidNodes.size() - 1; ind >= 0; --ind) {
        const unsigned int i = fluidNodes[ind]->coord;
         if (nodesToErase.find(i) != nodesToErase.end()) {
             fluidNodes.erase(fluidNodes.begin() + ind);
         }
    }
    
    // measure time for performance check (end)
    endSmoothenInterfaceStep_2 = std::chrono::steady_clock::now();

}

void LB::averageFromNeighbors(node* linkNode, double& newNodeN, tVect& newNodeU, double& newNodeVisc, double newNodeDistributions[]) {

    double fluidNeighbors(0.0);
    double wallNeighbors(0.0);
    //
    //unsigned int neighborCoord2[lbmDirec];
    //findNeighbors(neighborCoord2, linkNode);
    //                bool is_in(false);
    for (int j = 1; j < lbmDirec; ++j) {
        const node* linkNode2 = linkNode->d[j];
        //const unsigned int linkIndex2 = neighborCoord2[j2];
        if (linkNode2 != nullptr) {
            if (linkNode2->isWall()) {
                //newNodeN += 1.0;
                //newNodeVisc += fluidMaterial.lbMinVisc;
                wallNeighbors += coeff[j];
                for (int j2 = 0; j2 < lbmDirec; ++j2) {
                    newNodeDistributions[j2] += coeff[j] * coeff[j2];
                }
            } else if (linkNode2->isActive() && linkNode2->n > 0) {
                newNodeN += coeff[j] * linkNode2->n;
                newNodeU += coeff[j] * linkNode2->u;
                newNodeVisc += coeff[j] * linkNode2->visc;
                for (int j2 = 0; j2 < lbmDirec; ++j2) {
                    newNodeDistributions[j2] += coeff[j] * linkNode2->f[j2];
                }
                fluidNeighbors += coeff[j];
            } else if (linkNode2 != linkNode) {
                cout << "Big problem at generating node " << linkNode->coord << " (" << getX(linkNode->coord) << " " << getY(linkNode->coord) << " " << getZ(linkNode->coord) << ")" << endl;
                cout << "From properties of node " << linkNode2->coord << " (" << getX(linkNode2->coord) << " " << getY(linkNode2->coord) << " " << getZ(linkNode2->coord) << ") of type " << typeString(linkNode2->type) << " n= " << linkNode2->n << endl;
                exit(0);
            }
        }
    }
    if (fluidNeighbors > 0) {
        //cout<<"new neigh"<<newNeighbors<<endl;
        newNodeN /= fluidNeighbors;
        newNodeU /= (fluidNeighbors + wallNeighbors);
        newNodeVisc /= fluidNeighbors;
        for (int j = 0; j < lbmDirec; ++j) {
            newNodeDistributions[j] /= (fluidNeighbors + wallNeighbors);
        }
        //cout<<" new node n="<<newNodeU.dot(Xp)<<" number="<<newNeighbors<<endl;
    } else {
        cout << "Big problem:" << endl;
        //<< nodeHere->coord << "(type " << typeString(nodeHere->type) << ") pointing to " << linkNode->coord << "or " << linkIndex << " (type " << typeString(linkNode->type) << ") pointing back to " << linkNode->d[opp[j]]->coord << "or " << neighborCoord2[opp[j]] << " (type " << typeString(linkNode->d[opp[j]]->type) << ")" << endl;
    }
}

void LB::updateMutants(double& massSurplus) {

    // measure time for performance check (begin)
    startUpdateMutantsStep = std::chrono::steady_clock::now();

    for (int it = 0; it < emptiedNodes.size(); ++it) {
        node* nodeHere = emptiedNodes[it];
    }
    // resetting new gas macroscopic quantities
    for (int it = 0; it < emptiedNodes.size(); ++it) {
        node* nodeHere = emptiedNodes[it];

        bool isNotInterface = true;
        for (int it1 = 0; it1 < newInterfaceNodes.size(); ++it1) {
            if (nodeHere == newInterfaceNodes[it1]) {
                isNotInterface = false;
            }
        }
        if (isNotInterface) {

            // updating mass surplus
            massSurplus += nodeHere->mass; // maybe problems
            // deleting node
            eraseNode(nodeHere->coord);
        }
    }

    // resetting new fluid macroscopic quantities
    for (int it = 0; it < filledNodes.size(); ++it) {
        node* nodeHere = filledNodes[it];
        bool isNotInterface = true;
        for (int it1 = 0; it1 < newInterfaceNodes.size(); ++it1) {
            if (nodeHere == newInterfaceNodes[it1]) {
                isNotInterface = false;
            }
        }
        if (isNotInterface) {
            // updating mass surplus
            massSurplus += nodeHere->mass - nodeHere->n;
            // setting liquid fraction for new fluid cell (other macroscopic characteristics stay the same)
            nodeHere->mass = nodeHere->n;
        }
    }

    // resetting neighbors for newly created interface cells
    for (int it = 0; it < newInterfaceNodes.size(); ++it) {
        interfaceNodes.push_back(newInterfaceNodes[it]);
    }

    // measure time for performance check (end)
    endUpdateMutantsStep = std::chrono::steady_clock::now();


}

void LB::removeIsolated(double& massSurplus) {
    // remove isolated interface cells (surrounded either by only fluid or only solid cells)

    // measure time for performance check (start)
    startRemoveIsolatedStep = std::chrono::steady_clock::now();

    // checking if it is surrounded by fluid (in that case is converted to fluid). Solid is an exception
    // reverse cycle is needed because of deletion function
    for (int ind = interfaceNodes.size() - 1; ind >= 0; --ind) {
        node* nodeHere = interfaceNodes[ind];
        bool surroundedFluid = true;
        for (int j = 1; j < lbmDirec; ++j) {
            const node* linkNode = nodeHere->d[j];
            if (linkNode == 0) {
                surroundedFluid = false;
                break;
            }
        }
        if (surroundedFluid) {
            //            cout<<nodes[i].x<<" "<<nodes[i].y<<" "<<nodes[i].z<<" NEW FLUID NODE from surrounding\n";
            // update mass storage for balance
            massSurplus += nodeHere->mass - nodeHere->n;
            // update characteristics (inherited from the gas node)
            nodeHere->mass = nodeHere->n;
            nodeHere->setFluid();
            interfaceNodes.erase(interfaceNodes.begin() + ind);
            fluidNodes.push_back(nodeHere);
        }
    }

    // checking if it is surrounded by gas (in that case is converted to gas)
    // or, better, if it is not connected to fluid (could be connected to walls or particles)
    for (int ind = interfaceNodes.size() - 1; ind >= 0; --ind) {
        const node* nodeHere = interfaceNodes[ind];
        bool surroundedGas = true;
        for (int j = 1; j < lbmDirec; ++j) {
            const node* linkNode = nodeHere->d[j];
            if (linkNode != 0) {
                if (linkNode->isFluid()) {
                    surroundedGas = false;
                    break;
                }
            }
        }
        // updating mass surplus
        if (surroundedGas) {
            //            cout<<nodes[i].x<<" "<<nodes[i].y<<" "<<nodes[i].z<<" NEW GAS NODE from surrounding\n";
            // update mass
            massSurplus += nodeHere->mass;
            interfaceNodes.erase(interfaceNodes.begin() + ind);
            eraseNode(nodeHere->coord);
        }
    }

    // measure time for performance check (end)
    endRemoveIsolatedStep = std::chrono::steady_clock::now();

}

void LB::redistributeMass(const double& massSurplus) {
    // redistribute the mass surplus among interface cells

    // measure time for performance check (begin)
    startRedistributeMassStep = std::chrono::steady_clock::now();

    const double addMass = massSurplus / double(interfaceNodes.size());

#pragma omp parallel for
    for (int it = 0; it < interfaceNodes.size(); ++it) {
        interfaceNodes[it]->mass += addMass;
    }

    // measure time for performance check (end)
    endRedistributeMassStep = std::chrono::steady_clock::now();
}

void LB::computeSurfaceNormal() {
    // compute the surface normal, as the gradient of the liquidFraction function

    // measure time for performance check (begin)
    startComputeNormalsStep = std::chrono::steady_clock::now();

#pragma omp parallel for
    for (int it = 0; it < interfaceNodes.size(); ++it) {
        const node* nodeHere = interfaceNodes[it];
        tVect surfaceNormalHere = Zero;
        for (int j = 1; j < lbmDirec; ++j) {
            const node* linkNode = nodeHere->d[j];
            if (linkNode != nullptr) {
                if (linkNode->isActive()) {
                    surfaceNormalHere -= vDirec[j] * linkNode->mass;
                } else if (linkNode->isWall()) {
                    surfaceNormalHere -= vDirec[j] * 1.0;
                }
            }
        }
        surfaceNormalHere /= surfaceNormalHere.norm();
        interfaceNodes[it]->surfaceNormal = surfaceNormalHere;
    }
#pragma omp parallel for
    for (int it = 0; it < fluidNodes.size(); ++it) {
        fluidNodes[it]->surfaceNormal = Zero;
    }

    // measure time for performance check (end)
    endComputeNormalsStep = std::chrono::steady_clock::now();
}

void LB::enforceMassConservation() {
    // calculate total mass
    double thisMass = 0.0;
    for (int ia = 0; ia < activeNodes.size(); ++ia) {
        if (!activeNodes[ia]->isInsideParticle()) {
            thisMass += activeNodes[ia]->mass;
        }
    }

    // mass deficit
    const double massDeficit = (thisMass - totalMass);
    // cout<<endl<<"This="<<thisMass<<" tot="<<totalMass;

    // fix it
    redistributeMass(-0.01 * massDeficit);

}

//// particle coupling functions

void LB::computeHydroForces(node* nodeHere, elmtList& elmts, particleList& particles) {
    // resetting hydrodynamic forces on nodes
    nodeHere->hydroForce.reset();
    if (nodeHere->isInsideParticle()) { // && types[index].isFluid()
        // getting the index of the particle to compute force in the right object
        const unsigned int index = nodeHere->coord;
        const unsigned int particleIndex = nodeHere->getSolidIndex();
        const unsigned int clusterIndex = particles[particleIndex].clusterIndex;
        // calculating velocity of the solid boundary at the node (due to rotation of particles)
        // vectorized radius (real units)
        const tVect radius = getPosition(index) - particles[particleIndex].x0 / unit.Length + particles[particleIndex].radiusVec / unit.Length;
        // update velocity of the particle node (u=v_center+omega x radius) (real units)
        const tVect localVel = elmts[clusterIndex].x1 / unit.Speed + (elmts[clusterIndex].wGlobal.cross(radius)) / unit.AngVel;


            // calculate differential velocity
            const tVect diffVel = nodeHere->age*nodeHere->age*nodeHere->liquidFraction()*(nodeHere->u - localVel);

            // force on fluid
             nodeHere->hydroForce += -1.0 * diffVel; // -1.0*nodes[index]->liquidFraction()*diffVel;

        //             // force on fluid
        //            # pragma omp critical
        //            {
        //                for (int j=0; j<lbmDirec; ++j) {
        //                    const unsigned int link=neighbors[index].d[j];
        //                    if (types[link].isActive()) {
        //                        nodes[link]->hydroForce+=-1.0/nodes[index]->liquidFraction()*coeff[j]*diffVel;
        //                    }
        //                }
        //            }

        // force on particle
#pragma omp critical
        {
            elmts[clusterIndex].fluidVolume += nodeHere->mass;
            elmts[clusterIndex].FHydro += 1.0 * diffVel;
            elmts[clusterIndex].MHydro += 1.0 * radius.cross(diffVel);
        }
    }
}

void LB::findNewActive(nodeList& newPopUpNodes, elmtList& elmts, particleList& particles) {

    // SOLID TO ACTIVE CHECK
    // cycling through particle nodes
    //    #pragma omp parallel for ordered
    for (int it = 0; it < activeNodes.size(); ++it) {
        node* nodeHere = activeNodes[it];
        if (nodeHere->isInsideParticle()) {
            const unsigned int nodeIndexHere = nodeHere->coord;
            const tVect nodePosition = getPosition(nodeIndexHere);
            // solid index to identify cluster
            const unsigned int particleIndex = nodeHere->getSolidIndex();
            const unsigned int clusterIndex = particles[particleIndex].clusterIndex;
            // in this case check if it has been uncovered (must be out of all particles of the cluster) - we start with a true hypothesis
            bool newActive = true;
            // cycling through component particles
            unsigned int componentIndex = 0;
            for (int j = 0; j < elmts[clusterIndex].components.size(); ++j) {
                // getting indexes from particle composing the cluster
                componentIndex = elmts[clusterIndex].components[j];
                // checking if it has been uncovered in component j of the cluster
                // radius need to be increased by half a lattice unit
                // this is because solid boundaries are located halfway between solid and fluid nodes
                if (nodePosition.insideSphere(particles[componentIndex].x0 / unit.Length, hydrodynamicRadius * particles[componentIndex].r / unit.Length)) { //-0.5?
                    // if the node is still inside the element, the hypothesis of new active is not true anymore
                    newActive = false;
                    // and we can get out of the cycle
                    break;
                }
            }
            if (newActive) {
                // turning up the cell
                //            #pragma omp ordered
                nodeHere->setOutsideParticle();
                newPopUpNodes.push_back(nodeHere);
                //            cout<<"new active\n";
            }
        }
    }

    //    for (int it = 0; it < newPopUpNodes.size(); ++it) {
    //        //        cout<<"New active\n";
    //        for (int ind = particleNodes.size() - 1; ind >= 0; --ind) {
    //            if (newPopUpNodes[it] == particleNodes[ind]) {
    //                // deleting the cell from particle list
    //                particleNodes.erase(particleNodes.begin() + ind);
    //            }
    //        }
    //    }
}

void LB::findNewSolid(nodeList& newSolidNodes, elmtList& elmts, particleList& particles) {

    // ACTIVE TO SOLID CHECK
    // we check only first order neighbors of particle nodes. This is done in order to avoid cycling through all active cells

    //const unsigned int particleNodesBeforeCheck=particleNodes.size();

    //    #pragma omp parallel for ordered
    for (int it = 0; it < activeNodes.size(); ++it) {
        node* nodeHere = activeNodes[it];
        if (nodeHere->isInsideParticle()) {
            //const unsigned int nodeIndexHere = nodeHere->coord;
            //    for (intlist::iterator it=particleNodes.begin(); it<particleNodes.end(); ++it) {
            // solid index to identify cluster
            const unsigned int particleIndex = nodeHere->getSolidIndex();
            const unsigned int clusterIndex = particles[particleIndex].clusterIndex;
            // cycle through first neighbors
            for (int k = 1; k < lbmMainDirec; ++k) {
                node* linkNode = nodeHere->d[k];
                if (linkNode != 0) {
                    const unsigned int linkIndex = linkNode->coord;
                    // checking if solid particle is close to an active one -> we have an active node to check
                    if (!linkNode->isInsideParticle() && linkNode->isActive()) {
                        const tVect linkPosition = getPosition(linkIndex);
                        // check if neighbors has been covered (by any of the particles of the cluster) - we start with a false hypothesis
                        bool newSolid = false;
                        unsigned int componentIndex = 0;
                        // cycling through all components of the cluster
                        for (int j = 0; j < elmts[clusterIndex].components.size(); ++j) {
                            // getting component particle index
                            componentIndex = elmts[clusterIndex].components[j];
                            // check if it getting inside
                            // radius need to be increased by half a lattice unit
                            // this is because solid boundaries are located halfway between soli and fluid nodes
                            if (linkPosition.insideSphere(particles[componentIndex].x0 / unit.Length, hydrodynamicRadius * particles[componentIndex].r / unit.Length)) { //-0.5?
                                // if so, then the false hypothesis does not hold true anymore
                                newSolid = true;
                                // and we exit the cycle
                                break;
                            }
                        }
                        // an active cell can become part of multiple particle at the same time
                        // we need to check if it is already on the list
                        if (newSolid) {
                            // check if we are creating a duplicate - we start with a false hypothesis
                            bool alreadyInside = false;
                            for (int it1 = 0; it1 < newSolidNodes.size(); ++it1) {
                                //                    for (intlist::iterator it1=newSolidNodes.begin(); it1!=newSolidNodes.end(); ++it1) {
                                if (linkNode == newSolidNodes[it1]) {
                                    // if it is already in the list than the true hypothesis does not hold true anymore
                                    alreadyInside = true;
                                    // and we exit the cycle
                                    break;
                                }
                            }
                            // at this point, if it is both newSolid and not alreadyInside we can add it to particle nodes
                            if (!alreadyInside) {
                                // solid index is assigned here to avoid sending this information to the switching functions
                                //                        # pragma omp ordered
                                {
                                    linkNode->setSolidIndex(componentIndex);
                                    linkNode->setInsideParticle();
                                    // sending node index for updating
                                    newSolidNodes.push_back(linkNode);
                                    //particleNodes.push_back(linkNode);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void LB::updateEnergy(double& totalKineticEnergy) {

    //resetting energy
    fluidEnergy.reset();
    fluidImmersedEnergy.reset();

    // kinetic energy
    double tKin(0.0), tpKin(0.0), mass(0.0), massp(0.0);

    for (int it = 0; it < activeNodes.size(); ++it) {
        node* nodeHere = activeNodes[it];
        if (nodeHere->isInsideParticle()) {
            tpKin += 0.5 * nodeHere->mass * nodeHere->u.norm2();
            massp += nodeHere->mass;
        } else {
            tKin += 0.5 * nodeHere->mass * nodeHere->u.norm2();
            mass += nodeHere->mass;
        }
    }

    fluidEnergy.rotKin = 0.0;
    fluidEnergy.trKin = tKin;
    fluidEnergy.mass = mass;

    fluidImmersedEnergy.rotKin = 0.0;
    fluidImmersedEnergy.trKin = tpKin;
    fluidImmersedEnergy.mass = massp;

    // potential energy
    // defining reference plane
    wall zeroWall;
    double g(0.0), gp(0.0);
    const double gravityNorm = lbF.norm();
    if (gravityNorm != 0.0) {
        zeroWall.n = -1.0 * lbF / gravityNorm;
        zeroWall.p = tVect(0.0, 0.0, 0.0);
        if (zeroWall.n.dot(Xp) < 0) {
            zeroWall.p += tVect(double(lbSize[0]), 0.0, 0.0);
        }
        if (zeroWall.n.dot(Yp) < 0) {
            zeroWall.p += tVect(0.0, double(lbSize[1]), 0.0);
        }
        if (zeroWall.n.dot(Zp) < 0) {
            zeroWall.p += tVect(0.0, 0.0, double(lbSize[2]));
        }
        for (int it = 0; it < activeNodes.size(); ++it) {
            node* nodeHere = activeNodes[it];
            if (nodeHere->isInsideParticle()) {
                const double heigth = zeroWall.dist(getPosition(nodeHere->coord));
                gp += nodeHere->mass * heigth*gravityNorm;
            } else {
                const double heigth = zeroWall.dist(getPosition(nodeHere->coord));
                g += nodeHere->mass * heigth*gravityNorm;
            }
        }
    }
    fluidEnergy.grav = g;
    fluidImmersedEnergy.grav = gp;

    // elastic is missing

    // total
    fluidEnergy.updateTotal();
    fluidImmersedEnergy.updateTotal();
    
    totalKineticEnergy=totalKineticEnergy+fluidEnergy.trKin+fluidImmersedEnergy.trKin;

}

//void LB::solidToActive(unsIntList& newPopUpNodes, elmtList& elmts, double& massSurplus) {
//
//    //    for (int it=0; it<newPopUpNodes.size(); ++it) {
//    ////        cout<<"New active\n";
//    //        for (int ind=particleNodes.size()-1; ind>=0; --ind) {
//    //            if (newPopUpNodes[it]==particleNodes[ind]) {
//    //                // deleting the cell from particle list
//    //                particleNodes.erase(particleNodes.begin()+ind);
//    //            }
//    //        }
//    //    }
//    //
//    //    // PRELIMINARY STEPS  ////////////////////////////////////////////////////////////
//    //    // Initializing types for new pop-up
//    //    for (int it=0; it<newPopUpNodes.size(); ++it) {
//    //        // index of the popping up node
//    //        unsigned int index=newPopUpNodes[it];
//    //        // routine to check if it should be interface, fluid or gas
//    //        bool outOfFluid=true;
//    //        bool closeToGas=false;
//    //        for (int j=1; j<lbmDirec; ++j) {
//    //            unsigned int link=neighbors[index].d[j];
//    //            if (types[link].isFluid()) {
//    //                outOfFluid=false;
//    //            }
//    //            else if(types[link].isGas()) {
//    //                closeToGas=true;
//    //            }
//    //        }
//    //        if (outOfFluid) {
//    //            //then node is fluid
//    //            types[index].setGas();
//    //        }
//    //        if ((!outOfFluid)&&(closeToGas)) {
//    //            //then node is interface
//    //            types[index].setInterface();
//    //        }
//    //        if ((!outOfFluid)&&(!closeToGas)) {
//    //            // then node is fluid
//    //            types[index].setFluid();
//    //        }
//    //    }
//    //
//    //    // NEW POP-UP NODES - Initialization  ////////////////////////////////////////////////////////////
//    //   for (int it=0; it<newPopUpNodes.size(); ++it) {
//    //       unsigned int index=newPopUpNodes[it];
//    //       // check if cell is active (it could be gas) - the type has already been assigned in the lines over here
//    //        if (types[index].isActive()) {
//    ////            cout<<nodes[index].x<<" "<<nodes[index].y<<" "<<nodes[index].z<<" NEW ACTIVE NODE from pop-up\n";
//    //            //  macroscopic velocity and density for the new cell (average of neighbors)
//    //            tVect uAverage;
//    //            double nAverage=0.0;
//    //            double viscAverage=0.0;
//    //            double mass=0.0;
//    //
//    //            // Routine for average properties of new nodes. CHECK!!
//    //            // number of neighbor active cells
//    //            unsigned int ku=0;
//    //            unsigned int kn=0;
//    //            uAverage.reset();
//    //            // calculating average
//    //            for (int j=0; j<lbmDirec; j++) {
//    //                unsigned int link=neighbors[index].d[j];
//    //                // checking if neighbor is active
//    //                if (types[link].isActive()) {
//    //                    //we need to check if this cell is a new cell, in that case its properties are not good to use
//    //                    bool isNew=false;
//    //                    for (int it1=0; it1<newPopUpNodes.size(); ++it1) {
//    ////                    for (intlist::iterator it1=newPopUpNodes.begin(); it1<newPopUpNodes.end(); ++it1) {
//    //                        if (link==newPopUpNodes[it1]) {
//    //                            isNew=true;
//    //                            break;
//    //                        }
//    //                    }
//    //                    if (!isNew) {
//    //                         // incrementing velocity average
//    //                        uAverage=uAverage+nodes[link]->u;
//    //                        nAverage+=nodes[link]->n;
//    //                        viscAverage+=nodes[link]->visc;
//    //                        // incrementing number active neighbors
//    //                        ++ku;
//    //                        ++kn;
//    //                    }
//    //                }
//    ////                if (nodes[link].isParticle()) {
//    ////                    // incrementing velocity average
//    ////                    // getting local velocity
//    ////                    tVect radius=tVect(((double)nodes[link].x), ((double)nodes[link].y), ((double)nodes[link].z))-elmts[nodes[link].solidIndex].x0/unit.Length;
//    ////                    // local velocity (v_local=v_center+omega x radius) (lattice units)
//    ////                    tVect localVel=(elmts[nodes[link].solidIndex].x1/unit.Speed+(elmts[nodes[link].solidIndex].w.cross(radius))/unit.AngVel);
//    ////                    uAverage=uAverage+localVel;
//    ////                    // incrementing number of active neighbors
//    ////                    ++ku;
//    ////                }
//    //            }
//    //            if (ku==0) {
//    //                uAverage.reset();
//    //            }
//    //            else {
//    //                // dividing to get true average
//    //                uAverage=uAverage/ku;
//    //            }
//    //            if (kn==0) {
//    //                viscAverage=initDynVisc;
//    //                nAverage=1.0;
//    //            }
//    //            else {
//    //                // dividing to get true average
//    //                nAverage/=kn;
//    //                viscAverage/=kn;
//    //            }
//    //
//    //            if (types[index].isInterface()) {
//    ////                activeNodes.push_back(nodes[index].index);
//    //                interfaceNodes.push_back(index);
//    //                mass=0.01*nAverage;
//    //            }
//    //            else if (types[index].isFluid()) {
//    //                fluidNodes.push_back(index);
//    ////                activeNodes.push_back(nodes[index].index);
//    //                mass=nAverage;
//    //            }
//    //            nodes[index]->initialize(nAverage, uAverage, mass, viscAverage, lbF);
//    //            massSurplus-=mass;
//    ////            initializeNode(nodes[index], nAverage, uAverage, mass, viscAverage);
//    //            unsigned int zero=0;
//    //            types[index].setSolidIndex(zero);
//    //        }
//    //        else {
//    //            // remove node
//    //            delete nodes[index];
//    //            nodes[index]=0;
//    //        }
//    //    }
//
//}
//
//void LB::activeToSolid(unsIntList& newSolidNodes, elmtList& elmts, double& massSurplus) {
//
//    //    for (int it=0; it<newSolidNodes.size(); ++it) {
//    ////    for (intlist::iterator it=newSolidNodes.begin(); it<newSolidNodes.end(); ++it) {
//    ////        cout<<"New solid\n";
//    //        //adding cell to solid list
//    //        particleNodes.push_back(newSolidNodes[it]);
//    //        // deleting the cell from active lists (fluid or interface, depending on the type)
//    //        if (types[newSolidNodes[it]].isFluid()) {
//    //            for (int ind=fluidNodes.size()-1; ind>=0; --ind) {
//    //                if (newSolidNodes[it]==fluidNodes[ind]) {
//    //                    fluidNodes.erase(fluidNodes.begin()+ind);
//    //                    break;
//    //                }
//    //            }
//    //        }
//    //        else if (types[newSolidNodes[it]].isInterface()) {
//    //            for (int ind=interfaceNodes.size()-1; ind>=0; --ind) {
//    //                if (newSolidNodes[it]==interfaceNodes[ind]) {
//    //                    interfaceNodes.erase(interfaceNodes.begin()+ind);
//    //                    break;
//    //                }
//    //            }
//    //        }
//    //    }
//    //    // initializing type for new solid nodes and redistributing mass in case of interface
//    //    for (int it=0; it<newSolidNodes.size(); ++it) {
//    //        // index of new solid node
//    //        unsigned int index=newSolidNodes[it];
//    //        massSurplus+=nodes[index]->mass;
//    //        // initializing type
//    //        types[index].setParticle();
//    //    }
//    //
//    //    // NEW SOLID NODES - Initialization ///////////////////////////////////////////////////////////
//    //    for (int it=0; it<newSolidNodes.size(); ++it) {
//    //        unsigned int index=newSolidNodes[it];
//    ////        cout<<nodes[index].x<<" "<<nodes[index].y<<" "<<nodes[index].z<<" NEW SOLID NODE from absorption\n";
//    //        // reset velocity and mass (useful for plotting)
//    //        nodes[index]->n=0.0;
//    //        nodes[index]->visc=0.0;
//    //        nodes[index]->u.reset();
//    //        nodes[index]->mass=0.0;
//    //        nodes[index]->shearRate=0.0;
//    ////        nodes[index].force.reset();
//    //
//    //        // solid index is added in the movement function
//    //    }
//
//}
//
//// functions for index management
//

double LB::maxHeight(const unsigned int& dir) const {
    // dir can be 0, 1 2 (x, y, z)

    unsigned int maxHeight = 0;
    // find maximum height in a specified direction

    double lastFraction = 0.0;

    if (dir == 0) {
        for (int it = 0; it < activeNodes.size(); ++it) {
            const unsigned int index = activeNodes[it]->coord;
            const unsigned int xHere = getPositionX(index);
            if (xHere > maxHeight) {
                maxHeight = xHere;
                lastFraction = activeNodes[it]->mass;
            }
        }
    } else if (dir == 1) {

        for (int it = 0; it < activeNodes.size(); ++it) {
            const unsigned int index = activeNodes[it]->coord;
            const unsigned int yHere = getPositionY(index);
            if (yHere > maxHeight) {
                maxHeight = yHere;
                lastFraction = activeNodes[it]->mass;
            }
        }
    } else if (dir == 2) {

        for (int it = 0; it < activeNodes.size(); ++it) {
            const unsigned int index = activeNodes[it]->coord;
            const unsigned int zHere = getPositionZ(index);
            if (zHere > maxHeight) {
                maxHeight = zHere;
                lastFraction = activeNodes[it]->mass;
            }
        }
    }

    return double(maxHeight) - 1.0 + lastFraction;
}

unsigned int LB::getIndex(const unsigned int& x, const unsigned int& y, const unsigned int& z) {
    return x + y * lbSize[0] + z * lbSize[0] * lbSize[1];
}

tVect LB::getPosition(const unsigned int& index) const {
    unsigned int x, y, z;

    // index is calculated in this fashion:
    // index = x + y*X + z*X*Y
    // where X and Y are sizes of the lattice in x and y direction

    // from this stems that
    // x + y*X = index MOD X*Y
    // x = x + y*X MOD X
    // y = x + y*X DIV X
    // z = index DIV X*Y

    // see online documentation for class div_t (stdlib.h)
    div_t firstDiv, secondDiv;

    firstDiv = div(int(index), int(lbSize[0] * lbSize[1]));
    secondDiv = div(firstDiv.rem, int(lbSize[0]));

    x = secondDiv.rem;
    y = secondDiv.quot;
    z = firstDiv.quot;

    return tVect(double(x) - 0.5, double(y) - 0.5, double(z) - 0.5);
}

unsigned int LB::getX(const unsigned int& index) const {

    // it=x+y*X+z*X*Y
    static const unsigned int XY = double(lbSize[0] * lbSize[1]);
    static const unsigned int X = double(lbSize[0]);
    const double firstQuotient = floor(double(index) / XY);
    const double firstRemainder = index - firstQuotient*XY;

    const double secondQuotient = floor(double(firstRemainder) / X);
    const double secondRemainder = firstRemainder - secondQuotient*X;

    return int(secondRemainder);

    // ALTERNATIVE
    // see function getPosition for documentation
    //    div_t firstDiv, secondDiv;
    //
    //    firstDiv = div(int(index), int(lbSize[0] * lbSize[1]));
    //    //cout<<"firstDiv "<<firstDiv<<" int(lbSize[0] * lbSize[1]) "<<int(lbSize[0] * lbSize[1])<<endl;
    //    secondDiv = div(firstDiv.rem, int(lbSize[0]));
    //    //cout<<"secondDiv "<<secondDiv<<" int(lbSize[0]) "<<int(lbSize[0])<<endl;
    //    //cout<<"secondDiv.rem" <<secondDiv.rem<<endl;
    //    return secondDiv.rem;
}

unsigned int LB::getY(const unsigned int& index) const {

    // it=x+y*X+z*X*Y
    static const unsigned int XY = double(lbSize[0] * lbSize[1]);
    static const unsigned int X = double(lbSize[0]);
    const double firstQuotient = floor(double(index) / XY);
    const double firstRemainder = index - firstQuotient*XY;

    const double secondQuotient = floor(double(firstRemainder) / X);

    return int(secondQuotient);

    // ALTERNATIVE
    //    // see function getPosition for documentation
    //    div_t firstDiv, secondDiv;
    //
    //    firstDiv = div(int(index), int(lbSize[0] * lbSize[1]));
    //    secondDiv = div(firstDiv.rem, int(lbSize[0]));
    //
    //    return secondDiv.quot;
}

unsigned int LB::getZ(const unsigned int& index) const {

    // it=x+y*X+z*X*Y
    static const unsigned int XY = double(lbSize[0] * lbSize[1]);

    const double firstQuotient = floor(double(index) / XY);


    return int(firstQuotient);

    // ALTERNATIVE
    //    // see function getPosition for documentation
    //    div_t firstDiv;
    //
    //    firstDiv = div(int(index), int(lbSize[0] * lbSize[1]));
    //
    //    return firstDiv.quot;
}

double LB::getPositionX(const unsigned int& index) const {

    return double(getX(index)) - 0.5;

}

double LB::getPositionY(const unsigned int& index) const {

    return double(getY(index)) - 0.5;

}

double LB::getPositionZ(const unsigned int& index) const {

    return double(getZ(index)) - 0.5;

}




#include "node.h"

// basic functions

//tVect node::position(){
//    return tVect(double(x), double(y), double(z));
//}

// liquid fraction

double node::liquidFraction() const {
    return mass / n;
}

// functions working on distribution functions

void node::showPopulations() const {
    cout << "(";
    for (int j = 0; j < lbmDirec; ++j) {
        cout << j << ":" << f[j]-n*coeff[j] << ", ";
    }
    cout << ")";
}

void node::initialize(const double& density, const tVect& velocity, const double& massFunction, const double& viscosity, const tVect& F, const double& ageHere, const tVect& rotationSpeed) {

    //macroscopic properties of the initial condition
    n = density;
    mass = massFunction;
    u = velocity; //FORCE PART  lbF / 2.0 / n+
    visc = viscosity;
    friction = 0.0;
    age = static_cast<float>(ageHere);

    
    const tVect coriolisAcceleration = computeCoriolis(u, rotationSpeed);
    
    const tVect force = F + centrifugalForce + coriolisAcceleration;

    setEquilibrium(Zero, density, velocity);
    shiftVelocity(Zero);
    reconstruct();

    // add force term to new distributions
    addForce(force);

}

void node::copy(const node*& copyNode) {

    //macroscopic properties of the initial condition
    n = copyNode->n;
    mass = 0.1*copyNode->mass;
    u =   copyNode->u;
    visc = copyNode->visc;
    friction = copyNode->friction;
    centrifugalForce=copyNode->centrifugalForce;
    
    for (int j = 0; j < lbmDirec; ++j) {
        f[j] = copyNode->f[j];
        fs[j] = copyNode->fs[j];
    }
}

void node::scatterMass(double& _extraMass) {

    nodeList scatterNodes;
    scatterNodes.clear();
    for (int j = 1; j < lbmDirec; ++j) {
        node* linkNode = d[j];
        if (linkNode != nullptr) {
            if (linkNode->type == INTERFACE) {
                scatterNodes.push_back(linkNode);
            }
        }
    }
    const unsigned int totNodes = static_cast<unsigned int>(scatterNodes.size());
    if (totNodes > 0) {
        for (unsigned int i = 0; i < totNodes; ++i) {
            scatterNodes[i]->extraMass+=_extraMass/double(totNodes);
        }
        _extraMass=0.0;
    }
    // otherwise mass is not redistributed, but rather given to the global redistribution function
}

void node::restart(const double& restDensity, const tVect& restVelocity, const double& restMassFunction, const double& restViscosity, const double restF[]) {

    //macroscopic properties of the restart condition
    n = restDensity;
    mass = restMassFunction;
    u = restVelocity;
    visc = restViscosity;
    //    force=tVect(0.0,0.0,0.0);

    for (int j = 0; j < lbmDirec; ++j) {
        f[j] = fs[j] = restF[j];
    }
}

void node::setEquilibrium(const tVect& /*F*/, const double& nHere, const tVect& velHere) {
    static double C1 = 3.0;
    static double C2 = 4.5;
    static double C3 = 1.5;
    // static double F1 = 3.0;
    // static double F2 = 9.0;

    // const double tauF = 0.5 + visc*3.0; // = tau-0.5
    // cout<<visc<<endl;
    // const double omegaF = 1.0 - 0.5 / tauF;

    const double usq = velHere.norm2();

    for (int j = 0; j < lbmDirec; ++j) {
        // the following lines initialize f to be the local equilibrium values
        const double vu = velHere.dot(v[j]);
        // const tVect vmu = v[j]-velHere; // FORCE PART
        // const tVect forcePart = F1 * vmu + F2 * vu * v[j]; //FORCE PART
        //                dummyNode.f[j]=dummyNode.fs[j]=coeff[j]*dummyNode.n*(1.0+3.0*vu+4.5*vu*vu-1.5*usq)-2.0*(tau-0.5)*coeff[j]*forcePart.dot(F); //FORCE PART ((tau-0.5)=omegaf/omega)
        //                f[j]=fs[j]=coeff[j]*n*(1.0+C1*vu+C2*vu*vu-C3*usq)-coeff[j]*forcePart.dot(F); //FORCE PART ((tau-0.5)=omegaf/omega)
        f[j] = fs[j] = coeff[j] * nHere * (1.0 + C1 * vu + C2 * vu * vu - C3 * usq);// - omegaF * coeff[j] * forcePart.dot(F); //FORCE PART ((tau-0.5)=omegaf/omega)
    }
}

void node::reconstruct() {
    // reconstruction of macroscopical physical variables

    n = 0.0;
    u.reset();
    for (int j = 0; j < lbmDirec; ++j) {
        // density
        n += f[j];
        // momentum
        u += f[j] * v[j];
    }

    // velocity
    u /= n;
}


void node::shiftVelocity(const tVect& F) {
    const tVect totalForce = F + hydroForce;
    u += 0.5 * mass*totalForce;
}

void node::computeEquilibriumTRT(double feqp[], double feqm[]) {
    constexpr double C1 = 3.0;
    constexpr double C2 = 4.5;
    constexpr double C3 = 1.5;

    const double usq = u.norm2();
    double posSum=0.0;
    for (int j = 1; j < lbmDirec; ++j) {
        const double vu = u.dot(v[j]);
        feqp[j] = coeff[j] * n * (1.0 +  C2 * vu * vu - C3 * usq);
        feqm[j] = coeff[j] * n * (C1 * vu);
        posSum+=feqp[j];
    }
    feqp[0]=n-posSum;
    feqm[0]=0.0;
}

void node::computeEquilibrium(double feq[]) {
    constexpr double C1 = 3.0;
    constexpr double C2 = 4.5;
    constexpr double C3 = 1.5;

    const double usq = u.norm2();

    for (int j = 0; j < lbmDirec; ++j) {
        const double vu = u.dot(v[j]);
        feq[j] = coeff[j] * n * (1.0 + C1 * vu + C2 * vu * vu - C3 * usq);
    }
}

void node::computeApparentViscosity(const double feq[], const FluidMaterial& fluidMaterial) {

    // minimum and maximum viscosity
    const double tau = 0.5 + 3.0 * visc;

    tMat gamma(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    for (int j = 1; j < lbmDirec; ++j) {
        gamma += (f[j] - feq[j]) * vv[j];
    }
    gamma *= -1.5 / (tau * n);

    // shear rate (second invariant)
    const double shearRate = 2.0 * gamma.magnitude();

    // Bingham model
    double nuApp = 0.0;
    switch (fluidMaterial.rheologyModel) {
        case NEWTONIAN:
        {
            nuApp = fluidMaterial.initDynVisc;
            break;
        }
        case BINGHAM:
        {
            nuApp = fluidMaterial.plasticVisc + fluidMaterial.yieldStress / shearRate;
            break;
        }
        case FRICTIONAL:
        {
            // p=c_s^2 * n, scaled by atmospheric pressure (n_atm=1.0)
            const double pressure = std::max(fluidMaterial.minimumPressure, 0.33333333 * (n - 1.0));
            if (this->basal) {
                friction = fluidMaterial.basalFrictionCoefFluid;
            } else {
                 friction = fluidMaterial.frictionCoefFluid;
            }
            nuApp = fluidMaterial.initDynVisc + friction * pressure / shearRate;
            break;
        }
        case VOELLMY:
        {
            // p=c_s^2 * n, scaled by atmospheric pressure (n_atm=1.0)
            const double pressure = std::max(fluidMaterial.minimumPressure, 0.33333333 * (n - 1.0));
            if (this->basal) {
                friction = fluidMaterial.basalFrictionCoefFluid;
            } else {
                friction = fluidMaterial.frictionCoefFluid;
            }
            nuApp = friction * pressure / shearRate + fluidMaterial.rhod2 * shearRate;
            //nuApp = ( fluidMaterial.frictionCoefFluid * pressure + u.norm2()/fluidMaterial.voellmyCoefficient) / shearRate ;
            break;
        }
        case BAGNOLD:
        {
            // p=c_s^2 * n, scaled by atmospheric pressure (n_atm=1.0)
            nuApp = fluidMaterial.rhod2 * shearRate ;
            break;
        }
        case MUI:
        {

//            double smoothedPressureHere = 0.0;
//            double neigh(0.0);
//            for (int j = 0; j < lbmDirec; ++j) {
//                if (d[j] != nullptr) {
//                    if (d[j]->isFluid()) {
//                        smoothedPressureHere += coeff[j] * 0.33333333 * (d[j]->n - 1.0);
//                        neigh += coeff[j];
//                    }
//                }
//            }
//            smoothedPressureHere/=neigh;
//            smoothedPressure=smoothedPressureHere;//(9.0*smoothedPressure+smoothedPressureHere)/10.0;
            // p=c_s^2 * n, scaled by atmospheric pressure (n_atm=1.0)
            const double pressure = std::max(fluidMaterial.minimumPressure, 0.33333333*(n-1.0) ); //smoothedPressure
            const double inertialNumber = fluidMaterial.particleDiameter * shearRate / sqrt(pressure/fluidMaterial.particleDensity);
            const double regularizationFactor=1.0;//-exp(-shearRate/0.00005);
            if (this->basal) {
                friction = fluidMaterial.basalFrictionCoefFluid*regularizationFactor + fluidMaterial.deltaFriction / (fluidMaterial.baseInertial/ inertialNumber  + 1.0);
            } else {
                 friction = fluidMaterial.frictionCoefFluid*regularizationFactor + fluidMaterial.deltaFriction / (fluidMaterial.baseInertial/ inertialNumber  + 1.0);
            }
            
            nuApp = friction * pressure / shearRate;
//            // p=c_s^2 * n, scaled by atmospheric pressure (n_atm=1.0)
            //cout<<fluidMaterial.particleDensity<<endl;
//            const double pressureHere=max(0.0,0.33333333 * (n - 1.0));//max(smoothedPressure,0.0);
//            const double inertialNumber = fluidMaterial.particleDiameter * shearRate / sqrt(pressureHere / fluidMaterial.particleDensity);
//            //const double pressure = std::max(fluidMaterial.minimumPressure, smoothedPressure);
            if (inertialNumber > 0.05) {
                //nuApp = fluidMaterial.lbMinVisc; //fluidMaterial.turbConst*shearRate;//
            } else {
                //friction = fluidMaterial.frictionCoefFluid + fluidMaterial.deltaFriction / (fluidMaterial.baseInertial / inertialNumber + 1.0);
                //nuApp = friction * pressure / shearRate;
            }
            //            }

            if (this->isInterface()) { //(pressure<1.5*fluidMaterial.minimumPressure) {
            nuApp=fluidMaterial.lbMinVisc;
          }
            //cout<<nuApp<<endl;
            break;
        }

    }

    // Smagorinsky turbulence model
    if (fluidMaterial.turbulenceOn) {
        const double nuTurb = fluidMaterial.turbConst*shearRate;
        nuApp += nuTurb;
    }

    // limiting for stability
    visc = std::max(fluidMaterial.lbMinVisc, std::min(fluidMaterial.lbMaxVisc, nuApp));
    
//    if (this->isInterface()) {
//                cout<< visc<<" "<<fluidMaterial.lbMaxVisc<<endl;
//            }

}

void node::resetViscosity(const double& initViscosity) {

    //    else {
    // limiting for stability
    visc = initViscosity;
    //        if (flag) {
    //            visc=plasticVisc;
    //        }
    //    }
}

void node::solveCollisionTRT(const double feqp[], const double feqm[], const double& magicNumber) {
    // relaxation frequency
    const double omegap = 1.0 / (0.5 + 3.0 * visc);
    const double omegam = 6.0 * visc / (2.0*magicNumber + 3.0 * visc);
    
    // symmetric probability density functions (one per lattice direction)
    double fp[lbmDirec];
    // antisymmetric probability density functions (one per lattice direction)
    double fm[lbmDirec];
    
    // symmetric part
    fp[0]=f[0];
    for (int j = 1; j < lbmDirec; ++j) {
        fp[j] = 0.5*(f[j]+f[opp[j]]);
    }
    
    // antisymmetric part
    fm[0]=0;
    for (int j = 1; j < lbmDirec; ++j) {
        fm[j] = 0.5*(f[j]-f[opp[j]]);
    }
    
    for (int j = 0; j < lbmDirec; ++j) {
        f[j] += omegam * (feqm[j] - fm[j]) + omegap * (feqp[j] - fp[j]);
    }
}

void node::solveCollision(const double feq[]) {
    // relaxation frequency
    const double omega = 1.0 / (0.5 + 3.0 * visc);
    for (int j = 0; j < lbmDirec; j++) {
        f[j] += omega * (feq[j] - f[j]);
    }
}

void node::addForce(const tVect& F) {
    static double F1 = 3.0;
    static double F2 = 9.0;

    //tVect vmu, forcePart;

    const double omegaf = 1.0 - 1.0 / (1.0 + 6.0 * visc);
    const tVect totalForce = F + hydroForce;

    for (int j = 0; j < lbmDirec; j++) {
        const double vu = u.dot(v[j]);
        const tVect vmu = v[j] - u;
        //        vmu-=u;
        const tVect forcePart = F2 * vu * v[j] + F1*vmu;
        //        forcePart+=F1*vmu;
        f[j] +=  omegaf * coeff[j] * forcePart.dot(mass*totalForce);
    }
}

void node::addForce(double feq[], const tVect& F) {
    static double F1 = 3.0;

    //tVect vmu, forcePart;

    const double omegaf = 1.0 - 1.0 / (1.0 + 6.0 * visc);
    const tVect totalForce = F + hydroForce;

    for (int j = 0; j < lbmDirec; j++) {
        const tVect vmu = v[j] - u;
        //        vmu-=u;
        const tVect forcePart = feq[j]*F1*vmu/n;
        //        forcePart+=F1*vmu;
        f[j] +=  omegaf * forcePart.dot(mass*totalForce);
    }
}

void node::addForceTRT(const tVect& F) {
    static double F1 = 3.0;
    static double F2 = 9.0;

    //tVect vmu, forcePart;

    static const double Lambda=0.25;//9*visc*visc;
    const double omegam = 6.0 * visc / (2.0*Lambda + 3.0 * visc);
    const double omegaf = 1.0 - 0.5 * omegam;
    
    const tVect totalForce = F + hydroForce;

    for (int j = 0; j < lbmDirec; j++) {
        const double vu = u.dot(v[j]);
        const tVect vmu = v[j] - u;
        //        vmu-=u;
        const tVect forcePart = F2 * vu * v[j] + F1*vmu;
        //        forcePart+=F1*vmu;
        f[j] +=  omegaf * coeff[j] * forcePart.dot(totalForce);
    }
}

/*
void node::collideFNN(const double& turbConst, const tVect& F, const double& plasticVisc, const double& yieldStress) {

    // equilibrium distributions
    double feq[lbmDirec];

    // shift velocity field to F/2
    shiftVelocity(F);
    // compute equilibrium distributions
    computeEquilibrium(feq);
    // compute shear rate tensor, find invariant calculate viscosity (Bingham)
    computeShearRate(feq, turbConst, plasticVisc, yieldStress);
    // compute new distributions
    solveCollision(feq);
    // add force term to new distributions
    addForce(F);
}

void node::collideNN(const double& turbConst, const double& plasticVisc, const double& yieldStress) {

    // equilibrium distributions
    double feq[lbmDirec];

//    // shift velocity field to F/2
//    shiftVelocity(F);
    // compute equilibrium distributions
    computeEquilibrium(feq);
    // compute shear rate tensor, find invariant calculate viscosity (Bingham)
    computeShearRate(feq, turbConst, plasticVisc, yieldStress);
    // compute new distributions
    solveCollision(feq);
//    // add force term to new distributions
//    addForce(F);

}

void node::collideF(const tVect& F, const double& initVisc) {

    // equilibrium distributions
    double feq[lbmDirec];

    // shift velocity field to F/2
    shiftVelocity(F);
    // compute equilibrium distributions
    computeEquilibrium(feq);
    // compute shear rate tensor, find invariant calculate viscosity (Bingham)
//    computeShearRate(feq, plasticVisc, yieldStress, minVisc, maxVisc);
    visc=initVisc; // for NEWTONIAN - viscosity part
    // compute new distributions
    solveCollision(feq);
    // add force term to new distributions
    addForce(F);

}

void node::collide(const double& initVisc) {

    // equilibrium distributions
    double feq[lbmDirec];

    // shift velocity field to F/2
//    shiftVelocity(F);
    // compute equilibrium distributions
    computeEquilibrium(feq);
    // compute shear rate tensor, find invariant calculate viscosity (Bingham)
//    computeShearRate(feq, plasticVisc, yieldStress, minVisc, maxVisc);

//    double tau=0.5+3.0*visc*dt;
//    tMat gamma(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
//    for (int j=0; j<direc; ++j) {
//        gamma+=(f[j]-feq[j])*vv[j];
//    }
//    gamma*=1.5*dt/(tau*n);
//    // shear rate (second invariant)
//    shearRate=gamma.magnitude();

    visc=initVisc; // for NEWTONIAN - viscosity part
    // compute new distributions
    solveCollision(feq);
    // add force term to new distributions
//    addForce(F);

}
 */

void node::store() {
    for (int j = 0; j < lbmDirec; ++j) {
        fs[j] = f[j];
        //if (j>0) f[j]=0.0;
    }
}

//void node::stream(unsigned int& fdir, node& linknode, unsigned int& linkdir, double& prod, double& sum){
//    // this is a general expression for the streaming, where
//    // a population can change direction and size (Bounce-Back & Free Surface)
//    f[fdir]=prod*linknode.fs[linkdir]+sum;
//}

tVect node::bounceBackForce(const unsigned int& j, const double staticPres[], const double& BBi, const double& earthPressureCoeff, const double& maxVisc) const {
    //    tVect force=(2.0*fs[j]-2.0*coeff[j])*dt*v[j];
    //    cout<<"AAAA "<<force.norm()<<" "<<n<<" "<<n*force.norm()<<"\n";
    //    return (2.0*fs[j]-BBi)*dt*v[j];//*(n-1.0)/n;
    if (visc==maxVisc) {
        return (2.0 * (fs[j] + coeff[j]*(n-1.0)*(earthPressureCoeff-1.0)- staticPres[j]) - BBi) * v[j];
    } else {
        return (2.0 * (fs[j] - staticPres[j]) - BBi) * v[j];
    }
    //    return (2.0*fs[j]-BBi)*mass*mass/n*dt*v[j];
    //    return force;
    //    return -staticPres[j]*v[j];
}


double node::massStream(const unsigned int& sdir) const {
    return (f[opp[sdir]] - fs[sdir]);
}

// functions for identification of type

bool node::isActive() const {
    if (type == LIQUID) {
        return true;
    } else if (type == INTERFACE) {
        return true;
    } else
        return false;
}

bool node::isWall() const {
    switch (type) {
        case SLIP_STAT_WALL:
        case SLIP_DYN_WALL:
        case STAT_WALL:
        case DYN_WALL:
        case CYL:
        case OBJ:
        case TOPO:
            return true;
            break;
        default:
            return false;
            break;
    }
}

bool node::isCurvedWall() const {
    switch (type) {
        case CYL:
        case TOPO:
            return true;
            break;
        default:
            return false;
            break;
    }
}

// functions for change of type

void node::setSlipStatWall() {
    type = SLIP_STAT_WALL;
}

void node::setSlipDynWall() {
    type = SLIP_DYN_WALL;
}

void node::setStatWall() {
    type = STAT_WALL;
}

void node::setDynWall() {
    type = DYN_WALL;
}

void node::setCylinderWall() {
    type = CYL;
}

void node::setObjectWall() {
    type = OBJ;
}

void node::setTopography() {
    type = TOPO;
}

void node::setOutlet() {
    type = OUTLET;
}

void node::setInterface() {
    type = INTERFACE;
}

void node::setFluid() {
    type = LIQUID;
}

void node::setGas() {
    type = GAS;
}

void node::setPeriodic() {
    type = PERIODIC;
}

void node::setType(const types& typ) {
    // setting type according to identification number
    type=typ;
}

// particle flag

void node::setInsideParticle() {
    p = true;
}

void node::setOutsideParticle() {
    p = false;
}

// get type

//unsigned int node::getType() const {
//    switch (type) {
//        case LIQUID: return 0;
//        case GAS: return 2;
//        case INTERFACE: return 3;
//        case PERIODIC: return 4;
//        case SLIP_STAT_WALL: return 5;
//        case SLIP_DYN_WALL: return 6;
//        case STAT_WALL: return 7;
//        case DYN_WALL: return 8;
//        case CYL: return 9;
//        case OBJ: return 10;
//        case TOPO: return 11;
//        case OUTLET: return 12;
//        default:
//        {
//            cout << "Error, invalid node type" << endl;
//            exit(0);
//            break;
//        }
//    }
//    return type;
//}

// solid index functions

unsigned short int node::getSolidIndex() const {
    return solidIndex;
    
}

void node::setSolidIndex(const short unsigned int& ind) {
    solidIndex = ind;
}

types boundary2Type(const string& typeString) {
    // converts a string to a type for the boundary 
    if (typeString == "wall" || typeString == "stat_wall" || typeString == "STAT_WALL" || typeString == "7") {
        return STAT_WALL;
    }
    else if (typeString == "dyn_wall" || typeString == "DYN_WALL" || typeString == "8") {
        return DYN_WALL;
    }
    else if (typeString == "slip_stat_wall" || typeString == "SLIP_STAT_WALL" || typeString == "5") {
        return SLIP_STAT_WALL;
    }
    else if (typeString == "slip_dyn_wall" || typeString == "SLIP_DYN_WALL" || typeString == "6") {
        return SLIP_DYN_WALL;
    }
    else if (typeString == "periodic" || typeString == "PERIODIC" || typeString == "4") {
        return PERIODIC;
    }
    else if (typeString == "outlet" || typeString == "OUTLET" || typeString == "12") {
        return OUTLET;
    } else {
        cout << "Error, unvalid boundary conditon" << endl;
        exit(0);
    }
}

//types int2Type(const unsigned int& typeID) {
//    // converts a string to a type for the boundary 
//    switch (typeID) {
//        case 0: return LIQUID;
//        case 2: return GAS;
//        case 3: return INTERFACE;
//        case 4: return PERIODIC;
//        case 5: return SLIP_STAT_WALL;
//        case 6: return SLIP_DYN_WALL;
//        case 7: return STAT_WALL;
//        case 8: return DYN_WALL;
//        case 9: return CYL;
//        case 10: return OBJ;
//        case 11: return TOPO;
//        case 12: return OUTLET;
//        default:
//        {
//            cout << "Error, invalid node type" << endl;
//            exit(0);
//            break;
//        }
//    }
//}


// CURVED ////////////////////////////////

void curve::computeCoefficients() {
    for (int j = 1; j < lbmDirec; ++j) {
        m1[j] = (delta[j] - 1.0) / delta[j];
        m2[j] = 1.0 / delta[j];
        m3[j] = (delta[j] - 1.0) / (1.0 + delta[j]);
    }
}

double curve::getChi(const unsigned int& j, const double& tau) const {
    if (delta[j] >= 0.5) {
        return (2.0 * delta[j] - 1.0) / tau;
    } else {
        return (2.0 * delta[j] - 1.0) / (tau - 2.0);
    }
}

// UNITS ///////////////////////////////////

void measureUnits::setComposite() {
    Volume = Length * Length*Length;
    Speed = Length / Time; //+2
    Accel = Length / Time / Time; // +6
    AngVel = 1.0 / Time;
    KinVisc = Length * Length / Time; // 0
    DynVisc = Density * Length * Length / Time; // 0
    Force = Density * Length * Length * Length * Length / Time / Time; // +3
    Torque = Density * Length * Length * Length * Length * Length / Time / Time; // +1
    Mass = Density * Length * Length * Length; // -3
    Stress = Density * Length * Length / Time / Time;
    Pressure = Density * Length * Length / Time / Time;
    FlowRate = Density * Length * Length * Length / Time;
    Energy = Density * Length * Length * Length * Length * Length / Time / Time;
    invLength=1.0/Length;
}
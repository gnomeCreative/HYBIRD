#include "LBOpenMP.h"

#include "DEM.h"
#include "LB2.h"
#include "Problem.h"
#include "LBCommon.h"

void LBOpenMP::init(Problem& problem, cylinderList& cylinders, wallList& walls, particleList& particles, objectList& objects, bool externalSolveCoriolis, bool externalSolveCentrifugal) {
    // Convert from AoS format to SoA
    syncCylindersIn(cylinders);
    syncWallsIn(walls);
    syncParticlesIn(particles);
    syncObjectsIn(objects);

    // Lattice Boltzmann initialization steps

    // switchers for apparent accelerations
    h_PARAMS.solveCoriolis = externalSolveCoriolis;
    h_PARAMS.solveCentrifugal = externalSolveCentrifugal;

    // first comes the initialization of the data structures
    cout << "Initializing nodes containers and types" << endl;
    // total number of nodes
    h_PARAMS.totPossibleNodes = h_PARAMS.lbSize[0] * h_PARAMS.lbSize[1] * h_PARAMS.lbSize[2];

    // Allocate nodes (dense matrix, even gas nodes are represented)
    // This initialises them all as GAS
    allocateNodes(h_PARAMS.totPossibleNodes);

    // application of lattice boundaries
    initializeLatticeBoundaries();
    // then the initial node type must be identified for every node (if not specified, it is already Fluid)
    initializeTypes(walls, cylinders, objects);

    ifstream fluidFileID;
    if (h_PARAMS.lbRestart) {
        // open fluid restart file
        fluidFileID.open(init_params.lbRestartFile.c_str(), ios::in);
        ASSERT(fluidFileID.is_open());
        // check if the restart file size is ok
        unsigned int restartX, restartY, restartZ, restartNodes;
        fluidFileID >> restartX;
        fluidFileID >> restartY;
        fluidFileID >> restartZ;
        fluidFileID >> restartNodes;
        ASSERT(restartX == h_PARAMS.lbSize[0]);
        ASSERT(restartY == h_PARAMS.lbSize[1]);
        ASSERT(restartZ == h_PARAMS.lbSize[2]);
        // read active nodes from file and generate
        // @todo
        fprintf(stderr, "lbRestart is not yet supported.\n");
        throw std::exception();
        // restartInterface(fluidFileID, restartNodes);
    } else {
        // initialize interface
        initializeInterface(problem);
        // initialize variables for active nodes
        initializeVariables();
    }

    // initialize variables for wall nodes
    initializeWalls();

    // initialize h_nodes's fluid, interface and active lists
    initializeLists();

    // Setup hd_nodes, copy it to d_nodes
    /// initDeviceNodes(); // @todo v3 in subclass

    // application of particle initial position
    const double inside_mass = initializeParticleBoundaries(); // @todo v3 in subclass virtual?

    // in case mass needs to be kept constant, compute it here
    if (PARAMS.imposeFluidVolume) {
        // volume and mass is the same in lattice units
        PARAMS.totalMass = PARAMS.imposedFluidVolume / PARAMS.unit.Volume;
    } else {
        if (problemName == DRUM) {
            PARAMS.totalMass = PARAMS.fluidMass / PARAMS.unit.Mass;
        } else if (problemName == STAVA) {
            PARAMS.totalMass = 200000.0 / PARAMS.unit.Volume;
        } else {
            PARAMS.totalMass = inside_mass;
        }
    }
    if (PARAMS.increaseVolume) {
        PARAMS.deltaVolume /= PARAMS.unit.Volume;
        PARAMS.deltaTime /= PARAMS.unit.Time;
    }

    // syncParams(); // @todo v3 subclass

    printf("%u/%u nodes\n", h_nodes.activeCount, h_nodes.count);
    cout << "Done with initialization" << endl;
}
void LBOpenMP::step(DEM& dem, bool io_demSolver) {
    //@todo v3 All methods here virtual for CUDA sub class?
    if (io_demSolver) {
        this->syncDEMIn(dem.elmts, dem.particles, dem.walls, dem.objects);
        this->latticeBoltzmannCouplingStep(dem.newNeighborList);
    }

    if (dem.demTime >= dem.demInitialRepeat && h_nodes.activeCount) {
        this->latticeBoltzmannStep();

        if (io_demSolver) {
            // Shift element/wall/object forces and torques to physical units
            this->shiftToPhysical();
        }

        // Lattice Boltzmann core steps
        if (PARAMS.freeSurface) {
            this->latticeBoltzmannFreeSurfaceStep();
        }
    }
    if (io_demSolver) {
        this->syncDEMOut(dem.elmts, dem.particles, dem.walls, dem.objects);
    }
}
void LBOpenMP::latticeBoltzmannCouplingStep(bool& newNeighbourList) {
    // identifies which nodes need to have an update due to particle movement
    // the complexity arises from trying to keep the scaling as close to linear as possible
    // maybe the best idea is to do this in two step:
    // 1) the first is to check for new active nodes and initialise them
    // 2) the second is checking for new solid nodes.
    // this automatically check also for the hideous case of particle to particle double transition

    /**
     * @todo The parallelisation of each of these methods should be reviewed
     *       Most are 2D loops, the range of each being (1k to 100k)
     */

     // first we check if a new neighbour table has been defined. In that case, the indexing needs to be reinitialised
    if (newNeighbourList) {
        cout << endl << "New neighbour list" << endl;
        this->initializeParticleBoundaries();
        newNeighbourList = false;
    } else {
        // SOLID TO ACTIVE CHECK
        // @note Calling this directly after initializeParticleBoundaries() is redundant, hence else
        this->findNewActive();
    }

    // ACTIVE TO SOLID CHECK
    this->findNewSolid();

    if (PARAMS.freeSurface) {
        this->checkNewInterfaceParticles();
    }
}
void LBOpenMP::latticeBoltzmannStep() {
    // Reconstruct active list
    h_nodes.cleanLists<IMPL>();

    // Initializing the elements forces (lattice units)
    h_elements.initElements<IMPL>();

    // Initialise lattice boltzmann force vector
    if (!h_PARAMS.forceField) {
        h_PARAMS.lbF.reset();
        // syncParams(); // @todo v3 sync for device here
    }

    // reconstruct(), computeHydroForces(), collision()
    // Reconstruct macroscopic variables from microscopic distribution
    // Compute interaction forces with DEM elmts
    // Collision step
    this->reconstructHydroCollide();

    // Streaming operator
    this->streaming();
}
extern ProblemName problemName;
void LBOpenMP::latticeBoltzmannFreeSurfaceStep() {
    // in case mass needs to be kept constant, call enforcing function here
    if (PARAMS.imposeFluidVolume) {
        this->enforceMassConservation();
    } else if (PARAMS.increaseVolume) {
        if (PARAMS.time < PARAMS.deltaTime) {
            this->redistributeMass(PARAMS.deltaVolume / PARAMS.deltaTime);
        }
    } else {
        if (problemName == DRUM ||
            problemName == STAVA) {
            this->enforceMassConservation();
        }
    }

    // mass and free surface update
    this->updateMass();
    this->updateInterface();
    h_nodes.cleanLists<CPU>();
}
//
// Init chain
//
void LBOpenMP::allocateNodes(const unsigned int count) {
    // Allocate enough memory for these nodes
    assert(h_nodes.count == 0);  // No nodes should exist at the time this is called
    h_nodes.count = count;
    // Allocate host buffers
    //h_nodes.coord = static_cast<unsigned int*>(malloc(h_nodes.count * sizeof(unsigned int))); // TODO nolonger required
    h_nodes.f = static_cast<double*>(malloc(h_nodes.count * lbmDirec * sizeof(double)));
    h_nodes.fs = static_cast<double*>(malloc(h_nodes.count * lbmDirec * sizeof(double)));
    h_nodes.n = static_cast<double*>(malloc(h_nodes.count * sizeof(double)));
    h_nodes.u = static_cast<tVect*>(malloc(h_nodes.count * sizeof(tVect)));
    h_nodes.hydroForce = static_cast<tVect*>(malloc(h_nodes.count * sizeof(tVect)));
    h_nodes.centrifugalForce = static_cast<tVect*>(malloc(h_nodes.count * sizeof(tVect)));
    h_nodes.mass = static_cast<double*>(malloc(h_nodes.count * sizeof(double)));
    h_nodes.newMass = static_cast<double*>(malloc(h_nodes.count * sizeof(double)));
    h_nodes.visc = static_cast<double*>(malloc(h_nodes.count * sizeof(double)));
    h_nodes.basal = static_cast<bool*>(malloc(h_nodes.count * sizeof(bool)));
    h_nodes.friction = static_cast<double*>(malloc(h_nodes.count * sizeof(double)));
    h_nodes.age = static_cast<float*>(malloc(h_nodes.count * sizeof(float)));
    h_nodes.solidIndex = static_cast<unsigned int*>(malloc(h_nodes.count * sizeof(unsigned int)));
    h_nodes.d = static_cast<unsigned int*>(malloc(h_nodes.count * lbmDirec * sizeof(unsigned int)));
    h_nodes.type = static_cast<types*>(malloc(h_nodes.count * sizeof(types)));
    h_nodes.p = static_cast<bool*>(malloc(h_nodes.count * sizeof(bool)));
    // Zero initialisation
    memset(h_nodes.f, 0, h_nodes.count * lbmDirec * sizeof(double));
    memset(h_nodes.fs, 0, h_nodes.count * lbmDirec * sizeof(double));
    memset(h_nodes.n, 0, h_nodes.count * sizeof(double));
    memset(h_nodes.u, 0, h_nodes.count * sizeof(tVect));
    memset(h_nodes.hydroForce, 0, h_nodes.count * sizeof(tVect));
    memset(h_nodes.centrifugalForce, 0, h_nodes.count * sizeof(tVect));
    memset(h_nodes.mass, 0, h_nodes.count * sizeof(double));
    memset(h_nodes.newMass, 0, h_nodes.count * sizeof(double));
    std::fill(h_nodes.visc, h_nodes.visc + h_nodes.count, 1.0);
    memset(h_nodes.basal, 0, h_nodes.count * sizeof(bool));
    memset(h_nodes.friction, 0, h_nodes.count * sizeof(double));
    memset(h_nodes.age, 0, h_nodes.count * sizeof(float));
    memset(h_nodes.solidIndex, UINT_MAX, h_nodes.count * sizeof(unsigned int));
    memset(h_nodes.d, UINT_MAX, h_nodes.count * lbmDirec * sizeof(unsigned int));
    std::fill(h_nodes.type, h_nodes.type + h_nodes.count, GAS);
    memset(h_nodes.p, 0, h_nodes.count * sizeof(bool));
}
void LBOpenMP::initializeLatticeBoundaries() {
    // assign boundary characteristic to nodes (see class)
    // if not differently defined, type is 0 (fluid)

    // BOUNDARY CONDITIONS ///////////////////////////
    // solid boundary wins over all in corners, where more than 1 bc is defined
    cout << "Initializing boundaries" << endl;

    unsigned int indexHere = 0;
    // XY
    for (unsigned int x = 0; x < h_PARAMS.lbSize[0]; ++x) {
        for (unsigned int y = 0; y < h_PARAMS.lbSize[1]; ++y) {
            // bottom
            indexHere = h_PARAMS.getIndex(x, y, 0);
            if (h_nodes.type[indexHere] == GAS) {
                generateNode(indexHere, h_PARAMS.boundary[4]);
            }
            // top
            indexHere = h_PARAMS.getIndex(x, y, h_PARAMS.lbSize[2] - 1);
            if (h_nodes.type[indexHere] == GAS) {
                generateNode(indexHere, h_PARAMS.boundary[5]);
            }
        }
    }

    // YZ
    for (unsigned int y = 0; y < h_PARAMS.lbSize[1]; ++y) {
        for (unsigned int z = 0; z < h_PARAMS.lbSize[2]; ++z) {
            // bottom
            indexHere = h_PARAMS.getIndex(0, y, z);
            if (h_nodes.type[indexHere] == GAS) {
                generateNode(indexHere, h_PARAMS.boundary[0]);
            }
            // top
            indexHere = h_PARAMS.getIndex(h_PARAMS.lbSize[0] - 1, y, z);
            if (h_nodes.type[indexHere] == GAS) {
                generateNode(indexHere, h_PARAMS.boundary[1]);
            }
        }
    }

    // ZX
    for (unsigned int z = 0; z < h_PARAMS.lbSize[2]; ++z) {
        for (unsigned int x = 0; x < h_PARAMS.lbSize[0]; ++x) {
            // bottom
            indexHere = h_PARAMS.getIndex(x, 0, z);
            if (h_nodes.type[indexHere] == GAS) {
                generateNode(indexHere, h_PARAMS.boundary[2]);
            }
            // top
            indexHere = h_PARAMS.getIndex(x, h_PARAMS.lbSize[1] - 1, z);
            if (h_nodes.type[indexHere] == GAS) {
                generateNode(indexHere, h_PARAMS.boundary[3]);
            }
        }
    }
}
void LBOpenMP::initializeTypes(const wallList& walls, const cylinderList& cylinders, const objectList& objects) {
    initializeWallBoundaries(walls);
    // application of solid cylinders
    initializeCylinderBoundaries(cylinders);
    // application of objects
    initializeObjectBoundaries(objects);
    // initializing topography if one is present
    initializeTopography();
}
void LBOpenMP::initializeWallBoundaries(const wallList& walls) {
    // const double wallThickness = 2.0 * h_PARAMS.unit.Length;
    // SOLID WALLS ////////////////////////
    for (unsigned int iw = 0; iw < walls.size(); ++iw) {
        const tVect convertedWallp = walls[iw].p / h_PARAMS.unit.Length;
        const tVect normHere = walls[iw].n;
        const unsigned int indexHere = walls[iw].index;
        const bool slipHere = walls[iw].slip;
        const bool movingHere = walls[iw].moving;
        // @todo This was previously OpenMP parallel, but could be race condition in generateNode?
        for (unsigned int it = 0; it < h_PARAMS.totPossibleNodes; ++it) {
            // check if the node is solid
            // all walls have max thickness 2 nodes
            const tVect pos = h_PARAMS.getPosition(it);
            const double wallDistance = pos.distance2Plane(convertedWallp, normHere);
            if (wallDistance > -2.0 && wallDistance < 0.0) {
                //check for borders in limted walls
                if (walls[iw].limited) {
                    const double xHere = pos.x * h_PARAMS.unit.Length;
                    const double yHere = pos.y * h_PARAMS.unit.Length;
                    const double zHere = pos.z * h_PARAMS.unit.Length;
                    // check if beyond limits
                    if (xHere < walls[iw].xMin || xHere > walls[iw].xMax ||
                        yHere < walls[iw].yMin || yHere > walls[iw].yMax ||
                        zHere < walls[iw].zMin || zHere > walls[iw].zMax) {
                        continue;
                    }
                }
                // Node is inside a wall
                // generate node (tentatively as static wall)
                generateNode(it, STAT_WALL);
                // setting solidIndex
                h_nodes.solidIndex[it] = indexHere; // TODO indexHere is redundant, use iw?
                // setting type: 5-6=slip, 7-8=no-slip
                if (slipHere) {
                    // setting type for slip: 5=static, 6=moving
                    if (movingHere) {
                        h_nodes.type[it] = SLIP_DYN_WALL;
                    } else {
                        h_nodes.type[it] = SLIP_STAT_WALL;
                    }
                } else {
                    // setting type for no-slip: 7=static, 8=moving
                    if (movingHere) {
                        h_nodes.type[it] = DYN_WALL;
                    } else {
                        h_nodes.type[it] = STAT_WALL;
                    }
                }
            }
        }
    }

}
void LBOpenMP::initializeObjectBoundaries(const objectList& objects) {
    // SOLID WALLS ////////////////////////
    for (int io = 0; io < objects.size(); ++io) {
        const tVect convertedPosition = objects[io].x0 / h_PARAMS.unit.Length;
        const double convertedRadius = objects[io].r / h_PARAMS.unit.Length;
        const unsigned int indexHere = objects[io].index;
        // @todo This was previously OpenMP parallel, but could be race condition in generateNode?
        for (unsigned int it = 0; it < h_PARAMS.totPossibleNodes; ++it) {
            const tVect nodePosition = h_PARAMS.getPosition(it);
            if (nodePosition.insideSphere(convertedPosition, convertedRadius)) {
                generateNode(it, OBJ);
                h_nodes.solidIndex[it] = indexHere; // TODO indexHere is redundant, use io?
            }
        }
    }
}
void LBOpenMP::initializeCylinderBoundaries(const cylinderList& cylinders) {
    // SOLID CYLINDERS ////////////////////////
    for (int ic = 0; ic < cylinders.size(); ++ic) {
        const tVect convertedCylinderp1 = cylinders[ic].p1 / h_PARAMS.unit.Length;
        const tVect naxesHere = cylinders[ic].naxes;
        const double convertedRadius = cylinders[ic].R / h_PARAMS.unit.Length;
        const unsigned int indexHere = cylinders[ic].index;
        const bool slipHere = cylinders[ic].slip;
        const bool movingHere = cylinders[ic].moving;
        // @todo This was previously OpenMP parallel, but could be race condition in generateNode?
        for (unsigned int it = 0; it < h_PARAMS.totPossibleNodes; ++it) {
            // creating solid cells
            const bool isOutside = h_PARAMS.getPosition(it).insideCylinder(convertedCylinderp1, naxesHere, convertedRadius, convertedRadius + 3.0);
            const bool isInside = h_PARAMS.getPosition(it).insideCylinder(convertedCylinderp1, naxesHere, max(convertedRadius - 3.0, 0.0), convertedRadius);
            if ((cylinders[ic].type == FULL && isInside) ||
                (cylinders[ic].type == EMPTY && isOutside)) {
                //check for borders in limted walls
                if (cylinders[ic].limited) {
                    const tVect here = h_PARAMS.getPosition(it) * h_PARAMS.unit.Length;
                    // check if beyond limits
                    if (here.x < cylinders[ic].xMin || here.x > cylinders[ic].xMax ||
                        here.y < cylinders[ic].yMin || here.y > cylinders[ic].yMax ||
                        here.z < cylinders[ic].zMin || here.z > cylinders[ic].zMax) {
                        continue;
                    }
                }
                // Node is inside a cylinder
                // tentatively static
                generateNode(it, STAT_WALL);
                // setting solidIndex
                h_nodes.solidIndex[it] = indexHere;  // TODO indexHere is redundant, use ic?
                // setting type: 5-6=slip, 7-8=no-slip
                if (slipHere) {
                    // setting type for slip: 5=static, 6=moving
                    if (movingHere) {
                        h_nodes.type[it] = SLIP_DYN_WALL;
                    } else {
                        h_nodes.type[it] = SLIP_STAT_WALL;
                    }
                } else {
                    // setting type for no-slip: 7=static, 8=moving
                    if (movingHere) {
                        h_nodes.type[it] = DYN_WALL;
                    } else {
                        h_nodes.type[it] = STAT_WALL;
                    }
                }
            }
        }
    }
}
void LBOpenMP::initializeTopography() {
    // Based on initializeTopography()
    
    const double surfaceThickness = 1.75 * h_PARAMS.unit.Length;

    // TOPOGRAPHY ////////////////////////
    if (h_PARAMS.lbTopography) {
        lbTop.readFromFile(init_params.lbTopographyFile, h_PARAMS.translateTopographyX, h_PARAMS.translateTopographyY, h_PARAMS.translateTopographyZ);
        lbTop.show();
        // check if topography grid contains the fluid domain
        ASSERT(lbTop.coordX[0] < h_PARAMS.unit.Length);
        ASSERT(lbTop.coordY[0] < h_PARAMS.unit.Length);

        cout << "lbTop.coordX[lbTop.sizeX - 1]=" << lbTop.coordX[lbTop.sizeX - 1] << endl;
        cout << "lbSize[0]) * unit.Length=" << h_PARAMS.lbSize[0] * h_PARAMS.unit.Length << endl;
        ASSERT(lbTop.coordX[lbTop.sizeX - 1] > h_PARAMS.lbSize[0] * h_PARAMS.unit.Length);
        cout << "lbTop.coordY[lbTop.sizeY - 1]=" << lbTop.coordY[lbTop.sizeY - 1] << endl;
        cout << "lbSize[1]) * unit.Length=" << h_PARAMS.lbSize[1] * h_PARAMS.unit.Length << endl;
        ASSERT(lbTop.coordY[lbTop.sizeY - 1] > h_PARAMS.lbSize[1] * h_PARAMS.unit.Length);

        // @todo This was previously OpenMP parallel, critical section around generateNode()
        for (unsigned int ix = 1; ix < h_PARAMS.lbSize[0] - 1; ++ix) {
            for (unsigned int iy = 1; iy < h_PARAMS.lbSize[1] - 1; ++iy) {
                for (unsigned int iz = 1; iz < h_PARAMS.lbSize[2] - 1; ++iz) {
                    const tVect nodePosition = tVect(ix, iy, iz) * h_PARAMS.unit.Length;
                    const double distanceFromTopography = lbTop.distance(nodePosition);

                    if (distanceFromTopography < 0.0 && distanceFromTopography>-1.0 * surfaceThickness) {
                        const unsigned int it = ix + iy * h_PARAMS.lbSize[0] + iz * h_PARAMS.lbSize[0] * h_PARAMS.lbSize[1];
                        generateNode(it, STAT_WALL);
                        h_nodes.type[it] = TOPO;
                    }
                }
            }
        }
    }
}
void LBOpenMP::initializeInterface(const Problem &problem) {
    // TODO old problem switch is no longer supported!
    // creates an interface electing interface cells from active cells
    if (h_PARAMS.lbTopographySurface) {
        // Formerly setTopographySurface()
        // @todo This was previously OpenMP parallel, critical section around generateNode()
        for (unsigned int it = 0; it < h_PARAMS.totPossibleNodes; ++it) {
            if (h_nodes.type[it] == GAS) {
                // control is done in real coordinates
                const tVect nodePosition = h_PARAMS.getPosition(it) * h_PARAMS.unit.Length;
                const double surfaceIsoparameterHere = lbTop.surfaceIsoparameter(nodePosition);
                if (surfaceIsoparameterHere > 0.0 && surfaceIsoparameterHere <= 1.0) {// setting solidIndex
                    generateNode(it, LIQUID);
                }
            }
        }
    } else if(!problem.file.empty()) { // Problem file was loaded
        cout << "Initializing from problem file:" << endl;
        for (const auto &f : problem.fluids_basic) {
            cout << "BOX=min(" << f.min.x << ", " << f.min.y << ", " << f.min.z << ")" << endl;
            cout << "    max(" << f.max.x << ", " << f.max.y << ", " << f.max.z << ")" << endl;
        }
        for (const auto& f : problem.fluids_complex_str) {
            cout << "EXPRESSION=if (" << f <<") > 0" << endl;
        }
        unsigned int ct = 0;
        for (unsigned int it = 0; it < h_PARAMS.totPossibleNodes; ++it) {
            if (h_nodes.type[it] == GAS) {
                // creating fluid cells
                if (problem.isFluid(h_PARAMS.getPosition(it) * h_PARAMS.unit.Length)) {
                    generateNode(it, LIQUID);
                    ++ct;
                }
            }
        }
        cout << ct << " of " << h_PARAMS.totPossibleNodes << " nodes were init as liquid." << endl;
    } else {
        switch (problemName) {
        case SHEARCELL:
        case AVALANCHE:
        case DRUM:
        case NET:
        case BARRIER:
        case ZHOU:
        case OPENBARRIER:
        case HONGKONG:
        case STVINCENT:
        case STAVA:
        case NIGRO:
        case CAROLINE:
        case DAMBREAK:
        case GRAY_DAMBREAK:
        case GRAY_DAMBREAK_2D:
        case INCLINEFLOW:
        case HOURGLASS:
        case IERVOLINO:
        case IERVOLINO_2D:
        case IERVOLINO_CYLINDERTEST:
        case HEAP:
        case TRIAXIAL:
        case JOP:
        case WILL:
        case WILL_SETTLING:
        case MANGENEY:
        case GRAY:
        case ESERCITAZIONE:
        case FILIPPO_SILOS:
        case HK_SMALL:
        case HK_LARGE:
        case KELVIN:
        case SHEARCELL2023:
        case INTRUDER:
        case OBJMOVING:
            cerr << "Error: problem '" << problemName << "' is not supported by LB2::initializeInterface()" << endl;
            cerr << "config file should be upgraded to use a problem file free surface definition" << endl;
            std::abort();
        case NONE:
        default:
            {
                cout << "Initializing from problem file:" << endl;
                cout << "X=(" << double(h_PARAMS.freeSurfaceBorders[0]) * h_PARAMS.unit.Length << ", " << double(h_PARAMS.freeSurfaceBorders[1]) * h_PARAMS.unit.Length << ")" << endl;
                cout << "Y=(" << double(h_PARAMS.freeSurfaceBorders[2]) * h_PARAMS.unit.Length << ", " << double(h_PARAMS.freeSurfaceBorders[3]) * h_PARAMS.unit.Length << ")" << endl;
                cout << "Z=(" << double(h_PARAMS.freeSurfaceBorders[4]) * h_PARAMS.unit.Length << ", " << double(h_PARAMS.freeSurfaceBorders[5]) * h_PARAMS.unit.Length << ")" << endl;
                for (unsigned int it = 0; it < h_PARAMS.totPossibleNodes; ++it) {
                    if (h_nodes.type[it] == GAS) {
                        // creating fluid cells
                        const tVect pos = h_PARAMS.getPosition(it);
                        if ((pos.x > h_PARAMS.freeSurfaceBorders[0]) &&
                            (pos.x < h_PARAMS.freeSurfaceBorders[1]) &&
                            (pos.y > h_PARAMS.freeSurfaceBorders[2]) &&
                            (pos.y < h_PARAMS.freeSurfaceBorders[3]) &&
                            (pos.z > h_PARAMS.freeSurfaceBorders[4]) &&
                            (pos.z < h_PARAMS.freeSurfaceBorders[5])) {
                            generateNode(it, LIQUID);
                        }
                    }
                }
            }
            break;
        }
    }
}
void LBOpenMP::generateNode(unsigned int coord, types typeHere) {
    // set type
    h_nodes.type[coord] = typeHere;
    h_nodes.p[coord] = false;  // setOutsideParticle()
    h_nodes.age[coord] = 0.0;

    // TODO Add it to list of known non-gas nodes?

    // find neighbor indices
    const std::array<unsigned int, lbmDirec> neighborCoord = h_nodes.findNeighbors(coord);

    h_nodes.basal[coord] = false;

    // set centrifugal acceleration
    h_nodes.centrifugalForce[coord] = computeCentrifugal(h_nodes.getPosition(coord), PARAMS.rotationCenter, PARAMS.rotationSpeed);

    // assign neighbor nodes
    for (unsigned int j = 1; j < lbmDirec; ++j) {
        // linearized coordinate of neighbor nodes
        const unsigned int link = neighborCoord[j];
        // check if node at that location exists
        if (link < h_nodes.count && h_nodes.type[link] != GAS) {
            // assign neighbor for local node
            h_nodes.d[j * h_nodes.count + coord] = link;
            // if neighbor node is also active, link it to local node
            if (h_nodes.isActive(coord)) {
                h_nodes.d[opp[j] * h_nodes.count + link] = coord;
                if (h_nodes.isWall(link)) {
                    h_nodes.basal[coord] = true;
                }
            }
        } else {
            h_nodes.d[j * h_nodes.count + coord] = std::numeric_limits<unsigned int>::max();
        }
    }
}
void LBOpenMP::initializeVariables() {
    cout << "Initializing variables" << endl;
    // note that interface is not defined here. All fluid, interface and gas cells are uninitialized at the moment
    // calculate maximum height of the fluid

    // find "taller" and "deepest" points
    double minProjection = std::numeric_limits<double>::max();
    double maxProjection = -std::numeric_limits<double>::max();
        
    if (!PARAMS.solveCentrifugal) {
        // TODO openmp reduction?
        for (unsigned int i = 0; i < h_nodes.count; ++i) {
            if (h_nodes.isActive(i)) {
                const tVect position = h_nodes.getPosition(i);
                const double projection = position.dot(PARAMS.lbF);
                minProjection = std::min(minProjection, projection);
                maxProjection = std::max(maxProjection, projection);
            }
        }
        cout << "minProjection = " << minProjection << endl;
    } else {
        // TODO openmp reduction?
        for (unsigned int i = 0; i < h_nodes.count; ++i) {
            if (h_nodes.isActive(i)) {
                const tVect position = h_nodes.getPosition(i);
                const double projection = position.dot(h_nodes.centrifugalForce[i]);
                minProjection = std::min(minProjection, projection);
                maxProjection = std::max(maxProjection, projection);
            }
        }
        cout << "minProjection = " << minProjection << endl;
    }

    // checking for boundary between gas and fluid and assigning interface properties
    // at this point fluid cells contain actual fluid cells and potential interface cells, so we create the node anyway
    double massFluid = 0.0;
    double massInterface = 0.0;
    // TODO openmp?
    for (unsigned int i = 0; i < h_nodes.count; ++i) {
        if (h_nodes.type[i] == LIQUID) {
            // check if it is interface
            for (int j = 1; j < lbmDirec; ++j) {
                unsigned int linkNode = h_nodes.d[j * h_nodes.count + i];
                if (linkNode == std::numeric_limits<unsigned int>::max()) {
                    h_nodes.type[i] = INTERFACE;
                    break;
                }
            }
        }
        // now assign macroscopic quantities accordingly
        // FLUID NODES ////
        if (h_nodes.type[i] == LIQUID) {
            massFluid += 1.0;
            // setting macroscopic variables
            // density is calculated using hydrostatic profile
            const tVect position = h_nodes.getPosition(i);
            if (!PARAMS.solveCentrifugal) {
                const double projection = position.dot(PARAMS.lbF);
                h_nodes.initialize(i, PARAMS.fluidMaterial.initDensity + 3.0 * PARAMS.fluidMaterial.initDensity * (projection-minProjection), PARAMS.initVelocity, PARAMS.fluidMaterial.initDensity, PARAMS.fluidMaterial.initDynVisc, PARAMS.lbF, 1.0, Zero);
            } else {
                const double projection = position.dot(h_nodes.centrifugalForce[i]);
                h_nodes.initialize(i, PARAMS.fluidMaterial.initDensity + 3.0 * PARAMS.fluidMaterial.initDensity * (projection-minProjection), PARAMS.initVelocity, PARAMS.fluidMaterial.initDensity, PARAMS.fluidMaterial.initDynVisc, PARAMS.lbF, 1.0, PARAMS.rotationSpeed);
            }
        }// INTERFACE NODES ////
        else if (h_nodes.type[i] == INTERFACE) {
            massInterface += 0.5;
            // setting macroscopic variables
            h_nodes.initialize(i, PARAMS.fluidMaterial.initDensity, PARAMS.initVelocity, 0.5 * PARAMS.fluidMaterial.initDensity, PARAMS.fluidMaterial.initDynVisc, PARAMS.lbF, 1.0, PARAMS.rotationSpeed);
        }

    }
    cout << "Approximate volume = " << massFluid * PARAMS.unit.Volume << " (fluid body), " << massInterface * PARAMS.unit.Volume << " (interface), " << (massFluid + massInterface) * PARAMS.unit.Volume << " (tot), " << endl;
}
void LBOpenMP::initializeWalls() {
    cout << "Initializing wall nodes" << endl;
    const double zero = 0.0;

    std::vector<unsigned int> wallNodes;

    // initializing wall nodes
    // note that, in the hypothesis that these walls are not evolving, only nodes at the interface need creation
    // TODO openmp?
    for (unsigned int i = 0; i < h_nodes.count; ++i) {
        if (h_nodes.isWall(i)) {
            // initialize node
            // STATIC WALL NODES ////
            if (h_nodes.type[i] == STAT_WALL ||
                h_nodes.type[i] == SLIP_STAT_WALL ||
                h_nodes.type[i] == OBJ ||
                h_nodes.type[i] == TOPO) {
                // reset velocity and mass (useful for plotting)
                // density=0.0; velocity=(0.0,0.0,0.0), mass=0.0; viscosity=0.0; force=(0.0,0.0,0.0)
                h_nodes.initialize(i, PARAMS.fluidMaterial.initDensity, Zero, zero, zero, Zero, 1.0, Zero);
            }// DYNAMIC WALL NODES ////
            else if (h_nodes.type[i] == DYN_WALL || 
                     h_nodes.type[i] == SLIP_DYN_WALL || 
                     h_nodes.type[i] == CYL) {
                // need to define velocity. It could be part of a cylinder or wall, we check both
                tVect solidVelocity;
                const tVect nodePosition = h_nodes.getPosition(i);
                unsigned int solidIndex = h_nodes.solidIndex[i];
                // wall
                if (solidIndex < h_walls.count && nodePosition.insidePlane(h_walls.p[solidIndex] / PARAMS.unit.Length, h_walls.n[solidIndex])) {
                    solidVelocity = h_walls.getSpeed(solidIndex, nodePosition * PARAMS.unit.Length) / PARAMS.unit.Speed;
                }// cylinder
                else if (solidIndex < h_cylinders.count && !nodePosition.insideCylinder(h_cylinders.p1[solidIndex] / PARAMS.unit.Length, h_cylinders.naxes[solidIndex], 0.0, h_cylinders.R[solidIndex] / PARAMS.unit.Length)) {
                    solidVelocity = h_cylinders.getSpeed(solidIndex, nodePosition * PARAMS.unit.Length) / PARAMS.unit.Speed;
                }// objects
                else if (solidIndex < h_objects.count && nodePosition.insideSphere(h_objects.x0[solidIndex] / PARAMS.unit.Length, h_objects.r[solidIndex] / PARAMS.unit.Length)) {
                    solidVelocity = h_objects.x1[solidIndex] / PARAMS.unit.Speed;
                }
                // reset velocity and mass (useful for plotting)
                // density=0.0; velocity=solidVelocity, mass=0.0; viscosity=0.0; force=(0.0,0.0,0.0)
                h_nodes.initialize(i, PARAMS.fluidMaterial.initDensity, solidVelocity, zero, zero, Zero, 1.0, PARAMS.rotationSpeed);
            }
            // add node to list
            wallNodes.push_back(i);
        }
    }
    // Allocate the wall nodes storage
    h_nodes.wallCount = static_cast<unsigned int>(wallNodes.size());
    h_nodes.wallI = static_cast<unsigned int*>(malloc(h_nodes.wallCount * sizeof(unsigned int)));
    memcpy(h_nodes.wallI, wallNodes.data(), h_nodes.wallCount * sizeof(unsigned int));    
}
void LBOpenMP::initializeLists() {
    cout << "Resetting lists ...";

    // note that interface is not defined here. All fluid, interface and gas cells are 0 at the moment
    std::vector<unsigned int> fluidNodes;
    std::vector<unsigned int> interfaceNodes;

    // creating list and initialize macroscopic variables for all nodes except walls
    // TODO OpenMP?
    for (unsigned int i = 0; i < h_nodes.count; ++i) {
        if (h_nodes.type[i] == LIQUID) {
            fluidNodes.push_back(i);
        } else if (h_nodes.type[i] == INTERFACE) {
            interfaceNodes.push_back(i);
        }
    }

    // Array to Buffer
    assert(!h_nodes.fluidI);
    h_nodes.fluidCount = static_cast<unsigned int>(fluidNodes.size());
    h_nodes.fluidI = static_cast<unsigned int*>(malloc(h_nodes.fluidCount * sizeof(unsigned int)));
    h_nodes.fluidAlloc = h_nodes.fluidCount;
    memcpy(h_nodes.fluidI, fluidNodes.data(), h_nodes.fluidCount * sizeof(unsigned int));

    assert(!h_nodes.interfaceI);
    h_nodes.interfaceCount = static_cast<unsigned int>(interfaceNodes.size());
    h_nodes.interfaceI = static_cast<unsigned int*>(malloc(h_nodes.interfaceCount * sizeof(unsigned int)));
    h_nodes.interfaceAlloc = h_nodes.interfaceCount;
    memcpy(h_nodes.interfaceI, interfaceNodes.data(), h_nodes.interfaceCount * sizeof(unsigned int));

    // Build a sorted active nodes list
    fluidNodes.insert(fluidNodes.end(), interfaceNodes.begin(), interfaceNodes.end());
    std::sort(fluidNodes.begin(), fluidNodes.end());

    // Array to buffer
    assert(!h_nodes.activeI);
    h_nodes.activeCount = static_cast<unsigned int>(fluidNodes.size());
    h_nodes.activeI = static_cast<unsigned int*>(malloc(h_nodes.activeCount * sizeof(unsigned int)));
    memcpy(h_nodes.activeI, fluidNodes.data(), h_nodes.activeCount * sizeof(unsigned int));
    
    cout << " done" << endl;
}

//
// Data ingress/egress from/to legacy structures
//
void LBOpenMP::syncDEMIn(const elmtList& elmts, const particleList& particles, const wallList& walls, const objectList& objects) {
    // Sync DEM data to structure of arrays format (and device memory)
    syncElementsIn(elmts);
    syncParticlesIn(particles);
    syncWallsIn(walls);
    syncObjectsIn(objects);
    // @todo is cylinder missing?
}
void LBOpenMP::syncDEMOut(elmtList& elmts, particleList& particles, wallList& walls, objectList& objects) const {
    // Sync DEM data from structure of arrays format (and device memory)
    // @todo which of this sync is redundant?
    syncElementsOut(elmts);
    syncParticlesOut(particles);
    syncWallsOut(walls);
    syncObjectsOut(objects);
    // @todo is cylinder missing?
}
bool LBOpenMP::syncElementsIn(const elmtList &elements) {
    bool componentsHasGrown = false;
    if (h_elements.alloc < elements.size()) {
        // Grow host buffers
         if (h_elements.x1) {
             free(h_elements.x1);
             free(h_elements.wGlobal);
             free(h_elements.FHydro);
             free(h_elements.MHydro);
             free(h_elements.fluidVolume);
         }
         h_elements.alloc = (unsigned int)elements.size();
         h_elements.x1 = (tVect*)malloc(elements.size() * sizeof(tVect));
         h_elements.wGlobal = (tVect*)malloc(elements.size() * sizeof(tVect));
         h_elements.FHydro = (tVect*)malloc(elements.size() * sizeof(tVect));
         h_elements.MHydro = (tVect*)malloc(elements.size() * sizeof(tVect));
         h_elements.fluidVolume = (double*)malloc(elements.size() * sizeof(double));
    }
    // Update size
    h_elements.count = (unsigned int)elements.size();
    // Repackage host particle data from array of structures, to structure of arrays
     for (unsigned int i = 0; i < h_elements.count; ++i) {
         h_elements.x1[i] = elements[i].x1;
         h_elements.wGlobal[i] = elements[i].wGlobal;
         h_elements.FHydro[i] = elements[i].FHydro;
         // h_elements.MHydro[i] = elements[i].MHydro; // This is zero'd before use in latticeBoltzmannStep()
         // h_elements.fluidVolume[i] = elements[i].fluidVolume; // This is zero'd before use in latticeBoltzmannStep()
     }
    // Construct the components storage
    {
        // Allocate memory for componentsData
        unsigned int totalComponents = 0;
        for (const auto& e : elements)
            totalComponents += static_cast<unsigned int>(e.components.size());
        if (!h_elements.componentsIndex || totalComponents >= h_elements.componentsIndex[elements.size()]) {
            if (h_elements.componentsData)
                free(h_elements.componentsData);
            h_elements.componentsData = (unsigned int*)malloc(totalComponents * sizeof(unsigned int));
            componentsHasGrown = true;
        }
        // Allocate componentsIndex if first pass
        if (!h_elements.componentsIndex)
            h_elements.componentsIndex = (unsigned int*)malloc((elements.size() + 1) * sizeof(unsigned int));
        // Fill componentsIndex and componentsData
        totalComponents = 0;
        for (int i = 0; i < elements.size(); ++i) {
            h_elements.componentsIndex[i] = totalComponents;
            if (!elements[i].components.empty()) {
                memcpy(h_elements.componentsData + totalComponents, elements[i].components.data(), elements[i].components.size() * sizeof(unsigned int));
                totalComponents += static_cast<unsigned int>(elements[i].components.size());
            }
        }
        h_elements.componentsIndex[elements.size()] = totalComponents;
    }
    return componentsHasGrown;
}
void LBOpenMP::syncParticlesIn(const particleList &particles) {
    if (h_particles.alloc < particles.size()) {
        // Grow host buffers
        if (h_particles.clusterIndex) {
            free(h_particles.clusterIndex);
            free(h_particles.r);
            free(h_particles.x0);
            free(h_particles.radiusVec);
        }
        h_particles.alloc = (unsigned int)particles.size();
        h_particles.clusterIndex = (unsigned int*)malloc(particles.size() * sizeof(unsigned int));
        h_particles.r = (double*)malloc(particles.size() * sizeof(double));
        h_particles.x0 = (tVect*)malloc(particles.size() * sizeof(tVect));
        h_particles.radiusVec = (tVect*)malloc(particles.size() * sizeof(tVect));
    }
    // Update size
    h_particles.count = static_cast<unsigned int>(particles.size());
    // Repackage host particle data from array of structures, to structure of arrays
    for (unsigned int i = 0; i < h_particles.count; ++i) {
        h_particles.clusterIndex[i] = particles[i].clusterIndex;
        h_particles.r[i] = particles[i].r;
        h_particles.x0[i] = particles[i].x0;
        h_particles.radiusVec[i] = particles[i].radiusVec;
    }
}
void LBOpenMP::syncCylindersIn(const cylinderList &cylinders) {
    if (h_cylinders.count < cylinders.size()) {
        // Grow host buffers
        if (h_cylinders.p1) {
            free(h_cylinders.p1);
            free(h_cylinders.p2);
            free(h_cylinders.R);
            free(h_cylinders.naxes);
            free(h_cylinders.omega);
            free(h_cylinders.moving);
        }
        h_cylinders.p1 = (tVect*)malloc(cylinders.size() * sizeof(tVect));
        h_cylinders.p2 = (tVect*)malloc(cylinders.size() * sizeof(tVect));
        h_cylinders.R = (double*)malloc(cylinders.size() * sizeof(double));
        h_cylinders.naxes = (tVect*)malloc(cylinders.size() * sizeof(tVect));
        h_cylinders.omega = (tVect*)malloc(cylinders.size() * sizeof(tVect));
        h_cylinders.moving = (bool*)malloc(cylinders.size() * sizeof(bool));
    }
    // Update size
    h_cylinders.count = (unsigned int)cylinders.size();
    // Repackage host particle data from array of structures, to structure of arrays
    for (unsigned int i = 0; i < h_cylinders.count; ++i) {
        h_cylinders.p1[i] = cylinders[i].p1;
        h_cylinders.p2[i] = cylinders[i].p2;
        h_cylinders.R[i] = cylinders[i].R;
        h_cylinders.naxes[i] = cylinders[i].naxes;
        h_cylinders.omega[i] = cylinders[i].omega;
        h_cylinders.moving[i] = cylinders[i].moving;
    }
}
void LBOpenMP::syncWallsIn(const wallList &walls) {
    if (h_walls.count < walls.size()) {
        // Grow host buffers
        if (h_walls.n) {
            free(h_walls.n);
            free(h_walls.p);
            free(h_walls.rotCenter);
            free(h_walls.omega);
            free(h_walls.vel);
            free(h_walls.FHydro);
        }
        h_walls.n = (tVect*)malloc(walls.size() * sizeof(tVect));
        h_walls.p = (tVect*)malloc(walls.size() * sizeof(tVect));
        h_walls.rotCenter = (tVect*)malloc(walls.size() * sizeof(tVect));
        h_walls.omega = (tVect*)malloc(walls.size() * sizeof(tVect));
        h_walls.vel = (tVect*)malloc(walls.size() * sizeof(tVect));
        h_walls.FHydro = (tVect*)malloc(walls.size() * sizeof(tVect));
    }
    // Update size
    h_walls.count = static_cast<unsigned int>(walls.size());
    // Repackage host particle data from array of structures, to structure of arrays
    for (unsigned int i = 0; i < h_walls.count; ++i) {
        h_walls.n[i] = walls[i].n;
        h_walls.p[i] = walls[i].p;
        h_walls.rotCenter[i] = walls[i].rotCenter;
        h_walls.omega[i] = walls[i].omega;
        h_walls.vel[i] = walls[i].vel;
        // h_walls.FHydro[i] = walls[i].FHydro; // Zero'd before use in streaming()
    }
}
void LBOpenMP::syncObjectsIn(const objectList &objects) {
    if (h_objects.count < objects.size()) {
        // Grow host buffers
        if (h_objects.r) {
            free(h_objects.r);
            free(h_objects.x0);
            free(h_objects.x1);
            free(h_objects.FHydro);
        }
        h_objects.r = (double*)malloc(objects.size() * sizeof(double));
        h_objects.x0 = (tVect*)malloc(objects.size() * sizeof(tVect));
        h_objects.x1 = (tVect*)malloc(objects.size() * sizeof(tVect));
        h_objects.FHydro = (tVect*)malloc(objects.size() * sizeof(tVect));
    }
    // Update size
    h_objects.count = static_cast<unsigned int>(objects.size());
    // Repackage host particle data from array of structures, to structure of arrays
    for (unsigned int i = 0; i < h_objects.count; ++i) {
        h_objects.r[i] = objects[i].r;
        h_objects.x0[i] = objects[i].x0;
        h_objects.x1[i] = objects[i].x1;
        // h_objects.FHydro[i] = objects[i].FHydro; // Zero'd before use in streaming()
    }
}
void LBOpenMP::syncElementsOut(elmtList &elements) const {
    // Assumes memory is already allocated
    // Repackage device particle data from structure of arrays to array of structures 
    for (unsigned int i = 0; i < h_elements.count; ++i) {
        elements[i].x1 = h_elements.x1[i];
        elements[i].wGlobal = h_elements.wGlobal[i];
        elements[i].FHydro = h_elements.FHydro[i];
        elements[i].MHydro = h_elements.MHydro[i];
        elements[i].fluidVolume = h_elements.fluidVolume[i];
    }
}
void LBOpenMP::syncParticlesOut(particleList &particles) const {
    // Assumes memory is already allocated
    // Repackage device particle data from structure of arrays to array of structures 
    for (unsigned int i = 0; i < h_particles.count; ++i) {
        particles[i].clusterIndex = h_particles.clusterIndex[i];
        particles[i].r = h_particles.r[i];
        particles[i].x0 = h_particles.x0[i];
        particles[i].radiusVec = h_particles.radiusVec[i];
    }
}
void LBOpenMP::syncCylindersOut(cylinderList &cylinders) const {
    // Assumes memory is already allocated
    // Repackage device particle data from structure of arrays to array of structures 
    for (unsigned int i = 0; i < h_cylinders.count; ++i) {
        cylinders[i].p1 = h_cylinders.p1[i];
        cylinders[i].p2 = h_cylinders.p2[i];
        cylinders[i].R = h_cylinders.R[i];
        cylinders[i].naxes = h_cylinders.naxes[i];
        cylinders[i].omega = h_cylinders.omega[i];
        cylinders[i].moving = h_cylinders.moving[i];
    }
}
void LBOpenMP::syncWallsOut(wallList &walls) const {
    // Assumes memory is already allocated
    // Repackage device particle data from structure of arrays to array of structures 
    for (unsigned int i = 0; i < h_walls.count; ++i) {
        walls[i].n = h_walls.n[i];
        walls[i].p = h_walls.p[i];
        walls[i].rotCenter = h_walls.rotCenter[i];
        walls[i].omega = h_walls.omega[i];
        walls[i].vel = h_walls.vel[i];
        walls[i].FHydro = h_walls.FHydro[i];
    }
}
void LBOpenMP::syncObjectsOut(objectList &objects) const {
    // Assumes memory is already allocated
    // Repackage device particle data from structure of arrays to array of structures 
    for (unsigned int i = 0; i < h_objects.count; ++i) {
        objects[i].r = h_objects.r[i];
        objects[i].x0 = h_objects.x0[i];
        objects[i].x1 = h_objects.x1[i];
        objects[i].FHydro = h_objects.FHydro[i];
    }
}

//
// latticeBoltzmannCouplingStep() subroutines
//
double LBOpenMP::initializeParticleBoundaries() {
    // Reset all nodes to outside
    memset(h_nodes.p, 0, h_nodes.count * sizeof(bool));

    double totalParticleMass = 0;
#pragma omp parallel for reduction(+:totalParticleMass) 
    for (unsigned int i = 0; i < h_nodes.activeCount; ++i) {
        // Pass the active node index to the common implementation
        totalParticleMass += common_initializeParticleBoundaries(i, &h_nodes, &h_particles);
    }
    return totalParticleMass;
}
void LBOpenMP::findNewActive() {
#pragma omp parallel for
    for (unsigned int i = 0; i < h_nodes.activeCount; ++i) {
        // Pass the active node index to the common implementation
        common_findNewActive(i, &h_nodes, &h_particles, &h_elements);
    }
}
void LBOpenMP::findNewSolid() {
#pragma omp parallel for
    for (unsigned int i = 0; i < h_nodes.activeCount; ++i) {
        // Pass the active node index to the common implementation
        common_findNewSolid(i, &h_nodes, &h_particles, &h_elements);
    }
}
void LBOpenMP::checkNewInterfaceParticles() {
#pragma omp parallel for
    for (unsigned int e_i = 0; e_i < h_elements.count; ++e_i) {
        common_checkNewInterfaceParticles(e_i, &h_nodes, &h_particles, &h_elements);
    }
}

//
// latticeBoltzmannStep() subroutines
//
void LBOpenMP::reconstructHydroCollide() {
    /**
     * reconstruct()
     * computeHydroForces()
     * collision()
     */
    // @todo the inside of this loop overlaps heavily with gpu subclass
#pragma omp parallel for
    for (unsigned int i = 0; i < h_nodes.activeCount; ++i) {
        // Convert index to active node index
        const unsigned int an_i = h_nodes.activeI[i];

        // reconstruction of macroscopic variables from microscopic distribution
        // this step is necessary to proceed to the collision step
        h_nodes.reconstruct(an_i);

        // compute interaction forces
        if (h_nodes.count) {
            common_computeHydroForces(an_i, &h_nodes, &h_particles, &h_elements);
        }

        //collision operator
        h_nodes.collision(an_i);
    }
}
void LBOpenMP::streaming() {
    // STREAMING STEP
    // Init forces to zero
    h_walls.initForces<CPU>();
    h_objects.initForces<CPU>();
    // Init streaming support vector
    h_nodes.store<CPU>();

#pragma omp parallel for // @note extraMass reduction is not currently implemented
    for (unsigned int i = 0; i < h_nodes.activeCount; ++i) {
        common_streaming(i, &h_nodes, &h_walls);
    }

    // redistributing extra mass due to bounce back to interface cells
    // redistributeMass(extraMass);  // @todo extraMass hasn't been implemented properly
}
void LBOpenMP::shiftToPhysical() {
    for (unsigned int i = 0; i < h_elements.count; ++i) {
        h_elements.FHydro[i] *= PARAMS.unit.Force;
        h_elements.MHydro[i] *= PARAMS.unit.Torque;
        h_elements.fluidVolume[i] *= PARAMS.unit.Volume;
    }
    for (unsigned int i = 0; i < h_walls.count; ++i) {
        h_walls.FHydro[i] *= PARAMS.unit.Force;
    }
    for (unsigned int i = 0; i < h_objects.count; ++i) {
        h_objects.FHydro[i] *= PARAMS.unit.Force;
    }
}

//
// latticeBoltzmannFreeSurfaceStep() subroutines
//
void LBOpenMP::enforceMassConservation() {
    // calculate total mass of active nodes
    double thisMass = 0.0;
    for (unsigned int i = 0; i < h_nodes.activeCount; ++i) {
        const unsigned int an_i = h_nodes.activeI[i];
        if (!h_nodes.isInsideParticle(an_i)) {
            thisMass += h_nodes.mass[an_i];
        }
    }

    // mass deficit
    const double massDeficit = (thisMass - PARAMS.totalMass);

    // fix it
    redistributeMass(-0.01 * massDeficit);
}
void LBOpenMP::redistributeMass(const double& massSurplus) {
    const double addMass = massSurplus / h_nodes.interfaceCount;

#pragma omp parallel for
    for (unsigned int i = 0; i < h_nodes.interfaceCount; ++i) {
        const unsigned int in_i = h_nodes.interfaceI[i];
        h_nodes.mass[in_i] += addMass;
    }
}
void LBOpenMP::updateMass() {
#pragma omp parallel for
    for (unsigned int i = 0; i < h_nodes.interfaceCount; ++i) {
        // Convert index to active node index
        const unsigned int in_i = h_nodes.interfaceI[i];
        common_updateMassInterface(in_i, &h_nodes);
    }
#pragma omp parallel for
    for (unsigned int i = 0; i < h_nodes.fluidCount; ++i) {
        // Convert index to active node index
        const unsigned int fn_i = h_nodes.fluidI[i];
        common_updateMassFluid(fn_i, &h_nodes);
    }
}
void LBOpenMP::updateInterface() {
    // Initialise reduction variable
    double h_massSurplus = 0.0;
#pragma omp parallel for
    for (unsigned int i = 0; i < h_nodes.interfaceCount; ++i) {
        // Convert index to interface node index
        const unsigned int in_i = h_nodes.interfaceI[i];
        // filling lists of mutant nodes and changing their type
        common_findInterfaceMutants(in_i, &h_nodes);
    }
#pragma omp parallel for
    for (unsigned int i = 0; i < h_nodes.interfaceCount; ++i) {
        // Convert index to interface node index
        const unsigned int in_i = h_nodes.interfaceI[i];
        // fixing the interface (always one interface between fluid and gas)
        common_smoothenInterface_find(in_i, &h_nodes);
    }
    // Build temporary list of new/interface/new_gas nodes
    // Returns a buffer, where first element is length of the list
    const std::vector<unsigned int> tempList = buildTempNewList(lbmDirec * h_nodes.interfaceCount);
    for (const unsigned int &in_i : tempList) {
        // Convert index to interface node index
        // fixing the interface (always one interface between fluid and gas)
        common_smoothenInterface_update(in_i, &h_nodes);
    }
#pragma omp parallel for
    for (unsigned int i = 0; i < h_nodes.interfaceCount; ++i) {
        // Convert index to interface node index
        const unsigned int in_i = h_nodes.interfaceI[i];
        // updating characteristics of mutant nodes
        common_updateMutants(in_i, &h_nodes, &h_massSurplus);
    }
    // Rebuild interface list
    buildInterfaceList(lbmDirec * h_nodes.activeCount);
#pragma omp parallel for
    for (unsigned int i = 0; i < h_nodes.interfaceCount; ++i) {
        // Convert index to interface node index
        const unsigned int in_i = h_nodes.interfaceI[i];
        // remove isolated interface cells (both surrounded by gas and by fluid)
        common_removeIsolated(in_i, &h_nodes, &h_massSurplus);
    }
    // Rebuild all lists
    buildAllLists(h_nodes.interfaceCount, h_nodes.fluidCount + lbmDirec * h_nodes.interfaceCount);
    // distributing surplus to interface cells
    this->redistributeMass(h_massSurplus);
#ifdef DEBUG
    // computeSurfaceNormal()
#endif
}
std::vector<unsigned int> LBOpenMP::buildTempNewList(const unsigned int& _max_len) {
    const unsigned int max_len = min(_max_len, h_nodes.count);
    std::vector<unsigned int> t_list;
    t_list.reserve(max_len);
    for (unsigned int i = 0; i < count; ++i) {
        if (h_nodes.type[i] == GAS_TO_INTERFACE || h_nodes.type[i] == FLUID_TO_INTERFACE) {
            t_list.push_back(i);
        }
    }
    return t_list;
}
void LBOpenMP::buildInterfaceList(const unsigned int& _max_len) {
    const unsigned int max_len = min(_max_len, h_nodes.count);
    std::vector<unsigned int> t_interfaceList;
    t_interfaceList.reserve(max_len);
    for (unsigned int i = 0; i < count; ++i) {
        if (h_nodes.type[i] == INTERFACE) {
            t_interfaceList.push_back(i);
        }
    }
    if (h_nodes.interfaceAlloc < t_interfaceList.size()) {
        free(h_nodes.interfaceI);
        h_nodes.interfaceI = static_cast<unsigned int*>(malloc(t_interfaceList.size() * sizeof(unsigned int)));
    }
    h_nodes.interfaceCount = static_cast<unsigned int>(t_interfaceList.size());
    memcpy(h_nodes.interfaceI, t_interfaceList.data(), t_interfaceList.size() * sizeof(unsigned int));
}
void LBOpenMP::buildAllLists(const unsigned int& _max_interface_len, const unsigned int& _max_fluid_len) {
    const unsigned int max_interface_len = min(_max_interface_len, h_nodes.count);
    const unsigned int max_fluid_len = min(_max_fluid_len, h_nodes.count);
    const unsigned int max_len = min(_max_interface_len + _max_fluid_len, h_nodes.count);
    std::vector<unsigned int> t_interfaceList;
    t_interfaceList.reserve(max_interface_len);
    std::vector<unsigned int> t_fluidList;
    t_fluidList.reserve(max_fluid_len);
    std::vector<unsigned int> t_activeList;
    t_activeList.reserve(max_len);
    for (unsigned int i = 0; i < count; ++i) {
        if (h_nodes.type[i] == LIQUID) {
            t_fluidList.push_back(i);
            t_activeList.push_back(i);
        } else if (h_nodes.type[i] == INTERFACE) {
            t_interfaceList.push_back(i);
            t_activeList.push_back(i);
        }
    }
    if (h_nodes.interfaceAlloc < t_interfaceList.size()) {
        free(h_nodes.interfaceI);
        h_nodes.interfaceI = static_cast<unsigned int*>(malloc(t_interfaceList.size() * sizeof(unsigned int)));
    }
    if (h_nodes.fluidAlloc < t_fluidList.size()) {
        free(h_nodes.fluidI);
        h_nodes.fluidI = static_cast<unsigned int*>(malloc(t_fluidList.size() * sizeof(unsigned int)));
    }
    if (h_nodes.activeAlloc < t_activeList.size()) {
        free(h_nodes.activeI);
        h_nodes.activeI = static_cast<unsigned int*>(malloc(t_activeList.size() * sizeof(unsigned int)));
    }
    h_nodes.interfaceCount = static_cast<unsigned int>(t_interfaceList.size());
    h_nodes.fluidCount = static_cast<unsigned int>(t_fluidList.size());
    h_nodes.activeCount = static_cast<unsigned int>(t_activeList.size());
    memcpy(h_nodes.interfaceI, t_interfaceList.data(), t_interfaceList.size() * sizeof(unsigned int));
    memcpy(h_nodes.fluidI, t_fluidList.data(), t_fluidList.size() * sizeof(unsigned int));
    memcpy(h_nodes.activeI, t_activeList.data(), t_activeList.size() * sizeof(unsigned int))
}
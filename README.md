# HYBIRD

HYBIRD is a numerical code that combines a Lattice-Boltzmann Method (LBM) solver for fluid dynamics with a Discrete Element Method (DEM) solver for the solution of the dynamics of a particle ensemble. Both DEM and LBM are parallelized using OpenMP loops, and therefore optimized for a few processors, sharing memory.


## Building

HYBIRD uses CMake, as a cross-platform process, for configuring and generating build directives, e.g. Makefile.

### Requirements

Building HYBIRD has the following requirements.

* [CMake](https://cmake.org/) `>=3.22`
* C++17 capable C++ compiler
  * Windows compilation is only supported with [MinGW](https://www.mingw-w64.org/), Visual Studio is not supported.
  * Linux/WSL2: GCC >= 8.0

### Building with CMake

Most modern IDEs have CMake integration, the below guide explains the usage at a high level with command-line examples.

Building via CMake is a 3 step process.

1. Create a build directory, the build files and compiled executable will all be stored here.
2. Configure CMake, configure any build options (such as debug vs release) and it will generate the Makefile.
3. Compilation, this can be triggered via CMake, make, or your IDE.

The following example begins in the root of the cloned repository.

```bash
# Create the build directory and change into it
mkdir -p build && cd build

# Configure CMake from the command line passing configure-time options.
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build the required targets, in this case all of them, using 8 threads.
cmake --build . --parallel 8

# Alternatively make can be invoked directly
make -j8
```

### CMake Configuration Options

| Option                               | Value                       | Description                                                                                                |
| -------------------------------------| --------------------------- | ---------------------------------------------------------------------------------------------------------- |
| `CMAKE_BUILD_TYPE`                   | `Release` / `Debug` / `MinSizeRel` / `RelWithDebInfo` | Select the build configuration for single-target generators such as `make`       |

## User Guide

The user guide can be found [here](https://github.com/gnomeCreative/HYBIRD/wiki)
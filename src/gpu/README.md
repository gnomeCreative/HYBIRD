# HYBIRD - GPU

HYBIRD's GPU (CUDA) implementation was developed with the intention of moving the full combined LBM-DEM model to utilise GPU acceleration.

The original intent of the codebase was to have a hybrid codebase, which could be compiled either with or without CUDA support. Such that models could be developed/tested without CUDA availability, and later executed at full-scale on HPC or similar. Due to time-constraints, the CPU support is incomplete relative to GPU support and has not been tested. As such, the CUDA implementation should remain in a distinct branch (or even be moved to a separate repository) so that the original CPU/OpenMP implementation is still easily available.

The CUDA LBM has been validated through visual testing of the existing tutorial models, it is assumed this leaves a large amount of the code not covered by tests.

Several minor issues (#13, #14, #15, #16, 19, #21) were identified during development, for consistency with the CPU implementation they were not addressed in the CUDA implementation.

## Supported Features

Due to time constraints 3/4 of the model has been ported to execute using CUDA, meaning that the GPU LBM model (including coupling and free surface) is able to execute, whilst the DEM model remains executing using the original CPU code. CUDA DEM was ported, available in `cuda_dem` branch, which runs with both LBM and DEM using CUDA, however it was found due to the small scales of the DEM CUDA execution was slower than CPU and it would be difficult to reach parity without increasing scale towards 10k particles if not higher.

Only currently Paraview output files have been updated to be exported by the GPU implementation (within `IO2.cu`). Likewise, simulation checkpoint/resume has not been updated. Updating remaining IO methods to work with the CUDA implementation should not require significant technical understanding, all that has changed is the memory layout.

LBM curves were not implemented, due to them being created and destroyed with nodes during runtime the existing implementation was not suitable for CUDA. It may be possible to create all possible curves at init, including for GAS nodes (where they would sit idle), however it's unclear to me whether this would be suitable for moving cylinders (as presumably the curve should move with the cylinder). Further thought would be required to map curves to cylinders.

## Implementation Philosophy

The goal of the implementation was to retain as much of the original model code as possible, although changes were required due to the difference of GPU architectures.

### Memory Layout

**Structure of Arrays**

The original implementation can be described as Array of Structures (AoS), such that memory is laid out with arrays/lists/vectors of structures.

```c++
struct Foo{
 double a;
 double b;
};
std::vector<Foo> fooList;
...
fooList[0].a = 12.0;
```

As vectors are not practical in GPU code, memory has been restructured to a Structure of Arrays (SoA) format.

```c++
struct Foo{
  double *a;
  double *b;
};
Foo fooList;
...
fooList.a[0] = 12.0;
```

This does make accessing variables more verbose, however it greatly improves the performance of memory accesses (as neighbouring threads access neighbouring memory locations).

This does however become more technical for `Node` which has several member variables which are arrays (`d`, `f`, `fs`).

These are still represented as pointer arrays in the SoA format, however the length of this array is multiplied by the length of the array.

Within this array, `d` has been stored differently to `f` and `fs`.

* `d`: Is stored such that the first N items are each nodes element `[0]`
* `f`/`fs`: Are stored such that the first node's full array, is followed by the second node's full array.


For example they should be accessed

```c+
for (int i = 0; i<nodes->count; ++i) {
    for (int j = 1; j<lbmDirec; ++j) {
        nodes->d[i * lbmDirec + j]  = 0;
        nodes->f[j * nodes->count + i]  = 0;
        nodes->fs[j * nodes->count + i]  = 0;
    }
}
```

*This is an easy source of bugs, it may be worth creating accessors on `Node2` which take parameters `index` and `j`.*

`d` is the proper, and more performant approach, however `f` was implemented this way for convenience.

**CUDA Pointers**

All data structures have three copies defined within LB2, e.g.

```c++
    // Host copy of node buffers, may not always be current whilst in CUDA mode
    Node2 h_nodes;
    // Host copy of device node buffer pointers
    Node2 hd_nodes;
    // Pointer to device copy of device node buffer pointers
    // In CPU, this is a pointer to h_nodes
    Node2 *d_nodes = nullptr;
```

This is a typical CUDA pattern, where data is duplicated on host (CPU) and device (GPU).

`h_` prefix denotes a structure of the host mirror of data, it's buffer may be out of sync with what's on device.
`hd_` prefix denotes a host copy of the struct held on device, it's counts will typically be correct but the pointers are device pointers so can't be accessed on host without a host-device memcpy.
`d_` prefix denotes a device pointer to the device copy of the struct, this is frequently passed to CUDA kernels.

active/fluid/interface lists are held within the `Node2` struct, rather than as distinct lists. This is mostly for convenience of passing them back and forth.

**Dense Matrix**

The original HYBIRD implementation could be described as sparse matrix, such that empty (`GAS`) nodes are not represented. This approach greatly increases the complexity of a GPU parallel implementation, and therefore the GPU implementation moved to a dense matrix whereby the environment is represented by a number of nodes which is fixed at the time of initialisation, many of these are simply marked `GAS` and do not participate in the model. This comes at a cost of additional memory utilisation, as in most scenarios there are more GAS nodes than other types combined, however some rough maths suggested this would remain well within the available 80GB VRAM on current HPC GPUs (e.g. A100/H100).

In order to access fluid/interface/active nodes, there are still corresponding lists, these merely hold indexes into the main node array.

Likewise, as it's a dense matrix, nodes no longer track their coordinate as it matches their index in the array.

### Model Implementation

## hybird.cpp

This is the original file, with minor changes to redirect the configuration/execution to the LB2 implementation.

The method `goCycle()` was moved into `LB2` as `LB2::step()`. 

## LB2

This class is a CUDA reimplementation of `LB`, the core implementation should be familiar however additional scaffolding has been added to support the hybrid CPU-CUDA implementation and for synchronising memory between host and device.

For example, the method `LB::streaming()` now has several components.

The main entrypoint streaming is now templated.

```c++
class LB2 {
...
    template<int impl>
    void streaming();
```

Within `LB2.cu` it has two implementations

```c++
/**
 * This is the main implementation of streaming() within the active node loop
 * It has been copied from LB.cpp with memory accesses and reductions updated
 * It is used by both the CPU and CUDA implementations
 */
__host__ __device__ __forceinline__ void common_streaming(const unsigned int i, Node2* nodes, Wall2* walls) {
...
}
/**
 * This shell of a CUDA kernel, acts like the active node for loop
 * Each unique CUDA thread calls common_streaming() with it's index
 */
__global__ void d_streaming(Node2* d_nodes, Wall2* d_walls) {
    // Get unique CUDA thread index, which corresponds to active node 
    const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    // Kill excess threads early
    if (i >= d_nodes->activeCount) return;

    common_streaming(i, d_nodes, d_walls);
}
/**
 * The CUDA implementation launches the kernel d_streaming
 * First it must calculate how many threads to launch
 * The <<<gridSize, blockSize>>> launches the CUDA kernel on the GPU
 */
template<>
void LB2::streaming<CUDA>() {
...
    // Launch cuda kernel to update
    int blockSize = 0;  // The launch configurator returned block size
    int minGridSize = 0;  // The minimum grid size needed to achieve the // maximum occupancy for a full device // launch
    int gridSize = 0;  // The actual grid size needed, based on input size
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, d_streaming, 0, hd_nodes.activeCount);
    // Round up to accommodate required threads
    gridSize = (hd_nodes.activeCount + blockSize - 1) / blockSize;
    d_streaming<<<gridSize, blockSize>>>(d_nodes, d_walls);
    CUDA_CHECK();
...
}
/**
 * The CPU implementation wraps common_streaming() in a for loop
 */
template<>
void LB2::streaming<CPU>() {
...
#pragma omp parallel for 
    for (unsigned int i = 0; i < d_nodes->activeCount; ++i) {
        common_streaming(i, d_nodes, d_walls);
    }
}
```

All trivially parallelisable methods (e.g. no race conditions) have been implemented in this way.

For methods with minor reductions, they have also been implemented in this way with their reductions implemented via parallel-safe atomic operations.

Atomic operations are handled differently between CUDA and OpenMP, therefore pre-processor macros have been used so each compiler only sees their respective atomic operations. For example within `common_streaming()` you will find.

```
#ifdef __CUDA_ARCH__
                // CUDA atomics
                atomicAdd(&walls->FHydro[solidIndex].x, BBforce.x);
                atomicAdd(&walls->FHydro[solidIndex].y, BBforce.y);
                atomicAdd(&walls->FHydro[solidIndex].z, BBforce.z);
#else
                // CPU atomics
#pragma omp atomic update
                walls->FHydro[solidIndex] += BBforce;
#endif
```

Other methods, particularly those within the free surface step, which are not trivially parallelisable have required a modified implementation. For most of these methods, they have not been implemented within the CPU implementation, instead a `@todo` note will likely be found.

**Reconstruct/computeHydroForces/Collision**

These three methods have been grouped into a `reconstructHydroCollide()` method.
Collision was converted to a member function of Node.

**shiftToPhysical**

The shifting of forces inside elements/walls/objects to physical units has been merged into a common method which now runs after streaming.

**enforceMassConservation**

This method reduces over all active nodes that are not inside particles, to gain a little bit more performance than performing an atomic reduction. Thrust, a CUDA parallel primitives, library is used. Due to the indexes of active notes being stored in a secondary array, a utility struct `unmapper` has been implemented to allow this reduction to map between active node indexes and implement the condition over inside particle nodes allowing the desired data to be reduced.

**updateInterface**

This method is rather complicated so had several changes to make it safer in parallel, and was not implemented in CPU due to time-constraints/prioritisation.

* Nodes to do not instantly change between type, nor do filled/empty lists exists, instead we have introduced several intermediate types `INTERFACE_EMPTY`, `INTERFACE_FILLED`, `GAS_TO_INTERFACE`, `FLUID_TO_INTERFACE`. Nodes temporarily change to these types denoting their state, whilst remaining in their original lists. Before the end of `updateInterface()` no nodes should remain in these states.
* For parity with the original model, the `Node::isActive()` check has been updated to consider most of these states as active.
* `smoothenInterface()` has been split into two methods. 
    * `smoothenInterface_find()`: Within this method nodes to be smoothed have their type updated to `GAS_TO_INTERFACE` or `FLUID_TO_INTERFACE`. Additionally the source node is set inside `d[0]` (this is otherwise unused but exists). There's a small potential for race condition here if multiple nodes attempt to update the same node simultaneously, however it's expected that this would be harmless. It could be addressed with an atomic compare and swap, but this would require type to be backed by an `int` rather than `char`.
    * `smoothenInterface_update()`: This second method detects nodes in the `GAS_TO_INTERFACE` or `FLUID_TO_INTERFACE` state, and completes their transformation to their new type. This avoids the larger race condition, where multiple nodes would attempt to generate/convert a node simultaneously. A list of all nodes in this state is built and used by `buildTempNewList()` between the two functions.
* The rebuilding of lists, in CUDA, occurs with the `buildFluidList()` and `buildInterfaceList()` methods, these poll all nodes to report their index if they are of the correct type. `buildActiveList()` then concatenates these two lists to produce the active node list.
    * Note, this feels like `cleanLists()` is therefore redundant.

**cleanLists**

This has been updated to combine the interface and fluid lists and sort it. Sort them in order of index which should improve memory locality.

Again, cleanLists() was converted to a member function of Node.

**coupling**

This implementation is naive, each active node tests each particle for intersection. If the number of particle's scales up this could become a bottleneck, it would make more sense for nodes to check a small number of particles according to the neighbour grid.

## IO2

This class extends the original `IO.cpp`, to inherit from it's initialisation and member variables. Due to changes in memory layout, all IO methods require updating. Currently only Paraview methods have been updated.

### Redundant Files

`LB.h`/`LB.cpp` and `Node.h`/`Node.cpp` which both wholy concern the LBM are now theoretically redundant, they are only required due to `IO2.cu` inheriting from `IO.CU`. If this inheritance were removed, they can be removed from the project without any harm.

Other structures, e.g. `Element`, `Particle`, `Wall`, `Cylinder` exist twice. They have an array of structure format within the original CPU DEM code, and an array of structure format in the CUDA LBM code.

### Support

Contact robert.chisholm@sheffield.ac.uk at the Research Software Engineering team who wrote the CUDA implementation.
Or book a code clinic appointment, where another member of the team may be able to provide support.

https://rse.shef.ac.uk/


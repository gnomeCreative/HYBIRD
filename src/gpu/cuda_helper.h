#ifndef CUDAHELPERH
#define CUDAHELPERH

#include <cstdio>
#include <cstdlib>

// For use in C++ templates to specify version
enum Implementation {
#ifdef USE_CUDA
    CUDA,
#endif
    CPU
};

#ifdef USE_CUDA
#ifdef _MSC_VER
#pragma warning( disable: 4505 )    // unreferenced function has been removed (under visual studio)
#endif
#include <cuda_runtime.h>
#include <thrust/transform_reduce.h>

/**
 * Error check function for safe CUDA API calling
 * Wrap all calls to CUDA API functions with CUDA_CALL() to catch errors on failure
 * e.g. CUDA_CALL(cudaFree(myPtr));
 * CUDA_CHECk() can also be used to perform error checking after kernel launches and async methods
 * e.g. CUDA_CHECK()
 */
#if defined(_DEBUG) || defined(D_DEBUG)
#define CUDA_CALL(ans) { gpuAssert((ans), __FILE__, __LINE__); }
#define CUDA_CHECK() { gpuAssert(cudaDeviceSynchronize(), __FILE__, __LINE__); }
#else
#define CUDA_CALL(ans) { gpuAssert((ans), __FILE__, __LINE__); }
#define CUDA_CHECK() { gpuAssert(cudaPeekAtLastError(), __FILE__, __LINE__); }
#endif
inline void gpuAssert(cudaError_t code, const char *file, int line) {
    if (code != cudaSuccess) {
        if (line >= 0) {
            fprintf(stderr, "\nCUDA Error: %s(%d): %s\n", file, line, cudaGetErrorString(code));
        } else {
            fprintf(stderr, "\nCUDA Error: %s(%d): %s\n", file, line, cudaGetErrorString(code));
        }
        exit(EXIT_FAILURE);
    }
}

#define IMPL Implementation::CUDA
#else
#define IMPL Implementation::CPU
// Define these CUDA Symbols empty when building for CPU
#ifndef __host__
#define __host__
#define __device__
#define __forceinline__
#endif
#endif  // USE_CUDA

// Constexpr arrays in lattice.h require marking as __host__device__ to build available to both host and device code
#ifdef __CUDA_ARCH__
#define __host__device__ __device__
#else
#define __host__device__ 
#endif
#endif  // CUDAHELPERH

#ifndef CUDAHELPERH
#define CUDAHELPERH

#include <cstdio>
#include <cstdlib>

#ifdef USE_CUDA
#include <cuda_runtime.h>

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
            fprintf(stderr, "CUDA Error: %s(%d): %s", file, line, cudaGetErrorString(code));
        } else {
            fprintf(stderr, "CUDA Error: %s(%d): %s", file, line, cudaGetErrorString(code));
        }
        exit(EXIT_FAILURE);
    }
}

#define IMPL Implementation::CUDA
#else
#define IMPL Implementation::CPU
// Define these CUDA Symbols empty when building for CPU
#define __host__
#define __device__
#endif  // USE_CUDA

#endif  // CUDAHELPERH

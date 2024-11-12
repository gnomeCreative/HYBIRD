#ifndef CUB_TEMP_MEM_H
#define CUB_TEMP_MEM_H

/**
 * Based on CubTempMem from https://github.com/FLAMEGPU/FLAMEGPU2
 * /include/flamegpu/simulation/detail/CubTempMem.cuh
 * It has been extended into a singleton pattern
 */

#include <unordered_map>
#include <utility>
#include <memory>

#include "cuda_helper.h"

/**
 * cub is a CUDA library, that provides highly optimised GPU parallel primitive
 * operations such as sorting and reduction. Most, if not all, of these parallel
 * primitives require a small temporary memory allocation.
 *
 * The purpose of this class is to manage a shared, singleton, for that buffer
 */
class CubTempMem {
 public:
    /**
     * Singletons should not be cloneable.
     */
    CubTempMem(CubTempMem &other) = delete;
    /**
     * Singletons should not be assignable.
     */
    void operator=(const CubTempMem &) = delete;
    /**
     * Singleton accessor (intended for cub temp storage)
     */
    static CubTempMem &GetTempSingleton() {
        if(!_singletonT)
            _singletonT = std::make_unique<CubTempMem>();
        return *_singletonT;
    }
    /**
     * Singleton accessor (intended for swap buffers)
     */
    static CubTempMem &GetBufferSingleton() {
        if (!_singletonB)
            _singletonB = std::make_unique<CubTempMem>();
        return *_singletonB;
    }
    /**
     * Singleton destructor
     */
    static void FreeSingletons() {
        _singletonT.reset();
        _singletonB.reset();
    }
    /**
     * Release any allocated memory
     */
    ~CubTempMem() { CUDA_CALL(cudaFree(d_cub_temp)); }
    /**
     * Grow the size of the allocated buffer
     */
    void resize(size_t newSize) {
        if (newSize > d_cub_temp_size) {
            if (d_cub_temp) {
                CUDA_CALL(cudaFree(d_cub_temp));
            }
            CUDA_CALL(cudaMalloc(&d_cub_temp, newSize));
            d_cub_temp_size = newSize;
        }
    }
    /**
     * Return the buffer's pointer
     */
    void *getPtr() const { return d_cub_temp; }
    /**
     * Return the buffer's size
     */
    size_t &getSize() const { d_cub_temp_size_rtn = d_cub_temp_size; return d_cub_temp_size_rtn; }
    /**
     * Swap the stored pointer with the provided
     */
    template<typename T>
    size_t swapPtr(T *&_ptr, size_t size) {
        void* t_ptr = static_cast<void*>(_ptr);
        _ptr = static_cast<T*>(d_cub_temp);
        d_cub_temp = t_ptr;
        size_t t_sz = d_cub_temp_size;
        d_cub_temp_size = size;
        d_cub_temp_size_rtn = size;
        return t_sz;
    }

 protected:
    CubTempMem() = default;
    /**
     * Singleton storage
     */
    static std::unique_ptr<CubTempMem> _singletonT;
    static std::unique_ptr<CubTempMem> _singletonB;

 private:
    void* d_cub_temp = nullptr;
    size_t d_cub_temp_size = 0;
    // We have this version, so it can be passed directly to cub (requires non-const reference)
    // It is simply overwritten every time it is requested, incase it gets accidentally changed
    mutable size_t d_cub_temp_size_rtn = 0;
};

#endif  // CUB_TEMP_MEM_H

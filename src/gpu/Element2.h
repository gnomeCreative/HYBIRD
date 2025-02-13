#ifndef ELEMENT2_H
#define ELEMENT2_H
#include "cuda_helper.h"

/**
 * \brief Representation of objects, such as walls, within the DEM.
 * @todo/note original member vars wsolver, index, active are fixed at runtime
 */
struct Element2 {
    // The total number of elements
    // @note This is currently fixed at runtime, but may become variable in future
    unsigned int count = 0;
    // Size of allocated buffer, as it never shrinks
    unsigned int alloc = 0;

    // velocity of the center of the element
    tVect *x1 = nullptr;
    // angular velocity in global and local reference frame
    tVect *wGlobal = nullptr;
    // forces
    tVect *FHydro = nullptr;
    // moments
    tVect *MHydro = nullptr;
    // fluid mass entrained by the particle
    double *fluidVolume = nullptr;
    /**
     * elmt::components is a std::vector, whereby each elmt may have a different length
     *
     * As such to represent it in CUDA we store all of these vectors compactly in order
     * inside Element2::componentsData
     *
     * In order to access a particular element's components that element's componentsIndex is read
     * This value is the first index within componentsData containing component data for the element.
     * The following elements componentsIndex points to the item after this elements last component
     * e.g.
     * auto e_i = 1  // element index
     * for (auto i = componentsIndex[e_i]; i < componentsIndex[e_i+1]; ++i)
     *     auto component = componentsData[i];
     */
    unsigned int *componentsIndex = nullptr; // length == count+1
    unsigned int *componentsData = nullptr;  // length == componentsIndex[count]

    template<int impl>
    void initElements();
};

template <>
inline void Element2::initElements<CPU>() {
    //initializing this time step hydrodynamic force
    memset(FHydro, 0, sizeof(tVect) * count);
    memset(MHydro, 0, sizeof(tVect) * count);
    // initializing the fluid mass for buoyancy
    memset(fluidVolume, 0, sizeof(double) * count);
}
#ifdef USE_CUDA
template <>
inline void Element2::initElements<CUDA>() {
    //initializing this time step hydrodynamic force
    CUDA_CALL(cudaMemset(FHydro, 0, sizeof(tVect) * count));
    CUDA_CALL(cudaMemset(MHydro, 0, sizeof(tVect) * count));
    // initializing the fluid mass for buoyancy
    CUDA_CALL(cudaMemset(fluidVolume, 0, sizeof(double) * count));
}
#endif

#endif  // ELEMENT2_H

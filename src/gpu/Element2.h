#ifndef ELEMENT2_H
#define ELEMENT2_H

/**
 * \brief Representation of objects, such as walls, within the DEM.
 * @todo/note original member vars wsolver, index, active are fixed at runtime
 */
struct Element2 {
    // The total number of elements
    // @note This is currently fixed at runtime, but may become variable in future
    unsigned int count = 0;

    tVect *FHydro = nullptr;
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
};
#endif  // ELEMENT2_H

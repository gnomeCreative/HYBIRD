#ifndef PROBLEM_H
#define PROBLEM_H
#include <string>
#include <vector>

#include <exprtk.hpp>

#include "utils.h"

/**
 * Represents a problem input file, and supports loading them from .yml format
 * Problems contain fluid volumes, walls and cylinders
 */
class Problem {
 public:
    struct MinMaxPair {
        tVect min, max;
    };
    std::string file;
    std::string name;
    std::vector<MinMaxPair> fluids_basic;
    std::vector<exprtk::expression<double>> fluids_complex;
    // String representation of loaded expressions
    std::vector<std::string> fluids_complex_str;
    std::vector<wall> walls;
    std::vector<cylinder> cylinders;
    /**
     * Init exprtk components
     */
    Problem();
    /**
     * These objects are uncopyable
     * Copying breaks the linkages within the symbol table
     * @note It should be possible to write a manual copy constructor to rebuild complex fluids on copy
     */
    Problem(const Problem&) = delete;
    /**
     * Attempt to load the specified .yml problem file into the current problem
     * This method may only be called once
     * @param filePath Path to .yml containing problem
     * @note Static "factory" because constructors are not supposed to throw exceptions (e.g. if file not exist)
     */
    void loadFile(const std::string &filePath);

    /**
     * @param pos Coordinate to test
     * @return True if the node is fluid
     * @note Defined here so that it is inlined as this may be called millions of times during init
     */
    bool isFluid(const tVect& pos, const double& unitLength) const {
        // scaled position
        const tVect scaledPos = pos * unitLength;
        // Check basic fluid volumes
        if (std::any_of(fluids_basic.cbegin(), fluids_basic.cend(),
            // Lambda function
            [&pos, &scaledPos](const MinMaxPair& minmax) {
                return (scaledPos.x >= minmax.min.x && scaledPos.x <= minmax.max.x
                     && scaledPos.y >= minmax.min.y && scaledPos.y <= minmax.max.y
                     && scaledPos.z >= minmax.min.z && scaledPos.z <= minmax.max.z);
            })) {
                return true;
        }
        // Check complex fluid volumes
        exprtk_pos = scaledPos;
        return std::any_of(fluids_complex.cbegin(), fluids_complex.cend(),
            // Lambda function
            [](const exprtk::expression<double>& expr) {return expr.value() > 0; });
    }
private:
    mutable tVect exprtk_pos;
    exprtk::symbol_table<double> symbol_table;
    bool is_loaded = false;
};

#endif // PROBLEM_H

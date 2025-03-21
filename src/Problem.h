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
     * Attempt to load the specified .yml problem file and return it as a Problem object
     * @param filePath Path to .yml containing problem
     * @note Static "factory" because constructors are not supposed to throw exceptions (e.g. if file not exist)
     */
    static Problem loadFile(const std::string &filePath);

    /**
     * @param pos Coordinate to test
     * @return True if the node is fluid
     * @note Defined here so that it is inlined as this may be called millions of times during init
     */
    bool isFluid(const tVect& pos) const {
        // Check basic fluid volumes
        if (std::any_of(fluids_basic.cbegin(), fluids_basic.cend(),
            // Lambda function
            [&pos](const MinMaxPair& minmax) {
                return (pos.x >= minmax.min.x && pos.x <= minmax.max.x
                     && pos.y >= minmax.min.y && pos.y <= minmax.max.y
                     && pos.z >= minmax.min.z && pos.z <= minmax.max.z);
            })) {
                return true;
        }
        // Check complex fluid volumes
        exprtk_pos = pos;
        return std::any_of(fluids_complex.cbegin(), fluids_complex.cend(),
            // Lambda function
            [](const exprtk::expression<double>& expr) { return expr.value() > 0; });
    }
private:
    mutable tVect exprtk_pos;
    exprtk::symbol_table<double> symbol_table;
};

#endif // PROBLEM_H

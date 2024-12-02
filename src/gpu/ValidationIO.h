#ifndef VALIDATIONIO_H
#define VALIDATIONIO_H

#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <numeric>

#include "gpu/LB2.h"
#include "gpu/Node2.h"

/**
 * Very basic class for outputting validation CSVs to compare two version of the model
 *
 * Each run outputs it's own validation file
 * Each column corresponds to a different validation variable
 * Each row corresponds to a different step of the model run
 *
 * If there is stochastic behaviour within the model
 * it may be necessary to combine the average of multiple runs before comparing versions
 */
class ValidationIO {
    std::filesystem::path output_path;
    std::fstream fs;
public:
    ValidationIO(const std::string& output_directory) {
        output_path = output_directory;
        output_path /= "validation.csv";
    }
    ~ValidationIO() {
        if (fs.is_open()) {
            fs.close();
        }
    }
    void output(LB2& lb) {
        if (!fs.is_open()) {
            fs.open(output_path, std::ios::in | std::ios::binary);
            if (!fs.is_open()) {
                throw std::exception(("Unable to open validation file for writing: " + output_path.generic_string()).c_str());
            }
            // Write the header
            fs << "node_count"                      << ",";
            fs << "coord_average"                   << ",";
            fs << "n_average"                       << ",";
            fs << "uX_average"                      << ",";
            fs << "uY_average"                      << ",";
            fs << "uZ_average"                      << ",";
            fs << "hydroForceX_average"             << ",";
            fs << "hydroForceY_average"             << ",";
            fs << "hydroForceZ_average"             << ",";
            fs << "mass_average"                    << ",";
            fs << "visc_average"                    << ",";
            fs << "type_average"                    << ",";
            fs << "p_average"                       << ",";
            fs << "active_count"                    << ",";
            fs << "active_coord_average"            << ",";
            fs << "interface_count"                 << ",";
            fs << "interface_coord_average"         << ",";
            fs << "fluid_count"                     << ",";
            fs << "fluid_coord_average"             << ",";
            fs << "wall_count"                      << ",";
            fs << "wall_coord_average"              << "\n";
        }
        const Node2 nodes = lb.getNodes();
        // Write a new data row
        // node_count
        fs << nodes.count << ",";
        // coord_average
        fs << std::accumulate(nodes.coord, nodes.coord + nodes.count, static_cast<uint64_t>(0))/static_cast<float>(nodes.count) << ",";
        // n_average
        fs << std::accumulate(nodes.n, nodes.n + nodes.count, 0.0) / static_cast<float>(nodes.count) << ",";
        // uX_average, uY_average, uZ_average
        tVect t = std::accumulate(nodes.u, nodes.u + nodes.count, Zero);
        fs << t.x / static_cast<float>(nodes.count) << ",";
        fs << t.y / static_cast<float>(nodes.count) << ",";
        fs << t.z / static_cast<float>(nodes.count) << ",";
        // hydroForceX_average, hydroForceY_average, hydroForceZ_average
        t = std::accumulate(nodes.hydroForce, nodes.hydroForce + nodes.count, Zero);
        fs << t.x / static_cast<float>(nodes.count) << ",";
        fs << t.y / static_cast<float>(nodes.count) << ",";
        fs << t.z / static_cast<float>(nodes.count) << ",";
        // mass_average"
        fs << std::accumulate(nodes.mass, nodes.mass + nodes.count, 0.0) / static_cast<float>(nodes.count) << ",";
        // visc_average
        fs << std::accumulate(nodes.visc, nodes.visc + nodes.count, 0.0) / static_cast<float>(nodes.count) << ",";
        // type_average
        fs << std::accumulate(nodes.type, nodes.type + nodes.count, static_cast<uint64_t>(0)) / static_cast<float>(nodes.count) << ",";
        // p_average
        fs << std::count(nodes.p, nodes.p + nodes.count, true) / static_cast<float>(nodes.count) << ",";
        // active_count
        fs << nodes.activeCount << ",";
        // active_coord_average
        uint64_t sum = 0;
        for (int i = 0; i < nodes.activeCount; ++i)
            sum += nodes.coord[nodes.activeI[i]];
        fs << sum / static_cast<float>(nodes.activeCount) << ",";
        // interface_count
        fs << nodes.interfaceCount << ",";
        // interface_coord_average
        sum = 0;
        for (int i = 0; i < nodes.interfaceCount; ++i)
            sum += nodes.coord[nodes.interfaceI[i]];
        fs << sum / static_cast<float>(nodes.interfaceCount) << ",";
        // fluid_count
        fs << nodes.fluidCount << ",";
        // fluid_coord_average
        sum = 0;
        for (int i = 0; i < nodes.fluidCount; ++i)
            sum += nodes.coord[nodes.fluidI[i]];
        fs << sum / static_cast<float>(nodes.fluidCount) << ",";
        // wall_count
        fs << nodes.wallCount << ",";
        // wall_coord_average
        sum = 0;
        for (int i = 0; i < nodes.wallCount; ++i)
            sum += nodes.coord[nodes.wallI[i]];
        fs << sum / static_cast<float>(nodes.wallCount) << "\n";
    }
    void output(LB& lb) {
        if (!fs.is_open()) {
            fs.open(output_path, std::ios::in | std::ios::binary);
            if (!fs.is_open()) {
                throw std::exception(("Unable to open validation file for writing: " + output_path.generic_string()).c_str());
            }
            // Write the header
            fs << "node_count"                      << ",";
            fs << "coord_average"                   << ",";
            fs << "n_average"                       << ",";
            fs << "uX_average"                      << ",";
            fs << "uY_average"                      << ",";
            fs << "uZ_average"                      << ",";
            fs << "hydroForceX_average"             << ",";
            fs << "hydroForceY_average"             << ",";
            fs << "hydroForceZ_average"             << ",";
            fs << "mass_average"                    << ",";
            fs << "visc_average"                    << ",";
            fs << "type_average"                    << ",";
            fs << "p_average"                       << ",";
            fs << "active_count"                    << ",";
            fs << "active_coord_average"            << ",";
            fs << "interface_count"                 << ",";
            fs << "interface_coord_average"         << ",";
            fs << "fluid_count"                     << ",";
            fs << "fluid_coord_average"             << ",";
            fs << "wall_count"                      << ",";
            fs << "wall_coord_average"              << "\n";
        }
        lb.nodes;
        // Write a new data row
        // node_count
        fs << lb.nodes.size() << ",";
        // coord_average
        uint64_t su = 0;
        for (const auto& [_, nd] : lb.nodes)
            su += nd.coord;
        fs << su / static_cast<float>(lb.nodes.size()) << ",";
        // n_average
        double sd = 0;
        for (const auto& [_, nd] : lb.nodes)
            sd += nd.n;
        fs << sd / static_cast<float>(lb.nodes.size()) << ",";
        // uX_average, uY_average, uZ_average
        tVect sv = Zero;
        for (const auto& [_, nd] : lb.nodes)
            sv += nd.u;
        fs << sv.x / static_cast<float>(lb.nodes.size()) << ",";
        fs << sv.y / static_cast<float>(lb.nodes.size()) << ",";
        fs << sv.z / static_cast<float>(lb.nodes.size()) << ",";
        // hydroForceX_average, hydroForceY_average, hydroForceZ_average
        sv = Zero;
        for (const auto& [_, nd] : lb.nodes)
            sv += nd.hydroForce;
        fs << sv.x / static_cast<float>(lb.nodes.size()) << ",";
        fs << sv.y / static_cast<float>(lb.nodes.size()) << ",";
        fs << sv.z / static_cast<float>(lb.nodes.size()) << ",";
        // mass_average
        sd = 0;
        for (const auto& [_, nd] : lb.nodes)
            sd += nd.mass;
        fs << sd / static_cast<float>(lb.nodes.size()) << ",";
        // visc_average
        sd = 0;
        for (const auto& [_, nd] : lb.nodes)
            sd += nd.visc;
        fs << sd / static_cast<float>(lb.nodes.size()) << ",";
        // type_average
        su = 0;
        for (const auto& [_, nd] : lb.nodes)
            su += nd.type;
        fs << su / static_cast<float>(lb.nodes.size()) << ",";
        // p_average
        su = 0;
        for (const auto& [_, nd] : lb.nodes)
            su += nd.p ? 1: 0;
        fs << su / static_cast<float>(lb.nodes.size()) << ",";
        // active_count
        fs << lb.activeNodes.size() << ",";
        // active_coord_average
        su = 0;
        for (const auto &nd : lb.activeNodes)
            su += nd->coord;
        fs << su / static_cast<float>(lb.nodes.size()) << ",";
        // interface_count
        fs << lb.interfaceNodes.size() << ",";
        // interface_coord_average
        su = 0;
        for (const auto& nd : lb.interfaceNodes)
            su += nd->coord;
        fs << su / static_cast<float>(lb.nodes.size()) << ",";
        // fluid_count
        fs << lb.fluidNodes.size() << ",";
        // fluid_coord_average
        su = 0;
        for (const auto& nd : lb.fluidNodes)
            su += nd->coord;
        fs << su / static_cast<float>(lb.nodes.size()) << ",";
        // wall_count
        fs << lb.wallNodes.size() << ",";
        // wall_coord_average
        su = 0;
        for (const auto& nd : lb.wallNodes)
            su += nd->coord;
        fs << su / static_cast<float>(lb.nodes.size()) << ",";
    }
};

#endif  // VALIDATIONIO_H

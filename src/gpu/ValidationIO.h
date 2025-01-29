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
            fs.open(output_path, std::ios::out | std::ios::binary | std::ios::trunc);
            if (!fs.is_open()) {
                fprintf(stderr, "Unable to open validation file for writing: %s\n", output_path.generic_string().c_str());
                throw std::exception();
            }
            // Write the header
            fs << "node_count"                      << ",";
            fs << "n_average"                       << ",";
            fs << "n_minimum"                       << ",";
            fs << "n_maximum"                       << ",";
            fs << "uX_average"                      << ",";
            fs << "uX_minimum"                      << ",";
            fs << "uX_maximum"                      << ",";
            fs << "uY_average"                      << ",";
            fs << "uY_minimum"                      << ",";
            fs << "uY_maximum"                      << ",";
            fs << "uZ_average"                      << ",";
            fs << "uZ_minimum"                      << ",";
            fs << "uZ_maximum"                      << ",";
            fs << "hydroForceX_average"             << ",";
            fs << "hydroForceY_average"             << ",";
            fs << "hydroForceZ_average"             << ",";
            fs << "centrifugalForceX_average"       << ",";
            fs << "centrifugalForceY_average"       << ",";
            fs << "centrifugalForceZ_average"       << ",";
            fs << "mass_average"                    << ",";
            fs << "visc_average"                    << ",";
            fs << "visc_minimum"                    << ",";
            fs << "visc_maximum"                    << ",";
            fs << "type_average"                    << ",";
            fs << "p_average"                       << ",";
            fs << "active_count"                    << ",";
            fs << "active_coord_average"            << ",";
            fs << "interface_count"                 << ",";
            fs << "interface_coord_average"         << ",";
            fs << "fluid_count"                     << ",";
            fs << "fluid_coord_average"             << ",";
            fs << "wall_count"                      << ",";
            fs << "wall_coord_average"              << ",";
            for (int i = 0; i < lbmDirec; ++i) {
                fs << "f[" << i << "]_average"      << ",";
                fs << "f[" << i << "]_minimum"      << ",";
                fs << "f[" << i << "]_maximum"      << ",";
            }
            for (int i = 0; i < lbmDirec; ++i)
                fs << "d[" << i << "]" << ",";
            fs << "lbFX" << ",";
            fs << "lbFY" << ",";
            fs << "lbFZ" << ",";
            fs << "activenode_type_liquid_count" << ",";
            fs << "coord_average" << "\n";
        }
        const Node2 nodes = lb.getNodes();
        // Write a new data row
        // node_count
        unsigned int sd_count = 0;
        for (unsigned int i = 0; i < nodes.count; ++i) {
            if (nodes.type[i] != GAS) {
                ++sd_count;
            }
        }
        fs << sd_count << ",";
        // n_average
        double sd_avg = 0;
        double sd_min = std::numeric_limits<double>::max();
        double sd_max = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < nodes.count; ++i) {
            if (nodes.type[i] != GAS) {
                ++sd_count;
                sd_avg += nodes.n[i];
                sd_min = std::min(sd_min, nodes.n[i]);
                sd_max = std::max(sd_max, nodes.n[i]);
            }
        }
        fs << sd_avg / sd_count << ",";
        fs << sd_min << ",";
        fs << sd_max << ",";
        // uX_average, uY_average, uZ_average
        tVect t_avg = Zero;
        tVect t_min = tVect(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
        tVect t_max = tVect(-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());
        for (unsigned int i = 0; i < nodes.count; ++i) {
            if (nodes.type[i] != GAS) {
                t_avg.x = nodes.u[i].x;
                t_avg.y = nodes.u[i].y;
                t_avg.z = nodes.u[i].z;
                t_min.x = std::min(t_min.x, nodes.u[i].x);
                t_min.y = std::min(t_min.y, nodes.u[i].y);
                t_min.z = std::min(t_min.z, nodes.u[i].z);
                t_max.x = std::max(t_max.x, nodes.u[i].x);
                t_max.y = std::max(t_max.y, nodes.u[i].y);
                t_max.z = std::max(t_max.z, nodes.u[i].z);
            }
        }
        fs << t_avg.x / sd_count << ",";
        fs << t_min.x << ",";
        fs << t_max.x << ",";
        fs << t_avg.y / sd_count << ",";
        fs << t_min.y << ",";
        fs << t_max.y << ",";
        fs << t_avg.z / sd_count << ",";
        fs << t_min.z << ",";
        fs << t_max.z << ",";
        // hydroForceX_average, hydroForceY_average, hydroForceZ_average
        t_avg = Zero;
        for (unsigned int i = 0; i < nodes.count; ++i) {
            if (nodes.type[i] != GAS) {
                t_avg.x = nodes.hydroForce[i].x;
                t_avg.y = nodes.hydroForce[i].y;
                t_avg.z = nodes.hydroForce[i].z;
            }
        }
        fs << t_avg.x / sd_count << ",";
        fs << t_avg.y / sd_count << ",";
        fs << t_avg.z / sd_count << ",";
        // centrifugalForceX_average, centrifugalForceY_average, centrifugalForceZ_average
        t_avg = Zero;
        for (unsigned int i = 0; i < nodes.count; ++i) {
            if (nodes.type[i] != GAS) {
                t_avg.x = nodes.centrifugalForce[i].x;
                t_avg.y = nodes.centrifugalForce[i].y;
                t_avg.z = nodes.centrifugalForce[i].z;
            }
        }
        fs << t_avg.x / sd_count << ",";
        fs << t_avg.y / sd_count << ",";
        fs << t_avg.z / sd_count << ",";
        // mass_average
        sd_avg = 0;
        for (unsigned int i = 0; i < nodes.count; ++i) {
            if (nodes.type[i] != GAS) {
                sd_avg += nodes.mass[i];
            }
        }
        fs << sd_avg / sd_count << ",";
        // visc_average
        sd_avg = 0;
        sd_min = std::numeric_limits<double>::max();
        sd_max = -std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < nodes.count; ++i) {
            if (nodes.type[i] != GAS) {
                sd_avg += nodes.visc[i];
                sd_min = std::min(sd_min, nodes.visc[i]);
                sd_max = std::max(sd_max, nodes.visc[i]);
            }
        }
        fs << sd_avg / sd_count << ",";
        fs << sd_min << ",";
        fs << sd_max << ",";
        // type_average
        sd_avg = 0;
        for (unsigned int i = 0; i < nodes.count; ++i) {
            if (nodes.type[i] != GAS) {
                sd_avg += nodes.type[i];
            }
        }
        fs << sd_avg / sd_count << ",";
        // p_average
        sd_avg = 0;
        for (unsigned int i = 0; i < nodes.count; ++i) {
            if (nodes.type[i] != GAS && nodes.p[i]) {
                sd_avg += 1;
            }
        }
        fs << sd_avg << ",";
        // active_count
        fs << nodes.activeCount << ",";
        // active_coord_average
        uint64_t sum = 0;
        for (int i = 0; i < nodes.activeCount; ++i)
            sum += nodes.activeI[i];
        fs << sum / static_cast<float>(nodes.activeCount) << ",";
        // interface_count
        fs << nodes.interfaceCount << ",";
        // interface_coord_average
        sum = 0;
        for (int i = 0; i < nodes.interfaceCount; ++i)
            sum += nodes.interfaceI[i];
        fs << sum / static_cast<float>(nodes.interfaceCount) << ",";
        // fluid_count
        fs << nodes.fluidCount << ",";
        // fluid_coord_average
        sum = 0;
        for (int i = 0; i < nodes.fluidCount; ++i)
            sum += nodes.fluidI[i];
        fs << sum / static_cast<float>(nodes.fluidCount) << ",";
        // wall_count
        fs << nodes.wallCount << ",";
        // wall_coord_average
        sum = 0;
        for (int i = 0; i < nodes.wallCount; ++i)
            sum += nodes.wallI[i];
        fs << sum / static_cast<float>(nodes.wallCount) << ",";
        // f
        for (int j = 0; j < lbmDirec; ++j) {
            sd_avg = 0;
            sd_min = std::numeric_limits<double>::max();
            sd_max = -std::numeric_limits<double>::max();
            for (int i = 0; i < nodes.count; ++i) {
                if (nodes.type[i] != GAS) {
                    sd_avg += nodes.f[i * lbmDirec + j];
                    sd_min = std::min(sd_min, nodes.f[i * lbmDirec + j]);
                    sd_max = std::max(sd_max, nodes.f[i * lbmDirec + j]);
                }
            }
            fs << sd_avg / sd_count << ",";
            fs << sd_min << ",";
            fs << sd_max << ",";
        }
        // d (we must convert index to coord)
        for (int j = 0; j < lbmDirec; ++j) {
            sum = 0;
            unsigned int active_count = 0;
            for (int i = 0; i < nodes.count; ++i)
                if (nodes.isActive(i)) {
                    ++active_count;
                    if (nodes.d[j * nodes.count + i] != std::numeric_limits<unsigned int>::max())
                        sum += nodes.d[j * nodes.count + i];                    
                }
            fs << sum / static_cast<float>(active_count) << ",";
        }
        // lbf
        fs << h_PARAMS.lbF.x << ",";
        fs << h_PARAMS.lbF.y << ",";
        fs << h_PARAMS.lbF.z << ",";
        // type==LIQUID
        sum = 0;
        for (int i = 0; i < nodes.activeCount; ++i)
            if (nodes.type[nodes.activeI[i]] == LIQUID)
                ++sum;
        fs << sum << ",";
        // coord_average
        //fs << std::accumulate(nodes.coord, nodes.coord + nodes.count, static_cast<uint64_t>(0)) / static_cast<float>(nodes.count) << "\n";
        fs << 0 << "\n";
    }
};

#endif  // VALIDATIONIO_H

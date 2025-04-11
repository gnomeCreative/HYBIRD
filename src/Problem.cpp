#include "Problem.h"

#include <yaml-cpp/yaml.h>
#include <exprtk.hpp>

Problem::Problem() {
    symbol_table.add_variable("x", exprtk_pos.x);
    symbol_table.add_variable("y", exprtk_pos.y);
    symbol_table.add_variable("z", exprtk_pos.z);
}

inline bool parse_bool(const YAML::Node& node, const std::string &parent, const std::string &child) {
    if (node[child].IsSequence()) {
        std::stringstream err;
        err << "'" << parent << "::" << child << "', if present, should be parseable as type boolean.";
        std::cout << err.str();
        throw std::runtime_error(err.str().c_str());
    }
    return node[child].as<bool>();
}
inline tVect parse_tVect(const YAML::Node& node, const std::string& parent, const std::string& child) {
    if (!node[child].IsSequence() || node[child].size() != 3) {
        std::stringstream err;
        err << "'" << parent << "::" << child << "', if present, should be parseable as type double[3].";
        std::cout << err.str();
        throw std::runtime_error(err.str().c_str());
    }
    return { node[child][0].as<double>(), node[child][1].as<double>() , node[child][2].as<double>() };
}
inline exprtk::expression<double> parse_expression(exprtk::parser<double> &parser, exprtk::symbol_table<double> &symbol_table, const std::string &expression_string) {
    exprtk::expression<double> expression;
    expression.register_symbol_table(symbol_table);
    if (!parser.compile(expression_string, expression)) {
        std::stringstream err;
        err << "Unable to parse expression '" << expression_string << "'.";
        std::cout << err.str();
        throw std::runtime_error(err.str().c_str());
    }
    return expression;
}
Problem::MinMaxPair readBasicFluid(const YAML::Node &basic_fluid) {
    // Validate it contains expected components
    if (!basic_fluid["min"] || !basic_fluid["max"]
        || !basic_fluid["min"].IsSequence() || !basic_fluid["max"].IsSequence()
        || basic_fluid["min"].size() != 3 || basic_fluid["max"].size() != 3) {
        const std::string err = "'fluids_basic' nodes should contain both a 'min' and 'max' type double[3].";
        std::cout << err;
        throw std::runtime_error(err.c_str());
    }
    Problem::MinMaxPair rtn;
    rtn.min.x = basic_fluid["min"][0].as<double>();
    rtn.min.y = basic_fluid["min"][1].as<double>();
    rtn.min.z = basic_fluid["min"][2].as<double>();
    rtn.max.x = basic_fluid["max"][0].as<double>();
    rtn.max.y = basic_fluid["max"][1].as<double>();
    rtn.max.z = basic_fluid["max"][2].as<double>();
    return rtn;
}
wall readWall(const YAML::Node& wall_node) {
    // Validate it contains mandatory components
    if (!wall_node["normal"] || !wall_node["point"]
        || !wall_node["normal"].IsSequence() || !wall_node["point"].IsSequence()
        || wall_node["normal"].size() != 3 || wall_node["point"].size() != 3) {
        const std::string err = "'walls' nodes should contain both a 'normal' and 'point' as type double[3].";
        std::cout << err;
        throw std::runtime_error(err.c_str());
    }
    wall rtn;
    // Mandatory components
    rtn.n.x = wall_node["normal"][0].as<double>();
    rtn.n.y = wall_node["normal"][1].as<double>();
    rtn.n.z = wall_node["normal"][2].as<double>();
    rtn.p.x = wall_node["point"][0].as<double>();
    rtn.p.y = wall_node["point"][1].as<double>();
    rtn.p.z = wall_node["point"][2].as<double>();
    // Optional components
    if (wall_node["moving"]) {
        rtn.moving = parse_bool(wall_node, "wall", "moving");
    }
    if (wall_node["rotational_center"]) {
        rtn.rotCenter = parse_tVect(wall_node, "wall", "rotational_center");
    }
    if (wall_node["rotational_speed"]) {
        rtn.omega = parse_tVect(wall_node, "wall", "rotational_speed");
    }
    if (wall_node["translation_speed"]) {
        rtn.vel = parse_tVect(wall_node, "wall", "translation_speed");
    }
    if (wall_node["slipping"]) {
        rtn.slip = parse_bool(wall_node, "wall", "slipping");
    }
    if (wall_node["translating"]) {
        rtn.translating = parse_bool(wall_node, "wall", "translating");
    }
    if (wall_node["translation_vector"]) {
        rtn.trans = parse_tVect(wall_node, "wall", "translation_vector");
    }
    if (wall_node["limited"]) {
        rtn.limited = parse_bool(wall_node, "wall", "limited");
    }
    if (wall_node["min"]) {
        const tVect t = parse_tVect(wall_node, "wall", "min");
        rtn.xMin = t.x;
        rtn.yMin = t.y;
        rtn.zMin = t.z;
    }
    if (wall_node["max"]) {
        const tVect t = parse_tVect(wall_node, "wall", "max");
        rtn.xMax = t.x;
        rtn.yMax = t.y;
        rtn.zMax = t.z;
    }
    return rtn;
}
cylinder readCylinder(const YAML::Node& cylinder_node) {
    // Validate it contains mandatory components
    if (!cylinder_node["point1"] || !cylinder_node["point2"]
        || !cylinder_node["point1"].IsSequence() || !cylinder_node["point2"].IsSequence()
        || cylinder_node["point1"].size() != 3 || cylinder_node["point2"].size() != 3) {
        const std::string err = "'cylinders' nodes should contain both a 'point1' and 'point2' as type double[3].";
        std::cout << err;
        throw std::runtime_error(err.c_str());
    }
    if (!cylinder_node["radius"] || cylinder_node["radius"].IsSequence()) {
        const std::string err = "'cylinders' nodes should contain a 'radius' as type double.";
        std::cout << err;
        throw std::runtime_error(err.c_str());
    }
    cylinder rtn;
    // Mandatory components
    rtn.p1.x = cylinder_node["point1"][0].as<double>();
    rtn.p1.y = cylinder_node["point1"][1].as<double>();
    rtn.p1.z = cylinder_node["point1"][2].as<double>();
    rtn.p2.x = cylinder_node["point2"][0].as<double>();
    rtn.p2.y = cylinder_node["point2"][1].as<double>();
    rtn.p2.z = cylinder_node["point2"][2].as<double>();
    rtn.R = cylinder_node["radius"].as<double>();
    // Optional components
    if (cylinder_node["moving"]) {
        rtn.moving = parse_bool(cylinder_node, "cylinder", "moving");
    }
    if (cylinder_node["rotational_speed"]) {
        rtn.omega = parse_tVect(cylinder_node, "cylinder", "rotational_speed");
    }
    if (cylinder_node["slipping"]) {
        rtn.slip = parse_bool(cylinder_node, "cylinder", "slipping");
    }
    if (cylinder_node["empty"]) {
        rtn.type = parse_bool(cylinder_node, "cylinder", "empty") ? EMPTY : FULL;
    }
    if (cylinder_node["translating"]) {
        rtn.translating = parse_bool(cylinder_node, "cylinder", "translating");
    }
    if (cylinder_node["translation_vector"]) {
        rtn.trans = parse_tVect(cylinder_node, "cylinder", "translation_vector");
    }
    if (cylinder_node["limited"]) {
        rtn.limited = parse_bool(cylinder_node, "cylinder", "limited");
    }
    if (cylinder_node["min"]) {
        const tVect t = parse_tVect(cylinder_node, "cylinder", "min");
        rtn.xMin = t.x;
        rtn.yMin = t.y;
        rtn.zMin = t.z;
    }
    if (cylinder_node["max"]) {
        const tVect t = parse_tVect(cylinder_node, "cylinder", "max");
        rtn.xMax = t.x;
        rtn.yMax = t.y;
        rtn.zMax = t.z;
    }
    return rtn;
}


void Problem::loadFile(const std::string& filePath) {
    if (is_loaded) {
        throw std::runtime_error("Problem::loadFile() called twice on same object.");
    }
    YAML::Node problem;
    try {
        problem = YAML::LoadFile(filePath);
    } catch (YAML::BadFile&) {
        std::stringstream err;
        err << "Unable to open file '" << filePath << "' for reading, does it exist?";
        std::cout << err.str();
        throw std::runtime_error(err.str().c_str());
    }
    this->file = filePath;
    // Name
    if (problem["name"]) {
        this->name = problem["name"].as<std::string>();
    }
    // Cuboid fluid volumes
    if (problem["fluids_basic"]) {
        if (problem["fluids_basic"].IsSequence()) {
            for (const auto &fluid : problem["fluids_basic"]) {
                this->fluids_basic.push_back(readBasicFluid(fluid));
            }
        } else {
            // They didn't provide it as a list, parse anyway
            readBasicFluid(problem["fluids_basic"]);
        }
    }
    // Mathematical expression fluid volumes
    if (problem["fluids_complex"]) {
        exprtk::parser<double> parser;
        if (problem["fluids_complex"].IsSequence()) {
            for (const auto &fluid : problem["fluids_complex"]) {
                const std::string t = fluid.as<std::string>();
                this->fluids_complex_str.push_back(t);
                this->fluids_complex.emplace_back(parse_expression(parser, this->symbol_table, t));
            }
        } else {
            // They didn't provide it as a list, parse anyway
            const std::string t = problem["fluids_complex"].as<std::string>();
            this->fluids_complex_str.push_back(t);
            this->fluids_complex.emplace_back(parse_expression(parser, this->symbol_table, t));
        }        
    }
    // Walls
    if (problem["walls"]) {
        if (problem["walls"].IsSequence()) {
            for (const auto &wall : problem["walls"]) {
                this->walls.push_back(readWall(wall));
            }
        } else {
            // They didn't provide it as a list, parse anyway
            this->walls.push_back(readWall(problem["walls"]));
        }
        // Correct indices
        for (size_t i = 0; i < this->walls.size(); ++i) {
            this->walls[i].index = static_cast<unsigned int>(i);
        }
    }
    // Cylinders
    if (problem["cylinders"]) {
        if (problem["cylinders"].IsSequence()) {
            for (const auto& cylinder : problem["cylinders"]) {
                this->cylinders.push_back(readCylinder(cylinder));
            }
        } else {
            // They didn't provide it as a list, parse anyway
            this->cylinders.push_back(readCylinder(problem["cylinders"]));
        }
        // Correct indices
        for (size_t i = 0; i < this->cylinders.size(); ++i) {
            this->cylinders[i].index = static_cast<unsigned int>(i);
        }
    }
    is_loaded = true;
}
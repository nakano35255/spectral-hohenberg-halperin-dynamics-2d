#include "param_parser.h"

#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

// ---------------------------------------------------------------------- //
namespace {
    // ---------------------------------------------------------------------- //
    std::string join_tokens(const std::vector<std::string>& tokens, std::size_t first) {
        if (first >= tokens.size()) {
            return "";
        }

        std::string result = tokens[first];
        for (std::size_t i = first + 1; i < tokens.size(); ++i) {
            result += " ";
            result += tokens[i];
        }
        return result;
    }
    // ---------------------------------------------------------------------- //
    Variable::Type parse_variable_type(const std::string& token, const std::string& context) {
        if (token == "const") {
            return Variable::Type::Constant;
        }
        if (token == "expr") {
            return Variable::Type::Expression;
        }
        throw std::runtime_error(context + " expects const|expr");
    }
    // ---------------------------------------------------------------------- //
    Variable parse_variable_from_tail(const std::vector<std::string>& tokens, std::size_t type_index, const std::string& context) {
        if (tokens.size() <= type_index + 1) {
            throw std::runtime_error(context + " needs: const <value> or expr <expression>");
        }

        const Variable::Type type = parse_variable_type(tokens[type_index], context);
        const std::string text = join_tokens(tokens, type_index + 1);
        if (text.empty()) {
            throw std::runtime_error(context + " has empty value");
        }
        return Variable(type, text);
    }
    // ---------------------------------------------------------------------- //
    bool parse_on_off(const std::string& value, const std::string& context) {
        if (value == "on") return true;
        if (value == "off") return false;
        throw std::runtime_error(context + " expects on|off");
    }
    // ---------------------------------------------------------------------- //
} // namespace
// ---------------------------------------------------------------------- //
std::vector<std::string> ParamParser::tokenize(const std::string& line) const {
    std::vector<std::string> tokens;
    std::string clean_line = line.substr(0, line.find('#'));
    
    std::istringstream iss(clean_line);
    std::string word;
    while (iss >> word) {
        tokens.push_back(word);
    }
    return tokens;
}
// ---------------------------------------------------------------------- //
void ParamParser::execute_command(const std::vector<std::string>& tokens) {
    if (tokens.empty()) {
        return;
    }

    const std::string& cmd = tokens[0];

    // --------------------------------------------------
    // 1. System Setup
    // --------------------------------------------------
    if (cmd == "dimension") {
        if (tokens.size() < 2) throw std::runtime_error("dimension syntax: dimension <num>");
        if (std::stoi(tokens[1]) != 2) {
            throw std::runtime_error("This solver only supports dimension 2.");
        }
    } 
    else if (cmd == "components") {
        if (tokens.size() < 2) throw std::runtime_error("components syntax: components <num>");
        params.physics.num_components = std::stoi(tokens[1]);
        if (params.physics.num_components <= 0) {
            throw std::runtime_error("components must be positive.");
        }
        params.physics.resize(params.physics.num_components);
    }
    else if (cmd == "boundary") {
        if (tokens.size() != 3) throw std::runtime_error("boundary syntax: boundary <x> <y>");
        if ((tokens[1] != "p" && tokens[1] != "periodic") || (tokens[2] != "p" && tokens[2] != "periodic")) {
            throw std::runtime_error("This solver only supports periodic (p) boundary.");
        }
    }
    else if (cmd == "grid") {
        if (tokens.size() != 3) throw std::runtime_error("grid syntax: grid <Nx> <Ny>");
        params.grid.num_grid[0] = std::stoi(tokens[1]);
        params.grid.num_grid[1] = std::stoi(tokens[2]);
        if (params.grid.num_grid[0] <= 0 || params.grid.num_grid[1] <= 0) {
            throw std::runtime_error("grid sizes must be positive");
        }
        if (params.grid.num_grid[0] % 2 != 0 || params.grid.num_grid[1] % 2 != 0) {
            throw std::runtime_error("grid must be even.");
        }
    }
    else if (cmd == "length") {
        if (tokens.size() != 3) throw std::runtime_error("length syntax: length <Lx> <Ly>");
        params.grid.length[0] = std::stod(tokens[1]);
        params.grid.length[1] = std::stod(tokens[2]);
        if (params.grid.length[0] <= 0.0 || params.grid.length[1] <= 0.0) {
            throw std::runtime_error("domain lengths must be positive");
        }
    }

    // --------------------------------------------------
    // 2. Integrator
    // --------------------------------------------------
    else if (cmd == "timestep") {
        if (tokens.size() != 2) throw std::runtime_error("timestep syntax: timestep <dt>");
        params.runtime.dt = std::stod(tokens[1]);
        if (params.runtime.dt <= 0.0) {
            throw std::runtime_error("timestep must be positive");
        }
    }
    else if (cmd == "integrator") {
        if (tokens.size() != 2) throw std::runtime_error("integrator syntax: integrator <name>");
        params.runtime.integrator_type = tokens[1];
    }
    
    // --------------------------------------------------
    // 3. Physics
    // --------------------------------------------------
    else if (cmd == "variable") {
        parse_variable_command(tokens);
    }
    else if (cmd == "free_energy") {
        params.physics.free_energy_entry = parse_variable_from_tail(tokens, 1, "free_energy");
    }
    else if (cmd == "L_coeff") {
        parse_l_coeff_command(tokens);
    }
    else if (cmd == "eta") {
        params.physics.eta_entry = parse_variable_from_tail(tokens, 1, "eta");
    }
    else if (cmd == "zeta") {
        params.physics.zeta_entry = parse_variable_from_tail(tokens, 1, "zeta");
    }

    // --------------------------------------------------
    // 4. 初期条件 (set) と 出力 (dump)
    // --------------------------------------------------
    else if (cmd == "set") {
        parse_set_command(tokens);
    }
    else if (cmd == "thermo") {
        parse_thermo_command(tokens);
    }
    else if (cmd == "restart") {
        parse_restart_command(tokens);
    }

    // --------------------------------------------------
    // 6. フラグ管理と実行 (run)
    // --------------------------------------------------
    else if (cmd == "fix") {
        parse_fix_command(tokens);
    }
    else if (cmd == "measure") {
        parse_measure_command(tokens);
    }
    else if (cmd == "run") {
        if (tokens.size() != 2) throw std::runtime_error("run syntax: run <steps>");
        const int steps = std::stoi(tokens[1]);
        if (steps <= 0) {
            throw std::runtime_error("run steps must be positive");
        }

        Command command;
        command.type = Command::Type::Run;
        command.run.steps = steps;
        params.commands.push_back(command);
    }
    
    else {
        throw std::runtime_error("Unknown command: " + cmd);
    }
}
// ---------------------------------------------------------------------- //
void ParamParser::parse_variable_command(const std::vector<std::string>& tokens) {
    if (tokens.size() < 4) throw std::runtime_error("variable syntax: variable <name> const <value> | variable <name> expr <expression>");
    Variable variable = parse_variable_from_tail(tokens, 2, "variable");

    for (auto& entry : params.physics.variables) {
        if (entry.first == tokens[1]) {
            entry.second = variable;
            return;
        }
    }
    params.physics.variables.push_back({tokens[1], variable});
}
// ---------------------------------------------------------------------- //
void ParamParser::parse_l_coeff_command(const std::vector<std::string>& tokens) {
    if (tokens.size() < 5) throw std::runtime_error("L_coeff syntax: L_coeff <i> <j> const <value> | L_coeff <i> <j> expr <expression>");
    const int i = std::stoi(tokens[1]);
    const int j = std::stoi(tokens[2]);
    
    check_component_index(i, "L_coeff");
    check_component_index(j, "L_coeff");
    params.physics.mobility(i, j) = parse_variable_from_tail(tokens, 3, "L_coeff");
}
// ---------------------------------------------------------------------- //
void ParamParser::parse_fix_command(const std::vector<std::string>& tokens) {
    if (tokens.size() < 5) {
        throw std::runtime_error("fix syntax: fix <ID> all <type> <on|off> ...");
    }

    if (!params.commands.empty()) {
        throw std::runtime_error("fix must be specified before run or measure commands");
    }
    if (!params.initial.density_commands.empty() || !params.initial.momentum_commands.empty()) {
        throw std::runtime_error("fix must be specified before initial condition commands");
    }

    const std::string& group = tokens[2];
    const std::string& type = tokens[3];

    if (group != "all") {
        throw std::runtime_error("only group 'all' is supported for now");
    }

    const FixSpec& spec = find_fix_spec(type);
    const bool enabled = parse_on_off(tokens[4], "fix");

    if (spec.arg_kind == FixArgKind::None) {
        if (tokens.size() != 5) {
            throw std::runtime_error(std::string("fix ") + spec.name + " syntax: fix <ID> all " + spec.name + " <on|off>");
        }
    }
    else if (spec.arg_kind == FixArgKind::Seed) {
        if (enabled) {
            if (tokens.size() != 7 || tokens[5] != "seed") {
                throw std::runtime_error("fix noise syntax: fix <ID> all noise on seed <integer>");
            }
            params.fix.noise.seed = std::stoi(tokens[6]);
        } else if (tokens.size() != 5) {
            throw std::runtime_error("fix noise off syntax: fix <ID> all noise off");
        }
    }
    else if (spec.arg_kind == FixArgKind::Rate) {
        if (enabled) {
            if (tokens.size() != 7 || tokens[5] != "rate") {
                throw std::runtime_error("fix shear syntax: fix <ID> all shear on rate <value>");
            }
            params.fix.shear.rate = std::stod(tokens[6]);
        } else if (tokens.size() != 5) {
            throw std::runtime_error("fix shear off syntax: fix <ID> all shear off");
        }
    }

    params.fix.set(spec.flag, enabled);
}
// ---------------------------------------------------------------------- //
void ParamParser::parse_set_command(const std::vector<std::string>& tokens) {
    if (tokens.size() < 3) throw std::runtime_error("set syntax: set density <component|all> <type> <args...> | set momentum <type> <args...>");

    if (tokens[1] == "density") {
        if (tokens.size() < 4) {
            throw std::runtime_error("set density syntax: set density <component|all> <type> <args...>");
        }

        const std::string& target = tokens[2];
        const std::string& type = tokens[3];

        InitialConditionArgs args;
        for (std::size_t i = 4; i < tokens.size(); i += 2) {
            if (i + 1 >= tokens.size()) {
                throw std::runtime_error("set density arguments must be key-value pairs");
            }
            args.entries.emplace_back(tokens[i], tokens[i + 1]);
        }

        const DensityInitialConditionStyle& style = initial_condition_registry.get_density(type);

        if (target == "all") {
            if (params.physics.num_components <= 0) {
                throw std::runtime_error("set density all requires components to be specified first");
            }

            for (int component = 0; component < params.physics.num_components; ++component) {
                std::shared_ptr<DensityInitialConditionCommandBase> command = style.parse_command(component, args, params);

                params.initial.density_commands.push_back(std::move(command));
            }

            return;
        }

        const int component = std::stoi(target);
        check_component_index(component, "set density");

        std::shared_ptr<DensityInitialConditionCommandBase> command = style.parse_command(component, args, params);

        params.initial.density_commands.push_back(std::move(command));
        return;
    }

    if (tokens[1] == "momentum") {
        if (tokens.size() < 4) {
            throw std::runtime_error("set momentum syntax: set momentum <component|x|y|0|1> <type> <args...>");
        }

        const std::string& target = tokens[2];
        const std::string& type = tokens[3];

        InitialConditionArgs args;
        for (std::size_t i = 4; i < tokens.size(); i += 2) {
            if (i + 1 >= tokens.size()) {
                throw std::runtime_error("set momentum arguments must be key-value pairs");
            }
            args.entries.emplace_back(tokens[i], tokens[i + 1]);
        }

        const MomentumInitialConditionStyle& style = initial_condition_registry.get_momentum(type);

        if (target == "all") {
            for (int component = 0; component < 2; ++component) {
                std::shared_ptr<MomentumInitialConditionCommandBase> command = style.parse_command(component, args, params);

                params.initial.momentum_commands.push_back(std::move(command));
            }

            return;
        }

        int component = -1;
        if (target == "x" || target == "0") {
            component = 0;
        } else if (target == "y" || target == "1") {
            component = 1;
        } else {
            throw std::runtime_error(
                "set momentum target must be all, x, y, 0, or 1: " + target
            );
        }
        check_momentum_component_index(component, "set momentum");

        std::shared_ptr<MomentumInitialConditionCommandBase> command = style.parse_command(component, args, params);

        params.initial.momentum_commands.push_back(std::move(command));
        return;
    }

    throw std::runtime_error("unknown set target: " + tokens[1]);

}
// ---------------------------------------------------------------------- //
void ParamParser::parse_thermo_command(const std::vector<std::string>& tokens) {
    if (tokens.size() < 2) throw std::runtime_error("thermo syntax: thermo off | thermo observe <on|off> progress <on|off> nevery <integer>");

    if (tokens[1] == "off") {
        if (tokens.size() != 2) {
            throw std::runtime_error("thermo off syntax: thermo off");
        }
        params.thermo.observe = false;
        params.thermo.progress = false;
        return;
    }

    for (std::size_t i = 1; i < tokens.size(); i += 2) {
        if (i + 1 >= tokens.size()) {
            throw std::runtime_error("thermo arguments must be key-value pairs");
        }

        const std::string& key = tokens[i];
        const std::string& value = tokens[i + 1];

        if (key == "observe") {
            params.thermo.observe = parse_on_off(value, "thermo observe");
        } else if (key == "progress") {
            params.thermo.progress = parse_on_off(value, "thermo progress");
        } else if (key == "nevery") {
            params.thermo.nevery = std::stoi(value);
            if (params.thermo.nevery <= 0) {
                throw std::runtime_error("thermo nevery must be positive");
            }
        } else {
            throw std::runtime_error("unknown thermo argument: " + key);
        }
    }
}
// ---------------------------------------------------------------------- //
void ParamParser::parse_restart_command(const std::vector<std::string>& tokens) {
    if (tokens.size() < 2) throw std::runtime_error("restart syntax: restart off | restart on file <filename>");

    if (tokens[1] == "off") {
        if (tokens.size() != 2) throw std::runtime_error("restart off syntax: restart off");
        params.restart_output.enabled = false;
        params.restart_output.file.clear();
        return;
    }

    if (tokens[1] == "on") {
        if (tokens.size() != 4 || tokens[2] != "file") throw std::runtime_error("restart on syntax: restart on file <filename>");
        params.restart_output.enabled = true;
        params.restart_output.file = tokens[3];
        return;
    }

    throw std::runtime_error("restart expects on|off");
}
// ---------------------------------------------------------------------- //
 void ParamParser::parse_measure_command(const std::vector<std::string>& tokens) {
    if (tokens.size() < 4) {
        throw std::runtime_error("measure syntax: measure <ID> <type> <on|off> <args...>");
    }

    const std::string& id = tokens[1];
    const std::string& type = tokens[2];
    const std::string& enabled_token = tokens[3];

    bool enabled = false;
    if (enabled_token == "on") {
        enabled = true;
    }
    else if (enabled_token == "off") {
        enabled = false;
    }
    else {
        throw std::runtime_error("measure: expected 'on' or 'off', got '" + enabled_token + "'");
    }

    MeasureArgs args;
    for (std::size_t i = 4; i < tokens.size(); i += 2) {
        if (i + 1 >= tokens.size()) {
            throw std::runtime_error("measure arguments must be key-value pairs.");
        }
        args.entries.emplace_back(tokens[i], tokens[i + 1]);
    }

    const MeasureStyle& style = measure_registry.get(type);
    std::shared_ptr<MeasureCommandBase> command = style.parse_command(id, enabled, args, params);

    Command cmd;
    cmd.type = Command::Type::Measure;
    cmd.measure = std::move(command);
    params.commands.push_back(std::move(cmd));

}
// ---------------------------------------------------------------------- //
void ParamParser::validate_configuration() const {
    if (params.physics.num_components <= 0) {
        throw std::runtime_error("components must be specified");
    }

    if (!params.physics.free_energy_entry.specified()) {
        throw std::runtime_error("free_energy must be specified");
    }

    if (!params.physics.eta_entry.specified()) {
        throw std::runtime_error("eta must be specified");
    }

    if (!params.physics.zeta_entry.specified()) {
        throw std::runtime_error("zeta must be specified");
    }

    for (int i = 0; i < params.physics.num_components; ++i) {
        for (int j = 0; j < params.physics.num_components; ++j) {
            if (!params.physics.mobility(i, j).specified()) {
                throw std::runtime_error("missing L_coeff " + std::to_string(i) + " " + std::to_string(j));
            }
        }
    }

}
// ---------------------------------------------------------------------- //
void ParamParser::check_component_index(int index, const std::string& command_name) const {
    if (params.physics.num_components <= 0) {
        throw std::runtime_error(command_name + " requires components to be specified first");
    }
    if (index < 0 || index >= params.physics.num_components) {
        throw std::runtime_error(
            command_name + ": component index out of range: " + std::to_string(index)
        );
    }
}
// ---------------------------------------------------------------------- //
void ParamParser::check_momentum_component_index(int index, const std::string& command_name) const {
    if (index < 0 || index >= 2) {
        throw std::runtime_error(
            command_name + ": momentum component index must be 0 or 1: " + std::to_string(index)
        );
    }
}

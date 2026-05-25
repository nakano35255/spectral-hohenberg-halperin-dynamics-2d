#include "param_parser.h"

#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

// ---------------------------------------------------------------------- //
namespace {
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
    else if (cmd == "order_parameters" || cmd == "order_parameter") {
        if (tokens.size() != 2) throw std::runtime_error("order_parameters syntax: order_parameters <num>");
        params.physics.num_order_parameters = std::stoi(tokens[1]);
        if (params.physics.num_order_parameters < 0) {
            throw std::runtime_error("order_parameters must be nonnegative.");
        }
    }
    else if (cmd == "boundary") {
        if (tokens.size() != 3) throw std::runtime_error("boundary syntax: boundary <x> <y>");
        if ((tokens[1] != "p" && tokens[1] != "periodic") || (tokens[2] != "p" && tokens[2] != "periodic")) {
            throw std::runtime_error("This solver only supports periodic (p) boundary.");
        }
    }
    else if (cmd == "grid") {
        if (tokens.size() != 3) throw std::runtime_error("grid syntax: grid <Nx> <Ny>");
        params.grid.active_num_grid[0] = std::stoi(tokens[1]);
        params.grid.active_num_grid[1] = std::stoi(tokens[2]);
        if (params.grid.active_num_grid[0] <= 0 || params.grid.active_num_grid[1] <= 0) {
            throw std::runtime_error("grid sizes must be positive");
        }
        if (params.grid.active_num_grid[0] % 2 != 0 || params.grid.active_num_grid[1] % 2 != 0) {
            throw std::runtime_error("grid must be even.");
        }
        params.grid.update_compute_grid();
    }
    else if (cmd == "length") {
        if (tokens.size() != 3) throw std::runtime_error("length syntax: length <Lx> <Ly>");
        params.grid.length[0] = std::stod(tokens[1]);
        params.grid.length[1] = std::stod(tokens[2]);
        if (params.grid.length[0] <= 0.0 || params.grid.length[1] <= 0.0) {
            throw std::runtime_error("domain lengths must be positive");
        }
    }
    else if (cmd == "dealias") {
        if (tokens.size() != 2) throw std::runtime_error("dealias syntax: dealias <none|three_halves|3/2|two|2>");
        const std::string& rule = tokens[1];

        if (rule == "none" || rule == "off") params.grid.dealias_rule = DealiasRule::None;
        else if (rule == "three_halves" || rule == "3/2") params.grid.dealias_rule = DealiasRule::ThreeHalves;
        else if (rule == "two" || rule == "2") params.grid.dealias_rule = DealiasRule::Two;
        else throw std::runtime_error("unknown dealias rule: " + rule);

        params.grid.update_compute_grid();
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
    else if (cmd == "time_evolution") {
        if (tokens.size() != 2) throw std::runtime_error("time_evolution syntax: time_evolution <type>");
        params.runtime.time_evolution_type = tokens[1];
    }
    
    // --------------------------------------------------
    // 3. Physics
    // --------------------------------------------------
    else if (cmd == "model") {
        parse_model_command(tokens);
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
void ParamParser::parse_model_command(const std::vector<std::string>& tokens) {
    if (tokens.size() < 3) {
        throw std::runtime_error(
            "model syntax: model <thermo|transport> <type> [args...] | model <thermo|transport> coeff <args...>"
        );
    }

    const std::string& category = tokens[1];
    const std::string& action = tokens[2];

    if (category == "thermo") {
        if (action == "coeff") {
            if (!params.physics.thermo) {
                throw std::runtime_error("model thermo coeff requires model thermo <type> first");
            }

            ThermodynamicsArgs args;
            for (std::size_t i = 3; i < tokens.size(); i += 2) {
                if (i + 1 >= tokens.size()) {
                    throw std::runtime_error("model thermo coeff arguments must be key-value pairs");
                }
                args.entries.emplace_back(tokens[i], tokens[i + 1]);
            }

            const auto& style = thermodynamics_registry.get_thermo(params.physics.thermo->type);
            style.update_command(*params.physics.thermo, args, params);
            return;
        }

        const std::string& type = action;
        const auto& style = thermodynamics_registry.get_thermo(type);

        params.physics.thermo = style.create_default_command(params);

        ThermodynamicsArgs args;
        for (std::size_t i = 3; i < tokens.size(); i += 2) {
            if (i + 1 >= tokens.size()) {
                throw std::runtime_error("model thermo " + type + " arguments must be key-value pairs");
            }
            args.entries.emplace_back(tokens[i], tokens[i + 1]);
        }

        style.update_command(*params.physics.thermo, args, params);
        return;
    }

    if (category == "free_energy") {
        if (action == "coeff") {
            if (!params.physics.free_energy) {
                throw std::runtime_error("model free_energy coeff requires model free_energy <type> first");
            }

            FreeEnergyArgs args;
            for (std::size_t i = 3; i < tokens.size(); i += 2) {
                if (i + 1 >= tokens.size()) {
                    throw std::runtime_error("model free_energy coeff arguments must be key-value pairs");
                }
                args.entries.emplace_back(tokens[i], tokens[i + 1]);
            }

            const auto& style = free_energy_registry.get_free_energy(params.physics.free_energy->type);
            style.update_command(*params.physics.free_energy, args, params);
            return;
        }

        const std::string& type = action;
        const auto& style = free_energy_registry.get_free_energy(type);

        params.physics.free_energy = style.create_default_command(params);

        FreeEnergyArgs args;
        for (std::size_t i = 3; i < tokens.size(); i += 2) {
            if (i + 1 >= tokens.size()) {
                throw std::runtime_error("model free_energy " + type + " arguments must be key-value pairs");
            }
            args.entries.emplace_back(tokens[i], tokens[i + 1]);
        }

        style.update_command(*params.physics.free_energy, args, params);
        return;
    }

    if (category == "transport") {
        if (action == "coeff") {
            if (!params.physics.transport) {
                throw std::runtime_error("model transport coeff requires model transport <type> first");
            }

            TransportCoefficientArgs args;
            for (std::size_t i = 3; i < tokens.size(); i += 2) {
                if (i + 1 >= tokens.size()) {
                    throw std::runtime_error("model transport coeff arguments must be key-value pairs");
                }
                args.entries.emplace_back(tokens[i], tokens[i + 1]);
            }

            const auto& style = transport_coefficient_registry.get_transport(params.physics.transport->type);
            style.update_command(*params.physics.transport, args, params);
            return;
        }

        const std::string& type = action;
        const auto& style =  transport_coefficient_registry.get_transport(type);

        params.physics.transport = style.create_default_command(params);

        TransportCoefficientArgs args;
        for (std::size_t i = 3; i < tokens.size(); i += 2) {
            if (i + 1 >= tokens.size()) {
                throw std::runtime_error("model transport " + type + " arguments must be key-value pairs");
            }
            args.entries.emplace_back(tokens[i], tokens[i + 1]);
        }

        style.update_command(*params.physics.transport, args, params);
        return;
    }

    throw std::runtime_error("unknown model category: " + category);
}
// ---------------------------------------------------------------------- //
void ParamParser::parse_fix_command(const std::vector<std::string>& tokens) {
    if (tokens.size() < 5) {
        throw std::runtime_error("fix syntax: fix <ID> all <type> <on|off> ...");
    }

    if (!params.commands.empty()) {
        throw std::runtime_error("fix must be specified before run or measure commands");
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
            if ((tokens.size() - 5) % 2 != 0) {
                throw std::runtime_error("fix noise syntax: fix <ID> all noise on [seed <integer>] [kBT <value>]");
            }

            for (std::size_t i = 5; i < tokens.size(); i += 2) {
                const std::string& key = tokens[i];
                const std::string& value = tokens[i + 1];

                if (key == "seed") {
                    params.fix.noise.seed = std::stoi(value);
                    continue;
                }

                if (key == "kBT" || key == "kbt" || key == "temperature") {
                    params.fix.noise.kBT = std::stod(value);
                    if (params.fix.noise.kBT < 0.0) {
                        throw std::runtime_error("fix noise requires nonnegative kBT.");
                    }
                    continue;
                }

                throw std::runtime_error("unknown fix noise argument: " + key);
            }
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
    if (tokens.size() < 3) {
        throw std::runtime_error(
            "set syntax: set density <type> <args...> | "
            "set momentum <x|y|all> <type> <args...> | "
            "set order_parameter <component|all> <type> <args...>"
        );
    }

    if (tokens[1] == "density" || tokens[1] == "rho") {
        if (tokens.size() < 3) {throw std::runtime_error("set density syntax: set density <type> <args...>");}

        const std::string& type = tokens[2];

        InitialConditionArgs args;
        for (std::size_t i = 3; i < tokens.size(); i += 2) {
            if (i + 1 >= tokens.size()) {
                throw std::runtime_error("set density arguments must be key-value pairs");
            }
            args.entries.emplace_back(tokens[i], tokens[i + 1]);
        }

        const DensityInitialConditionStyle& style = initial_condition_registry.get_density(type);
        std::shared_ptr<DensityInitialConditionCommandBase> command = style.parse_command(args, params);

        params.initial.density_commands.push_back(std::move(command));
        return;
    }

    if (tokens[1] == "momentum" || tokens[1] == "j") {
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

        int direction = -1;
        if (target == "x" || target == "0") {
            direction = 0;
        } else if (target == "y" || target == "1") {
            direction = 1;
        } else {
            throw std::runtime_error("set momentum target must be all, x, y, 0, or 1: " + target);
        }
        check_momentum_direction_index(direction, "set momentum");

        std::shared_ptr<MomentumInitialConditionCommandBase> command = style.parse_command(direction, args, params);

        params.initial.momentum_commands.push_back(std::move(command));
        return;
    }

    if (tokens[1] == "order_parameter" || tokens[1] == "psi") {
        if (tokens.size() < 4) {
            throw std::runtime_error("set order_parameter syntax: set order_parameter <component|all> <type> <args...>");
        }

        const std::string& target = tokens[2];
        const std::string& type = tokens[3];

        InitialConditionArgs args;
        for (std::size_t i = 4; i < tokens.size(); i += 2) {
            if (i + 1 >= tokens.size()) {
                throw std::runtime_error("set order_parameter arguments must be key-value pairs");
            }
            args.entries.emplace_back(tokens[i], tokens[i + 1]);
        }

        const OrderParameterInitialConditionStyle& style = initial_condition_registry.get_order_parameter(type);

        if (target == "all") {
            if (params.physics.num_order_parameters <= 0) {
                throw std::runtime_error("set order_parameter all requires at least one order parameter.");
            }

            for (int component = 0; component < params.physics.num_order_parameters; ++component) {
                std::shared_ptr<OrderParameterInitialConditionCommandBase> command = style.parse_command(component, args, params);
                params.initial.order_parameter_commands.push_back(std::move(command));
            }

            return;
        }

        const int component = std::stoi(target);
        check_order_parameter_index(component, "set order_parameter");

        std::shared_ptr<OrderParameterInitialConditionCommandBase> command = style.parse_command(component, args, params);

        params.initial.order_parameter_commands.push_back(std::move(command));
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
    if (params.physics.num_order_parameters < 0) {
        throw std::runtime_error("order_parameters must be nonnegative");
    }

    if (!params.physics.transport) {
        throw std::runtime_error("model transport must be specified");
    }

}
// ---------------------------------------------------------------------- //
void ParamParser::check_order_parameter_index(int index, const std::string& command_name) const {
    if (index < 0 || index >= params.physics.num_order_parameters) {
        throw std::runtime_error(
            command_name + ": order-parameter index out of range: " + std::to_string(index)
        );
    }
}
// ---------------------------------------------------------------------- //
void ParamParser::check_momentum_direction_index(int index, const std::string& command_name) const {
    if (index < 0 || index >= 2) {
        throw std::runtime_error(
            command_name + ": momentum direction index must be 0 or 1: " + std::to_string(index)
        );
    }
}

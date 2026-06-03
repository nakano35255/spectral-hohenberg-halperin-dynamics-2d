#include "param_parser.h"

#include <algorithm>
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
    int parse_xy_token(const std::string& value, const std::string& context) {
        if (value == "x" || value == "0") return 0;
        if (value == "y" || value == "1") return 1;
        throw std::runtime_error(context + " must be x|y|0|1.");
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
    else if (cmd == "fix") {
        parse_fix_command(tokens);
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
    // 6. 実行 (run)
    // --------------------------------------------------
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
// ---------------------------------------------------------------------- //
void ParamParser::parse_fix_command(const std::vector<std::string>& tokens) {
    if (tokens.size() < 5) {
        throw std::runtime_error("fix syntax: fix <ID> <target> <nonlinear|noise|force/sine|force/gradient> <on|off> ...");
    }

    if (!params.commands.empty()) {
        throw std::runtime_error("fix must be specified before run or measure commands");
    }

    const std::string& target_token = tokens[2];
    const std::string& type = tokens[3];
    const bool enabled = parse_on_off(tokens[4], "fix " + type);

    std::string target;

    if (target_token == "momentum" || target_token == "j") {
        target = "momentum";
    }
    else if (target_token == "order_parameter" || target_token == "psi") {
        if (params.physics.num_order_parameters <= 0) {
            throw std::runtime_error("fix " + type + " " + target_token + " requires order_parameters > 0.");
        }
        target = "order_parameter";
    }
    else if (target_token == "all") {
        target = "all";
    }
    else {
        throw std::runtime_error("unknown fix target: " + target_token);
    }

    if (type == "nonlinear") {
        parse_fix_nonlinear_command(tokens, tokens[1], target, enabled);
        return;
    }

    if (type == "noise") {
        parse_fix_noise_command(tokens, tokens[1], target, enabled);
        return;
    }

    if (type == "force/sine") {
        parse_fix_sine_force_command(tokens, tokens[1], target, enabled);
        return;
    }

    if (type == "force/gradient") {
        parse_fix_gradient_force_command(tokens, tokens[1], target, enabled);
        return;
    }

    throw std::runtime_error("unknown fix style: " + type);
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
void ParamParser::parse_fix_nonlinear_command(const std::vector<std::string>& tokens, const std::string& /*id*/, const std::string& target, bool enabled) {
    if (tokens.size() != 5) {
        throw std::runtime_error("fix nonlinear syntax: fix <ID> <target> nonlinear <on|off>");
    }

    if (target == "momentum" || target == "all") {
        params.fix.momentum_advection = enabled;
    }

    if (target == "order_parameter" || target == "all") {
        params.fix.order_parameter_advection = enabled && (params.physics.num_order_parameters > 0);
    }
}
// ---------------------------------------------------------------------- //
void ParamParser::parse_fix_noise_command(const std::vector<std::string>& tokens, const std::string& /*id*/, const std::string& target, bool enabled) {
    if (!enabled && tokens.size() != 5) {
        throw std::runtime_error("fix noise off syntax: fix <ID> <target> noise off");
    }

    if (enabled) {
        if ((tokens.size() - 5) % 2 != 0) {
            throw std::runtime_error( "fix noise syntax: fix <ID> <target> noise on [seed <integer>] [kBT <value>]"
            );
        }

        for (std::size_t i = 5; i < tokens.size(); i += 2) {
            const std::string& key = tokens[i];
            const std::string& value = tokens[i + 1];

            if (key == "seed") {
                params.fix.noise.seed = std::stoi(value);
                continue;
            }

            if (key == "kBT" || key == "temperature") {
                params.fix.noise.kBT = std::stod(value);
                if (params.fix.noise.kBT < 0.0) {
                    throw std::runtime_error("fix noise requires nonnegative kBT.");
                }
                continue;
            }

            if (key == "chi" || key == "order_parameter_chi" || key == "psi_chi") {
                if (target == "momentum") {
                    throw std::runtime_error("fix momentum noise does not use chi.");
                }
                if (params.physics.num_order_parameters <= 0) {
                    throw std::runtime_error("fix noise chi requires order_parameters > 0.");
                }

                params.fix.noise.order_parameter_noise_chi = std::stod(value);
                if (params.fix.noise.order_parameter_noise_chi < 0.0) {
                    throw std::runtime_error("fix noise requires nonnegative chi.");
                }
                continue;
            }

            throw std::runtime_error("unknown fix noise argument: " + key);
        }
    }

    if (target == "momentum" || target == "all") {
        params.fix.noise.momentum_enabled = enabled;
    }
    if (target == "order_parameter" || target == "all") {
        params.fix.noise.order_parameter_enabled = enabled && (params.physics.num_order_parameters > 0);
    }

    params.fix.noise.enabled = params.fix.noise.momentum_enabled || params.fix.noise.order_parameter_enabled;
}
// ---------------------------------------------------------------------- //
void ParamParser::parse_fix_sine_force_command(const std::vector<std::string>& tokens, const std::string& id, const std::string& target, bool enabled) {
    auto& forces = params.fix.sine_forces;
    forces.erase(
        std::remove_if(forces.begin(), forces.end(),
            [&](const SineForceFixConfig& cfg) { return cfg.id == id; }),
        forces.end()
    );

    if (!enabled) {
        if (tokens.size() != 5) {
            throw std::runtime_error("fix force/sine off syntax: fix <ID> <target> force/sine off");
        }
        return;
    }

    if (target == "all") {
        throw std::runtime_error("fix force/sine does not support target all.");
    }
    if ((tokens.size() - 5) % 2 != 0) {
        throw std::runtime_error("fix force/sine arguments must be key-value pairs.");
    }

    SineForceFixConfig cfg;
    cfg.id = id;
    cfg.enabled = true;
    cfg.momentum_enabled = (target == "momentum");
    cfg.order_parameter_enabled = (target == "order_parameter");

    bool has_component = false;
    bool has_axis = false;
    bool has_nk = false;
    bool has_amplitude = false;

    for (std::size_t i = 5; i < tokens.size(); i += 2) {
        const std::string& key = tokens[i];
        const std::string& value = tokens[i + 1];

        if (key == "component") {
            if (cfg.momentum_enabled) {
                cfg.component = parse_xy_token(value, "fix momentum force/sine component");
            } else {
                cfg.component = std::stoi(value);
                check_order_parameter_index(cfg.component, "fix order_parameter force/sine");
            }
            has_component = true;
        } else if (key == "axis") {
            cfg.axis = parse_xy_token(value, "fix force/sine axis");
            has_axis = true;
        } else if (key == "nk") {
            cfg.nk = std::stoi(value);
            has_nk = true;
        } else if (key == "amplitude") {
            cfg.amplitude = std::stod(value);
            has_amplitude = true;
        } else {
            throw std::runtime_error("unknown fix force/sine argument: " + key);
        }
    }

    if (!has_component || !has_axis || !has_nk || !has_amplitude) {
        throw std::runtime_error("fix force/sine requires component, axis, nk, and amplitude.");
    }
    if (cfg.nk <= 0) {
        throw std::runtime_error("fix force/sine requires nk > 0.");
    }

    const int active_size = (cfg.axis == 0) ? params.grid.active_num_grid[0] : params.grid.active_num_grid[1];
    if (cfg.nk >= active_size / 2) {
        throw std::runtime_error("fix force/sine requires nk < active_N_axis/2.");
    }

    forces.push_back(cfg);
}
// ---------------------------------------------------------------------- //
void ParamParser::parse_fix_gradient_force_command(const std::vector<std::string>& tokens, const std::string& id, const std::string& target, bool enabled) {
    auto& forces = params.fix.gradient_forces;
    forces.erase(
        std::remove_if(forces.begin(), forces.end(),
            [&](const GradientForceFixConfig& cfg) { return cfg.id == id; }),
        forces.end()
    );

    if (!enabled) {
        if (tokens.size() != 5) {
            throw std::runtime_error("fix force/gradient off syntax: fix <ID> <target> force/gradient off");
        }
        return;
    }

    if (target != "order_parameter") {
        throw std::runtime_error("fix force/gradient only supports target order_parameter.");
    }
    if ((tokens.size() - 5) % 2 != 0) {
        throw std::runtime_error("fix force/gradient arguments must be key-value pairs.");
    }

    GradientForceFixConfig cfg;
    cfg.id = id;
    cfg.enabled = true;

    bool has_component = false;
    bool has_direction = false;
    bool has_amplitude = false;

    for (std::size_t i = 5; i < tokens.size(); i += 2) {
        const std::string& key = tokens[i];
        const std::string& value = tokens[i + 1];

        if (key == "component") {
            cfg.component = std::stoi(value);
            check_order_parameter_index(cfg.component, "fix order_parameter force/gradient");
            has_component = true;
        } else if (key == "direction") {
            cfg.direction = parse_xy_token(value, "fix force/gradient direction");
            has_direction = true;
        } else if (key == "amplitude") {
            cfg.amplitude = std::stod(value);
            has_amplitude = true;
        } else {
            throw std::runtime_error("unknown fix force/gradient argument: " + key);
        }
    }

    if (!has_component || !has_direction || !has_amplitude) {
        throw std::runtime_error("fix force/gradient requires component, direction, and amplitude.");
    }

    forces.push_back(cfg);
}
// ---------------------------------------------------------------------- //
void ParamParser::validate_configuration() const {
    if (params.physics.num_order_parameters < 0) {
        throw std::runtime_error("order_parameters must be nonnegative");
    }

    if (!params.physics.transport) {
        throw std::runtime_error("model transport must be specified");
    }

    const bool quiescent = params.runtime.time_evolution_type == "euler/quiescent" || params.runtime.time_evolution_type == "srk3/quiescent";

    if (params.fix.order_parameter_advection && params.physics.num_order_parameters == 0) {
        throw std::runtime_error("fix nonlinear order_parameter requires at least one order parameter.");
    }

    if (quiescent && params.fix.momentum_advection) {
        throw std::runtime_error("quiescent time evolution cannot use momentum nonlinear advection.");
    }

    if (quiescent && params.fix.order_parameter_advection) {
        throw std::runtime_error("quiescent time evolution cannot use order_parameter nonlinear advection.");
    }

    for (const auto& force : params.fix.sine_forces) {
        if (quiescent && force.momentum_enabled) {
            throw std::runtime_error("quiescent time evolution cannot use momentum force/sine.");
        }
    }

    for (const auto& force : params.fix.gradient_forces) {
        if (quiescent && force.enabled) {
            throw std::runtime_error("quiescent time evolution cannot use force/gradient.");
        }
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

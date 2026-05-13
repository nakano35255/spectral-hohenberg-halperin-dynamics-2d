#include "param_parser.h"

#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace {

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

Variable::Type parse_variable_type(const std::string& token, const std::string& context) {
    if (token == "const") {
        return Variable::Type::Constant;
    }
    if (token == "expr") {
        return Variable::Type::Expression;
    }
    throw std::runtime_error(context + " expects const|expr");
}

Variable parse_variable_from_tail(const std::vector<std::string>& tokens,
                                  std::size_t type_index,
                                  const std::string& context) {
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

struct ParsedMeasureCommand : MeasureCommandBase {
    std::vector<std::pair<std::string, std::string>> args;

    void print(std::ostream& os) const override {
        os << id << " " << type << " " << (enabled ? "on" : "off");
        for (const auto& arg : args) {
            os << " " << arg.first << " " << arg.second;
        }
    }
};

} // namespace

ParamParser::ParamParser(Params& params_in)
    : params(params_in) {}

void ParamParser::parse_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open input file: " + filename);
    }

    std::string line;
    int line_number = 0;
    while (std::getline(file, line)) {
        ++line_number;
        try {
            execute_command(tokenize(line));
        } catch (const std::exception& e) {
            throw std::runtime_error(
                filename + ":" + std::to_string(line_number) + ": " + e.what()
            );
        }
    }

    validate_configuration();
}

std::vector<std::string> ParamParser::tokenize(const std::string& line) const {
    std::vector<std::string> tokens;
    std::string current;
    bool in_quote = false;

    for (char c : line) {
        if (!in_quote && c == '#') {
            break;
        }

        if (c == '"') {
            if (in_quote) {
                tokens.push_back(current);
                current.clear();
                in_quote = false;
            } else {
                if (!current.empty()) {
                    tokens.push_back(current);
                    current.clear();
                }
                in_quote = true;
            }
            continue;
        }

        if (!in_quote && std::isspace(static_cast<unsigned char>(c))) {
            if (!current.empty()) {
                tokens.push_back(current);
                current.clear();
            }
            continue;
        }

        current.push_back(c);
    }

    if (in_quote) {
        throw std::runtime_error("unterminated quoted string");
    }
    if (!current.empty()) {
        tokens.push_back(current);
    }
    return tokens;
}

void ParamParser::execute_command(const std::vector<std::string>& tokens) {
    if (tokens.empty()) {
        return;
    }

    const std::string& cmd = tokens[0];

    if (cmd == "dimension") {
        if (tokens.size() != 2) {
            throw std::runtime_error("dimension syntax: dimension 2");
        }
        params.grid.dimension = std::stoi(tokens[1]);
        if (params.grid.dimension != 2) {
            throw std::runtime_error("this solver supports only dimension 2");
        }
    } else if (cmd == "components") {
        if (tokens.size() != 2) {
            throw std::runtime_error("components syntax: components <N>");
        }
        params.physics.num_components = std::stoi(tokens[1]);
        if (params.physics.num_components <= 0) {
            throw std::runtime_error("components must be positive");
        }
        params.physics.resize(params.physics.num_components);
        params.initial.resize_densities(params.physics.num_components);
    } else if (cmd == "boundary") {
        if (tokens.size() != 3) {
            throw std::runtime_error("boundary syntax: boundary p p");
        }
        params.grid.boundary[0] = tokens[1];
        params.grid.boundary[1] = tokens[2];
        if ((params.grid.boundary[0] != "p" && params.grid.boundary[0] != "periodic") ||
            (params.grid.boundary[1] != "p" && params.grid.boundary[1] != "periodic")) {
            throw std::runtime_error("this solver supports only periodic boundaries");
        }
    } else if (cmd == "grid") {
        if (tokens.size() != 3) {
            throw std::runtime_error("grid syntax: grid <Nx> <Ny>");
        }
        params.grid.num_grid[0] = std::stoi(tokens[1]);
        params.grid.num_grid[1] = std::stoi(tokens[2]);
        if (params.grid.num_grid[0] <= 0 || params.grid.num_grid[1] <= 0) {
            throw std::runtime_error("grid sizes must be positive");
        }
    } else if (cmd == "length") {
        if (tokens.size() != 3) {
            throw std::runtime_error("length syntax: length <Lx> <Ly>");
        }
        params.grid.length[0] = std::stod(tokens[1]);
        params.grid.length[1] = std::stod(tokens[2]);
        if (params.grid.length[0] <= 0.0 || params.grid.length[1] <= 0.0) {
            throw std::runtime_error("domain lengths must be positive");
        }
    } else if (cmd == "integrator") {
        if (tokens.size() != 2) {
            throw std::runtime_error("integrator syntax: integrator <name>");
        }
        params.runtime.integrator = tokens[1];
    } else if (cmd == "timestep") {
        if (tokens.size() != 2) {
            throw std::runtime_error("timestep syntax: timestep <dt>");
        }
        params.runtime.dt = std::stod(tokens[1]);
        if (params.runtime.dt <= 0.0) {
            throw std::runtime_error("timestep must be positive");
        }
    } else if (cmd == "variable") {
        parse_variable_command(tokens);
    } else if (cmd == "free_energy") {
        params.physics.free_energy = parse_variable_from_tail(tokens, 1, "free_energy");
    } else if (cmd == "L_coeff") {
        parse_l_coeff_command(tokens);
    } else if (cmd == "eta") {
        params.physics.eta = parse_variable_from_tail(tokens, 1, "eta");
    } else if (cmd == "zeta") {
        params.physics.zeta = parse_variable_from_tail(tokens, 1, "zeta");
    } else if (cmd == "fix") {
        parse_fix_command(tokens);
    } else if (cmd == "set") {
        parse_set_command(tokens);
    } else if (cmd == "thermo") {
        parse_thermo_command(tokens);
    } else if (cmd == "measure") {
        parse_measure_command(tokens);
    } else if (cmd == "restart") {
        parse_restart_command(tokens);
    } else if (cmd == "run") {
        if (tokens.size() != 2) {
            throw std::runtime_error("run syntax: run <steps>");
        }
        const int steps = std::stoi(tokens[1]);
        if (steps <= 0) {
            throw std::runtime_error("run steps must be positive");
        }
        Command command;
        command.type = Command::Type::Run;
        command.run.steps = steps;
        params.commands.push_back(command);
    } else {
        throw std::runtime_error("Unknown command: " + cmd);
    }
}

void ParamParser::parse_variable_command(const std::vector<std::string>& tokens) {
    if (tokens.size() < 4) {
        throw std::runtime_error("variable syntax: variable <name> const <value> | variable <name> expr <expression>");
    }

    Variable variable = parse_variable_from_tail(tokens, 2, "variable");

    for (auto& entry : params.physics.variables) {
        if (entry.first == tokens[1]) {
            entry.second = variable;
            return;
        }
    }
    params.physics.variables.push_back({tokens[1], variable});
}

void ParamParser::parse_l_coeff_command(const std::vector<std::string>& tokens) {
    if (tokens.size() < 5) {
        throw std::runtime_error("L_coeff syntax: L_coeff <i> <j> const <value> | L_coeff <i> <j> expr <expression>");
    }

    const int i = std::stoi(tokens[1]);
    const int j = std::stoi(tokens[2]);
    check_component_index(i, "L_coeff");
    check_component_index(j, "L_coeff");
    params.physics.mobility(i, j) = parse_variable_from_tail(tokens, 3, "L_coeff");
}

void ParamParser::parse_fix_command(const std::vector<std::string>& tokens) {
    if (tokens.size() < 5) {
        throw std::runtime_error("fix syntax: fix <ID> all noise <on|off> ...");
    }

    const std::string& group = tokens[2];
    const std::string& style = tokens[3];
    const std::string& state = tokens[4];

    if (group != "all") {
        throw std::runtime_error("only group 'all' is supported for now");
    }
    if (style != "noise") {
        throw std::runtime_error("unknown fix style: " + style);
    }

    if (state == "on") {
        if (tokens.size() != 7 || tokens[5] != "seed") {
            throw std::runtime_error("fix noise syntax: fix <ID> all noise on seed <integer>");
        }
        params.physics.has_noise = true;
        params.physics.noise.enabled = true;
        params.physics.noise.seed = std::stoi(tokens[6]);
    } else if (state == "off") {
        if (tokens.size() != 5) {
            throw std::runtime_error("fix noise off syntax: fix <ID> all noise off");
        }
        params.physics.has_noise = false;
        params.physics.noise.enabled = false;
    } else {
        throw std::runtime_error("fix noise expects on|off");
    }
}

void ParamParser::parse_set_command(const std::vector<std::string>& tokens) {
    if (tokens.size() < 3) {
        throw std::runtime_error("set syntax: set density <component> <type> <args...> | set momentum <type> <args...>");
    }

    if (tokens[1] == "density") {
        if (tokens.size() < 4) {
            throw std::runtime_error("set density syntax: set density <component> <type> <args...>");
        }
        const int component = std::stoi(tokens[2]);
        check_component_index(component, "set density");

        auto& density = params.initial.densities.at(static_cast<std::size_t>(component));
        density.type = tokens[3];
        density.args.clear();
        for (std::size_t i = 4; i < tokens.size(); i += 2) {
            if (i + 1 >= tokens.size()) {
                throw std::runtime_error("set density arguments must be key-value pairs");
            }
            density.args.push_back({tokens[i], tokens[i + 1]});
        }
        return;
    }

    if (tokens[1] == "momentum") {
        params.initial.momentum.type = tokens[2];
        params.initial.momentum.args.clear();
        for (std::size_t i = 3; i < tokens.size(); i += 2) {
            if (i + 1 >= tokens.size()) {
                throw std::runtime_error("set momentum arguments must be key-value pairs");
            }
            params.initial.momentum.args.push_back({tokens[i], tokens[i + 1]});
        }
        return;
    }

    if (tokens[1] == "all" && tokens[2] == "restart") {
        if (tokens.size() != 5 || tokens[3] != "file") {
            throw std::runtime_error("set all restart syntax: set all restart file <filename>");
        }
        params.initial.use_restart = true;
        params.initial.use_stationary = false;
        params.initial.restart_file = tokens[4];
        return;
    }

    if (tokens[1] == "all" && tokens[2] == "stationary") {
        if (tokens.size() != 3 && tokens.size() != 5) {
            throw std::runtime_error("set all stationary syntax: set all stationary [seed <integer>]");
        }
        params.initial.use_restart = false;
        params.initial.use_stationary = true;
        params.initial.stationary_seed = 12345;
        if (tokens.size() == 5) {
            if (tokens[3] != "seed") {
                throw std::runtime_error("set all stationary syntax: set all stationary [seed <integer>]");
            }
            params.initial.stationary_seed = std::stoi(tokens[4]);
        }
        return;
    }

    throw std::runtime_error("unknown set target: " + tokens[1]);
}

void ParamParser::parse_measure_command(const std::vector<std::string>& tokens) {
    if (tokens.size() < 4) {
        throw std::runtime_error("measure syntax: measure <ID> <style> <on|off> ...");
    }

    auto measure = std::make_shared<ParsedMeasureCommand>();
    measure->id = tokens[1];
    measure->type = tokens[2];

    if (tokens[3] == "on") {
        measure->enabled = true;
    } else if (tokens[3] == "off") {
        measure->enabled = false;
    } else {
        throw std::runtime_error("measure expects on|off");
    }

    for (std::size_t i = 4; i < tokens.size(); i += 2) {
        if (i + 1 >= tokens.size()) {
            throw std::runtime_error("measure arguments must be key-value pairs");
        }
        measure->args.push_back({tokens[i], tokens[i + 1]});
    }

    Command command;
    command.type = Command::Type::Measure;
    command.measure = measure;
    params.commands.push_back(command);
}

void ParamParser::parse_thermo_command(const std::vector<std::string>& tokens) {
    if (tokens.size() < 2) {
        throw std::runtime_error("thermo syntax: thermo off | thermo observe <on|off> progress <on|off> nevery <integer>");
    }

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

void ParamParser::parse_restart_command(const std::vector<std::string>& tokens) {
    if (tokens.size() < 2) {
        throw std::runtime_error("restart syntax: restart off | restart on file <filename>");
    }

    if (tokens[1] == "off") {
        if (tokens.size() != 2) {
            throw std::runtime_error("restart off syntax: restart off");
        }
        params.restart.enabled = false;
        params.restart.file.clear();
        return;
    }

    if (tokens[1] == "on") {
        if (tokens.size() != 4 || tokens[2] != "file") {
            throw std::runtime_error("restart on syntax: restart on file <filename>");
        }
        params.restart.enabled = true;
        params.restart.file = tokens[3];
        return;
    }

    throw std::runtime_error("restart expects on|off");
}

void ParamParser::validate_configuration() const {
    if (params.grid.dimension != 2) {
        throw std::runtime_error("dimension 2 must be specified");
    }
    if (params.physics.num_components <= 0) {
        throw std::runtime_error("components must be specified before model coefficients");
    }
    if (!params.physics.free_energy.specified()) {
        throw std::runtime_error("free_energy must be specified");
    }
    if (!params.physics.eta.specified()) {
        throw std::runtime_error("eta must be specified");
    }
    if (!params.physics.zeta.specified()) {
        throw std::runtime_error("zeta must be specified");
    }

    for (int i = 0; i < params.physics.num_components; ++i) {
        for (int j = 0; j < params.physics.num_components; ++j) {
            if (!params.physics.mobility(i, j).specified()) {
                throw std::runtime_error(
                    "missing L_coeff " + std::to_string(i) + " " + std::to_string(j)
                );
            }
        }
    }

    if (params.initial.densities.size() != static_cast<std::size_t>(params.physics.num_components)) {
        throw std::runtime_error("initial density storage is not resized");
    }
    if (params.commands.empty()) {
        throw std::runtime_error("input script must contain at least one command such as run");
    }
}

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

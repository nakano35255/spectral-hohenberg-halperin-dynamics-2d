#include "simulationinfo.h"

#include <ostream>
#include <stdexcept>

namespace {

const char* variable_type_name(Variable::Type type) {
    switch (type) {
    case Variable::Type::Constant:
        return "const";
    case Variable::Type::Expression:
        return "expr";
    }
    return "unknown";
}

void print_variable(std::ostream& os, const std::string& label, const Variable& variable) {
    os << "  " << label << ": "
       << variable_type_name(variable.type) << " " << variable.text << '\n';
}

const char* command_type_name(Command::Type type) {
    switch (type) {
    case Command::Type::Run:
        return "run";
    case Command::Type::Measure:
        return "measure";
    }
    return "unknown";
}

} // namespace

void RuntimeConfig::print_config(std::ostream& os) const {
    os << "[Runtime]\n";
    os << "  integrator: " << integrator << '\n';
    os << "  timestep: " << dt << '\n';
}

void PhysicsConfig::resize(int n) {
    num_components = n;
    mobility_entries.assign(static_cast<std::size_t>(n * n), Variable());
}

Variable& PhysicsConfig::mobility(int i, int j) {
    return mobility_entries.at(static_cast<std::size_t>(i * num_components + j));
}

const Variable& PhysicsConfig::mobility(int i, int j) const {
    return mobility_entries.at(static_cast<std::size_t>(i * num_components + j));
}

void PhysicsConfig::print_config(std::ostream& os) const {
    os << "[Physics]\n";
    os << "  components: " << num_components << '\n';
    os << "  barodiffusion: " << (has_barodiffusion ? "on" : "off") << '\n';
    os << "  noise: " << (has_noise ? "on" : "off");
    if (has_noise) {
        os << " seed " << noise.seed;
    }
    os << '\n';

    os << "  variables:\n";
    for (const auto& variable : variables) {
        os << "    " << variable.first << " "
           << variable_type_name(variable.second.type) << " "
           << variable.second.text << '\n';
    }

    print_variable(os, "free_energy", free_energy);

    os << "  L_coeff:\n";
    for (int i = 0; i < num_components; ++i) {
        for (int j = 0; j < num_components; ++j) {
            const auto& value = mobility(i, j);
            os << "    " << i << " " << j << " "
               << variable_type_name(value.type) << " " << value.text << '\n';
        }
    }

    print_variable(os, "eta", eta);
    print_variable(os, "zeta", zeta);
}

void InitialConditionConfig::print_config(std::ostream& os, const int num_components) const {
    os << "[Initial]\n";
    for (int i = 0; i < num_components; ++i) {
        os << "  density " << i << ": " << densities.at(static_cast<std::size_t>(i)).type;
        for (const auto& arg : densities.at(static_cast<std::size_t>(i)).args) {
            os << " " << arg.first << " " << arg.second;
        }
        os << '\n';
    }

    os << "  momentum: " << momentum.type;
    for (const auto& arg : momentum.args) {
        os << " " << arg.first << " " << arg.second;
    }
    os << '\n';

    os << "  restart input: " << (use_restart ? "on" : "off");
    if (use_restart) {
        os << " file " << restart_file;
    }
    os << '\n';

    os << "  stationary: " << (use_stationary ? "on" : "off");
    if (use_stationary) {
        os << " seed " << stationary_seed;
    }
    os << '\n';
}

void ThermoConfig::print_config(std::ostream& os) const {
    os << "[Thermo]\n";
    os << "  observe: " << (observe ? "on" : "off") << '\n';
    os << "  progress: " << (progress ? "on" : "off") << '\n';
    os << "  nevery: " << nevery << '\n';
}

void RestartConfig::print_config(std::ostream& os) const {
    os << "[Restart]\n";
    os << "  enabled: " << (enabled ? "on" : "off") << '\n';
    if (enabled) {
        os << "  file: " << file << '\n';
    }
}

void GridConfig::print_config(std::ostream& os) const {
    os << "[Grid]\n";
    os << "  dimension: " << dimension << '\n';
    os << "  grid: " << num_grid[0] << " " << num_grid[1] << '\n';
    os << "  length: " << length[0] << " " << length[1] << '\n';
    os << "  boundary: " << boundary[0] << " " << boundary[1] << '\n';
    os << "  spacing: " << dx() << " " << dy() << '\n';
}

int Params::total_run_steps() const {
    int total = 0;
    for (const auto& command : commands) {
        if (command.type == Command::Type::Run) {
            total += command.run.steps;
        }
    }
    return total;
}

void Params::write_summary(std::ostream& os) const {
    os << "Spectral Fluctuating Isothermal Fluid input summary\n";
    os << "===================================================\n";
    grid.print_config(os);
    runtime.print_config(os);
    physics.print_config(os);
    initial.print_config(os, physics.num_components);
    thermo.print_config(os);
    restart.print_config(os);

    os << "[Commands]\n";
    for (std::size_t i = 0; i < commands.size(); ++i) {
        const auto& command = commands[i];
        os << "  " << i << ": " << command_type_name(command.type);
        if (command.type == Command::Type::Run) {
            os << " " << command.run.steps;
        } else if (command.type == Command::Type::Measure && command.measure) {
            os << " ";
            command.measure->print(os);
        }
        os << '\n';
    }
    os << "  total run steps: " << total_run_steps() << '\n';
}

bool parse_on_off(const std::string& value, const std::string& context) {
    if (value == "on") {
        return true;
    }
    if (value == "off") {
        return false;
    }
    throw std::runtime_error(context + " expects on|off");
}

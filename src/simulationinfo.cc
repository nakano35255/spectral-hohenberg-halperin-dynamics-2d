#include "simulationinfo.h"

#include <ostream>
#include <stdexcept>
#include <iomanip>
#include <sstream>

// ---------------------------------------------------------------------- //
namespace {
    // ---------------------------------------------------------------------- //
    constexpr int LABEL_WIDTH = 25;
    constexpr int RULE_WIDTH = 70;

    void print_rule(std::ostream& os) {
        os << std::string(RULE_WIDTH, '-') << "\n";
    }
    // ---------------------------------------------------------------------- //
    void print_section(std::ostream& os, const std::string& title) {
        os << "\n[" << title << "]\n";
    }
    // ---------------------------------------------------------------------- //
    template <typename T>
    void print_entry(std::ostream& os, const std::string& label, const T& value) {
        os << "  "
           << std::left << std::setw(LABEL_WIDTH)
           << label << ": " << value << '\n';
    }
    // ---------------------------------------------------------------------- //
    template <typename T, typename U>
    std::string pair_value(const T& first, const U& second) {
        std::ostringstream ss;
        ss << first << " " << second;
        return ss.str();
    }
    // ---------------------------------------------------------------------- //
    const char* dealias_rule_name(DealiasRule rule) {
        if (rule == DealiasRule::None) {
            return "none";
        }
        if (rule == DealiasRule::ThreeHalves) {
            return "three_halves";
        }
        if (rule == DealiasRule::Two) {
            return "two";
        }
        return "unknown";
    }
    // ---------------------------------------------------------------------- //
} // namespace
// ---------------------------------------------------------------------- //
void GridConfig::print_config(std::ostream& os) const {
    print_section(os, "Grid Setup");
    print_entry(os, "Grid points", pair_value(active_num_grid[0], active_num_grid[1]));
    print_entry(os, "Dealias rule", dealias_rule_name(dealias_rule));
    print_entry(os, "Domain Length", pair_value(length[0], length[1]));
    print_entry(os, "dx, dy", pair_value(active_dx(), active_dy()));
    os << '\n';
    print_rule(os);
}
// ---------------------------------------------------------------------- //
void RuntimeConfig::print_config(std::ostream& os) const {
    print_section(os, "Runtime Setup");
    print_entry(os, "Time step (dt)", dt);
    print_entry(os, "Time evolution", time_evolution_type);
    os << '\n';
    print_rule(os);
}
// ---------------------------------------------------------------------- //
void PhysicsConfig::print_config(std::ostream& os) const {
    print_section(os, "Physics");
    print_entry(os, "order_parameters", num_order_parameters);

    if (thermo) {
        thermo->print(os);
    } else {
        print_entry(os, "thermodynamics", "none");
    }

    if (free_energy) {
        free_energy->print(os);
    } else {
        print_entry(os, "free_energy", "none");
    }

    if (transport) {
        transport->print(os);
    } else {
        print_entry(os, "transport", "none");
    }

    os << "\n";
    print_rule(os);
}
// ---------------------------------------------------------------------- //
void FixConfig::print_config(std::ostream& os) const {
    print_section(os, "Fix");

    print_entry(os, "noise", noise.enabled ? "ON" : "OFF");
    if (noise.enabled) {
        print_entry(os, "  momentum", noise.momentum_enabled ? "ON" : "OFF");
        print_entry(os, "  order_parameter", noise.order_parameter_enabled ? "ON" : "OFF");
        print_entry(os, "  Seed", noise.seed);
        print_entry(os, "  kBT", noise.kBT);
        if (noise.order_parameter_enabled) {
            print_entry(os, "  order_parameter_chi", noise.order_parameter_noise_chi);
        }
    }

    print_entry(os, "momentum_advection", momentum_advection ? "ON" : "OFF");
    print_entry(os, "order_parameter_advection", order_parameter_advection ? "ON" : "OFF");

    os << "\n";
    print_rule(os);
}
// ---------------------------------------------------------------------- //
void InitialConditionConfig::print_config(std::ostream& os) const {
    print_section(os, "Initial Condition");

    for (const auto& command : density_commands) {
        command->print(os);
    }

    for (const auto& command : momentum_commands) {
        command->print(os);
    }

    for (const auto& command : order_parameter_commands) {
        command->print(os);
    }

    os << "\n";
    print_rule(os);
}
// ---------------------------------------------------------------------- //
void ThermoConfig::print_config(std::ostream& os) const {
    print_section(os, "Thermo");
    print_entry(os, "Observe", observe ? "ON" : "OFF");
    print_entry(os, "Progress", progress ? "ON" : "OFF");
    if (observe || progress) {
        print_entry(os, "Nevery", nevery);
    }
    os << "\n";
    print_rule(os);
}
// ---------------------------------------------------------------------- //
void RestartOutputConfig::print_config(std::ostream& os) const {
    print_section(os, "Restart Output");
    print_entry(os, "Enabled", enabled ? "ON" : "OFF");
    if (enabled) {
        print_entry(os, "File", file);
    }
    os << "\n";
    print_rule(os);
}
// ---------------------------------------------------------------------- //
int Params::total_run_steps() const {
    int total = 0;
    for (const auto& command : commands) {
        if (command.type == Command::Type::Run) {
            total += command.run.steps;
        }
    }
    return total;
}
// ---------------------------------------------------------------------- //
void Params::write_summary(std::ostream& os) const {
    os << std::string(70, '=') << "\n";
    os << "Spectral Fluctuating Isothermal Fluid Solver (v0.0)\n";
    os << std::string(70, '=') << "\n";
    grid.print_config(os);
    runtime.print_config(os);
    physics.print_config(os);
    fix.print_config(os);
    initial.print_config(os);
    thermo.print_config(os);
    restart_output.print_config(os);


}

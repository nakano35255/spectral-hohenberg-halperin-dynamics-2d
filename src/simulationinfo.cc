#include "simulationinfo.h"

#include <ostream>
#include <stdexcept>
#include <iomanip>
#include <sstream>

// ---------------------------------------------------------------------- //
namespace {
    // ---------------------------------------------------------------------- //
    constexpr int LABEL_WIDTH = 25;
    constexpr int VARIABLE_LABEL_WIDTH = 23;
    constexpr int MATRIX_COL_WIDTH = 17;
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
    void print_variable(std::ostream& os, const std::string& label, const Variable& variable) {
        print_entry(os, label, variable.text);
    }
    // ---------------------------------------------------------------------- //
    void print_variable_definition(std::ostream& os,
                                   const std::string& name,
                                   const Variable& variable) {
        os << "    "
           << std::left << std::setw(VARIABLE_LABEL_WIDTH)
           << name << ": " << variable.text << '\n';
    }
    // ---------------------------------------------------------------------- //
} // namespace
// ---------------------------------------------------------------------- //
void GridConfig::print_config(std::ostream& os) const {
    print_section(os, "Grid Setup");
    print_entry(os, "Grid points", pair_value(num_grid[0], num_grid[1]));
    print_entry(os, "Length", pair_value(length[0], length[1]));
    print_entry(os, "dx, dy", pair_value(dx(), dy()));
    os << '\n';
    print_rule(os);
}
// ---------------------------------------------------------------------- //
void RuntimeConfig::print_config(std::ostream& os) const {
    print_section(os, "Runtime Setup");
    print_entry(os, "Time step (dt)", dt);
    print_entry(os, "Integrator", integrator_type);
    os << '\n';
    print_rule(os);
}
// ---------------------------------------------------------------------- //
Variable& PhysicsConfig::mobility(int i, int j) {
    return L_entries.at(static_cast<std::size_t>(i * num_components + j));
}
const Variable& PhysicsConfig::mobility(int i, int j) const {
    return L_entries.at(static_cast<std::size_t>(i * num_components + j));
}
void PhysicsConfig::print_config(std::ostream& os) const {
    print_section(os, "Physics");
    print_entry(os, "components", num_components);

    os << "  variables:\n";
    for (const auto& variable : variables) {
        print_variable_definition(os, variable.first, variable.second);
    }

    os << '\n';
    print_variable(os, "free_energy", free_energy_entry);

    os << "  "
       << std::left << std::setw(LABEL_WIDTH)
       << "L_coeff" << ":\n";

    const std::string matrix_indent =
        "  " + std::string(static_cast<std::size_t>(LABEL_WIDTH), ' ') + "  ";
    for (int i = 0; i < num_components; ++i) {
        os << matrix_indent;
        for (int j = 0; j < num_components; ++j) {
            const auto& value = mobility(i, j);
            os << std::left << std::setw(MATRIX_COL_WIDTH) << value.text;
        }
        os << '\n';
    }

    print_variable(os, "eta", eta_entry);
    print_variable(os, "zeta", zeta_entry);
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
const char* FixCommand::style_name(FixCommand::Style style) {
    switch (style) {
    case FixCommand::Style::Noise:
        return "noise";
    case FixCommand::Style::Shear:
        return "shear";
    case FixCommand::Style::Nonlinear:
        return "nonlinear";
    case FixCommand::Style::Barodiffusion:
        return "barodiffusion";
    }
    return "unknown";
}
// ---------------------------------------------------------------------- //
void FixCommand::print(std::ostream& os) const {
    const std::string state =
        std::string(FixCommand::style_name(style)) + " " + (enabled ? "ON" : "OFF");

    print_entry(os, "Fix " + id, state);
    print_entry(os, "  Group", group);

    if (enabled && style == FixCommand::Style::Noise) {
        print_entry(os, "  Seed", noise.seed);
    }

    if (enabled && style == FixCommand::Style::Shear) {
        print_entry(os, "  Rate", shear.rate);
        print_entry(os, "  Flow direction", shear.flow_direction);
    }
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
    initial.print_config(os);
    thermo.print_config(os);
    restart_output.print_config(os);


}

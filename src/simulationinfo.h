#ifndef SFI_SIMULATIONINFO_H
#define SFI_SIMULATIONINFO_H

#include <iosfwd>
#include <string>
#include <complex>
#include <utility>
#include <vector>
#include "measure_registry.h"

using Complex = std::complex<double>;

// ---------------------------------------------------------------------- //
struct Variable {
    enum class Type {
        Constant,
        Expression
    };
    Type type = Type::Expression;
    std::string text;

    Variable() = default;
    Variable(Type type_in, std::string value) : type(type_in), text(std::move(value)) {}
    
    bool specified() const {
        return !text.empty();
    }
};
// ---------------------------------------------------------------------- //
struct RuntimeConfig {
    double dt = 0.005;
    std::string integrator = "srk3";
    void print_config(std::ostream& os) const;
};
// ---------------------------------------------------------------------- //
struct NoiseSpec {
    bool enabled = false;
    int seed = 12345;
};
struct PhysicsConfig {
    int num_components = 0;
    bool has_barodiffusion = true;
    bool has_noise = false;
    std::vector<std::pair<std::string, Variable>> variables;
    Variable free_energy;
    std::vector<Variable> mobility_entries;
    Variable eta;
    Variable zeta;
    NoiseSpec noise;

    void resize(int n);
    Variable& mobility(int i, int j);
    const Variable& mobility(int i, int j) const;
    void print_config(std::ostream& os) const;
};
// ---------------------------------------------------------------------- //
struct DensityICSpec {
    std::string type = "zero";
    std::vector<std::pair<std::string, std::string>> args;
};
struct MomentumICSpec {
    std::string type = "zero";
    std::vector<std::pair<std::string, std::string>> args;
};
struct InitialConditionConfig {
    std::vector<DensityICSpec> densities;
    MomentumICSpec momentum;
    bool use_restart = false;
    std::string restart_file = "";
    bool use_stationary = false;
    int stationary_seed = 12345;
    void resize_densities(int num_components) {
        if (num_components <= 0) return;
        densities.assign(num_components, DensityICSpec());
    }
    void print_config(std::ostream& os, const int num_components) const;
};
// ---------------------------------------------------------------------- //
struct ThermoConfig {
    bool observe = false;
    bool progress = false;
    int nevery = 1000;
    void print_config(std::ostream& os) const;
};
// ---------------------------------------------------------------------- //
struct RestartConfig {
    bool enabled = false;
    std::string file = "";
    void print_config(std::ostream& os) const;
};
// ---------------------------------------------------------------------- //
struct RunCommand {
    int steps = 0;
};
struct Command {
    enum class Type {
        Run,
        Measure
    };
    Type type;
    RunCommand run;
    std::shared_ptr<MeasureCommandBase> measure;
};
// ---------------------------------------------------------------------- //
struct GridConfig {
    static constexpr double PI = 3.14159265358979323846;

    int dimension = 2;
    std::string boundary[2] = {"p", "p"};
    int num_grid[2] = {128, 128};
    double length[2] = {2.0 * PI, 2.0 * PI};

    double dx() const {
        return length[0] / static_cast<double>(num_grid[0]); 
    }
    double dy() const {
        return length[1] / static_cast<double>(num_grid[1]); 
    }
    void print_config(std::ostream& os) const;
};
// ---------------------------------------------------------------------- //

struct Params {
    GridConfig grid;
    RuntimeConfig runtime;
    PhysicsConfig physics;
    InitialConditionConfig initial;
    ThermoConfig thermo;
    RestartConfig restart;
    std::vector<Command> commands;

    int total_run_steps() const;
    void write_summary(std::ostream& os) const;
};
// ---------------------------------------------------------------------- //

bool parse_on_off(const std::string& value, const std::string& context);

#endif

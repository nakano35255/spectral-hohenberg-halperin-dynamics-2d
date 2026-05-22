#ifndef SFI_SIMULATIONINFO_H
#define SFI_SIMULATIONINFO_H

#include <iosfwd>
#include <string>
#include <complex>
#include <utility>
#include <vector>
#include "fix_flag.h"
#include "model_thermodynamics_registry.h"
#include "model_transport_coefficient_registry.h"
#include "measure_registry.h"
#include "initial_condition_registry.h"

using Complex = std::complex<double>;

// ---------------------------------------------------------------------- //
struct GridConfig {
    static constexpr double PI = 3.14159265358979323846;

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
struct RuntimeConfig {
    double dt = 0.005;
    std::string time_evolution_type = "srk3/compressible";
    void print_config(std::ostream& os) const;
};
// ---------------------------------------------------------------------- //
struct PhysicsConfig {
    int num_components = 0;
    std::shared_ptr<ThermodynamicsCommandBase> thermo;
    std::shared_ptr<TransportCoefficientCommandBase> transport;

    void print_config(std::ostream& os) const;
};
// ---------------------------------------------------------------------- //
struct NoiseSpec {
    int seed = 12345;
};
struct ShearSpec {
    double rate = 0.0;
    std::string flow_direction = "x";
};
struct FixConfig {
    std::uint32_t flags = 0;
    NoiseSpec noise;
    ShearSpec shear;

    void set(FixFlag flag, bool value) {
        if (value) {
            flags |= fix_bit(flag);
        } else {
            flags &= ~fix_bit(flag);
        }
    }

    void print_config(std::ostream& os) const;
};
// ---------------------------------------------------------------------- //
struct InitialConditionConfig {
    std::vector<std::shared_ptr<DensityInitialConditionCommandBase>> density_commands;
    std::vector<std::shared_ptr<MomentumInitialConditionCommandBase>> momentum_commands;

    void print_config(std::ostream& os) const;
};
// ---------------------------------------------------------------------- //
struct ThermoConfig {
    bool observe = false;
    bool progress = false;
    int nevery = 1000;
    void print_config(std::ostream& os) const;
};
// ---------------------------------------------------------------------- //
struct RestartOutputConfig {
    bool enabled = false;
    std::string file = "";
    void print_config(std::ostream& os) const;
};
// ---------------------------------------------------------------------- //
struct RunCommand {
    int steps = 0;
};
// ---------------------------------------------------------------------- //
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
struct Params {
    GridConfig grid;
    RuntimeConfig runtime;
    PhysicsConfig physics;
    InitialConditionConfig initial;
    FixConfig fix;
    ThermoConfig thermo;
    RestartOutputConfig restart_output;
    std::vector<Command> commands;

    int total_run_steps() const;
    void write_summary(std::ostream& os) const;
};
// ---------------------------------------------------------------------- //
#endif

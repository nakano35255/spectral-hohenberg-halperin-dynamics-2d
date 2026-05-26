#ifndef SFI_SIMULATIONINFO_H
#define SFI_SIMULATIONINFO_H

#include <iosfwd>
#include <string>
#include <complex>
#include <utility>
#include <vector>
#include <stdexcept>

#include "model_free_energy_registry.h"
#include "model_thermodynamics_registry.h"
#include "model_transport_coefficient_registry.h"
#include "measure_registry.h"
#include "initial_condition_registry.h"

using Complex = std::complex<double>;

// ---------------------------------------------------------------------- //
enum class PhysicalFieldPlan {
    None,
    PsiOnly,
    AllFields
};
enum class DealiasRule {
    None,
    ThreeHalves,
    Two
};
struct RHSEvaluationPlan {
    PhysicalFieldPlan physical_fields = PhysicalFieldPlan::None;
    DealiasRule dealias_rule = DealiasRule::None;
};
// ---------------------------------------------------------------------- //
struct GridConfig {
    static constexpr double PI = 3.14159265358979323846;

    int active_num_grid[2] = {128, 128};
    int num_grid[2] = {128, 128};
    double length[2] = {2.0 * PI, 2.0 * PI};

    DealiasRule dealias_rule = DealiasRule::None;

    static int compute_grid_size(int active_size, DealiasRule rule) {
        if (active_size <= 0) throw std::runtime_error("grid size must be positive.");
        if (rule == DealiasRule::None) return active_size;
        if (rule == DealiasRule::ThreeHalves) {
            if (active_size % 2 != 0) {
                throw std::runtime_error("3/2 dealiasing requires even active grid size.");
            }
            return 3 * active_size / 2;
        }
        if (rule == DealiasRule::Two) return 2 * active_size;
        throw std::runtime_error("Unknown dealias rule.");
    }
    void update_compute_grid() {
        num_grid[0] = compute_grid_size(active_num_grid[0], dealias_rule);
        num_grid[1] = compute_grid_size(active_num_grid[1], dealias_rule);
    }
    double dx() const {
        return length[0] / static_cast<double>(num_grid[0]); 
    }
    double dy() const {
        return length[1] / static_cast<double>(num_grid[1]); 
    }
    double active_dx() const {
        return length[0] / static_cast<double>(active_num_grid[0]);
    }
    double active_dy() const {
        return length[1] / static_cast<double>(active_num_grid[1]);
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
    int num_order_parameters = 0;
    std::shared_ptr<ThermodynamicsCommandBase> thermo;
    std::shared_ptr<FreeEnergyCommandBase> free_energy;
    std::shared_ptr<TransportCoefficientCommandBase> transport;

    void print_config(std::ostream& os) const;
};
// ---------------------------------------------------------------------- //
struct NoiseFixConfig {
    bool enabled = false;
    int seed = 12345;
    double kBT = 1.0;
};
struct FixConfig {
    NoiseFixConfig noise;
    bool momentum_advection = false;
    bool order_parameter_advection = false;

    void print_config(std::ostream& os) const;
};
// ---------------------------------------------------------------------- //
struct InitialConditionConfig {
    std::vector<std::shared_ptr<DensityInitialConditionCommandBase>> density_commands;
    std::vector<std::shared_ptr<MomentumInitialConditionCommandBase>> momentum_commands;
    std::vector<std::shared_ptr<OrderParameterInitialConditionCommandBase>> order_parameter_commands;

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
    FixConfig fix;
    InitialConditionConfig initial;
    ThermoConfig thermo;
    RestartOutputConfig restart_output;
    std::vector<Command> commands;

    int total_run_steps() const;
    void write_summary(std::ostream& os) const;
};
// ---------------------------------------------------------------------- //
#endif

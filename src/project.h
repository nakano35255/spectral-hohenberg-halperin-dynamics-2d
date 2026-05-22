#ifndef PROJECT_H
#define PROJECT_H

#include "simulationinfo.h"
#include "domain.h"
#include "state.h"
#include "buffer_physical_state.h"
#include "fourier_transform.h"
#include "fluid_constraint.h"
#include "monitor.h"
#include "time_evolution_mask.h"
#include "model_thermodynamics.h"
#include "model_thermodynamics_registry_builtin.h"
#include "model_transport_coefficient.h"
#include "model_transport_coefficient_registry_builtin.h"
#include "measure_registry_builtin.h"
#include "measure_manager.h"
#include "initial_condition_registry_builtin.h"
#include "initial_condition.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <memory>

class Project {
private:
    Params params;
    ThermodynamicsRegistry thermodynamics_registry;
    TransportCoefficientRegistry transport_coefficient_registry;
    MeasureRegistry measure_registry;
    InitialConditionRegistry initial_registry;

    std::unique_ptr<Thermodynamics> thermodynamics;
    std::unique_ptr<TransportCoefficient> transport_coefficient;

    Domain2D domain;
    State state;
    PhysicalStateBuffer buf_physical_state;
    FourierTransform2D fourier;
    FluidConstraint constraint;
    TimeEvolutionMask time_evolution_mask;

    MeasurementManager measurements;
    SimulationMonitor monitor;
    int step;
    int run_index;
    double time;

    static std::unique_ptr<Thermodynamics> create_thermodynamics(const Params& params, const ThermodynamicsRegistry& registry) {
        const std::string& thermo_type = params.physics.thermo->type;
        const ThermodynamicsStyle& thermo_style = registry.get_thermo(thermo_type);

        return thermo_style.create(params, params.physics.thermo);
    }

    static std::unique_ptr<TransportCoefficient> create_transport_coefficient(const Params& params, const TransportCoefficientRegistry& registry) {
        const std::string& transport_type = params.physics.transport->type;
        const TransportCoefficientStyle& transport_style = registry.get_transport(transport_type);

        return transport_style.create(params, params.physics.transport);
    }

    void apply_initial_conditions() {
        state.clear();

        for (const auto& command : params.initial.density_commands) {
            const DensityInitialConditionStyle& style = initial_registry.get_density(command->type);
            std::unique_ptr<DensityInitialCondition> initial_condition = style.create_initial_condition(params, command);
            initial_condition->apply(state, domain);
        }

        for (const auto& command : params.initial.momentum_commands) {
            const MomentumInitialConditionStyle& style = initial_registry.get_momentum(command->type);
            std::unique_ptr<MomentumInitialCondition> initial_condition = style.create_initial_condition(params, command);
            initial_condition->apply(state, domain);
        }

        constraint.initialize(state);
        constraint.apply(state);
    }

    void execute_command(const Command& command) {
        if (command.type == Command::Type::Run) {
            run_steps(command.run.steps);
        }
        else if (command.type == Command::Type::Measure) {
            monitor.print_measure_command(*command.measure);
            measurements.apply_measure_command(command.measure);
        }
    }

    void run_steps(int nsteps) {
        ++run_index;

        for (int local_step = 1; local_step <= nsteps; ++local_step) {
            ++step;
            time += params.runtime.dt;
            monitor.print_progress(run_index, local_step, nsteps, step, time);
            measurements.observe(state, buf_physical_state, fourier, domain, step, time);
        }

        monitor.finish_run_segment(run_index, step, time);
    }

public:
    explicit Project(const Params& p)
        : params(p),
          thermodynamics_registry(build_thermodynamics_registry()),
          transport_coefficient_registry(build_transport_coefficient_registry()),
          measure_registry(build_measure_registry()),
          initial_registry(build_initial_condition_registry()),
          thermodynamics(create_thermodynamics(params, thermodynamics_registry)),
          transport_coefficient(create_transport_coefficient(params, transport_coefficient_registry)),
          domain(params),
          state(domain, params),
          buf_physical_state(domain, params),
          fourier(domain),
          constraint(domain, params),
          time_evolution_mask(params),
          measurements(params, measure_registry),
          monitor(params, "output.log", domain.rank() == 0),
          step(0),
          run_index(0),
          time(0.0) {}

    void run() {
        monitor.start();

        apply_initial_conditions();

        for (const auto& command : params.commands) {
            execute_command(command);
        }

        measurements.finalize();
        monitor.finish();
    }
};

#endif

#ifndef SHHD_PROJECT_H
#define SHHD_PROJECT_H

#include "simulationinfo.h"
#include "domain.h"
#include "state.h"
#include "spectral_mask.h"
#include "buffer_physical_state.h"
#include "fourier_transform.h"
#include "monitor.h"
#include "model_free_energy.h"
#include "model_free_energy_null.h"
#include "model_free_energy_registry_builtin.h"
#include "model_thermodynamics.h"
#include "model_thermodynamics_null.h"
#include "model_thermodynamics_registry_builtin.h"
#include "model_transport_coefficient.h"
#include "model_transport_coefficient_registry_builtin.h"
#include "fluid_solver.h"
#include "measure_registry_builtin.h"
#include "measure_manager.h"
#include "initial_condition_registry_builtin.h"
#include "initial_condition.h"
#include "fcalculator_dynamics_mode.h"

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
    MeasureRegistry measure_registry;
    InitialConditionRegistry initial_registry;

    Domain2D domain;
    SpectralMask2D spectral_mask;
    State state;
    FourierTransform2D fourier;

    std::unique_ptr<Thermodynamics> thermodynamics_;
    std::unique_ptr<FreeEnergy> free_energy_;
    std::unique_ptr<TransportCoefficient> transport_coefficient_;
    FluidSolver solver;

    MeasurementManager measurements;
    SimulationMonitor monitor;

    int step;
    int run_index;
    double time;

    static std::unique_ptr<Thermodynamics> create_thermodynamics(const Params& params) {
        if (!params.physics.thermo) {
            return std::make_unique<NullThermodynamics>();
        }

        ThermodynamicsRegistry registry = build_thermodynamics_registry();
        const std::string& thermo_type = params.physics.thermo->type;
        const ThermodynamicsStyle& style = registry.get_thermo(thermo_type);
        return style.create(params, params.physics.thermo);
    }

    static std::unique_ptr<FreeEnergy> create_free_energy(const Params& params) {
        if (!params.physics.free_energy) {
            return std::make_unique<NullFreeEnergy>(params.physics.num_order_parameters);
        }

        FreeEnergyRegistry registry = build_free_energy_registry();
        const std::string& free_energy_type = params.physics.free_energy->type;
        const FreeEnergyStyle& style = registry.get_free_energy(free_energy_type);
        return style.create(params, params.physics.free_energy);
    }

    static std::unique_ptr<TransportCoefficient> create_transport_coefficient(const Params& params) {
        if (!params.physics.transport) {
            throw std::runtime_error("Project requires model transport.");
        }

        TransportCoefficientRegistry registry = build_transport_coefficient_registry();
        const std::string& transport_type = params.physics.transport->type;
        const TransportCoefficientStyle& style = registry.get_transport(transport_type);
        return style.create(params, params.physics.transport);
    }

    void apply_initial_conditions() {
        state.clear();

        for (const auto& command : params.initial.density_commands) {
            const DensityInitialConditionStyle& style = initial_registry.get_density(command->type);
            std::unique_ptr<DensityInitialCondition> initial_condition = style.create_initial_condition(params, command);
            initial_condition->apply(state, domain, spectral_mask);
        }

        for (const auto& command : params.initial.momentum_commands) {
            const MomentumInitialConditionStyle& style = initial_registry.get_momentum(command->type);
            std::unique_ptr<MomentumInitialCondition> initial_condition = style.create_initial_condition(params, command);
            initial_condition->apply(state, domain, spectral_mask);
        }

        for (const auto& command : params.initial.order_parameter_commands) {
            const OrderParameterInitialConditionStyle& style = initial_registry.get_order_parameter(command->type);
            std::unique_ptr<OrderParameterInitialCondition> initial_condition = style.create_initial_condition(params, command);
            initial_condition->apply(state, domain, spectral_mask);
        }

        if (is_incompressible_mode(parse_dynamics_mode(params.runtime.time_evolution_type))) {
            validate_initial_momentum_is_transverse(state, domain, spectral_mask);
        }
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

        const FluxRequest flux_request = measurements.flux_request();
        solver.set_flux_request(flux_request);

        for (int local_step = 1; local_step <= nsteps; ++local_step) {
            solver.step(state, time);

            ++step;
            time += params.runtime.dt;

            monitor.print_progress(run_index, local_step, nsteps, step, time);
            measurements.observe(state, fourier, solver.flux_buffer(), step, time);
        }

        monitor.finish_run_segment(run_index, step, time);
    }

public:
    explicit Project(const Params& p)
        : params(p),
          measure_registry(build_measure_registry()),
          initial_registry(build_initial_condition_registry()),
          domain(params),
          spectral_mask(params, domain),
          state(params, domain),
          fourier(domain),
          thermodynamics_(create_thermodynamics(params)),
          free_energy_(create_free_energy(params)),
          transport_coefficient_(create_transport_coefficient(params)),
          solver(
              params,
              domain,
              *thermodynamics_,
              *free_energy_,
              *transport_coefficient_
            ),
          measurements(
              params,
              domain,
              measure_registry,
              *thermodynamics_,
              *free_energy_,
              *transport_coefficient_
          ),
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

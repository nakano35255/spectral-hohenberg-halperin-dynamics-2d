#ifndef PROJECT_H
#define PROJECT_H

#include <string>
#include "simulationinfo.h"
#include "monitor.h"
#include "measure_registry_builtin.h"

class Project {
private:
    Params params;
    MeasureRegistry registry;
    SimulationMonitor monitor;
    int step;
    int run_index;
    double time;

    void execute_command(const Command& command) {
        if (command.type == Command::Type::Run) {
            run_steps(command.run.steps);
        }
        else if (command.type == Command::Type::Fix) {
            monitor.print_fix_command(command.fix);
            // TODO: active_fixes.apply(command.fix);
        }
        else if (command.type == Command::Type::Measure) {
            monitor.print_measure_command(*command.measure);
            // TODO: measurements.apply_measure_command(...);
        }
    }

    void run_steps(int nsteps) {
        ++run_index;

        for (int local_step = 1; local_step <= nsteps; ++local_step) {
            ++step;
            monitor.print_progress(run_index, local_step, nsteps, step, time);
        }

        monitor.finish_run_segment(run_index, step, time);
    }

public:
    explicit Project(const Params& p)
        : params(p),
        registry(build_measure_registry()),
        monitor(params, "output.log"),
        step(0),
        run_index(0),
        time(0.0) {}

    void run() {
        monitor.start();

        for (const auto& command : params.commands) {
            execute_command(command);
        }

        monitor.finish();
    }
};

#endif

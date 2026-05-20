#include "monitor.h"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <stdexcept>

// ---------------------------------------------------------------------- //
SimulationMonitor::SimulationMonitor(const Params& params, const std::string& log_filename, bool active)
    : params(params),
      config(params.thermo),
      active(active) {
    if (active) {
        log_file.open(log_filename);
    }
}
// ---------------------------------------------------------------------- //
SimulationMonitor::~SimulationMonitor() {
    if (log_file.is_open()) {
        log_file.close();
    }
}
// ---------------------------------------------------------------------- //
void SimulationMonitor::write_thermo_header(std::ostream& os) {
    os << "\nTHERMO columns: step time";
    // TODO:
    // Add density, momentum, energy, or other diagnostic columns here
    // once State/FFT/field classes are introduced.
    os << "\n";
}
// ---------------------------------------------------------------------- //
void SimulationMonitor::observe_thermo(int step, double time) {
    if (!active) {
        return;
    }
    if (!config.observe) {
        return;
    }
    if (config.nevery <= 0 || step % config.nevery != 0) {
        return;
    }
    if (!log_file.is_open()) {
        return;
    }

    if (!thermo_header_written) {
        write_thermo_header(log_file);
        thermo_header_written = true;
    }

    log_file << "THERMO "
                << step << " "
                << std::scientific << std::setprecision(8) << time;
    // TODO:
    // Physical observables will be appended here.
    // For example:
    //   mean density of each component
    //   total momentum
    //   kinetic energy
    //   fluctuation statistics

    log_file << "\n";
}
// ---------------------------------------------------------------------- //
void SimulationMonitor::start() {
    if (!active) {
        return;
    }
    if (log_file.is_open()) {
        params.write_summary(log_file);
    }
    params.write_summary(std::cout);

    std::cout << "\nStarting Simulation..." << std::endl;
    if (log_file.is_open()) {
        log_file << "Starting Simulation...\n";
    }

    start_time = std::chrono::high_resolution_clock::now();
}
// ---------------------------------------------------------------------- //
void SimulationMonitor::print_measure_command(const MeasureCommandBase& command) {
    if (!active) {
        return;
    }
    if (!printed_first_command) {
        std::cout << "\n";
        if (log_file.is_open()) {
            log_file << "\n";
        }
        printed_first_command = true;
    }

    command.print(std::cout);
    std::cout << "\n";
    if (log_file.is_open()) {
        command.print(log_file);
        log_file << "\n";
    }
}
void SimulationMonitor::print_progress(int run_index, int local_step, int run_steps, int global_step, double time) {
    if (!active) {
        return;
    }
    observe_thermo(global_step, time);

    if (!config.progress) {
        return;
    }
    if (config.nevery <= 0) {
        return;
    }
    if (local_step % config.nevery != 0 && local_step != run_steps) {
        return;
    }
    std::cout << "Run " << run_index
                << ": local step " << local_step
                << " / " << run_steps
                << "  global step " << global_step
                << "\r" << std::flush;
}
// ---------------------------------------------------------------------- //
void SimulationMonitor::finish_run_segment(int run_index, int global_step, double time) {
    if (!active) {
        return;
    }
    completed_steps = global_step;

    std::cout << "\nRun " << run_index
              << " finished: global step " << global_step
              << " time " << time << "\n\n";
    std::cout << std::string(70, '-') << "\n\n";

    if (log_file.is_open()) {
        log_file << "\nRun " << run_index
                 << " finished: global step " << global_step
                 << " time " << time << "\n\n";
        log_file << std::string(70, '-') << "\n\n";
    }

}
// ---------------------------------------------------------------------- //
void SimulationMonitor::finish() {
    if (!active) {
        return;
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    std::cout << "\nSimulation Finished. Total time: "
              << std::fixed << std::setprecision(3) << elapsed.count() << " s" << std::endl;

    if (log_file.is_open()) {
        log_file << "\n" << std::string(60, '=') << "\n";
        log_file << "Simulation finished.\n";
        log_file << "Total execution time : " << std::fixed << std::setprecision(3)
                 << elapsed.count() << " seconds\n";
        if (elapsed.count() > 0) {
            double steps_per_sec = completed_steps / elapsed.count();
            log_file << "Performance          : " << steps_per_sec << " steps/sec\n";
        }
        log_file << std::string(70, '=') << std::endl;
    }
}
// ---------------------------------------------------------------------- //

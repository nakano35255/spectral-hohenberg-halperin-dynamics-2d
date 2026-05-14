#ifndef MONITOR_H
#define MONITOR_H

#include "simulationinfo.h"

#include <chrono>
#include <fstream>
#include <iosfwd>
#include <string>

class SimulationMonitor {
private:
    const Params& params;
    ThermoConfig config;

    std::ofstream log_file;
    std::chrono::high_resolution_clock::time_point start_time;

    int completed_steps = 0;
    bool printed_first_command = false;
    bool printed_fix_command_header = false;
    bool thermo_header_written = false;

    void write_thermo_header(std::ostream& os);
    void observe_thermo(int step, double time);

    

public:
    explicit SimulationMonitor(const Params& params, const std::string& log_filename = "output.log");
    ~SimulationMonitor();

    void start();
    void print_measure_command(const MeasureCommandBase& command);
    void print_fix_command(const FixCommand& command);
    void print_progress(int run_index, int local_step, int run_steps, int global_step, double time);
    void finish_run_segment(int run_index, int global_step, double time);
    void finish();
};

#endif

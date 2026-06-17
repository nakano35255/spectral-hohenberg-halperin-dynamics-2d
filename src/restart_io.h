#ifndef SHHD_RESTART_IO_H
#define SHHD_RESTART_IO_H

#include "domain.h"
#include "simulationinfo.h"
#include "state.h"

#include <string>

void write_restart_file(const std::string& file, const Params& params, const Domain2D& domain, const State& state, int step, double time);
void read_restart_file(const std::string& file, const Params& params, const Domain2D& domain, State& state, int& step, double& time);

#endif
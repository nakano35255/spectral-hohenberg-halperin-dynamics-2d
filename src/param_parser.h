#ifndef SFI_PARAM_PARSER_H
#define SFI_PARAM_PARSER_H

#include "simulationinfo.h"
#include "model_free_energy_registry.h"
#include "model_thermodynamics_registry.h"
#include "model_transport_coefficient_registry.h"
#include "measure_registry.h"
#include "initial_condition_registry.h"

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

class ParamParser {
private:
    Params& params;
    const ThermodynamicsRegistry& thermodynamics_registry;
    const FreeEnergyRegistry& free_energy_registry;
    const TransportCoefficientRegistry& transport_coefficient_registry;
    const MeasureRegistry& measure_registry;
    const InitialConditionRegistry &initial_condition_registry;

    std::vector<std::string> tokenize(const std::string& line) const;

    void execute_command(const std::vector<std::string>& tokens);
    void validate_configuration() const;

    void parse_model_command(const std::vector<std::string>& tokens);
    void parse_fix_command(const std::vector<std::string>& tokens);
    void parse_set_command(const std::vector<std::string>& tokens);
    void parse_measure_command(const std::vector<std::string>& tokens);
    void parse_thermo_command(const std::vector<std::string>& tokens);
    void parse_restart_command(const std::vector<std::string>& tokens);

    void check_order_parameter_index(int index, const std::string& command_name) const;
    void check_momentum_direction_index(int index, const std::string& command_name) const;

public:
    ParamParser(Params& p, const ThermodynamicsRegistry& tmr, const FreeEnergyRegistry& fer, const TransportCoefficientRegistry& tcr, const MeasureRegistry& mr, const InitialConditionRegistry& icr)
      : params(p),
        thermodynamics_registry(tmr),
        free_energy_registry(fer),
        transport_coefficient_registry(tcr),
        measure_registry(mr),
        initial_condition_registry(icr) {}

    ~ParamParser() {}

    void parse_file(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open input file: " + filename);
        }

        std::string line;
        while (std::getline(file, line)) {
            std::vector<std::string> tokens = tokenize(line);
            execute_command(tokens);
        }

        validate_configuration();
    }
};

#endif

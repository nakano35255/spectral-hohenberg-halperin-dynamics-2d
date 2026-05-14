#ifndef SFI_PARAM_PARSER_H
#define SFI_PARAM_PARSER_H

#include "simulationinfo.h"
#include "measure_registry.h"

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

class ParamParser {
private:
    Params& params;
    const MeasureRegistry& registry;

    std::vector<std::string> tokenize(const std::string& line) const;

    void execute_command(const std::vector<std::string>& tokens);
    void validate_configuration() const;

    void parse_variable_command(const std::vector<std::string>& tokens);
    void parse_l_coeff_command(const std::vector<std::string>& tokens);
    void parse_fix_command(const std::vector<std::string>& tokens);
    void parse_set_command(const std::vector<std::string>& tokens);
    void parse_measure_command(const std::vector<std::string>& tokens);
    void parse_thermo_command(const std::vector<std::string>& tokens);
    void parse_restart_command(const std::vector<std::string>& tokens);

    void check_component_index(int index, const std::string& command_name) const;

public:
    ParamParser(Params& p, const MeasureRegistry& r) : params(p), registry(r) {}
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

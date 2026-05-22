#include "simulationinfo.h"
#include "param_parser.h"
#include "monitor.h"
#include "project.h"
#include "model_thermodynamics_registry_builtin.h"
#include "model_transport_coefficient_registry_builtin.h"
#include "measure_registry_builtin.h"
#include "initial_condition_registry_builtin.h"

#include <exception>
#include <iostream>
#include <mpi.h>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int exit_code = 0;

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_script>\n";
        MPI_Finalize();
        return 1;
    }

    Params params;
    ThermodynamicsRegistry thermo_registry = build_thermodynamics_registry();
    TransportCoefficientRegistry transport_registry = build_transport_coefficient_registry();
    MeasureRegistry measure_registry = build_measure_registry();
    InitialConditionRegistry initial_registry = build_initial_condition_registry();
    ParamParser parser(params, thermo_registry, transport_registry, measure_registry, initial_registry);

    try {
        parser.parse_file(argv[1]);
        Project proj(params);

        proj.run();

    } catch (const std::exception& e) {
        std::cerr << "\n[Fatal Error] " << e.what() << '\n';
        exit_code = 1;
    }

    MPI_Finalize();
    return exit_code;
}

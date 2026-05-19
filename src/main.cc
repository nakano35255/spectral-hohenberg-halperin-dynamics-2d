#include "simulationinfo.h"
#include "param_parser.h"
#include "monitor.h"
#include "project.h"
#include "measure_registry_builtin.h"

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
    MeasureRegistry registry = build_measure_registry();
    ParamParser parser(params, registry);

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

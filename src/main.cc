#include "param_parser.h"
#include "simulationinfo.h"

#include <exception>
#include <iostream>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_script>\n";
        return 1;
    }

    Params params;
    ParamParser parser(params);

    try {
        parser.parse_file(argv[1]);
        params.write_summary(std::cout);
    } catch (const std::exception& e) {
        std::cerr << "\n[Fatal Error] " << e.what() << '\n';
        return 1;
    }

    return 0;
}

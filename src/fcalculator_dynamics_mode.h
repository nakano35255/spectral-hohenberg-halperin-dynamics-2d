#ifndef SFI_FCALCULATOR_DYNAMICS_MODE_H
#define SFI_FCALCULATOR_DYNAMICS_MODE_H

#include <stdexcept>
#include <string>

enum class DynamicsMode {
    Compressible,
    Incompressible,
    Quiescent
};

inline DynamicsMode parse_dynamics_mode(const std::string& time_evolution_type) {
    if (time_evolution_type == "euler/compressible" ||
        time_evolution_type == "srk3/compressible") {
        return DynamicsMode::Compressible;
    }

    if (time_evolution_type == "euler/incompressible" ||
        time_evolution_type == "srk3/incompressible") {
        return DynamicsMode::Incompressible;
    }

    if (time_evolution_type == "euler/quiescent" ||
        time_evolution_type == "srk3/quiescent") {
        return DynamicsMode::Quiescent;
    }

    throw std::runtime_error("unknown time_evolution type: " + time_evolution_type);
}

inline bool is_compressible_mode(DynamicsMode mode) {
    return mode == DynamicsMode::Compressible;
}

inline bool is_incompressible_mode(DynamicsMode mode) {
    return mode == DynamicsMode::Incompressible;
}

inline bool is_quiescent_mode(DynamicsMode mode) {
    return mode == DynamicsMode::Quiescent;
}

#endif
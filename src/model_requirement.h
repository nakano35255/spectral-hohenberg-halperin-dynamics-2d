#ifndef SFI_MODEL_REQUIREMENT_H
#define SFI_MODEL_REQUIREMENT_H

#include "fix_flag.h"

#include <cstdint>

struct ThermodynamicsRequirement {
     bool free_energy = true;
     bool chemical_potential = true;
     bool bulk_pressure = true;
     bool interfacial_force = false;
};

struct TransportCoefficientRequirement {
     bool shear_viscosity = true;
     bool bulk_viscosity = true;
     bool mobility = true;
     bool noise_factor = true;
};

struct ModelRequirement {
     ThermodynamicsRequirement thermodynamics;
     TransportCoefficientRequirement transport_coefficient;
};

inline ModelRequirement full_requirement() {
    ModelRequirement req;
    req.thermodynamics.free_energy = true; // The free energy must be specified whenever any thermodynamic quantity is required.
    req.thermodynamics.chemical_potential = true;
    req.thermodynamics.bulk_pressure = true;
    req.thermodynamics.interfacial_force = false;

    req.transport_coefficient.shear_viscosity = true;
    req.transport_coefficient.bulk_viscosity = true;
    req.transport_coefficient.mobility = true;
    req.transport_coefficient.noise_factor = true;
    return req;
}

inline ModelRequirement build_model_requirement(std::uint32_t fix_flags) {
     ModelRequirement req = full_requirement();

     if (fix_enabled(fix_flags, FixFlag::Quiescent)) {
          req.thermodynamics.bulk_pressure = false;
          req.transport_coefficient.shear_viscosity = false;
          req.transport_coefficient.bulk_viscosity = false;
     }

     if (!fix_enabled(fix_flags, FixFlag::Noise)) {
          req.transport_coefficient.noise_factor = false;
     }

     if (fix_enabled(fix_flags, FixFlag::Incompressible)) {
          req.thermodynamics.bulk_pressure = false;
          req.transport_coefficient.bulk_viscosity = false;
     }

     if (fix_enabled(fix_flags, FixFlag::Interface)) {
          req.thermodynamics.interfacial_force = true;
     }

     return req;
}

#endif

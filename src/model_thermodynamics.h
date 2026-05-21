#ifndef SFI_MODEL_THERMODYNAMICS_H
#define SFI_MODEL_THERMODYNAMICS_H

#include <vector>

struct ThermodynamicsRequest {
    bool free_energy = false;
    bool chemical_potential = true;
    bool pressure = true;
};

class ThermodynamicsModel {
public:
    virtual ~ThermodynamicsModel() = default;

    virtual double free_energy(const std::vector<double>& rho) const = 0;

    virtual void chemical_potential(const std::vector<double>& rho, std::vector<double>& mu) const = 0;

    virtual double pressure(const std::vector<double>& rho) const = 0;
};

#endif

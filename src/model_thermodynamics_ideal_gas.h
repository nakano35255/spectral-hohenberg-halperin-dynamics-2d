#ifndef SFI_MODEL_THERMODYNAMICS_IDEAL_GAS_H
#define SFI_MODEL_THERMODYNAMICS_IDEAL_GAS_H

#include "model_thermodynamics.h"

#include <vector>

class IdealGasThermodynamicsModel : public ThermodynamicsModel {
private:
    double kB_ = 1.0;
    double temperature_ = 1.0;
    std::vector<double> mass_;
    std::vector<double> reference_density_;

    void check_density_size(const std::vector<double>& rho) const;

public:
    IdealGasThermodynamicsModel(
        double kB,
        double temperature,
        std::vector<double> mass,
        std::vector<double> reference_density
    );

    double free_energy(const std::vector<double>& rho) const override;

    void chemical_potential(
        const std::vector<double>& rho,
        std::vector<double>& mu
    ) const override;

    double pressure(const std::vector<double>& rho) const override;
};

#endif

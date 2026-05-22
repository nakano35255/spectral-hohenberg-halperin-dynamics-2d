#ifndef SFI_MODEL_TRANSPORT_COEFFICIENT_CONSTANT_H
#define SFI_MODEL_TRANSPORT_COEFFICIENT_CONSTANT_H

#include "model_transport_coefficient.h"

#include <vector>

class ConstantTransportCoefficient : public TransportCoefficient {
private:
    int num_components_ = 0;
    double shear_viscosity_ = 0.0;
    double bulk_viscosity_ = 0.0;
    std::vector<double> mobility_;
    std::vector<double> noise_factor_;

    void check_density_size(const std::vector<double>& rho) const;
    void build_noise_factor();

public:
    ConstantTransportCoefficient(
        int num_components,
        double shear_viscosity,
        double bulk_viscosity,
        std::vector<double> mobility
    );

    double shear_viscosity(const std::vector<double>& rho) const override;
    double bulk_viscosity(const std::vector<double>& rho) const override;

    void mobility(const std::vector<double>& rho, std::vector<double>& L) const override;
    void noise_factor(const std::vector<double>& rho, std::vector<double>& B) const override;
};

#endif

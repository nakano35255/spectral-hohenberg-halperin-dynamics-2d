#ifndef SFI_MODEL_TRANSPORT_COEFFICIENT_H
#define SFI_MODEL_TRANSPORT_COEFFICIENT_H

#include <vector>

struct TransportCoefficientRequest {
    bool shear_viscosity = true;
    bool bulk_viscosity = true;
    bool mobility = true;
};

class TransportCoefficientModel {
public:
    virtual ~TransportCoefficientModel() = default;

    virtual double shear_viscosity(const std::vector<double>& rho) const = 0;
    virtual double bulk_viscosity(const std::vector<double>& rho) const = 0;

    virtual void mobility(const std::vector<double>& rho, std::vector<double>& L) const = 0;
};

#endif

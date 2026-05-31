#ifndef SHHD_MODEL_TRANSPORT_COEFFICIENT_H
#define SHHD_MODEL_TRANSPORT_COEFFICIENT_H

#include "simulationinfo.h"

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

class TransportCoefficient {
private:
    int num_order_parameters_ = 0;
    double shear_viscosity_ = 0.0;
    double bulk_viscosity_ = 0.0;
    std::vector<double> order_parameter_mobility_;

public:
    TransportCoefficient(
        int num_order_parameters,
        double shear_viscosity,
        double bulk_viscosity,
        std::vector<double> order_parameter_mobility
    ) : num_order_parameters_(num_order_parameters),
        shear_viscosity_(shear_viscosity),
        bulk_viscosity_(bulk_viscosity),
        order_parameter_mobility_(std::move(order_parameter_mobility))
    {
        if (num_order_parameters_ < 0) {
            throw std::runtime_error("transport coefficient requires nonnegative number of order parameters.");
        }
        if (shear_viscosity_ < 0.0) {
            throw std::runtime_error("transport coefficient requires nonnegative eta.");
        }
        if (bulk_viscosity_ < 0.0) {
            throw std::runtime_error("transport coefficient requires nonnegative zeta.");
        }
        if (order_parameter_mobility_.size() != static_cast<std::size_t>(num_order_parameters_)) {
            throw std::runtime_error("transport coefficient mobility size mismatch.");
        }
        for (double mobility : order_parameter_mobility_) {
            if (mobility < 0.0) {
                throw std::runtime_error("transport coefficient requires nonnegative order parameter mobility.");
            }
        }
    }

    virtual ~TransportCoefficient() = default;

    int num_order_parameters() const {
        return num_order_parameters_;
    }

    double shear_viscosity() const {
        return shear_viscosity_;
    }

    double bulk_viscosity() const {
        return bulk_viscosity_;
    }

    const std::vector<double>& order_parameter_mobility() const {
        return order_parameter_mobility_;
    }

    virtual bool has_physical_shear_viscosity() const {
        return false;
    }

    virtual bool has_physical_bulk_viscosity() const {
        return false;
    }

    virtual double physical_shear_viscosity(double rho, const double* psi) const {
        (void)rho;
        (void)psi;
        return shear_viscosity_;
    }

    virtual double physical_bulk_viscosity(double rho, const double* psi) const {
        (void)rho;
        (void)psi;
        return bulk_viscosity_;
    }
};

#endif

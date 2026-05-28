#ifndef SFI_MODEL_FREE_ENERGY_H
#define SFI_MODEL_FREE_ENERGY_H

#include "simulationinfo.h"

#include <stdexcept>

class FreeEnergy {
private:
    int num_order_parameters_ = 0;

public:
    explicit FreeEnergy(int num_order_parameters)
        : num_order_parameters_(num_order_parameters)
    {
        if (num_order_parameters_ < 0) {
            throw std::runtime_error("free energy requires nonnegative number of order parameters.");
        }
    }
    virtual ~FreeEnergy() = default;

    int num_order_parameters() const {
        return num_order_parameters_;
    }

    // for chemical potential
    virtual double chemical_potential_k0_coefficient(int order_parameter) const {
        check_order_parameter_index(order_parameter);
        return 0.0;
    }
    virtual double chemical_potential_k2_coefficient(int order_parameter) const {
        check_order_parameter_index(order_parameter);
        return 0.0;
    }
    virtual double chemical_potential_k4_coefficient(int order_parameter) const {
        check_order_parameter_index(order_parameter);
        return 0.0;
    }
    virtual bool has_physical_chemical_potential() const {
        return false;
    }
    virtual double physical_chemical_potential(int order_parameter, const double* psi) const {
        check_order_parameter_index(order_parameter);
        (void)psi;
        return 0.0;
    }

    // for free energy
    virtual double physical_free_energy_density(const double* psi) const {
        (void)psi;
        return 0.0;
    }

protected:
    void check_order_parameter_index(int order_parameter) const {
        if (order_parameter < 0 || order_parameter >= num_order_parameters_) {
            throw std::runtime_error("free energy order parameter index out of range.");
        }
    }
};

#endif

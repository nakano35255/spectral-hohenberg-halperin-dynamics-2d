#include "model_free_energy_phi4.h"

#include <cstddef>
#include <stdexcept>
#include <utility>

// ---------------------------------------------------------------------- //
Phi4FreeEnergy::Phi4FreeEnergy(
    int num_order_parameters,
    std::vector<double> k0,
    std::vector<double> k2,
    std::vector<double> phi4
) : FreeEnergy(num_order_parameters),
    k0_(std::move(k0)),
    k2_(std::move(k2)),
    phi4_(std::move(phi4))
{
    const std::size_t expected_size = static_cast<std::size_t>(num_order_parameters);
    if (k0_.size() != expected_size || k2_.size() != expected_size || phi4_.size() != expected_size) {
        throw std::runtime_error("phi4 free energy coefficient size mismatch.");
    }

    for (double coefficient : phi4_) {
        if (coefficient < 0.0) {
            throw std::runtime_error("phi4 free energy requires nonnegative phi4 coefficient.");
        }
        if (coefficient != 0.0) {
            has_physical_chemical_potential_ = true;
        }
    }
}
// ---------------------------------------------------------------------- //
double Phi4FreeEnergy::chemical_potential_k0_coefficient(int order_parameter) const {
    check_order_parameter_index(order_parameter);
    return k0_[static_cast<std::size_t>(order_parameter)];
}
// ---------------------------------------------------------------------- //
double Phi4FreeEnergy::chemical_potential_k2_coefficient(int order_parameter) const {
    check_order_parameter_index(order_parameter);
    return k2_[static_cast<std::size_t>(order_parameter)];
}
// ---------------------------------------------------------------------- //
bool Phi4FreeEnergy::has_physical_chemical_potential() const {
    return has_physical_chemical_potential_;
}
// ---------------------------------------------------------------------- //
double Phi4FreeEnergy::physical_chemical_potential(
    int order_parameter,
    const double* psi
) const {
    check_order_parameter_index(order_parameter);

    const double value = psi[order_parameter];
    return phi4_[static_cast<std::size_t>(order_parameter)] * value * value * value;
}
// ---------------------------------------------------------------------- //

#include "model_free_energy_phi4.h"

#include <cstddef>
#include <stdexcept>
#include <utility>

// ---------------------------------------------------------------------- //
Phi4FreeEnergy::Phi4FreeEnergy(
    int num_order_parameters,
    std::vector<double> a,
    std::vector<double> b,
    std::vector<double> u
) : FreeEnergy(num_order_parameters),
    a_(std::move(a)),
    b_(std::move(b)),
    u_(std::move(u))
{
    const std::size_t expected_size = static_cast<std::size_t>(num_order_parameters);
    if (a_.size() != expected_size || b_.size() != expected_size || u_.size() != expected_size) {
        throw std::runtime_error("phi4 free energy coefficient size mismatch.");
    }

    for (double coefficient : u_) {
        if (coefficient < 0.0) {
            throw std::runtime_error("phi4 free energy requires nonnegative u coefficient.");
        }
        if (coefficient != 0.0) {
            has_physical_chemical_potential_ = true;
        }
    }
}
// ---------------------------------------------------------------------- //
double Phi4FreeEnergy::chemical_potential_k0_coefficient(int order_parameter) const {
    check_order_parameter_index(order_parameter);
    return a_[static_cast<std::size_t>(order_parameter)];
}
// ---------------------------------------------------------------------- //
double Phi4FreeEnergy::chemical_potential_k2_coefficient(int order_parameter) const {
    check_order_parameter_index(order_parameter);
    return b_[static_cast<std::size_t>(order_parameter)];
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
    return u_[static_cast<std::size_t>(order_parameter)] * value * value * value;
}
// ---------------------------------------------------------------------- //

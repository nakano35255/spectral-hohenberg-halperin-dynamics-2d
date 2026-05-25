#include "model_free_energy_quadratic.h"

#include <cstddef>
#include <stdexcept>
#include <utility>

// ---------------------------------------------------------------------- //
QuadraticFreeEnergy::QuadraticFreeEnergy(
    int num_order_parameters,
    std::vector<double> coefficients
) : FreeEnergy(num_order_parameters),
    coefficients_(std::move(coefficients)) {
    if (coefficients_.size() != static_cast<std::size_t>(num_order_parameters)) {
        throw std::runtime_error("quadratic free energy coefficient size mismatch.");
    }
}
// ---------------------------------------------------------------------- //
double QuadraticFreeEnergy::chemical_potential_k0_coefficient(int order_parameter) const {
    check_order_parameter_index(order_parameter);
    return coefficients_[static_cast<std::size_t>(order_parameter)];
}
// ---------------------------------------------------------------------- //

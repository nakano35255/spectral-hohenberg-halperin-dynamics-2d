#include "model_thermodynamics_ideal_gas.h"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>

IdealGasThermodynamics::IdealGasThermodynamics(
    double kB,
    double temperature,
    std::vector<double> mass,
    std::vector<double> reference_density
) : kB_(kB),
    temperature_(temperature),
    mass_(std::move(mass)),
    reference_density_(std::move(reference_density)) {
    if (kB_ <= 0.0) {
        throw std::runtime_error("ideal_gas thermodynamics requires positive kB.");
    }
    if (temperature_ <= 0.0) {
        throw std::runtime_error("ideal_gas thermodynamics requires positive T.");
    }
    if (mass_.empty()) {
        throw std::runtime_error("ideal_gas thermodynamics requires at least one component.");
    }
    if (mass_.size() != reference_density_.size()) {
        throw std::runtime_error("ideal_gas thermodynamics mass/reference-density size mismatch.");
    }

    for (std::size_t k = 0; k < mass_.size(); ++k) {
        if (mass_[k] <= 0.0) {
            throw std::runtime_error("ideal_gas thermodynamics requires positive molecular mass.");
        }
        if (reference_density_[k] <= 0.0) {
            throw std::runtime_error("ideal_gas thermodynamics requires positive reference density.");
        }
    }
}

void IdealGasThermodynamics::check_density_size(const std::vector<double>& rho) const {
    if (rho.size() != mass_.size()) {
        throw std::runtime_error("ideal_gas thermodynamics expected " + std::to_string(mass_.size()) + " density components, got " + std::to_string(rho.size()));
    }
}

double IdealGasThermodynamics::free_energy(const std::vector<double>& rho) const {
    check_density_size(rho);

    double f = 0.0;
    for (std::size_t k = 0; k < rho.size(); ++k) {
        if (rho[k] <= 0.0) {
            throw std::runtime_error("ideal_gas thermodynamics requires positive density.");
        }

        const double prefactor = kB_ * temperature_ / mass_[k];
        f += prefactor * rho[k] * (std::log(rho[k] / reference_density_[k]) - 1.0);
    }

    return f;
}

void IdealGasThermodynamics::chemical_potential(const std::vector<double>& rho, std::vector<double>& mu) const {
    check_density_size(rho);

    mu.assign(rho.size(), 0.0);
    for (std::size_t k = 0; k < rho.size(); ++k) {
        if (rho[k] <= 0.0) {
            throw std::runtime_error("ideal_gas thermodynamics requires positive density.");
        }

        const double prefactor = kB_ * temperature_ / mass_[k];
        mu[k] = prefactor * std::log(rho[k] / reference_density_[k]);
    }
}

double IdealGasThermodynamics::pressure(const std::vector<double>& rho) const {
    check_density_size(rho);

    double p = 0.0;
    for (std::size_t k = 0; k < rho.size(); ++k) {
        if (rho[k] <= 0.0) {
            throw std::runtime_error("ideal_gas thermodynamics requires positive density.");
        }

        p += kB_ * temperature_ / mass_[k] * rho[k];
    }

    return p;
}

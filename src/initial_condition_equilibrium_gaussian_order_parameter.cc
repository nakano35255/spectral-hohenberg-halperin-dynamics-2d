#include "initial_condition_equilibrium_gaussian_order_parameter.h"
#include "initial_condition_equilibrium_gaussian_order_parameter_style.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace {
    constexpr double PI = 3.14159265358979323846;

    int signed_ky(int gy, int ny) {
        return (gy <= ny / 2) ? gy : gy - ny;
    }

    std::uint64_t splitmix64(std::uint64_t value) {
        value += 0x9e3779b97f4a7c15ULL;
        value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
        value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
        return value ^ (value >> 31);
    }

    double uniform01_from_mode(int seed, int channel, int gx, int ky, int draw) {
        std::uint64_t key = static_cast<std::uint64_t>(static_cast<std::uint32_t>(seed));
        key ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(channel)) << 48;
        key ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(gx)) << 32;
        key ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(ky)) << 8;
        key ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(draw));
        const std::uint64_t random = splitmix64(key);
        return static_cast<double>(random >> 11) * (1.0 / 9007199254740992.0);
    }

    Complex gaussian_mode(int seed, int channel, const SpectralMode2D& mode, int ny) {
        const int ky = signed_ky(mode.gy, ny);
        const int phase_ky = (mode.gx == 0) ? std::abs(ky) : ky;

        const double u1 = std::max(uniform01_from_mode(seed, channel, mode.gx, phase_ky, 0), 1.0e-300);
        const double u2 = uniform01_from_mode(seed, channel, mode.gx, phase_ky, 1);
        const double r = std::sqrt(-2.0 * std::log(u1));
        const double theta = 2.0 * PI * u2;

        Complex value(r * std::cos(theta), r * std::sin(theta));
        if (mode.gx == 0 && ky < 0) {
            value = std::conj(value);
        }
        return value;
    }
}


EquilibriumGaussianOrderParameterInitialCondition::EquilibriumGaussianOrderParameterInitialCondition(
    const Params& params,
    std::shared_ptr<const OrderParameterInitialConditionCommandBase> command
) : OrderParameterInitialCondition(params, command) {
    (void)params;

    const auto cfg = std::dynamic_pointer_cast<const EquilibriumGaussianOrderParameterInitialConditionCommand>(command);
    if (!cfg) {
        throw std::runtime_error("EquilibriumGaussianOrderParameterInitialCondition: invalid command type.");
    }

    seed_ = cfg->seed;
    kBT_ = cfg->kBT;
    psi0_ = cfg->psi0;
    k0_coefficient_ = cfg->k0_coefficient;
    k2_coefficient_ = cfg->k2_coefficient;
    k4_coefficient_ = cfg->k4_coefficient;
}

void EquilibriumGaussianOrderParameterInitialCondition::apply(
    State& state,
    const Domain2D& domain,
    const SpectralMask2D& spectral_mask
) const {
    Complex* psi = state.psi_hat_data(component_);
    std::fill(psi, psi + domain.spectral_size(), Complex(0.0, 0.0));

    const double grid_size =
        static_cast<double>(domain.nx_global()) *
        static_cast<double>(domain.ny_global());

    for (const SpectralMode2D& mode : spectral_mask.active_modes()) {
        if (mode.gx == 0 && mode.gy == 0) {
            psi[mode.index] = Complex(psi0_ * grid_size, 0.0);
            break;
        }
    }

    if (kBT_ == 0.0) {
        return;
    }

    const double volume = domain.lx() * domain.ly();
    const double cell_area = volume / grid_size;

    for (const SpectralMode2D& mode : spectral_mask.active_modes()) {
        if (mode.k2 == 0.0) continue;

        const double k2 = mode.k2;
        const double coefficient = k0_coefficient_ + k2_coefficient_ * k2 + k4_coefficient_ * k2 * k2;
        if (coefficient <= 0.0) {
            throw std::runtime_error("set order_parameter equilibrium/gaussian requires positive linear chemical-potential coefficient for every active nonzero mode.");
        }

        const double sigma = std::sqrt(kBT_ * grid_size / (2.0 * coefficient * cell_area));
        psi[mode.index] = sigma * gaussian_mode(seed_, component_, mode, domain.ny_global());
    }
}

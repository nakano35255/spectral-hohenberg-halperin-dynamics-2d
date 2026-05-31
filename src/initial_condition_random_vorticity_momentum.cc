#include "initial_condition_random_vorticity_momentum.h"
#include "initial_condition_random_vorticity_momentum_style.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <mpi.h>
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

    double uniform01_from_mode(int seed, int gx, int ky) {
        std::uint64_t key = static_cast<std::uint64_t>(static_cast<std::uint32_t>(seed));
        key ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(gx)) << 32;
        key ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(ky));
        const std::uint64_t random = splitmix64(key);
        return static_cast<double>(random >> 11) * (1.0 / 9007199254740992.0);
    }

    Complex random_phase(int seed, int gx, int ky) {
        const double theta = 2.0 * PI * uniform01_from_mode(seed, gx, ky);
        return Complex(std::cos(theta), std::sin(theta));
    }

    double reference_density(const State& state, const Domain2D& domain) {
        const double grid_size = static_cast<double>(domain.nx_global()) * static_cast<double>(domain.ny_global());

        double local_rho0 = 0.0;
        const Box2D& box = domain.spectral_box();
        if (box.low[0] <= 0 && 0 <= box.high[0] && box.low[1] <= 0 && 0 <= box.high[1]) {
            const std::size_t local_nkx = static_cast<std::size_t>(box.size[0]);
            const std::size_t index = static_cast<std::size_t>(0 - box.low[1]) * local_nkx + static_cast<std::size_t>(0 - box.low[0]);
            local_rho0 = state.rho_hat_data()[index].real() / grid_size;
        }

        double rho0 = 0.0;
        MPI_Allreduce(&local_rho0, &rho0, 1, MPI_DOUBLE, MPI_SUM, domain.comm());

        if (rho0 == 0.0) {
            throw std::runtime_error("set momentum random_vorticity requires nonzero reference density.");
        }

        return rho0;
    }

    Complex vorticity_hat(const SpectralMode2D& mode, int ny, double n0, double sigma, int seed) {
        if (mode.k2 == 0.0) {
            return Complex(0.0, 0.0);
        }

        const int ky = signed_ky(mode.gy, ny);
        const int phase_ky = (mode.gx == 0) ? std::abs(ky) : ky;
        const double mode_nx = static_cast<double>(mode.gx);
        const double mode_ny = static_cast<double>(ky);
        const double n = std::sqrt(mode_nx * mode_nx + mode_ny * mode_ny);
        const double normalized_distance = (n - n0) / sigma;
        const double envelope = std::exp(-0.5 * normalized_distance * normalized_distance);
        Complex value = envelope * random_phase(seed, mode.gx, phase_ky);

        if (mode.gx == 0 && ky < 0) {
            value = std::conj(value);
        }

        return value;
    }

    Complex velocity_hat_from_vorticity(const SpectralMode2D& mode, Complex omega_hat, int direction) {
        if (mode.k2 == 0.0) {
            return Complex(0.0, 0.0);
        }

        if (direction == 0) {
            return Complex(0.0, mode.ky / mode.k2) * omega_hat;
        }
        if (direction == 1) {
            return -Complex(0.0, mode.kx / mode.k2) * omega_hat;
        }

        throw std::runtime_error("RandomVorticityMomentumInitialCondition: direction must be 0 or 1.");
    }
}

RandomVorticityMomentumInitialCondition::RandomVorticityMomentumInitialCondition(
    const Params& params,
    std::shared_ptr<const MomentumInitialConditionCommandBase> command
) : MomentumInitialCondition(params, command) {
    (void)params;

    const auto cfg = std::dynamic_pointer_cast<const RandomVorticityMomentumInitialConditionCommand>(command);
    if (!cfg) {
        throw std::runtime_error("RandomVorticityMomentumInitialCondition: invalid command type.");
    }

    n0_ = cfg->n0;
    sigma_ = cfg->sigma;
    urms_ = cfg->urms;
    seed_ = cfg->seed;
}

void RandomVorticityMomentumInitialCondition::apply(State& state, const Domain2D& domain, const SpectralMask2D& spectral_mask) const {
    const int ny = domain.ny_global();
    const double grid_size = static_cast<double>(domain.nx_global()) * static_cast<double>(domain.ny_global());

    double local_velocity_norm = 0.0;

    for (const SpectralMode2D& mode : spectral_mask.active_modes()) {
        const Complex omega = vorticity_hat(mode, ny, n0_, sigma_, seed_);
        const Complex ux = velocity_hat_from_vorticity(mode, omega, 0);
        const Complex uy = velocity_hat_from_vorticity(mode, omega, 1);
        const double weight = mode.gx == 0 ? 1.0 : 2.0;
        local_velocity_norm += weight * (std::norm(ux) + std::norm(uy));
    }

    double global_velocity_norm = 0.0;
    MPI_Allreduce(&local_velocity_norm, &global_velocity_norm, 1, MPI_DOUBLE, MPI_SUM, domain.comm());

    if (global_velocity_norm == 0.0) {
        throw std::runtime_error("set momentum random_vorticity generated zero velocity field.");
    }

    const double rho0 = reference_density(state, domain);
    const double scale = urms_ * grid_size / std::sqrt(global_velocity_norm);

    Complex* momentum = nullptr;
    if (direction_ == 0) momentum = state.jx_hat_data();
    else if (direction_ == 1) momentum = state.jy_hat_data();
    else { throw std::runtime_error("RandomVorticityMomentumInitialCondition: direction must be 0 or 1.");}

    std::fill(momentum, momentum + domain.spectral_size(), Complex(0.0, 0.0));

    for (const SpectralMode2D& mode : spectral_mask.active_modes()) {
        const Complex omega = vorticity_hat(mode, ny, n0_, sigma_, seed_);
        momentum[mode.index] = rho0 * scale * velocity_hat_from_vorticity(mode, omega, direction_);
    }
}

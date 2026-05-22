#ifndef SFI_FLUID_CONSTRAINT_H
#define SFI_FLUID_CONSTRAINT_H

#include "domain.h"
#include "fix_flag.h"
#include "simulationinfo.h"
#include "state.h"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

class FluidConstraint {
private:
    const Params& params_;
    int num_components_;
    std::size_t local_spectral_size_;
    std::vector<Complex> rho_total_hat_;
    bool initialized_ = false;

    void zero_field(Complex* field) const {
        std::fill(field, field + local_spectral_size_, Complex(0.0, 0.0));
    }

    void apply_quiescent(State& state) const {
        if (!initialized_) {
            throw std::runtime_error("FluidConstraint must be initialized before applying quiescent constraints.");
        }

        Complex* last_density = state.rho_hat_data(num_components_ - 1);
        for (std::size_t i = 0; i < local_spectral_size_; ++i) {
            Complex value = rho_total_hat_[i];
            for (int component = 0; component < num_components_ - 1; ++component) {
                value -= state.data()[static_cast<std::size_t>(component) * local_spectral_size_ + i];
            }
            last_density[i] = value;
        }

        zero_field(state.jx_hat_data());
        zero_field(state.jy_hat_data());
    }

public:
    FluidConstraint(const Domain2D& domain, const Params& params)
        : params_(params),
          num_components_(params.physics.num_components),
          local_spectral_size_(domain.spectral_size()) {}

    void initialize(const State& state) {
        if (!fix_enabled(params_.fix.flags, FixFlag::Quiescent)) {
            initialized_ = true;
            return;
        }

        rho_total_hat_.assign(local_spectral_size_, Complex(0.0, 0.0));
        for (int component = 0; component < num_components_; ++component) {
            const Complex* rho = state.data() + static_cast<std::size_t>(component) * local_spectral_size_;
            for (std::size_t i = 0; i < local_spectral_size_; ++i) {
                rho_total_hat_[i] += rho[i];
            }
        }
        initialized_ = true;
    }

    void apply(State& state) const {
        if (fix_enabled(params_.fix.flags, FixFlag::Quiescent)) {
            apply_quiescent(state);
        }
    }
};

#endif

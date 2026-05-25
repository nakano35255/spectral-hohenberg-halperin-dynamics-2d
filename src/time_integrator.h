#ifndef SFI_TIME_INTEGRATOR_H
#define SFI_TIME_INTEGRATOR_H

#include "simulationinfo.h"
#include "domain.h"
#include "state.h"
#include "spectral_mask.h"

#include <algorithm>
#include <functional>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

using DensityDetRHSFunc = std::function<void(const State&, Complex*, double)>;
using OrderParameterDetRHSFunc = std::function<void(int, const State&, Complex*, double)>;
using MomentumDetRHSFunc = std::function<void(const State&, Complex*, Complex*, double)>;

using DensityStoRHSFunc = std::function<void(const State&, Complex*)>;
using OrderParameterStoRHSFunc = std::function<void(int, const State&, Complex*)>;
using MomentumStoRHSFunc = std::function<void(const State&, Complex*, Complex*)>;


struct RHSOperators {
    DensityDetRHSFunc rho_det;
    OrderParameterDetRHSFunc psi_det;
    MomentumDetRHSFunc j_det;

    OrderParameterStoRHSFunc psi_sto;
    MomentumStoRHSFunc j_sto;
};

class TimeIntegrator {
protected:
    const Domain2D& domain_;
    const Params& params_;
    const SpectralMask2D& spectral_mask_;

    const int num_order_parameters_;
    const int num_fields_;
    const std::size_t local_spectral_size_;
    const std::size_t total_spectral_size_;
    const double dt_;
    const double sqrt_dt_;

    void clear_state(State& state) const {
        std::fill(state.data(), state.data() + total_spectral_size_, Complex(0.0, 0.0));
    }

    void copy_state(State& dst, const State& src) const {
        std::copy(src.data(), src.data() + total_spectral_size_, dst.data());
    }

    bool owns_spectral_mode(int gx, int gy) const {
        const Box2D& box = domain_.spectral_box();
        return !box.empty() && gx >= box.low[0] && gx <= box.high[0] && gy >= box.low[1] && gy <= box.high[1];
    }

    Complex& spectral_field(State& state, int field, int gx, int gy) const {
        return state.data()[static_cast<std::size_t>(field) * local_spectral_size_ + domain_.spectral_local_index(gx, gy)];
    }

    void enforce_self_conjugate_line(State& state, int kx) const {
        if (!owns_spectral_mode(kx, 0)) {
            return;
        }

        const int ny = domain_.ny_global();

        for (int field = 0; field < num_fields_; ++field) {
            for (int gy = 0; gy < ny; ++gy) {
                if (!spectral_mask_.active(kx, gy)) {
                    continue;
                }

                const int conjugate_gy = (ny - gy) % ny;

                if (conjugate_gy == gy) {
                    Complex& value = spectral_field(state, field, kx, gy);
                    value = Complex(value.real(), 0.0);
                    continue;
                }

                if (gy > conjugate_gy) {
                    continue;
                }

                Complex& a = spectral_field(state, field, kx, gy);
                Complex& b = spectral_field(state, field, kx, conjugate_gy);

                const Complex projected = 0.5 * (a + std::conj(b));
                a = projected;
                b = std::conj(projected);
            }
        }
    }

    void enforce_real_symmetry(State& state) const {
        spectral_mask_.apply_many(num_fields_, state.data());

        enforce_self_conjugate_line(state, 0);
    }

public:
    TimeIntegrator(
        const Domain2D& domain,
        const Params& params,
        const SpectralMask2D& spectral_mask
    )
        : domain_(domain),
          params_(params),
          spectral_mask_(spectral_mask),
          num_order_parameters_(params.physics.num_order_parameters),
          num_fields_(params.physics.num_order_parameters + 3),
          local_spectral_size_(domain.spectral_size()),
          total_spectral_size_(static_cast<std::size_t>(num_fields_) * local_spectral_size_),
          dt_(params.runtime.dt),
          sqrt_dt_(std::sqrt(dt_))
    {
        if (num_order_parameters_ < 0) {
            throw std::runtime_error("TimeIntegrator requires a nonnegative number of order parameters.");
        }
    }

    virtual ~TimeIntegrator() = default;

    virtual void step(State& u, double t, const RHSOperators& rhs) = 0;

};

#endif

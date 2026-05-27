#include "fcalculator_workspace.h"

#include <algorithm>
#include <stdexcept>

// ---------------------------------------------------------------------- //
FCalculatorWorkspace::FCalculatorWorkspace(
    const Domain2D& domain,
    const Params& params,
    DynamicsMode dynamics_mode,
    const PhysicalRHSFieldRequests& physical_field_requests,
    int temporary_field_capacity
)
    : domain_(domain),
      fourier_(domain),
      dynamics_mode_(dynamics_mode),
      field_layout_(PhysicalFieldLayout::from_requests(
          dynamics_mode,
          physical_field_requests,
          params.physics.num_order_parameters
      )),
      num_order_parameters_(params.physics.num_order_parameters),
      temporary_field_capacity_(temporary_field_capacity),
      local_physical_size_(domain.physical_size()),
      local_spectral_size_(domain.spectral_size())
{
    if (num_order_parameters_ < 0) {
        throw std::runtime_error("FCalculatorWorkspace requires nonnegative num_order_parameters.");
    }

    if (temporary_field_capacity_ < 0) {
        throw std::runtime_error("FCalculatorWorkspace requires nonnegative temporary field capacity.");
    }

    spectral_pack_.assign(static_cast<std::size_t>(field_layout_.ntransform) * local_spectral_size_, Complex(0.0, 0.0));
    physical_pack_.assign(static_cast<std::size_t>(field_layout_.nphysical) * local_physical_size_, 0.0);
    temporary_physical_.assign(static_cast<std::size_t>(temporary_field_capacity_) * local_physical_size_, 0.0);
    temporary_spectral_.assign(static_cast<std::size_t>(temporary_field_capacity_) * local_spectral_size_, Complex(0.0, 0.0));
}
// ---------------------------------------------------------------------- //
double FCalculatorWorkspace::reference_density(const State& state) {
    if (reference_density_ready_) {
        return reference_density_;
    }

    const double grid_size = static_cast<double>(domain_.nx_global()) * static_cast<double>(domain_.ny_global());

    double local_rho0 = 0.0;

    const Box2D& box = domain_.spectral_box();
    if (box.low[0] <= 0 && 0 <= box.high[0] && box.low[1] <= 0 && 0 <= box.high[1]) {
        const std::size_t local_nkx = static_cast<std::size_t>(box.size[0]);
        const std::size_t index = static_cast<std::size_t>(0 - box.low[1]) * local_nkx + static_cast<std::size_t>(0 - box.low[0]);

        local_rho0 = state.rho_hat_data()[index].real() / grid_size;
    }

    MPI_Allreduce(&local_rho0, &reference_density_, 1, MPI_DOUBLE, MPI_SUM, domain_.comm());

    if (reference_density_ == 0.0) {
        throw std::runtime_error("reference density must be nonzero.");
    }

    reference_density_ready_ = true;
    return reference_density_;
}
// ---------------------------------------------------------------------- //
void FCalculatorWorkspace::ensure_physical_fields(const State& state, double time) {
    if (cached_state_ == &state && cached_time_ == time && physical_fields_ready_) {
        return;
    }

    cached_state_ = &state;
    cached_time_ = time;
    physical_fields_ready_ = false;

    std::fill(spectral_pack_.begin(), spectral_pack_.end(), Complex(0.0, 0.0));
    std::fill(physical_pack_.begin(), physical_pack_.end(), 0.0);

    const PhysicalFieldOffsets& offsets = field_layout_.offsets;

    if (field_layout_.transform_psi) {
        std::copy(
            state.psi_hat_data(0),
            state.psi_hat_data(0) + static_cast<std::size_t>(num_order_parameters_) * local_spectral_size_,
            spectral_pack_.data() + static_cast<std::size_t>(offsets.psi) * local_spectral_size_
        );
    }

    if (field_layout_.transform_rho) {
        std::copy(
            state.rho_hat_data(),
            state.rho_hat_data() + local_spectral_size_,
            spectral_pack_.data() + static_cast<std::size_t>(offsets.rho) * local_spectral_size_
        );
    }

    if (field_layout_.transform_j) {
        std::copy(
            state.jx_hat_data(),
            state.jx_hat_data() + local_spectral_size_,
            spectral_pack_.data() + static_cast<std::size_t>(offsets.jx) * local_spectral_size_
        );

        std::copy(
            state.jy_hat_data(),
            state.jy_hat_data() + local_spectral_size_,
            spectral_pack_.data() + static_cast<std::size_t>(offsets.jy) * local_spectral_size_
        );
    }

    if (field_layout_.ntransform > 0) {
        fourier_.backward_many(field_layout_.ntransform, spectral_pack_.data(), physical_pack_.data());
    }

    if (field_layout_.make_velocity) {
        const double* jx = physical_pack_.data() + static_cast<std::size_t>(offsets.jx) * local_physical_size_;
        const double* jy = physical_pack_.data() + static_cast<std::size_t>(offsets.jy) * local_physical_size_;
        double* vx = physical_pack_.data() + static_cast<std::size_t>(offsets.vx) * local_physical_size_;
        double* vy = physical_pack_.data() + static_cast<std::size_t>(offsets.vy) * local_physical_size_;

        if (is_compressible_mode(dynamics_mode_)) {
            const double* rho = physical_pack_.data() + static_cast<std::size_t>(offsets.rho) * local_physical_size_;

            for (std::size_t i = 0; i < local_physical_size_; ++i) {
                if (rho[i] == 0.0) {
                    throw std::runtime_error("compressible velocity encountered zero density.");
                }

                vx[i] = jx[i] / rho[i];
                vy[i] = jy[i] / rho[i];
            }
        }
        else if (is_incompressible_mode(dynamics_mode_)) {
            const double rho0 = reference_density(state);

            for (std::size_t i = 0; i < local_physical_size_; ++i) {
                vx[i] = jx[i] / rho0;
                vy[i] = jy[i] / rho0;
            }
        }
    }

    physical_fields_ready_ = true;
}
// ---------------------------------------------------------------------- //
const double* FCalculatorWorkspace::physical_psi_data() const {
    if (!physical_fields_ready_ || field_layout_.offsets.psi < 0) {
        throw std::runtime_error("physical psi is not ready.");
    }

    return physical_pack_.data() + static_cast<std::size_t>(field_layout_.offsets.psi) * local_physical_size_;
}
// ---------------------------------------------------------------------- //
const double* FCalculatorWorkspace::physical_rho_data() const {
    if (!physical_fields_ready_ || field_layout_.offsets.rho < 0) {
        throw std::runtime_error("physical rho is not ready.");
    }

    return physical_pack_.data() + static_cast<std::size_t>(field_layout_.offsets.rho) * local_physical_size_;
}
// ---------------------------------------------------------------------- //
const double* FCalculatorWorkspace::physical_jx_data() const {
    if (!physical_fields_ready_ || field_layout_.offsets.jx < 0) {
        throw std::runtime_error("physical jx is not ready.");
    }

    return physical_pack_.data() + static_cast<std::size_t>(field_layout_.offsets.jx) * local_physical_size_;
}
// ---------------------------------------------------------------------- //
const double* FCalculatorWorkspace::physical_jy_data() const {
    if (!physical_fields_ready_ || field_layout_.offsets.jy < 0) {
        throw std::runtime_error("physical jy is not ready.");
    }

    return physical_pack_.data() +  static_cast<std::size_t>(field_layout_.offsets.jy) * local_physical_size_;
}
// ---------------------------------------------------------------------- //
const double* FCalculatorWorkspace::physical_vx_data() const {
    if (!physical_fields_ready_ || field_layout_.offsets.vx < 0) {
        throw std::runtime_error("physical vx is not ready.");
    }

    return physical_pack_.data() + static_cast<std::size_t>(field_layout_.offsets.vx) * local_physical_size_;
}
// ---------------------------------------------------------------------- //
const double* FCalculatorWorkspace::physical_vy_data() const {
    if (!physical_fields_ready_ || field_layout_.offsets.vy < 0) {
        throw std::runtime_error("physical vy is not ready.");
    }

    return physical_pack_.data() + static_cast<std::size_t>(field_layout_.offsets.vy) * local_physical_size_;
}
// ---------------------------------------------------------------------- //
double* FCalculatorWorkspace::temporary_physical_fields(int nfields) {
    if (nfields < 0 || nfields > temporary_field_capacity_) {
        throw std::runtime_error(
            "FCalculatorWorkspace physical temporary field request exceeds capacity."
        );
    }

    const std::size_t used_size = static_cast<std::size_t>(nfields) * local_physical_size_;
    std::fill(temporary_physical_.begin(), temporary_physical_.begin() + used_size, 0.0);

    return temporary_physical_.data();
}
// ---------------------------------------------------------------------- //
Complex* FCalculatorWorkspace::temporary_spectral_fields(int nfields) {
    if (nfields < 0 || nfields > temporary_field_capacity_) {
        throw std::runtime_error(
            "FCalculatorWorkspace spectral temporary field request exceeds capacity."
        );
    }

    const std::size_t used_size = static_cast<std::size_t>(nfields) * local_spectral_size_;
    std::fill(temporary_spectral_.begin(), temporary_spectral_.begin() + used_size, Complex(0.0, 0.0));

    return temporary_spectral_.data();
}
// ---------------------------------------------------------------------- //
void FCalculatorWorkspace::forward_many(int nfields, const double* physical, Complex* spectral) {
    fourier_.forward_many(nfields, physical, spectral);
}
// ---------------------------------------------------------------------- //
void FCalculatorWorkspace::backward_many(int nfields, const Complex* spectral, double* physical) {
    fourier_.backward_many(nfields, spectral, physical);
}
// ---------------------------------------------------------------------- //
std::size_t FCalculatorWorkspace::local_physical_size() const {
    return local_physical_size_;
}
// ---------------------------------------------------------------------- //
std::size_t FCalculatorWorkspace::local_spectral_size() const {
    return local_spectral_size_;
}
// ---------------------------------------------------------------------- //
int FCalculatorWorkspace::temporary_field_capacity() const {
    return temporary_field_capacity_;
}
// ---------------------------------------------------------------------- //
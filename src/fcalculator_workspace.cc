#include "fcalculator_workspace.h"

#include <algorithm>
#include <stdexcept>
#include <string>

// ---------------------------------------------------------------------- //
FCalculatorWorkspace::FCalculatorWorkspace(
    const Domain2D& domain,
    const Params& params,
    const PhysicalRHSPlan& physical_rhs_plan,
    int temporary_field_capacity
)
    : domain_(domain),
      params_(params),
      fourier_(domain),
      physical_rhs_plan_(physical_rhs_plan),
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

    const std::string& type = params_.runtime.time_evolution_type;

    is_quiescent_ = type == "euler/quiescent" || type == "srk3/quiescent";
    is_incompressible_ = type == "euler/incompressible" || type == "srk3/incompressible";
    is_compressible_ = type == "euler/compressible" || type == "srk3/compressible";

    if (!is_quiescent_ && !is_incompressible_ && !is_compressible_) {
        throw std::runtime_error("unknown time_evolution type: " + type);
    }

    const bool plan_need_psi = physical_rhs_plan_.psi_det.need_psi || physical_rhs_plan_.j_det.need_psi;
    const bool plan_need_rho = physical_rhs_plan_.psi_det.need_rho || physical_rhs_plan_.j_det.need_rho;
    const bool plan_need_j = physical_rhs_plan_.psi_det.need_j || physical_rhs_plan_.j_det.need_j;
    const bool plan_need_velocity = physical_rhs_plan_.psi_det.need_velocity || physical_rhs_plan_.j_det.need_velocity;

    if (is_quiescent_ && (plan_need_rho || plan_need_j || plan_need_velocity)) {
        throw std::runtime_error("quiescent workspace cannot request hydrodynamic physical fields.");
    }

    if (plan_need_velocity && !plan_need_j) {
        throw std::runtime_error("velocity request requires physical momentum.");
    }

    if (!is_compressible_ && plan_need_rho) {
        throw std::runtime_error("density IFFT is only supported for compressible time evolution.");
    }

    request_.need_psi = plan_need_psi;
    request_.need_j = plan_need_j;
    request_.need_velocity = plan_need_velocity;
    request_.need_rho = is_compressible_ && (plan_need_rho || plan_need_velocity);

    if (request_.need_psi) {
        if (num_order_parameters_ <= 0) {
            throw std::runtime_error("physical psi requested with no order parameters.");
        }

        offsets_.psi = ntransform_;
        ntransform_ += num_order_parameters_;
    }

    if (request_.need_rho) {
        offsets_.rho = ntransform_;
        ntransform_ += 1;
    }

    if (request_.need_j) {
        offsets_.jx = ntransform_;
        offsets_.jy = ntransform_ + 1;
        ntransform_ += 2;
    }

    nphysical_ = ntransform_;

    if (request_.need_velocity) {
        offsets_.vx = nphysical_;
        offsets_.vy = nphysical_ + 1;
        nphysical_ += 2;
    }

    reference_density_needed_ = is_incompressible_ && request_.need_velocity;
    spectral_pack_.assign(static_cast<std::size_t>(ntransform_) * local_spectral_size_, Complex(0.0, 0.0));
    physical_pack_.assign(static_cast<std::size_t>(nphysical_) * local_physical_size_, 0.0);
    temporary_physical_.assign(static_cast<std::size_t>(temporary_field_capacity_) * local_physical_size_, 0.0);
    temporary_spectral_.assign(static_cast<std::size_t>(temporary_field_capacity_) * local_spectral_size_, Complex(0.0, 0.0));
}
// ---------------------------------------------------------------------- //
double FCalculatorWorkspace::reference_density(const State& state) {
    if (reference_density_ready_) {
        return reference_density_;
    }

    const double grid_size =
        static_cast<double>(domain_.nx_global())
      * static_cast<double>(domain_.ny_global());

    double local_rho0 = 0.0;

    const Box2D& box = domain_.spectral_box();
    if (box.low[0] <= 0 && 0 <= box.high[0]
        && box.low[1] <= 0 && 0 <= box.high[1]) {
        const std::size_t local_nkx = static_cast<std::size_t>(box.size[0]);
        const std::size_t index =
            static_cast<std::size_t>(0 - box.low[1]) * local_nkx
          + static_cast<std::size_t>(0 - box.low[0]);

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

    if (offsets_.psi >= 0) {
        std::copy(
            state.psi_hat_data(0),
            state.psi_hat_data(0) + static_cast<std::size_t>(num_order_parameters_) * local_spectral_size_,
            spectral_pack_.data() + static_cast<std::size_t>(offsets_.psi) * local_spectral_size_
        );
    }

    if (offsets_.rho >= 0) {
        std::copy(
            state.rho_hat_data(),
            state.rho_hat_data() + local_spectral_size_,
            spectral_pack_.data() + static_cast<std::size_t>(offsets_.rho) * local_spectral_size_
        );
    }

    if (offsets_.jx >= 0) {
        std::copy(
            state.jx_hat_data(),
            state.jx_hat_data() + local_spectral_size_,
            spectral_pack_.data() + static_cast<std::size_t>(offsets_.jx) * local_spectral_size_
        );

        std::copy(
            state.jy_hat_data(),
            state.jy_hat_data() + local_spectral_size_,
            spectral_pack_.data() + static_cast<std::size_t>(offsets_.jy) * local_spectral_size_
        );
    }

    if (ntransform_ > 0) {
        fourier_.backward_many(ntransform_, spectral_pack_.data(), physical_pack_.data());
    }

    if (request_.need_velocity) {
        const double* jx = physical_pack_.data() + static_cast<std::size_t>(offsets_.jx) * local_physical_size_;
        const double* jy = physical_pack_.data() + static_cast<std::size_t>(offsets_.jy) * local_physical_size_;
        double* vx = physical_pack_.data() + static_cast<std::size_t>(offsets_.vx) * local_physical_size_;
        double* vy = physical_pack_.data() + static_cast<std::size_t>(offsets_.vy) * local_physical_size_;

        if (is_compressible_) {
            const double* rho = physical_pack_.data() + static_cast<std::size_t>(offsets_.rho) * local_physical_size_;

            for (std::size_t i = 0; i < local_physical_size_; ++i) {
                if (rho[i] == 0.0) {
                    throw std::runtime_error("compressible velocity encountered zero density.");
                }

                vx[i] = jx[i] / rho[i];
                vy[i] = jy[i] / rho[i];
            }
        }
        else if (is_incompressible_) {
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
    if (!physical_fields_ready_ || offsets_.psi < 0) {
        throw std::runtime_error("physical psi is not ready.");
    }

    return physical_pack_.data() + static_cast<std::size_t>(offsets_.psi) * local_physical_size_;
}
// ---------------------------------------------------------------------- //
const double* FCalculatorWorkspace::physical_rho_data() const {
    if (!physical_fields_ready_ || offsets_.rho < 0) {
        throw std::runtime_error("physical rho is not ready.");
    }

    return physical_pack_.data() + static_cast<std::size_t>(offsets_.rho) * local_physical_size_;
}
// ---------------------------------------------------------------------- //
const double* FCalculatorWorkspace::physical_jx_data() const {
    if (!physical_fields_ready_ || offsets_.jx < 0) {
        throw std::runtime_error("physical jx is not ready.");
    }

    return physical_pack_.data() + static_cast<std::size_t>(offsets_.jx) * local_physical_size_;
}
// ---------------------------------------------------------------------- //
const double* FCalculatorWorkspace::physical_jy_data() const {
    if (!physical_fields_ready_ || offsets_.jy < 0) {
        throw std::runtime_error("physical jy is not ready.");
    }

    return physical_pack_.data() + static_cast<std::size_t>(offsets_.jy) * local_physical_size_;
}
// ---------------------------------------------------------------------- //
const double* FCalculatorWorkspace::physical_vx_data() const {
    if (!physical_fields_ready_ || offsets_.vx < 0) {
        throw std::runtime_error("physical vx is not ready.");
    }

    return physical_pack_.data() + static_cast<std::size_t>(offsets_.vx) * local_physical_size_;
}
// ---------------------------------------------------------------------- //
const double* FCalculatorWorkspace::physical_vy_data() const {
    if (!physical_fields_ready_ || offsets_.vy < 0) {
        throw std::runtime_error("physical vy is not ready.");
    }

    return physical_pack_.data() + static_cast<std::size_t>(offsets_.vy) * local_physical_size_;
}
// ---------------------------------------------------------------------- //
double* FCalculatorWorkspace::temporary_physical_fields(int nfields) {
    if (nfields < 0 || nfields > temporary_field_capacity_) {
        throw std::runtime_error("FCalculatorWorkspace physical temporary field request exceeds capacity.");
    }

    const std::size_t used_size = static_cast<std::size_t>(nfields) * local_physical_size_;

    std::fill(temporary_physical_.begin(), temporary_physical_.begin() + used_size, 0.0);
    return temporary_physical_.data();
}
// ---------------------------------------------------------------------- //
Complex* FCalculatorWorkspace::temporary_spectral_fields(int nfields) {
    if (nfields < 0 || nfields > temporary_field_capacity_) {
        throw std::runtime_error("FCalculatorWorkspace spectral temporary field request exceeds capacity.");
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
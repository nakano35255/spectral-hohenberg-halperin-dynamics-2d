#ifndef SFI_FCALCULATOR_H
#define SFI_FCALCULATOR_H

#include "domain.h"
#include "fourier_transform.h"
#include "model_thermodynamics.h"
#include "model_transport_coefficient.h"
#include "simulationinfo.h"
#include "state.h"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

class FCalculator {
private:
    const Params& params_;
    const Domain2D& domain_;
    const Thermodynamics& thermodynamics_;
    const TransportCoefficient& transport_coefficient_;

    FourierTransform2D fourier_;

    int num_components_ = 0;
    int num_fields_ = 0;
    std::size_t local_physical_size_ = 0;
    std::size_t local_spectral_size_ = 0;

    std::vector<double> rho_;
    std::vector<double> physical_fields_;
    std::vector<double> mu_;
    std::vector<double> pressure_;
    std::vector<double> grad_mu_x_;
    std::vector<double> grad_mu_y_;
    std::vector<double> diffusion_flux_x_;
    std::vector<double> diffusion_flux_y_;
    std::vector<double> velocity_x_;
    std::vector<double> velocity_y_;
    std::vector<double> grad_vx_x_;
    std::vector<double> grad_vx_y_;
    std::vector<double> grad_vy_x_;
    std::vector<double> grad_vy_y_;
    std::vector<double> stress_xx_;
    std::vector<double> stress_xy_;
    std::vector<double> stress_yy_;

    std::vector<Complex> mu_hat_;
    std::vector<Complex> pressure_hat_;
    std::vector<Complex> grad_mu_x_hat_;
    std::vector<Complex> grad_mu_y_hat_;
    std::vector<Complex> diffusion_flux_x_hat_;
    std::vector<Complex> diffusion_flux_y_hat_;
    std::vector<Complex> velocity_x_hat_;
    std::vector<Complex> velocity_y_hat_;
    std::vector<Complex> grad_vx_x_hat_;
    std::vector<Complex> grad_vx_y_hat_;
    std::vector<Complex> grad_vy_x_hat_;
    std::vector<Complex> grad_vy_y_hat_;
    std::vector<Complex> stress_xx_hat_;
    std::vector<Complex> stress_xy_hat_;
    std::vector<Complex> stress_yy_hat_;
    std::vector<Complex> density_rhs_;
    std::vector<Complex> momentum_rhs_x_;
    std::vector<Complex> momentum_rhs_y_;

    std::vector<double> rho_point_;
    std::vector<double> mu_point_;
    std::vector<double> mobility_;

    const Complex* cached_state_data_ = nullptr;
    double cached_time_ = 0.0;
    bool rho_cache_valid_ = false;
    bool momentum_physical_cache_valid_ = false;
    bool density_cache_valid_ = false;
    bool momentum_cache_valid_ = false;
    bool momentum_cache_uses_pressure_ = false;
    bool momentum_cache_uses_bulk_viscosity_ = false;

    std::size_t physical_offset(int field) const {
        return static_cast<std::size_t>(field) * local_physical_size_;
    }

    std::size_t spectral_offset(int field) const {
        return static_cast<std::size_t>(field) * local_spectral_size_;
    }

    void clear_density_output() {
        std::fill(density_rhs_.begin(), density_rhs_.end(), Complex(0.0, 0.0));
    }

    void clear_momentum_output() {
        std::fill(momentum_rhs_x_.begin(), momentum_rhs_x_.end(), Complex(0.0, 0.0));
        std::fill(momentum_rhs_y_.begin(), momentum_rhs_y_.end(), Complex(0.0, 0.0));
    }

    void multiply_gradient(const std::vector<Complex>& field_hat, std::vector<Complex>& grad_x_hat, std::vector<Complex>& grad_y_hat, int nfields) {
        std::fill(grad_x_hat.begin(), grad_x_hat.end(), Complex(0.0, 0.0));
        std::fill(grad_y_hat.begin(), grad_y_hat.end(), Complex(0.0, 0.0));

        const Box2D& box = domain_.spectral_box();
        for (int gy = box.low[1]; gy <= box.high[1]; ++gy) {
            const double ky = domain_.ky(gy);
            for (int gx = box.low[0]; gx <= box.high[0]; ++gx) {
                const double kx = domain_.kx(gx);
                const std::size_t index = domain_.spectral_local_index(gx, gy);
                const Complex ikx(0.0, kx);
                const Complex iky(0.0, ky);

                for (int field = 0; field < nfields; ++field) {
                    const std::size_t offset = spectral_offset(field);
                    grad_x_hat[offset + index] = ikx * field_hat[offset + index];
                    grad_y_hat[offset + index] = iky * field_hat[offset + index];
                }
            }
        }
    }

    void compute_chemical_potential() {
        for (std::size_t i = 0; i < local_physical_size_; ++i) {
            for (int component = 0; component < num_components_; ++component) {
                rho_point_[static_cast<std::size_t>(component)] = rho_[physical_offset(component) + i];
            }

            thermodynamics_.chemical_potential(rho_point_, mu_point_);

            for (int component = 0; component < num_components_; ++component) {
                mu_[physical_offset(component) + i] = mu_point_[static_cast<std::size_t>(component)];
            }
        }

        fourier_.forward_many(num_components_, mu_.data(), mu_hat_.data());
    }

    void compute_pressure() {
        for (std::size_t i = 0; i < local_physical_size_; ++i) {
            for (int component = 0; component < num_components_; ++component) {
                rho_point_[static_cast<std::size_t>(component)] = rho_[physical_offset(component) + i];
            }

            pressure_[i] = thermodynamics_.pressure(rho_point_);
        }

        fourier_.forward(pressure_.data(), pressure_hat_.data());
    }

    void compute_diffusion_rhs() {
        multiply_gradient(mu_hat_, grad_mu_x_hat_, grad_mu_y_hat_, num_components_);
        fourier_.backward_many(num_components_, grad_mu_x_hat_.data(), grad_mu_x_.data());
        fourier_.backward_many(num_components_, grad_mu_y_hat_.data(), grad_mu_y_.data());

        std::fill(diffusion_flux_x_.begin(), diffusion_flux_x_.end(), 0.0);
        std::fill(diffusion_flux_y_.begin(), diffusion_flux_y_.end(), 0.0);

        for (std::size_t i = 0; i < local_physical_size_; ++i) {
            for (int component = 0; component < num_components_; ++component) {
                rho_point_[static_cast<std::size_t>(component)] = rho_[physical_offset(component) + i];
            }

            transport_coefficient_.mobility(rho_point_, mobility_);

            for (int row = 0; row < num_components_; ++row) {
                double flux_x = 0.0;
                double flux_y = 0.0;
                for (int col = 0; col < num_components_; ++col) {
                    const double L = mobility_[static_cast<std::size_t>(row) * static_cast<std::size_t>(num_components_) + static_cast<std::size_t>(col)];
                    flux_x += L * grad_mu_x_[physical_offset(col) + i];
                    flux_y += L * grad_mu_y_[physical_offset(col) + i];
                }
                diffusion_flux_x_[physical_offset(row) + i] = flux_x;
                diffusion_flux_y_[physical_offset(row) + i] = flux_y;
            }
        }

        fourier_.forward_many(num_components_, diffusion_flux_x_.data(), diffusion_flux_x_hat_.data());
        fourier_.forward_many(num_components_, diffusion_flux_y_.data(), diffusion_flux_y_hat_.data());

        const Box2D& box = domain_.spectral_box();
        for (int gy = box.low[1]; gy <= box.high[1]; ++gy) {
            const double ky = domain_.ky(gy);
            for (int gx = box.low[0]; gx <= box.high[0]; ++gx) {
                const double kx = domain_.kx(gx);
                const std::size_t index = domain_.spectral_local_index(gx, gy);
                const Complex ikx(0.0, kx);
                const Complex iky(0.0, ky);

                for (int component = 0; component < num_components_; ++component) {
                    const std::size_t offset = spectral_offset(component);
                    density_rhs_[offset + index] =
                        ikx * diffusion_flux_x_hat_[offset + index]
                      + iky * diffusion_flux_y_hat_[offset + index];
                }
            }
        }
    }

    void compute_velocity() {
        const std::size_t jx_offset = physical_offset(num_components_);
        const std::size_t jy_offset = physical_offset(num_components_ + 1);

        for (std::size_t i = 0; i < local_physical_size_; ++i) {
            double rho_total = 0.0;
            for (int component = 0; component < num_components_; ++component) {
                rho_total += rho_[physical_offset(component) + i];
            }

            if (rho_total <= 0.0) {
                throw std::runtime_error("FCalculator requires positive total density.");
            }

            velocity_x_[i] = physical_fields_[jx_offset + i] / rho_total;
            velocity_y_[i] = physical_fields_[jy_offset + i] / rho_total;
        }

        fourier_.forward(velocity_x_.data(), velocity_x_hat_.data());
        fourier_.forward(velocity_y_.data(), velocity_y_hat_.data());
    }

    void compute_viscous_rhs(bool use_bulk_viscosity) {
        multiply_gradient(velocity_x_hat_, grad_vx_x_hat_, grad_vx_y_hat_, 1);
        multiply_gradient(velocity_y_hat_, grad_vy_x_hat_, grad_vy_y_hat_, 1);

        fourier_.backward(grad_vx_x_hat_.data(), grad_vx_x_.data());
        fourier_.backward(grad_vx_y_hat_.data(), grad_vx_y_.data());
        fourier_.backward(grad_vy_x_hat_.data(), grad_vy_x_.data());
        fourier_.backward(grad_vy_y_hat_.data(), grad_vy_y_.data());

        for (std::size_t i = 0; i < local_physical_size_; ++i) {
            for (int component = 0; component < num_components_; ++component) {
                rho_point_[static_cast<std::size_t>(component)] = rho_[physical_offset(component) + i];
            }

            const double eta = transport_coefficient_.shear_viscosity(rho_point_);
            const double zeta = use_bulk_viscosity ? transport_coefficient_.bulk_viscosity(rho_point_) : 0.0;
            const double div_v = grad_vx_x_[i] + grad_vy_y_[i];

            stress_xx_[i] = 2.0 * eta * grad_vx_x_[i] + (zeta - eta) * div_v;
            stress_yy_[i] = 2.0 * eta * grad_vy_y_[i] + (zeta - eta) * div_v;
            stress_xy_[i] = eta * (grad_vx_y_[i] + grad_vy_x_[i]);
        }

        fourier_.forward(stress_xx_.data(), stress_xx_hat_.data());
        fourier_.forward(stress_xy_.data(), stress_xy_hat_.data());
        fourier_.forward(stress_yy_.data(), stress_yy_hat_.data());

        const Box2D& box = domain_.spectral_box();
        for (int gy = box.low[1]; gy <= box.high[1]; ++gy) {
            const double ky = domain_.ky(gy);
            for (int gx = box.low[0]; gx <= box.high[0]; ++gx) {
                const double kx = domain_.kx(gx);
                const std::size_t index = domain_.spectral_local_index(gx, gy);
                const Complex ikx(0.0, kx);
                const Complex iky(0.0, ky);

                momentum_rhs_x_[index] += ikx * stress_xx_hat_[index] + iky * stress_xy_hat_[index];
                momentum_rhs_y_[index] += ikx * stress_xy_hat_[index] + iky * stress_yy_hat_[index];
            }
        }
    }

    void compute_pressure_rhs() {
        const Box2D& box = domain_.spectral_box();
        for (int gy = box.low[1]; gy <= box.high[1]; ++gy) {
            const double ky = domain_.ky(gy);
            for (int gx = box.low[0]; gx <= box.high[0]; ++gx) {
                const double kx = domain_.kx(gx);
                const std::size_t index = domain_.spectral_local_index(gx, gy);
                const Complex ikx(0.0, kx);
                const Complex iky(0.0, ky);

                momentum_rhs_x_[index] += -ikx * pressure_hat_[index];
                momentum_rhs_y_[index] += -iky * pressure_hat_[index];
            }
        }
    }

    void ensure_state_cache_identity(const State& current, double t) {
        if (cached_state_data_ == current.data() && cached_time_ == t) {
            return;
        }

        cached_state_data_ = current.data();
        cached_time_ = t;
        rho_cache_valid_ = false;
        momentum_physical_cache_valid_ = false;
        density_cache_valid_ = false;
        momentum_cache_valid_ = false;
    }

    void ensure_rho_cache(const State& current, double t) {
        ensure_state_cache_identity(current, t);
        if (rho_cache_valid_) {
            return;
        }

        fourier_.backward_many(num_components_, current.data(), rho_.data());
        rho_cache_valid_ = true;
    }

    void ensure_momentum_physical_cache(const State& current, double t) {
        ensure_rho_cache(current, t);
        if (momentum_physical_cache_valid_) {
            return;
        }

        fourier_.backward_many(
            2,
            current.data() + spectral_offset(num_components_),
            physical_fields_.data() + physical_offset(num_components_)
        );
        momentum_physical_cache_valid_ = true;
    }

    void ensure_density_cache(const State& current, double t) {
        ensure_rho_cache(current, t);
        if (density_cache_valid_) {
            return;
        }

        clear_density_output();
        compute_chemical_potential();
        compute_diffusion_rhs();

        density_cache_valid_ = true;
    }

    void ensure_momentum_cache(const State& current, double t, bool use_pressure, bool use_bulk_viscosity) {
        ensure_momentum_physical_cache(current, t);
        if (momentum_cache_valid_
            && momentum_cache_uses_pressure_ == use_pressure
            && momentum_cache_uses_bulk_viscosity_ == use_bulk_viscosity) {
            return;
        }

        clear_momentum_output();
        compute_velocity();
        if (use_pressure) {
            compute_pressure();
            compute_pressure_rhs();
        }
        compute_viscous_rhs(use_bulk_viscosity);

        momentum_cache_valid_ = true;
        momentum_cache_uses_pressure_ = use_pressure;
        momentum_cache_uses_bulk_viscosity_ = use_bulk_viscosity;
    }

public:
    FCalculator(
        const Params& params,
        const Domain2D& domain,
        const Thermodynamics& thermodynamics,
        const TransportCoefficient& transport_coefficient
    )
        : params_(params),
          domain_(domain),
          thermodynamics_(thermodynamics),
          transport_coefficient_(transport_coefficient),
          fourier_(domain),
          num_components_(params.physics.num_components),
          num_fields_(params.physics.num_components + 2),
          local_physical_size_(domain.physical_size()),
          local_spectral_size_(domain.spectral_size()),
          rho_(static_cast<std::size_t>(num_components_) * local_physical_size_, 0.0),
          physical_fields_(static_cast<std::size_t>(num_fields_) * local_physical_size_, 0.0),
          mu_(static_cast<std::size_t>(num_components_) * local_physical_size_, 0.0),
          pressure_(local_physical_size_, 0.0),
          grad_mu_x_(static_cast<std::size_t>(num_components_) * local_physical_size_, 0.0),
          grad_mu_y_(static_cast<std::size_t>(num_components_) * local_physical_size_, 0.0),
          diffusion_flux_x_(static_cast<std::size_t>(num_components_) * local_physical_size_, 0.0),
          diffusion_flux_y_(static_cast<std::size_t>(num_components_) * local_physical_size_, 0.0),
          velocity_x_(local_physical_size_, 0.0),
          velocity_y_(local_physical_size_, 0.0),
          grad_vx_x_(local_physical_size_, 0.0),
          grad_vx_y_(local_physical_size_, 0.0),
          grad_vy_x_(local_physical_size_, 0.0),
          grad_vy_y_(local_physical_size_, 0.0),
          stress_xx_(local_physical_size_, 0.0),
          stress_xy_(local_physical_size_, 0.0),
          stress_yy_(local_physical_size_, 0.0),
          mu_hat_(static_cast<std::size_t>(num_components_) * local_spectral_size_, Complex(0.0, 0.0)),
          pressure_hat_(local_spectral_size_, Complex(0.0, 0.0)),
          grad_mu_x_hat_(static_cast<std::size_t>(num_components_) * local_spectral_size_, Complex(0.0, 0.0)),
          grad_mu_y_hat_(static_cast<std::size_t>(num_components_) * local_spectral_size_, Complex(0.0, 0.0)),
          diffusion_flux_x_hat_(static_cast<std::size_t>(num_components_) * local_spectral_size_, Complex(0.0, 0.0)),
          diffusion_flux_y_hat_(static_cast<std::size_t>(num_components_) * local_spectral_size_, Complex(0.0, 0.0)),
          velocity_x_hat_(local_spectral_size_, Complex(0.0, 0.0)),
          velocity_y_hat_(local_spectral_size_, Complex(0.0, 0.0)),
          grad_vx_x_hat_(local_spectral_size_, Complex(0.0, 0.0)),
          grad_vx_y_hat_(local_spectral_size_, Complex(0.0, 0.0)),
          grad_vy_x_hat_(local_spectral_size_, Complex(0.0, 0.0)),
          grad_vy_y_hat_(local_spectral_size_, Complex(0.0, 0.0)),
          stress_xx_hat_(local_spectral_size_, Complex(0.0, 0.0)),
          stress_xy_hat_(local_spectral_size_, Complex(0.0, 0.0)),
          stress_yy_hat_(local_spectral_size_, Complex(0.0, 0.0)),
          density_rhs_(static_cast<std::size_t>(num_components_) * local_spectral_size_, Complex(0.0, 0.0)),
          momentum_rhs_x_(local_spectral_size_, Complex(0.0, 0.0)),
          momentum_rhs_y_(local_spectral_size_, Complex(0.0, 0.0)),
          rho_point_(static_cast<std::size_t>(num_components_), 0.0),
          mu_point_(static_cast<std::size_t>(num_components_), 0.0),
          mobility_(static_cast<std::size_t>(num_components_) * static_cast<std::size_t>(num_components_), 0.0) {
        if (num_components_ <= 0) {
            throw std::runtime_error("FCalculator requires a positive number of components.");
        }
    }

    void density_det(int component, const State& current, Complex* out, double t) {
        if (component < 0 || component >= num_components_) {
            throw std::out_of_range("FCalculator density component out of range.");
        }
        ensure_density_cache(current, t);

        const Complex* rhs = density_rhs_.data() + spectral_offset(component);
        std::copy(rhs, rhs + local_spectral_size_, out);
    }

    void momentum_det(const State& current, Complex* out_jx, Complex* out_jy, double t) {
        ensure_momentum_cache(current, t, true, true);

        std::copy(momentum_rhs_x_.begin(), momentum_rhs_x_.end(), out_jx);
        std::copy(momentum_rhs_y_.begin(), momentum_rhs_y_.end(), out_jy);
    }

    void momentum_det_incompressible(const State& current, Complex* out_jx, Complex* out_jy, double t) {
        ensure_momentum_cache(current, t, false, false);

        std::copy(momentum_rhs_x_.begin(), momentum_rhs_x_.end(), out_jx);
        std::copy(momentum_rhs_y_.begin(), momentum_rhs_y_.end(), out_jy);
    }
};

#endif

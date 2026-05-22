#include "model_transport_coefficient_constant.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
    std::size_t matrix_index(int row, int col, int num_cols) {
        return static_cast<std::size_t>(row) * static_cast<std::size_t>(num_cols) + static_cast<std::size_t>(col);
    }

    double matrix_scale(const std::vector<double>& matrix) {
        double scale = 0.0;
        for (double value : matrix) {
            scale = std::max(scale, std::abs(value));
        }
        return std::max(1.0, scale);
    }

    std::vector<double> make_composition_basis(int num_components) {
        const int reduced_size = num_components - 1;
        std::vector<double> basis(static_cast<std::size_t>(num_components) * static_cast<std::size_t>(reduced_size), 0.0);

        for (int col = 0; col < reduced_size; ++col) {
            const double denominator = std::sqrt(static_cast<double>((col + 1) * (col + 2)));
            for (int row = 0; row <= col; ++row) {
                basis[matrix_index(row, col, reduced_size)] = 1.0 / denominator;
            }
            basis[matrix_index(col + 1, col, reduced_size)] = -static_cast<double>(col + 1) / denominator;
        }

        return basis;
    }

    std::vector<double> project_to_composition_space(
        const std::vector<double>& matrix,
        const std::vector<double>& basis,
        int num_components
    ) {
        const int reduced_size = num_components - 1;
        std::vector<double> reduced(static_cast<std::size_t>(reduced_size) * static_cast<std::size_t>(reduced_size), 0.0);

        for (int a = 0; a < reduced_size; ++a) {
            for (int b = 0; b < reduced_size; ++b) {
                double value = 0.0;
                for (int i = 0; i < num_components; ++i) {
                    for (int j = 0; j < num_components; ++j) {
                        value += basis[matrix_index(i, a, reduced_size)]
                               * matrix[matrix_index(i, j, num_components)]
                               * basis[matrix_index(j, b, reduced_size)];
                    }
                }
                reduced[matrix_index(a, b, reduced_size)] = value;
            }
        }

        return reduced;
    }

    std::vector<double> cholesky_semidefinite(
        const std::vector<double>& matrix,
        int size,
        double tolerance
    ) {
        std::vector<double> factor(static_cast<std::size_t>(size) * static_cast<std::size_t>(size), 0.0);

        for (int col = 0; col < size; ++col) {
            double diagonal = matrix[matrix_index(col, col, size)];
            for (int k = 0; k < col; ++k) {
                const double value = factor[matrix_index(col, k, size)];
                diagonal -= value * value;
            }

            if (diagonal < -tolerance) {
                throw std::runtime_error("constant transport coefficient mobility matrix is not positive semidefinite.");
            }

            if (diagonal <= tolerance) {
                factor[matrix_index(col, col, size)] = 0.0;
                for (int row = col + 1; row < size; ++row) {
                    double residual = matrix[matrix_index(row, col, size)];
                    for (int k = 0; k < col; ++k) {
                        residual -= factor[matrix_index(row, k, size)] * factor[matrix_index(col, k, size)];
                    }
                    if (std::abs(residual) > tolerance) {
                        throw std::runtime_error("constant transport coefficient mobility matrix is not positive semidefinite.");
                    }
                }
                continue;
            }

            factor[matrix_index(col, col, size)] = std::sqrt(diagonal);
            for (int row = col + 1; row < size; ++row) {
                double value = matrix[matrix_index(row, col, size)];
                for (int k = 0; k < col; ++k) {
                    value -= factor[matrix_index(row, k, size)] * factor[matrix_index(col, k, size)];
                }
                factor[matrix_index(row, col, size)] = value / factor[matrix_index(col, col, size)];
            }
        }

        return factor;
    }

    std::vector<double> lift_noise_factor(
        const std::vector<double>& basis,
        const std::vector<double>& reduced_factor,
        int num_components
    ) {
        const int reduced_size = num_components - 1;
        std::vector<double> factor(static_cast<std::size_t>(num_components) * static_cast<std::size_t>(reduced_size), 0.0);

        for (int row = 0; row < num_components; ++row) {
            for (int col = 0; col < reduced_size; ++col) {
                double value = 0.0;
                for (int k = 0; k < reduced_size; ++k) {
                    value += basis[matrix_index(row, k, reduced_size)]
                           * reduced_factor[matrix_index(k, col, reduced_size)];
                }
                factor[matrix_index(row, col, reduced_size)] = value;
            }
        }

        return factor;
    }
}

ConstantTransportCoefficient::ConstantTransportCoefficient(
    int num_components,
    double shear_viscosity,
    double bulk_viscosity,
    std::vector<double> mobility
) : num_components_(num_components),
    shear_viscosity_(shear_viscosity),
    bulk_viscosity_(bulk_viscosity),
    mobility_(std::move(mobility)) {
    if (num_components_ <= 0) {
        throw std::runtime_error("constant transport coefficient requires positive number of components.");
    }
    if (shear_viscosity_ < 0.0) {
        throw std::runtime_error("constant transport coefficient requires nonnegative eta.");
    }
    if (bulk_viscosity_ < 0.0) {
        throw std::runtime_error("constant transport coefficient requires nonnegative zeta.");
    }

    const std::size_t expected_size = static_cast<std::size_t>(num_components_) * static_cast<std::size_t>(num_components_);
    if (mobility_.size() != expected_size) {
        throw std::runtime_error("constant transport coefficient mobility matrix size mismatch.");
    }

    build_noise_factor();
}

void ConstantTransportCoefficient::check_density_size(const std::vector<double>& rho) const {
    if (rho.size() != static_cast<std::size_t>(num_components_)) {
        throw std::runtime_error("constant transport coefficient expected " + std::to_string(num_components_) + " density components, got " + std::to_string(rho.size()));
    }
}

void ConstantTransportCoefficient::build_noise_factor() {
    const double scale = matrix_scale(mobility_);
    const double tolerance = 1.0e-10 * scale * static_cast<double>(num_components_);

    for (int col = 0; col < num_components_; ++col) {
        double column_sum = 0.0;
        for (int row = 0; row < num_components_; ++row) {
            column_sum += mobility_[matrix_index(row, col, num_components_)];
        }

        if (std::abs(column_sum) > tolerance) {
            throw std::runtime_error("constant transport coefficient mobility matrix must have zero column sums; column " + std::to_string(col) + " sum is " + std::to_string(column_sum));
        }
    }

    for (int row = 0; row < num_components_; ++row) {
        for (int col = row + 1; col < num_components_; ++col) {
            const double left = mobility_[matrix_index(row, col, num_components_)];
            const double right = mobility_[matrix_index(col, row, num_components_)];
            if (std::abs(left - right) > tolerance) {
                throw std::runtime_error("constant transport coefficient mobility matrix must be symmetric; L" + std::to_string(row) + std::to_string(col) + " != L" + std::to_string(col) + std::to_string(row));
            }
        }
    }

    const int reduced_size = num_components_ - 1;
    if (reduced_size == 0) {
        noise_factor_.clear();
    } else {
        const std::vector<double> basis = make_composition_basis(num_components_);
        const std::vector<double> reduced = project_to_composition_space(mobility_, basis, num_components_);
        const std::vector<double> reduced_factor = cholesky_semidefinite(reduced, reduced_size, tolerance);
        noise_factor_ = lift_noise_factor(basis, reduced_factor, num_components_);
    }

    for (int row = 0; row < num_components_; ++row) {
        for (int col = 0; col < num_components_; ++col) {
            double reconstructed = 0.0;
            for (int k = 0; k < reduced_size; ++k) {
                reconstructed += noise_factor_[matrix_index(row, k, reduced_size)]
                               * noise_factor_[matrix_index(col, k, reduced_size)];
            }

            const double expected = mobility_[matrix_index(row, col, num_components_)];
            if (std::abs(reconstructed - expected) > 10.0 * tolerance) {
                throw std::runtime_error("constant transport coefficient mobility noise factor does not reconstruct the mobility matrix.");
            }
        }
    }
}

double ConstantTransportCoefficient::shear_viscosity(const std::vector<double>& rho) const {
    check_density_size(rho);
    return shear_viscosity_;
}

double ConstantTransportCoefficient::bulk_viscosity(const std::vector<double>& rho) const {
    check_density_size(rho);
    return bulk_viscosity_;
}

void ConstantTransportCoefficient::mobility(const std::vector<double>& rho, std::vector<double>& L) const {
    check_density_size(rho);
    L = mobility_;
}

void ConstantTransportCoefficient::noise_factor(const std::vector<double>& rho, std::vector<double>& B) const {
    check_density_size(rho);
    B = noise_factor_;
}

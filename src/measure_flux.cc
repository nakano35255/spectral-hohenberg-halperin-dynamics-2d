#include "measure_flux.h"
#include "measure_flux_style.h"

#include <iomanip>
#include <mpi.h>
#include <stdexcept>
#include <vector>
#include <optional>

// ---------------------------------------------------------------------- //
FluxMeasure::FluxMeasure(
    const Params& params,
    const Domain2D& domain,
    std::shared_ptr<const MeasureCommandBase> command
) : Measure(params, domain, command),
    num_order_parameters_(params.physics.num_order_parameters)
{
    if (num_order_parameters_ < 0) {
        throw std::runtime_error("FluxMeasure requires a nonnegative number of order parameters.");
    }

    auto cfg = std::dynamic_pointer_cast<const FluxMeasureCommand>(command);
    if (!cfg) {
        throw std::runtime_error("FluxMeasure: invalid command type.");
    }

    nevery_ = cfg->nevery;
    file_ = cfg->file;
    request_ = cfg->request;
}
// ---------------------------------------------------------------------- //
FluxRequest FluxMeasure::flux_request() const {
    return request_;
}
// ---------------------------------------------------------------------- //
void FluxMeasure::open_output_if_needed() {
    int open_ok = 1;

    if (domain_.rank() == 0 && !output_.is_open()) {
        output_.open(file_);
        if (!output_) {
            open_ok = 0;
        }
        else {
            output_ << "# step time";
            for (const std::string& name : selected_field_names()) {
                output_ << ' ' << name;
            }
            output_ << "\n";
            output_ << std::scientific << std::setprecision(16);
        }
    }

    MPI_Bcast(&open_ok, 1, MPI_INT, 0, domain_.comm());
    if (open_ok == 0) {
        throw std::runtime_error("FluxMeasure: cannot open file: " + file_);
    }
}
// ---------------------------------------------------------------------- //
std::vector<std::string> FluxMeasure::selected_field_names() const {
    std::vector<std::string> names;

    if (request_.density) {
        names.push_back("density_flux_x");
        names.push_back("density_flux_y");
    }

    if (request_.order_parameter) {
        for (int op = 0; op < num_order_parameters_; ++op) {
            names.push_back("order_parameter" + std::to_string(op) + "_flux_x");
            names.push_back("order_parameter" + std::to_string(op) + "_flux_y");
        }
    }

    if (request_.momentum) {
        names.push_back("momentum_flux_xx");
        names.push_back("momentum_flux_xy");
        names.push_back("momentum_flux_yy");
    }

    return names;
}
// ---------------------------------------------------------------------- //
std::vector<FluxMeasure::FluxField> FluxMeasure::selected_fields(const FluxBuffer& flux) const {
    std::vector<FluxField> fields;

    if (request_.density) {
        fields.push_back({"density_flux_x", flux.density_flux_x_hat_data()});
        fields.push_back({"density_flux_y", flux.density_flux_y_hat_data()});
    }

    if (request_.order_parameter) {
        for (int op = 0; op < num_order_parameters_; ++op) {
            fields.push_back({"order_parameter" + std::to_string(op) + "_flux_x", flux.order_parameter_flux_x_hat_data(op)});
            fields.push_back({"order_parameter" + std::to_string(op) + "_flux_y", flux.order_parameter_flux_y_hat_data(op)});
        }
    }

    if (request_.momentum) {
        fields.push_back({"momentum_flux_xx", flux.momentum_flux_xx_hat_data()});
        fields.push_back({"momentum_flux_xy", flux.momentum_flux_xy_hat_data()});
        fields.push_back({"momentum_flux_yy", flux.momentum_flux_yy_hat_data()});
    }

    return fields;
}
// ---------------------------------------------------------------------- //
std::optional<std::size_t> FluxMeasure::zero_mode_index() const {
    const Box2D& local_box = domain_.spectral_box();

    if (local_box.low[0] <= 0 && 0 <= local_box.high[0] && local_box.low[1] <= 0 && 0 <= local_box.high[1]) {
        const std::size_t local_nkx = static_cast<std::size_t>(local_box.size[0]);
        return static_cast<std::size_t>(0 - local_box.low[1]) * local_nkx
             + static_cast<std::size_t>(0 - local_box.low[0]);
    }

    return std::nullopt;
}
// ---------------------------------------------------------------------- //
void FluxMeasure::observe(
    const State& /*state*/,
    FourierTransform2D& /*fft*/,
    MeasureWorkspace& /*workspace*/,
    const FluxBuffer& flux,
    int step,
    double time
) {
    if (nevery_ <= 0 || step % nevery_ != 0) {
        return;
    }

    open_output_if_needed();

    const std::vector<FluxField> fields = selected_fields(flux);
    const std::size_t nfields = fields.size();
    const double grid_size = static_cast<double>(domain_.nx_global()) * static_cast<double>(domain_.ny_global());

    std::vector<double> local_values(nfields, 0.0);
    const auto zero_index = zero_mode_index();

    if (zero_index) {
        for (std::size_t field = 0; field < nfields; ++field) {
            local_values[field] = fields[field].data[*zero_index].real() / grid_size;
        }
    }

    std::vector<double> global_values(nfields, 0.0);
    MPI_Reduce(
        local_values.data(), global_values.data(), static_cast<int>(nfields), MPI_DOUBLE, MPI_SUM,
        0, domain_.comm()
    );

    int write_ok = 1;

    if (domain_.rank() == 0) {
        output_ << step << ' ' << time;
        for (double value : global_values) {
            output_ << ' ' << value;
        }
        output_ << '\n';

        if (!output_) {
            write_ok = 0;
        }
    }

    MPI_Bcast(&write_ok, 1, MPI_INT, 0, domain_.comm());
    if (write_ok == 0) {
        throw std::runtime_error("FluxMeasure: failed to write file: " + file_);
    }
}
// ---------------------------------------------------------------------- //
void FluxMeasure::finalize() {
    if (output_.is_open()) {
        output_.close();
    }
}
// ---------------------------------------------------------------------- //

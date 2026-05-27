#include "initial_condition_uniform_density.h"
#include "initial_condition_uniform_density_style.h"

#include <algorithm>

UniformDensityInitialCondition::UniformDensityInitialCondition(const Params& params, std::shared_ptr<const DensityInitialConditionCommandBase> command) :
    DensityInitialCondition(params, command)
{
    auto cfg = std::dynamic_pointer_cast<const UniformDensityInitialConditionCommand>(command);
    if (!cfg) {
        throw std::runtime_error("UniformDensityInitialCondition: invalid command type.");
    }

    value_ = cfg->value;
}

void UniformDensityInitialCondition::apply(
    State& state,
    const Domain2D& domain,
    const SpectralMask2D& spectral_mask
) const {
    const double amplitude = value_ * static_cast<double>(domain.nx_global()) * static_cast<double>(domain.ny_global());

    Complex* rho = state.rho_hat_data();
    std::fill(rho, rho + domain.spectral_size(), Complex(0.0, 0.0));

    for (const SpectralMode2D& mode : spectral_mask.active_modes()) {
        if (mode.gx == 0 && mode.gy == 0) {
            rho[mode.index] = Complex(amplitude, 0.0);
            break;
        }
    }
}

#include "measure_manager.h"

#include <algorithm>
#include <stdexcept>

MeasurementManager::MeasurementManager(
    const Params& params,
    const Domain2D& domain,
    const MeasureRegistry& registry,
    const Thermodynamics& thermodynamics,
    const FreeEnergy& free_energy,
    const TransportCoefficient& transport_coefficient
) : params_(params),
    domain_(domain),
    thermodynamics_(thermodynamics),
    free_energy_(free_energy),
    transport_coefficient_(transport_coefficient),
    registry_(registry),
    workspace_(domain, params) {}

void MeasurementManager::apply_measure_command(std::shared_ptr<MeasureCommandBase> command) {
    if (!command) {
        throw std::runtime_error("MeasurementManager received a null measure command.");
    }

    remove_measure(command->id);

    if (!command->enabled) {
        return;
    }

    const MeasureStyle& style = registry_.get(command->type);
    measures_.push_back(style.create_measure(params_, domain_, thermodynamics_, free_energy_, transport_coefficient_, command));
}

void MeasurementManager::remove_measure(const std::string& id) {
    auto it = std::find_if(measures_.begin(), measures_.end(), [&](const std::unique_ptr<Measure>& measure) {
            return measure->id() == id;
        }
    );

    if (it == measures_.end()) {
        return;
    }

    (*it)->finalize();
    measures_.erase(it);
}

void MeasurementManager::observe(const State& state, FourierTransform2D& fft, int step, double time) {
    workspace_.begin_step(step);

    for (const auto& measure : measures_) {
        measure->observe(state, fft, workspace_, step, time);
    }
}

void MeasurementManager::finalize() {
    for (const auto& measure : measures_) {
        measure->finalize();
    }
    measures_.clear();
}

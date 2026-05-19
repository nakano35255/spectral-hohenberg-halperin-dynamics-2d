#include "measure_manager.h"

#include <algorithm>
#include <stdexcept>

MeasurementManager::MeasurementManager(const Params& params, const MeasureRegistry& registry)
    : params_(params),
      registry_(registry) {}

void MeasurementManager::apply_measure_command(std::shared_ptr<MeasureCommandBase> command) {
    if (!command) {
        throw std::runtime_error("MeasurementManager received a null measure command.");
    }

    remove_measure(command->id);

    if (!command->enabled) {
        return;
    }

    const MeasureStyle& style = registry_.get(command->type);
    measures_.push_back(style.create_measure(params_, command));
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

void MeasurementManager::observe(const State& state, PhysicalStateBuffer& physical, FourierTransform2D& fft, const Domain2D& domain, int step, double time
) {
    for (const auto& measure : measures_) {
        measure->observe(state, physical, fft, domain, step, time);
    }
}

void MeasurementManager::finalize() {
    for (const auto& measure : measures_) {
        measure->finalize();
    }
    measures_.clear();
}

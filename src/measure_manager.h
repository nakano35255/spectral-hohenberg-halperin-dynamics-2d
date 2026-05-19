#ifndef MEASURE_MANAGER_H
#define MEASURE_MANAGER_H

#include "buffer_physical_state.h"
#include "domain.h"
#include "fourier_transform.h"
#include "measure.h"
#include "measure_registry.h"
#include "simulationinfo.h"
#include "state.h"

#include <memory>
#include <string>
#include <vector>

class MeasurementManager {
private:
    const Params& params_;
    const MeasureRegistry& registry_;
    std::vector<std::unique_ptr<Measure>> measures_;

public:
    MeasurementManager(const Params& params, const MeasureRegistry& registry);

    void apply_measure_command(std::shared_ptr<MeasureCommandBase> command);
    void remove_measure(const std::string& id);

    void observe(
        const State& state,
        PhysicalStateBuffer& physical,
        FourierTransform2D& fft,
        const Domain2D& domain,
        int step,
        double time
    );

    void finalize();
};

#endif

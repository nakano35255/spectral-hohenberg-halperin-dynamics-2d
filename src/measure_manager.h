#ifndef SHHD_MEASURE_MANAGER_H
#define SHHD_MEASURE_MANAGER_H

#include "buffer_physical_state.h"
#include "domain.h"
#include "fourier_transform.h"
#include "measure.h"
#include "measure_registry.h"
#include "simulationinfo.h"
#include "state.h"
#include "model_free_energy.h"
#include "model_thermodynamics.h"
#include "model_transport_coefficient.h"
#include "measure_workspace.h"

#include <memory>
#include <string>
#include <vector>

class MeasurementManager {
private:
    const Params& params_;
    const Domain2D& domain_;
    const Thermodynamics& thermodynamics_;
    const FreeEnergy& free_energy_;
    const TransportCoefficient& transport_coefficient_;

    const MeasureRegistry& registry_;
    MeasureWorkspace workspace_;

    std::vector<std::unique_ptr<Measure>> measures_;

public:
    MeasurementManager(
        const Params& params,
        const Domain2D& domain,
        const MeasureRegistry& registry,
        const Thermodynamics& thermodynamics,
        const FreeEnergy& free_energy,
        const TransportCoefficient& transport_coefficient
    );

    void apply_measure_command(std::shared_ptr<MeasureCommandBase> command);
    void remove_measure(const std::string& id);

    FluxRequest flux_request() const;

    void observe(const State& state, FourierTransform2D& fft, const FluxBuffer& flux, int step, double time);

    void finalize();
};

#endif

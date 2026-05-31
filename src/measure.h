#ifndef SHHD_MEASURE_H
#define SHHD_MEASURE_H

#include "simulationinfo.h"
#include "domain.h"
#include "state.h"
#include "buffer_physical_state.h"
#include "buffer_flux.h"
#include "fourier_transform.h"
#include "measure_workspace.h"

#include <memory>
#include <utility>
#include <string>

class Measure {
protected:
    const Params& params_;
    const Domain2D& domain_;
    std::shared_ptr<const MeasureCommandBase> command_;

public:
    Measure(const Params& params, const Domain2D& domain, std::shared_ptr<const MeasureCommandBase> command)
      : params_(params),
        domain_(domain),
        command_(std::move(command)) {}

    virtual ~Measure() = default;

    virtual FluxRequest flux_request() const {
        return {};
    }

    virtual void observe(const State& state, FourierTransform2D& fft, MeasureWorkspace& workspace, const FluxBuffer& flux, int step, double time) = 0;
    virtual void finalize() {}

    const std::string& id() const {
        return command_->id;
    }

};

#endif

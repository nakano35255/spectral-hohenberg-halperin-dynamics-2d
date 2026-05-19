#ifndef MEASURE_H
#define MEASURE_H

#include "simulationinfo.h"
#include "domain.h"
#include "state.h"
#include "buffer_physical_state.h"
#include "fourier_transform.h"

#include <memory>
#include <utility>
#include <string>

class Measure {
protected:
    std::shared_ptr<const MeasureCommandBase> command_;

public:
    Measure(const Params& params, std::shared_ptr<const MeasureCommandBase> command)
        : command_(std::move(command)) {
            (void)params;
        }

    virtual ~Measure() = default;

    virtual void observe(const State& state, PhysicalStateBuffer& physical, FourierTransform2D& fft, const Domain2D& domain, int step, double time) = 0;
    virtual void finalize() {}

    const std::string& id() const {
        return command_->id;
    }

};

#endif

#ifndef SFI_MEASURE_WORKSPACE_H
#define SFI_MEASURE_WORKSPACE_H

#include "buffer_physical_state.h"
#include "fourier_transform.h"
#include "state.h"

class MeasureWorkspace {
private:
    PhysicalStateBuffer physical_;
    int cached_step_ = -1;
    bool physical_ready_ = false;
    int num_fields_ = 0;

public:
    MeasureWorkspace(const Domain2D& domain, const Params& params)
        : physical_(domain, params),
          num_fields_(params.physics.num_order_parameters + 3) {}

    void begin_step(int step) {
        if (cached_step_ != step) {
            cached_step_ = step;
            physical_ready_ = false;
        }
    }

    PhysicalStateBuffer& ensure_physical(const State& state, FourierTransform2D& fft) {
        if (!physical_ready_) {
            fft.backward_many(num_fields_, state.data(), physical_.data());
            physical_ready_ = true;
        }
        return physical_;
    }

    bool physical_ready() const {
        return physical_ready_;
    }
};

#endif
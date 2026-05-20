#ifndef SFI_TIME_EVOLUTION_MASK_H
#define SFI_TIME_EVOLUTION_MASK_H

#include "fix_flag.h"
#include "simulationinfo.h"

#include <vector>

class TimeEvolutionMask {
private:
    int num_components_;
    int num_fields_;
    std::vector<bool> active_;

public:
    explicit TimeEvolutionMask(const Params& params)
        : num_components_(params.physics.num_components),
          num_fields_(num_components_ + 2),
          active_(static_cast<std::size_t>(num_fields_), true) {
        if (params.fix.enabled(FixFlag::Quiescent)) {
            active_[static_cast<std::size_t>(last_density_field())] = false;
            active_[static_cast<std::size_t>(momentum_x_field())] = false;
            active_[static_cast<std::size_t>(momentum_y_field())] = false;
        }

    }

    int last_density_field() const {
        return num_components_ - 1;
    }

    int momentum_x_field() const {
        return num_components_;
    }

    int momentum_y_field() const {
        return num_components_ + 1;
    }

    bool active_field(int field) const {
        return active_[static_cast<std::size_t>(field)];
    }
};

#endif

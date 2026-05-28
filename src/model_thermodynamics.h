#ifndef SFI_MODEL_THERMODYNAMICS_H
#define SFI_MODEL_THERMODYNAMICS_H

#include "simulationinfo.h"

#include <stdexcept>

class Thermodynamics {
private:
    double sound_speed_ = 0.0;

public:
    explicit Thermodynamics(double sound_speed)
        : sound_speed_(sound_speed) {
        if (sound_speed_ < 0.0) {
            throw std::runtime_error("thermodynamics requires nonnegative cT.");
        }
    }
    virtual ~Thermodynamics() = default;

    double sound_speed() const {
        return sound_speed_;
    }

    virtual double linear_pressure_coefficient() const {
        return sound_speed_ * sound_speed_;
    }

    virtual bool has_physical_pressure() const {
        return false;
    }

    virtual double physical_pressure(double rho) const {
        (void)rho;
        return 0.0;
    }

    virtual double physical_thermodynamic_energy(double rho) const {
        (void)rho;
        return 0.0;
    }
};

#endif

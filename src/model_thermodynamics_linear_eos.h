#ifndef SFI_MODEL_THERMODYNAMICS_LINEAR_EOS_H
#define SFI_MODEL_THERMODYNAMICS_LINEAR_EOS_H

#include "model_thermodynamics.h"

class LinearEOSThermodynamics : public Thermodynamics {
public:
    explicit LinearEOSThermodynamics(
        double sound_speed
    );
};

#endif

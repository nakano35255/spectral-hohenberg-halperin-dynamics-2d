#include "measure_registry_builtin.h"
#include "measure_snapshot_style.h"

MeasureRegistry build_measure_registry() {
    MeasureRegistry registry;
    registry.register_style(std::make_unique<SnapshotMeasureStyle>());
    return registry;
}

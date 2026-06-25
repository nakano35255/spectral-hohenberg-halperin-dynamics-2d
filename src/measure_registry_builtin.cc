#include "measure_registry_builtin.h"
#include "measure_snapshot_style.h"
#include "measure_time_series_style.h"
#include "measure_yokota_green_kubo_style.h"
#include "measure_ave_profile_style.h"
#include "measure_long_time_tail_style.h"
#ifdef SHHD_ENABLE_PASSIVE_SCALAR
#include "PASSIVE_SCALAR/measure_budget_spectrum_style.h"
#endif

MeasureRegistry build_measure_registry() {
    MeasureRegistry registry;
    registry.register_style(std::make_unique<SnapshotMeasureStyle>());
    registry.register_style(std::make_unique<TimeSeriesMeasureStyle>());
    registry.register_style(std::make_unique<YokotaGreenKuboMeasureStyle>());
    registry.register_style(std::make_unique<AveProfileMeasureStyle>());
    registry.register_style(std::make_unique<LongTimeTailMeasureStyle>());

#ifdef SHHD_ENABLE_PASSIVE_SCALAR
    registry.register_style(std::make_unique<BudgetSpectrumMeasureStyle>());
#endif

    return registry;
}

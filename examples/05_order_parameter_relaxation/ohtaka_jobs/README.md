# Ohtaka i8cpu replicas for order-parameter relaxation

This directory contains helpers for running many independent order-parameter
relaxation replicas on Ohtaka.

Default layout:

- partition: `i8cpu`
- allocation: `1` node, `128` MPI ranks
- replicas per job: `16`
- MPI ranks per replica: `8`
- default active grid: `256 x 256`
- default length: `8192 x 8192`
- default grid spacing: `dx = dy = 32`
- default time step: `dt = 4.0`
- default run time: `T = 200000.0`

The job generates one input file per sample from the defaults in
`prepare_ohtaka_inputs.py` and the overrides passed by the jobscript. Seeds are
drawn randomly using Python OS entropy and recorded in `seeds.dat`.

Before submitting, build on Ohtaka:

```sh
make -f Makefile.ohtaka -j 8
```

Submit the default 16-sample job:

```sh
cd /home/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d
sbatch examples/05_order_parameter_relaxation/ohtaka_jobs/job_ohtaka_i8cpu_16samples.sh
```

To append more samples under the same `RUN_NAME`, shift the sample offset:

```sh
SAMPLE_OFFSET=16 sbatch examples/05_order_parameter_relaxation/ohtaka_jobs/job_ohtaka_i8cpu_16samples.sh
SAMPLE_OFFSET=32 sbatch examples/05_order_parameter_relaxation/ohtaka_jobs/job_ohtaka_i8cpu_16samples.sh
```

The default output root is:

```text
/work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/05_order_parameter_relaxation/raw_data/sine2d_grid256_L8192_dt4_T200000/
```

Its structure is:

```text
runs/input_000.script
logs/stdout_000.log
logs/stderr_000.log
samples/sample_000/snapshots/physical_step001000.snapshot
samples/sample_000/snapshots/physical_step002000.snapshot
...
seeds.dat
config.dat
```

Useful overrides:

```sh
RUN_NAME=my_run_name \
SAMPLES=16 \
TASKS_PER_SAMPLE=8 \
SAMPLE_OFFSET=0 \
GRID_X=256 \
GRID_Y=256 \
LENGTH_X=8192 \
LENGTH_Y=8192 \
DT=4.0 \
RUN_TIME=200000.0 \
sbatch examples/05_order_parameter_relaxation/ohtaka_jobs/job_ohtaka_i8cpu_16samples.sh
```

If you want the job to rebuild before running:

```sh
BUILD_BEFORE_RUN=1 sbatch examples/05_order_parameter_relaxation/ohtaka_jobs/job_ohtaka_i8cpu_16samples.sh
```

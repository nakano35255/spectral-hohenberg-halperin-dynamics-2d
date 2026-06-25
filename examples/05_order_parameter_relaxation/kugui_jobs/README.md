# Kugui Packed Relaxation

This directory contains PBS helpers for long `sine2d` relaxation runs on
Kugui. It uses the same offset-based batching convention as the Ohtaka job
helpers elsewhere in this repository:

- one PBS job reserves one `F1cpu` node
- the job script generates input files for 16 replicas by default
- the job script launches 16 background `mpiexec` runs by default
- `OFFSET` selects which block of 16 replica ids is used

Default parameters:

- queue: `F1cpu`
- allocation: `1` node, `128` MPI ranks
- placement: `place=excl`
- replicas per job: `SAMPLES=16`
- ranks per replica: `TASKS_PER_SAMPLE=8`
- replica ids: `OFFSET .. OFFSET + SAMPLES - 1`, with `OFFSET=0` by default
- active grid: `256 x 256`
- length: `8192 x 8192`
- dealias: `three_halves`
- timestep: `dt = 4.0`
- run time: `2.4e6`, i.e. `600000` steps
- snapshot interval: `48000`, i.e. `12000` steps
- snapshots per sample: `50`
- transport: `eta = 0.12`, `M[0,0] = 0.12`


## Build

Build before submitting jobs:

```sh
cd ~/spectral-hohenberg-halperin-dynamics-2d

module purge
module load intel/2022.2.1 intel-mpi/2021.7.1

export HEFFTE_ROOT=$HOME/local/heffte-oneapi-fftw
export FFTW_ROOT=$HOME/local/fftw-oneapi
export LD_LIBRARY_PATH=$HEFFTE_ROOT/lib:$FFTW_ROOT/lib:${LD_LIBRARY_PATH:-}

make -f Makefile.ohtaka clean
make -f Makefile.ohtaka -j 8 all
```

The job script also supports rebuilding inside the job:

```sh
qsub -v BUILD_BEFORE_RUN=1 \
  examples/05_order_parameter_relaxation/kugui_jobs/job_kugui_relax_sine2d_grid256_L8192_dt4.pbs
```


## Relaxation

Submit replicas 0 through 15 from the repository root on Kugui:

```sh
qsub examples/05_order_parameter_relaxation/kugui_jobs/job_kugui_relax_sine2d_grid256_L8192_dt4.pbs
```

Submit replicas 16 through 31:

```sh
qsub -v OFFSET=16 \
  examples/05_order_parameter_relaxation/kugui_jobs/job_kugui_relax_sine2d_grid256_L8192_dt4.pbs
```

Each replica writes to its own sample id:

```text
/work/k0565/k056500/spectral-hohenberg-halperin-dynamics-2d/examples/05_order_parameter_relaxation/raw_data/
  sine2d_grid256_L8192_dt4_T2400000/
    runs/input_000.script
    logs/stdout_000.log
    logs/stderr_000.log
    samples/sample_000/snapshots/physical_step*.snapshot
```


## Common Overrides

PBS variables should be passed with `qsub -v`.

Use fewer ranks while still reserving one `F1cpu` node:

```sh
qsub -v OFFSET=0,TASKS_PER_SAMPLE=4 \
  examples/05_order_parameter_relaxation/kugui_jobs/job_kugui_relax_sine2d_grid256_L8192_dt4.pbs
```

Longer run with the same number of output snapshots:

```sh
qsub -v OFFSET=0,RUN_TIME=4800000.0,SNAPSHOT_DTOUT=96000.0,RUN_NAME=sine2d_grid256_L8192_dt4_T4800000 \
  examples/05_order_parameter_relaxation/kugui_jobs/job_kugui_relax_sine2d_grid256_L8192_dt4.pbs
```

Short i2cpu probe with only replica 0:

```sh
qsub -q i2cpu -l walltime=00:30:00 \
  -v SAMPLES=1,OFFSET=0,RUN_TIME=4000.0,SNAPSHOT_DTOUT=4000.0,THERMO_DTOUT=400.0,RUN_NAME=probe_sine2d_grid256_L8192_dt4_T4000 \
  examples/05_order_parameter_relaxation/kugui_jobs/job_kugui_relax_sine2d_grid256_L8192_dt4.pbs
```

`SAMPLE_OFFSET` is accepted as an alias for `OFFSET`.


## Notes

- These scripts are PBS scripts; use `qsub`, not `sbatch`.
- The default queue is `F1cpu`. The job requests `12:00:00`.
- One job runs `SAMPLES` replicas and launches `SAMPLES` background `mpiexec`
  commands.
- The script intentionally does not split hostfiles and disables Intel MPI
  pinning with `I_MPI_PIN=0`; each replica is launched as
  `mpiexec -n "$TASKS_PER_SAMPLE" ... &`.
- Physical snapshots are text files and can be large, so outputs are written
  under `/work/.../examples/05_order_parameter_relaxation/raw_data/`.
- Override `WORK_BASE` if your Kugui work directory is different.
- The current domain decomposition is x-slab based, so do not request more MPI
  ranks per replica than the available spectral x size. For default
  `grid 256 256` with `three_halves`, the practical upper bound is `193`
  ranks.

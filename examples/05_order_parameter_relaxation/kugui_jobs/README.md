# Kugui Single-Replica Relaxation

This directory contains PBS helpers for long `sine2d` relaxation runs on
Kugui. It is the Kugui counterpart of
`examples/05_order_parameter_relaxation/ohtaka_jobs`, but the scheduler model
is different:

- Ohtaka version: one large Slurm allocation, many replicas launched together.
- Kugui version: one `F1cpu` PBS job runs one replica.

This avoids running many independent `mpirun` commands inside the same PBS
allocation, which can lead to poor CPU placement and very slow runs.

Default parameters:

- queue: `F1cpu`
- allocation: `1` node, `128` MPI ranks
- one job runs one replica only
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
qsub -v REPLICA_ID=0,BUILD_BEFORE_RUN=1 \
  examples/05_order_parameter_relaxation/kugui_jobs/job_kugui_relax_sine2d_grid256_L8192_dt4.pbs
```


## Relaxation

Submit replica 0 from the repository root on Kugui:

```sh
qsub -v REPLICA_ID=0 \
  examples/05_order_parameter_relaxation/kugui_jobs/job_kugui_relax_sine2d_grid256_L8192_dt4.pbs
```

Submit replicas 0 through 15:

```sh
for rid in $(seq 0 15); do
  qsub -v REPLICA_ID="$rid" \
    examples/05_order_parameter_relaxation/kugui_jobs/job_kugui_relax_sine2d_grid256_L8192_dt4.pbs
done
```

Each replica writes:

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
qsub -v REPLICA_ID=0,TASKS_PER_RUN=64 \
  examples/05_order_parameter_relaxation/kugui_jobs/job_kugui_relax_sine2d_grid256_L8192_dt4.pbs
```

Longer run with the same number of output snapshots:

```sh
qsub -v REPLICA_ID=0,RUN_TIME=4800000.0,SNAPSHOT_DTOUT=96000.0,RUN_NAME=sine2d_grid256_L8192_dt4_T4800000 \
  examples/05_order_parameter_relaxation/kugui_jobs/job_kugui_relax_sine2d_grid256_L8192_dt4.pbs
```

Short i2cpu probe:

```sh
qsub -q i2cpu -l walltime=00:30:00 \
  -v REPLICA_ID=0,RUN_TIME=4000.0,SNAPSHOT_DTOUT=4000.0,THERMO_DTOUT=400.0,RUN_NAME=probe_sine2d_grid256_L8192_dt4_T4000 \
  examples/05_order_parameter_relaxation/kugui_jobs/job_kugui_relax_sine2d_grid256_L8192_dt4.pbs
```

`SAMPLE_ID` is still accepted as an alias for `REPLICA_ID`, for compatibility
with older submit commands.


## Notes

- These scripts are PBS scripts; use `qsub`, not `sbatch`.
- The default queue is `F1cpu`. The job requests `12:00:00`.
- One job runs one replica and launches exactly one `mpirun`.
- Physical snapshots are text files and can be large, so outputs are written
  under `/work/.../examples/05_order_parameter_relaxation/raw_data/`.
- Override `WORK_BASE` if your Kugui work directory is different.
- The current domain decomposition is x-slab based, so do not request more MPI
  ranks than the available spectral x size. For default `grid 256 256` with
  `three_halves`, the practical upper bound is `193` ranks.

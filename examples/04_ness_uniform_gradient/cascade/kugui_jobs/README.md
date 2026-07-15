# Kugui Large Restart Run

This directory contains Kugui PBS job helpers for a large `examples/04_ness_uniform_gradient/cascade` style
run. It is the Kugui counterpart of `examples/04_ness_uniform_gradient/cascade/ohtaka_jobs`, but the
scheduler model is different:

- Ohtaka version: one large Slurm allocation, many replicas launched together.
- Kugui version: one `L1cpu` PBS job runs one replica and one segment.

This is intended for manually submitting a few replicas and extending each
replica by restart segments.

Default parameters:

- queue: `L1cpu`
- allocation: `1` node, `128` MPI ranks
- one job runs one replica only
- active grid: `512 x 512`
- length: `16384 x 16384`
- dealias: `three_halves`
- `dt = 16`
- `eta = M[0,0] = 0.004`
- force/gradient amplitude: `2 / 16384 = 0.0001220703125`
- relaxation segment time: `RELAX_TIME=400000000.0`
- budget measurement time: `BUDGET_TIME=50000000.0`

The default grid is `512 x 512`, because the current Kugui plan is to use
`L1cpu` one sample at a time. The scripts can still be overridden for a
`1024 x 1024` run.


## Build

Build with `PASSIVE_SCALAR` enabled before submitting jobs:

```sh
cd ~/spectral-hohenberg-halperin-dynamics-2d

module purge
module load intel/2022.2.1 intel-mpi/2021.7.1

export HEFFTE_ROOT=$HOME/local/heffte-oneapi-fftw
export FFTW_ROOT=$HOME/local/fftw-oneapi
export LD_LIBRARY_PATH=$HEFFTE_ROOT/lib:$FFTW_ROOT/lib:${LD_LIBRARY_PATH:-}

make -f Makefile.ohtaka yes-PASSIVE-SCALAR
make -f Makefile.ohtaka clean
make -f Makefile.ohtaka -j 8 all
```

The job scripts also support rebuilding inside the job:

```sh
qsub -v BUILD_BEFORE_RUN=1 examples/04_ness_uniform_gradient/cascade/kugui_jobs/job_kugui_relax_D0_0p004_grid512_dt16.pbs
```


## Relaxation

Submit the first relaxation segment of replica 0:

```sh
qsub -v REPLICA_ID=0,RELAX_SEGMENT=1 \
  examples/04_ness_uniform_gradient/cascade/kugui_jobs/job_kugui_relax_D0_0p004_grid512_dt16.pbs
```

Continue the same replica from the restart written by segment 1:

```sh
qsub -v REPLICA_ID=0,RELAX_SEGMENT=2 \
  examples/04_ness_uniform_gradient/cascade/kugui_jobs/job_kugui_relax_D0_0p004_grid512_dt16.pbs
```

Run another sample manually:

```sh
qsub -v REPLICA_ID=1,RELAX_SEGMENT=1 \
  examples/04_ness_uniform_gradient/cascade/kugui_jobs/job_kugui_relax_D0_0p004_grid512_dt16.pbs
```

Each relaxation segment writes:

```text
/work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/04_ness_uniform_gradient/cascade/kugui_jobs/
  relax_D0_0p004_grid512_L16384_dt16/
    replica_000/
      restarts/relax_001.restart
      restarts/relax_002.restart
      results/time_series_relax_001.dat
      results/time_series_relax_002.dat
      runs/input_relax_001.script
      runs/input_relax_002.script
      logs/stdout_relax_001.log
      logs/stderr_relax_001.log
```

If `relax_NNN.restart` already exists, the corresponding job exits without
rerunning it.


## Budget Spectrum

After checking the relaxation time-series, run budget/spectrum from a selected
restart. For example, to use `relax_002.restart` from replica 0:

```sh
qsub -v REPLICA_ID=0,RESTART_INDEX=2 \
  examples/04_ness_uniform_gradient/cascade/kugui_jobs/job_kugui_budget_D0_0p004_grid512_dt16.pbs
```

The budget job writes:

```text
replica_000/
  results/time_series_budget_from_002.dat
  results/budget_shell_from_002.dat
  restarts/budget_from_002.restart
  runs/input_budget_from_002.script
  logs/stdout_budget_from_002.log
  logs/stderr_budget_from_002.log
```


## Structure Factor Shell

After checking the relaxation time-series, run the order-parameter static
correlation shell from the same restart used for the budget/spectrum job.
For example, to use `relax_005.restart` from replica 0 for `D0=0.004`:

```sh
qsub -v REPLICA_ID=0,RESTART_INDEX=5 \
  examples/04_ness_uniform_gradient/cascade/kugui_jobs/job_kugui_structure_D0_0p004_grid512_dt16.pbs
```

For `D0=0.12` after two relaxation segments:

```sh
qsub -v REPLICA_ID=0,RESTART_INDEX=2 \
  examples/04_ness_uniform_gradient/cascade/kugui_jobs/job_kugui_structure_D0_0p12_grid512_dt04.pbs
```

For `D0=4.0` after one relaxation segment:

```sh
qsub -v REPLICA_ID=0,RESTART_INDEX=1 \
  examples/04_ness_uniform_gradient/cascade/kugui_jobs/job_kugui_structure_D0_4p00_grid512_dt04.pbs
```

The structure job writes:

```text
replica_000/
  results/time_series_structure_from_NNN.dat
  results/static_corr_shell_from_NNN.dat
  restarts/structure_from_NNN.restart
  runs/input_structure_from_NNN.script
  logs/stdout_structure_from_NNN.log
  logs/stderr_structure_from_NNN.log
```

The measured static-correlation command is:

```text
measure sc_shell correlation/static on nevery 20 nblock 200 file .../static_corr_shell_from_NNN.dat mode shell average running cross off target psi[0]
```


## Common Overrides

PBS variables should be passed with `qsub -v`.

Short probe:

```sh
qsub -v REPLICA_ID=0,RELAX_SEGMENT=1,RELAX_TIME=204800.0,TIME_SERIES_DTOUT=1024.0,RUN_NAME=f1_probe_D0_0p004_grid512_L16384_dt16 \
  examples/04_ness_uniform_gradient/cascade/kugui_jobs/job_kugui_relax_D0_0p004_grid512_dt16.pbs
```

Use fewer MPI ranks while still reserving one `L1cpu` node:

```sh
qsub -v REPLICA_ID=0,TASKS_PER_RUN=64 \
  examples/04_ness_uniform_gradient/cascade/kugui_jobs/job_kugui_relax_D0_0p004_grid512_dt16.pbs
```

Run the same scripts with `1024 x 1024` and `L = 32768`:

```sh
qsub -v REPLICA_ID=0,RELAX_SEGMENT=1,GRID_N=1024,LENGTH_N=32768,GRADIENT_AMPLITUDE=0.00006103515625,RUN_NAME=relax_D0_0p004_grid1024_L32768_dt16 \
  examples/04_ness_uniform_gradient/cascade/kugui_jobs/job_kugui_relax_D0_0p004_grid512_dt16.pbs
```


## Notes

- These scripts are PBS scripts; use `qsub`, not `sbatch`.
- The default queue is `L1cpu`. All jobs request `96:00:00`.
- One job runs one replica and writes one final restart.
- For multi-segment relaxation, increment `RELAX_SEGMENT` manually.
- The restart files are text files and can be large, so outputs are written
  under `/work/.../examples/04_ness_uniform_gradient/cascade/kugui_jobs/`.
- The current domain decomposition is x-slab based, so do not request more MPI
  ranks than the available spectral x size. For default `grid 512 512` with
  `three_halves`, the practical upper bound is `385` ranks.

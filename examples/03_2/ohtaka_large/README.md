# Ohtaka Large Restart Run

This directory contains Ohtaka job helpers for a large `examples/03_2` style run.

Default parameters:

- partition: `F72cpu`
- allocation: `4` nodes, `512` MPI ranks
- active grid: `1024 x 1024`
- length: `32768 x 32768`
- dealias: `three_halves`
- `dt = 16`
- `eta = M[0,0] = 0.004`
- force/gradient amplitude: `2 / 32768 = 0.00006103515625`

The workflow is intentionally split into two jobs.

1. `job_ohtaka_relax_D0_0p004_grid1024_dt16.sh`
   - runs several restart segments sequentially
   - writes one restart file per segment
   - writes one time-series file per segment for steady-state checks
   - skips completed restart segments, so the same job can be resubmitted

2. `job_ohtaka_budget_D0_0p004_grid1024_dt16.sh`
   - reads a chosen relaxation restart
   - measures `budget/spectrum` in `mode shell`
   - writes a final restart after the budget run

Before submitting budget jobs, build with `PASSIVE_SCALAR` enabled:

```sh
make -f Makefile.ohtaka yes-PASSIVE-SCALAR
make -f Makefile.ohtaka clean
make -f Makefile.ohtaka -j 8
```

The job scripts also support rebuilding inside the job:

```sh
BUILD_BEFORE_RUN=1 sbatch examples/03_2/ohtaka_large/job_ohtaka_relax_D0_0p004_grid1024_dt16.sh
```

After checking the relaxation time-series, choose the restart index for the
budget run:

```sh
RESTART_INDEX=8 sbatch examples/03_2/ohtaka_large/job_ohtaka_budget_D0_0p004_grid1024_dt16.sh
```

The restart files are text files. For this grid they can be large, so the jobs
write them under `/work/.../examples/03_2/ohtaka_large/`.

The MPI rank count is kept at `512` because the current domain decomposition is
an x-slab decomposition and cannot use more ranks than the spectral x size.

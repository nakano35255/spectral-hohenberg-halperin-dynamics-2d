# Ohtaka Large Restart Run

This directory contains Ohtaka job helpers for a large `examples/03_2` style run.

Default parameters:

- partition: `F72cpu`
- allocation: `72` nodes, `9216` MPI ranks
- parallel layout: `18` independent replicas, each using `4` nodes and `512` MPI ranks
- active grid: `1024 x 1024`
- length: `32768 x 32768`
- dealias: `three_halves`
- `dt = 16`
- `eta = M[0,0] = 0.004`
- force/gradient amplitude: `2 / 32768 = 0.00006103515625`

The workflow is intentionally split into two jobs.

1. `job_ohtaka_relax_D0_0p004_grid1024_dt16.sh`
   - runs `18` independent relaxation samples in parallel by default
   - writes one restart file per sample: `relax_001.restart`
   - writes one time-series file per sample for steady-state and wall-time checks
   - skips completed samples, so the same job can be resubmitted

2. `job_ohtaka_budget_D0_0p004_grid1024_dt16.sh`
   - reads `relax_001.restart` from each replica by default
   - measures `budget/spectrum` in `mode shell`
   - writes a final restart after each budget run

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

For a short F72cpu probe using one `4`-node, `512`-rank run:

```sh
RUN_NAME=f72_probe_D0_0p004_grid1024_L32768_dt16 \
REPLICAS=1 \
RELAX_TIME=51200.0 \
TIME_SERIES_DTOUT=1024.0 \
sbatch examples/03_2/ohtaka_large/job_ohtaka_relax_D0_0p004_grid1024_dt16.sh
```

Use a separate `RUN_NAME` for probes so a short `relax_001.restart` does not
make the production job skip `replica_000`.

The relaxation job runs for `RELAX_TIME=20000000.0` by default with a
`16:00:00` wall-time limit. A F72cpu probe with `RELAX_TIME=51200.0` took
`136.839 s` for `3200` steps, so `RELAX_TIME=20000000.0` corresponds to about
`1.25e6` steps and roughly `14.9 h` of solver time, leaving margin for the final
restart write.

For a shorter wall-time probe:

```sh
RELAX_TIME=10000000.0 sbatch examples/03_2/ohtaka_large/job_ohtaka_relax_D0_0p004_grid1024_dt16.sh
```

After checking the relaxation time-series, run the budget job from
`relax_001.restart`:

```sh
sbatch examples/03_2/ohtaka_large/job_ohtaka_budget_D0_0p004_grid1024_dt16.sh
```

To rerun or analyze only selected replicas:

```sh
REPLICA_IDS="0 3 7" sbatch examples/03_2/ohtaka_large/job_ohtaka_budget_D0_0p004_grid1024_dt16.sh
```

The restart files are text files. For this grid they can be large, so the jobs
write them under `/work/.../examples/03_2/ohtaka_large/`.

Output is organized by replica:

```text
/work/.../examples/03_2/ohtaka_large/relax_D0_0p004_grid1024_L32768_dt16/
  replica_000/
    results/
    restarts/
    runs/
    logs/
  replica_001/
  ...
```

Each individual run uses `512` ranks because the current domain decomposition is
an x-slab decomposition and cannot use more ranks than the spectral x size. The
full F72cpu allocation is used by launching `18` such runs concurrently.

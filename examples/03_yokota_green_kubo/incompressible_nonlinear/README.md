# Incompressible nonlinear Yokota Green-Kubo viscosity

This example estimates the wave-number-dependent shear viscosity from the
Yokota Green-Kubo stress integral.  The setup is chosen to match
`examples/02_kolomogorov_flow/incompressible/viscosity` as closely as possible:

- `grid 256 256`
- `length 256.0 256.0`
- `dealias three_halves`
- `timestep 0.01`
- `time_evolution srk3/incompressible`
- `kBT 1.0`
- `density 1.0`
- nonlinear momentum advection enabled

Unlike Example 02, no sinusoidal force is applied.  The viscosity is estimated
from equilibrium stress fluctuations.  The current incompressible setup uses
`mode diagonal`, i.e.

```math
(k_x,k_y)=\frac{2\pi n}{L}(1,1),
\qquad n=1,\ldots,N/2-1.
```

For this Yokota Green-Kubo estimator, the plotted scalar wave number is the
diagonal component `k=kx=ky=2*pi*n/L`, not the vector magnitude `sqrt(kx^2+ky^2)`.

## Layout

This directory keeps the incompressible nonlinear data in the same style as
Examples 02 and 03:

- `raw_data/`: downloaded production data used for the main analysis.
- `processed_data/`: compact CSV tables generated from `raw_data/`.
- `legacy/`: older or superseded runs kept for comparison.
- `tmp/`: incomplete downloads, probes, and scratch analysis.
- `figures/`: regenerated figures, typically written with
  `estimate_viscosity.py --figure-dir ...`.
- `ohtaka_jobs/`: Ohtaka job scripts and their input generator.
- `kugui_jobs/`: Kugui job scripts, when needed.

Raw run directories under `raw_data/` or `legacy/` should keep the generated
`runs/` and `results/` subdirectories.  For example:

```text
raw_data/eta0_0p1_grid256_L256_dt0p01_T25000_diag_n576/
  runs/
  results/
```

## Ohtaka production inputs

The Ohtaka helper writes generated inputs and outputs under `/work`.

```sh
sbatch examples/03_yokota_green_kubo/incompressible_nonlinear/ohtaka_jobs/job_ohtaka_ykgk_eta0_0p1.sh
sbatch examples/03_yokota_green_kubo/incompressible_nonlinear/ohtaka_jobs/job_ohtaka_ykgk_eta0_0p5.sh
```

Additional independent batches can be appended to the same run directory by
setting `SAMPLE_OFFSET`.  For example, after the default `0..575` batch:

```sh
SAMPLE_OFFSET=576 sbatch examples/03_yokota_green_kubo/incompressible_nonlinear/ohtaka_jobs/job_ohtaka_ykgk_eta0_0p1.sh
```

The jobs use `dt=0.01` and `T=25000.0`, with `ykgk_block_time=5000.0` so each
sample contributes 5 Green-Kubo blocks.  The production inputs only write the
Yokota Green-Kubo output, without a separate time-series measure.  Because one
Green-Kubo sample measures all diagonal wave numbers at once, the F72cpu jobs
fill the allocation by running many independent samples with 16 MPI ranks per
sample.

After downloading results locally, keep curated production data under
`raw_data/`. Use `tmp/` only for incomplete downloads, probes, and scratch
analysis; older superseded runs can later be moved to `legacy/`.

Analyze the sample files as a single ensemble:

```sh
python3 examples/03_yokota_green_kubo/incompressible_nonlinear/estimate_viscosity.py \
  --input-glob "examples/03_yokota_green_kubo/incompressible_nonlinear/raw_data/eta0_0p1_grid256_L256_dt0p01_T25000_diag_n576/results/yokota_green_kubo_*.dat" \
  --input-script "examples/03_yokota_green_kubo/incompressible_nonlinear/raw_data/eta0_0p1_grid256_L256_dt0p01_T25000_diag_n576/runs/input_000.script" \
  --output-dir "examples/03_yokota_green_kubo/incompressible_nonlinear/processed_data/eta0_0p1_grid256_L256_dt0p01_T25000_diag_n576" \
  --figure-dir "examples/03_yokota_green_kubo/incompressible_nonlinear/figures/eta0_0p1_grid256_L256_dt0p01_T25000_diag_n576"
```

The processed output contains:

```text
processed_data/eta0_0p1_grid256_L256_dt0p01_T25000_diag_n576/
  metadata.csv
  yokota_green_kubo_summary.csv
```

The Kolmogorov-flow comparison figure is written under `figures/`:

```sh
python3 examples/03_yokota_green_kubo/incompressible_nonlinear/plot_compare_kolmogorov.py
```

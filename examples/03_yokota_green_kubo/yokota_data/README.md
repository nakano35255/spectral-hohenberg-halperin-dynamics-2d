# Yokota comparison data

This directory collects the data and generated artifacts used to compare the
Yokota molecular-dynamics reference with the fluctuating-hydrodynamics runs.

- `md_data/`: digitized Yokota et al. MD data and FNS fit tables/figures.
- `raw_data/`: FHD production outputs used for the Yokota MD comparison.
  Each run directory keeps the generated `config*.dat`, `runs/`, and
  `results/` files. Per-sample seeds are recorded in `runs/input_*.script`, so
  separate seed manifest files are intentionally omitted.
- `processed_data/`: compact CSV summaries generated from `raw_data/`.
- `figures/`: generated MD comparison and diagnostic figures.
- `ohtaka_jobs/`: Ohtaka job scripts for the MD comparison parameter scan.
- `tmp/`: incomplete downloads, probes, and scratch staging for the comparison.

The older standalone `eta0=0.1` and `eta0=0.5` runs remain under
`../incompressible_nonlinear/`.

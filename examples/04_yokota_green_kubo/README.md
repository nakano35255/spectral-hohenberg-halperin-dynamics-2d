# Example 04: Yokota Green-Kubo viscosity estimate

This example estimates the wave-number-dependent shear viscosity from the
Yokota Green-Kubo stress integral.

The incompressible nonlinear setup is chosen to match
`examples/02_Kolomogorov_flow/incompressible/viscosity` as closely as possible:

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

## Single input

From the repository root:

```sh
mkdir -p examples/04_yokota_green_kubo/results
./src/out.exe examples/04_yokota_green_kubo/input_incompressible.script
python3 examples/04_yokota_green_kubo/estimate_viscosity.py \
  --input examples/04_yokota_green_kubo/results/yokota_green_kubo_incompressible.dat \
  --input-script examples/04_yokota_green_kubo/input_incompressible.script \
  --output-dir examples/04_yokota_green_kubo/results
```

## Ohtaka production inputs

The Ohtaka helper writes generated inputs and outputs under `/work`.

```sh
sbatch examples/04_yokota_green_kubo/incompressible_nonlinear/jobs/job_ohtaka_ykgk_eta0_0p1.sh
sbatch examples/04_yokota_green_kubo/incompressible_nonlinear/jobs/job_ohtaka_ykgk_eta0_0p5.sh
```

The jobs use the same `dt=0.01` and `T=50000.0` as the Example 02 viscosity
runs, with `ykgk_block_time=5000.0` so each sample contributes 10 Green-Kubo
blocks.  Because one Green-Kubo sample measures all diagonal wave numbers at
once, the F72cpu jobs fill the allocation by running many independent samples
with 16 MPI ranks per sample.

After a job finishes, analyze the sample files as a single ensemble:

```sh
python3 examples/04_yokota_green_kubo/estimate_viscosity.py \
  --input-glob "/work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/04_yokota_green_kubo/incompressible_nonlinear/eta0_0p1_grid256_L256_dt0p01_T50000_diag_n576/results/yokota_green_kubo_*.dat" \
  --input-script "/work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/04_yokota_green_kubo/incompressible_nonlinear/eta0_0p1_grid256_L256_dt0p01_T50000_diag_n576/runs/input_000.script" \
  --output-dir "/work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/04_yokota_green_kubo/incompressible_nonlinear/eta0_0p1_grid256_L256_dt0p01_T50000_diag_n576/analysis"
```

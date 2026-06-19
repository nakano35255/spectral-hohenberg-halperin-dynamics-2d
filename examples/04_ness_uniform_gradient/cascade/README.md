# Example 04 cascade: passive-scalar budget spectrum

This example is a small steady-state check for the `PASSIVE_SCALAR` package.
It uses the same basic setup as `examples/04_ness_uniform_gradient`:

- `srk3/incompressible`
- one quadratic order parameter
- nonlinear order-parameter advection
- momentum noise
- a uniform-gradient order-parameter force

The first `run` segment relaxes the system. The `budget/spectrum` measures are
created only after that segment, so `budget_shell.dat` and `budget_2d.dat`
contain measurements from the second segment only.

The grid is intentionally small (`32 x 32`) and the gradient amplitude is
larger than in the production `03` runs. This keeps the example quick and makes
the budget terms easier to see in a short run. For production data, increase the
box size, lower the forcing, and extend both the relaxation and averaging runs.

Large Ohtaka restart jobs for `grid 1024 1024`, `dt = 16`, and
`eta = M[0,0] = 0.004` are collected in `ohtaka_jobs/`.

Kugui PBS jobs for one-sample-at-a-time `F1cpu` runs are collected in
`kugui_jobs/`. The default Kugui setup uses `grid 512 512`, `dt = 16`, and
`eta = M[0,0] = 0.004`, with restart segments submitted manually per replica.

Build with the package enabled:

```sh
make yes-PASSIVE-SCALAR
make clean
make
```

Run from the repository root:

```sh
./src/out.exe examples/04_ness_uniform_gradient/cascade/input.script
python3 examples/04_ness_uniform_gradient/cascade/plot_budget.py
```

The main output files are:

```text
examples/04_ness_uniform_gradient/cascade/raw_data/budget_shell.dat
examples/04_ness_uniform_gradient/cascade/raw_data/budget_2d.dat
examples/04_ness_uniform_gradient/cascade/raw_data/budget_shell_summary.txt
examples/04_ness_uniform_gradient/cascade/raw_data/budget_shell.png
```

`budget_shell.dat` is the most compact correctness check. Each block contains
the shell-averaged transfer, dissipation, production, and total budget spectra.

# Example 04: sine-force viscosity estimate

This example drives a shear flow with a sinusoidal body force,

```math
\partial_t j_x = \cdots + F \sin(k y),
```

and measures the stationary profiles of `vx` and `pi_xy` with
`measure ave/profile`.

The incompressible input uses `dx = dy = 1`, `eta = 0.1`, momentum noise,
and nonlinear advection.  The box size is kept at `16 x 16` so that the
driven shear mode relaxes on the time scale of the example run.

In the linear steady response,

```math
v_x(y) \simeq A_v \sin(k y),
\qquad
\Pi_{xy}(y) \simeq A_\Pi \cos(k y),
```

with

```math
A_v = F / (\eta k^2),
\qquad
A_\Pi = -F / k.
```

The checker estimates

```math
\eta_v = F / (k^2 A_v),
\qquad
\eta_\Pi = -A_\Pi / (k A_v).
```

For the noisy nonlinear run, `eta_v` is the forced-response estimator.  The
stress-based value is a force-balance consistency check; it should approach
the same value only after sufficient time averaging of the total momentum
flux.

Run from the repository root:

```sh
./src/out.exe examples/04_sine_force_viscosity/input_incompressible.script
./src/out.exe examples/04_sine_force_viscosity/input_compressible.script
python3 examples/04_sine_force_viscosity/estimate_viscosity.py
```

The result summary is written to
`examples/04_sine_force_viscosity/results/viscosity_estimate.txt`.

The checker also writes pointwise profile diagnostics,

```text
examples/04_sine_force_viscosity/results/incompressible_profile_diagnostics.dat
examples/04_sine_force_viscosity/results/compressible_profile_diagnostics.dat
```

These files contain the reconstructed force profile `f_x`, the measured
`vx` and `pi_xy`, and the local ratios `vx/sin(ky)` and `pi_xy/cos(ky)`.
They also include pointwise estimates of `eta` from `f_x/(k^2 vx)` and
`-pi_xy/(d_y vx)`. Rows near zeros of the denominator are written as `nan`.

The corresponding figures are written to

```text
examples/04_sine_force_viscosity/results/incompressible_profile.png
examples/04_sine_force_viscosity/results/compressible_profile.png
```

Each figure shows the spatial profiles of `f_x`, `vx`, `pi_xy`,
`f_x/(k^2 vx)`, and the stress-based local estimate
`-pi_xy/(k A_v cos(ky))`, where `A_v` is obtained by projecting `vx` onto
`sin(ky)`.

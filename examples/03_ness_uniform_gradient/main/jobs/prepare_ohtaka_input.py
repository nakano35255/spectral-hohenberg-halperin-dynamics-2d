#!/usr/bin/env python3
import argparse
import os
import struct
from pathlib import Path


WORK_EXAMPLE = Path("/work/k0565/k056500/spectral-hohenberg-halperin-dynamics-2d/examples/03_ness_uniform_gradient")
OUTPUT_ROOT = Path(os.environ.get("SHHD_OUTPUT_ROOT", WORK_EXAMPLE))


def label(prefix, value):
    return f"{prefix}_{value}"


def steps_from_time(time_value, dt, name):
    steps = round(time_value / dt)
    if steps <= 0 or abs(steps * dt - time_value) > 1.0e-12:
        raise RuntimeError(f"{name} must be a positive integer multiple of dt")
    return steps


def generate_input(args, seed, time_series_path):
    dt = float(args.dt)
    eta = args.eta if args.eta is not None else f"{float(args.d0) * args.schmidt_number:.16g}"
    run_steps = steps_from_time(args.run_time, dt, "run-time")
    time_series_nevery = args.time_series_nevery
    if args.time_series_dtout is not None:
        time_series_nevery = steps_from_time(args.time_series_dtout, dt, "time-series-dtout")

    lines = [
        "dimension           2",
        "boundary            p p",
        "",
        f"grid                {args.grid[0]} {args.grid[1]}",
        f"length              {args.length[0]} {args.length[1]}",
        f"dealias             {args.dealias}",
        "",
        "order_parameters    1",
        "",
        f"timestep            {args.dt}",
        f"time_evolution      {args.time_evolution}",
        "",
        f"model free_energy   quadratic a[0] {args.free_energy_a}",
        f"model transport     constant eta {eta} zeta 0.0 M[0,0] {args.d0}",
        "",
        "fix                 1 order_parameter nonlinear on",
        f"fix                 2 momentum noise on seed {seed} kBT {args.kBT}",
        f"fix                 3 order_parameter force/gradient on component 0 direction x amplitude {args.gradient_amplitude}",
        "",
        f"set                 density uniform value {args.density}",
        "set                 momentum all uniform value 0.0",
        f"set                 order_parameter all uniform value {args.order_parameter}",
        "",
        f"thermo              observe on progress off nevery {args.thermo_nevery}",
        f"measure             1 time_series on nevery {time_series_nevery} file {time_series_path.as_posix()} target E_K |psi[0]|^2 |d_psi[0]|^2 Jpsi[0]_x",
        "restart             off",
        "",
        f"run                 {run_steps}",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=16)
    parser.add_argument("--D0", dest="d0", required=True)
    parser.add_argument("--eta", default=None)
    parser.add_argument("--schmidt-number", type=float, default=1.0)
    parser.add_argument("--dt", required=True)
    parser.add_argument("--run-time", type=float, default=10000.0)
    parser.add_argument("--time-series-dtout", type=float, default=None)
    parser.add_argument("--grid", nargs=2, type=int, default=[128, 128])
    parser.add_argument("--length", nargs=2, default=["4096", "4096"])
    parser.add_argument("--dealias", default="three_halves")
    parser.add_argument("--time-evolution", default="srk3/incompressible")
    parser.add_argument("--free-energy-a", default="1.0")
    parser.add_argument("--gradient-amplitude", default="0.000048828125")
    parser.add_argument("--kBT", default="1.0")
    parser.add_argument("--density", default="1.0")
    parser.add_argument("--order-parameter", default="0.0")
    parser.add_argument("--thermo-nevery", type=int, default=100)
    parser.add_argument("--time-series-nevery", type=int, default=50)
    args = parser.parse_args()

    if args.samples <= 0:
        raise RuntimeError("samples must be positive")
    if args.thermo_nevery <= 0 or args.time_series_nevery <= 0:
        raise RuntimeError("nevery values must be positive")

    case_dir = Path(label("D0", args.d0)) / label("dt", args.dt)
    run_dir = OUTPUT_ROOT / "runs" / case_dir
    result_dir = OUTPUT_ROOT / "results" / case_dir
    run_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)

    for sample in range(args.samples):
        seed = struct.unpack("<I", os.urandom(4))[0] % 2147483646 + 1
        sid = f"{sample:03d}"
        input_path = run_dir / f"input_{sid}.script"
        time_series = result_dir / f"time_series_{sid}.dat"
        input_path.write_text(generate_input(args, seed, time_series))


if __name__ == "__main__":
    main()

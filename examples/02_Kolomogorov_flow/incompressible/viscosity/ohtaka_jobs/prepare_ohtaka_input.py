#!/usr/bin/env python3
import argparse
import os
import struct
from pathlib import Path


CASE_NAME = Path(__file__).resolve().parents[1].name
WORK_EXAMPLE_BASE = Path("/work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/02_kolomogorov_flow/incompressible")
OUTPUT_ROOT = Path(os.environ.get("SHHD_OUTPUT_ROOT", WORK_EXAMPLE_BASE / CASE_NAME))


def steps_from_time(time_value, dt, name):
    if time_value is None:
        return None
    steps = round(time_value / dt)
    if steps <= 0 or abs(steps * dt - time_value) > 1.0e-12:
        raise RuntimeError(f"{name} must be a positive integer multiple of dt")
    return steps


def resolve_run_steps(args, dt):
    run_steps = steps_from_time(args.run_time, dt, "run-time")
    relax_steps = steps_from_time(args.relax_time, dt, "relax-time")
    measure_steps = steps_from_time(args.measure_time, dt, "measure-time")

    if run_steps is not None:
        if relax_steps is not None or measure_steps is not None:
            raise RuntimeError("--run-time cannot be combined with --relax-time or --measure-time")
        return run_steps
    if relax_steps is not None and measure_steps is not None:
        return relax_steps + measure_steps
    if relax_steps is not None:
        return relax_steps
    if measure_steps is not None:
        return measure_steps
    return args.run_steps


def nevery_from_time(time_value, dt, default_steps, name):
    if time_value is None:
        return default_steps
    return steps_from_time(time_value, dt, name)


def generate_input(args, seed, time_series_path, profile_path):
    dt = float(args.dt)
    grid_x, grid_y = args.grid
    length_x, length_y = args.length
    run_steps = resolve_run_steps(args, dt)
    thermo_nevery = args.thermo_nevery
    time_series_nevery = nevery_from_time(args.time_series_dtout, dt, args.time_series_nevery, "time-series-dtout")
    profile_nevery = nevery_from_time(args.profile_dtout, dt, args.profile_nevery, "profile-dtout")
    profile_nblock = nevery_from_time(args.profile_block_time, dt, args.profile_nblock, "profile-block-time")

    if profile_nblock % profile_nevery != 0:
        raise RuntimeError("profile-block-time must be an integer multiple of profile-dtout")

    lines = [
        "dimension           2",
        "boundary            p p",
        "",
        f"grid                {grid_x} {grid_y}",
        f"length              {length_x} {length_y}",
        f"dealias             {args.dealias}",
        "",
        "order_parameters    0",
        "",
        f"timestep            {args.dt}",
        f"time_evolution      {args.time_evolution}",
        "",
        f"model transport     constant eta {args.eta}",
        "",
        "fix                 1 momentum nonlinear on",
        f"fix                 2 momentum force/sine on component x axis y nk {args.force_nk} amplitude {args.force_amplitude}",
        f"fix                 3 momentum noise on seed {seed} kBT {args.kBT}",
        "",
        f"set                 density uniform value {args.density}",
        "set                 momentum all uniform value 0.0",
        "",
        f"thermo              observe on progress off nevery {thermo_nevery}",
        "restart             off",
        "",
        f"measure             1 time_series on nevery {time_series_nevery} file {time_series_path.as_posix()} target E_T E_K",
        f"measure             2 ave/profile on axis y nevery {profile_nevery} nblock {profile_nblock} file {profile_path.as_posix()} average block target vx pi_xy",
        "",
        f"run                 {run_steps}",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=16)
    parser.add_argument("--grid", nargs=2, type=int, default=[64, 64])
    parser.add_argument("--length", nargs=2, default=["64.0", "64.0"])
    parser.add_argument("--dealias", default="three_halves")
    parser.add_argument("--dt", default="0.01")
    parser.add_argument("--time-evolution", default="srk3/incompressible")
    parser.add_argument("--eta", default="0.1")
    parser.add_argument("--force-nk", type=int, default=1)
    parser.add_argument("--force-amplitude", default="0.02")
    parser.add_argument("--kBT", default="1.0")
    parser.add_argument("--density", default="1.0")
    parser.add_argument("--run-steps", type=int, default=200000)
    parser.add_argument("--run-time", type=float)
    parser.add_argument("--relax-time", type=float)
    parser.add_argument("--measure-time", type=float)
    parser.add_argument("--thermo-nevery", type=int, default=10000)
    parser.add_argument("--time-series-nevery", type=int, default=1000)
    parser.add_argument("--profile-nevery", type=int, default=1000)
    parser.add_argument("--profile-nblock", type=int, default=10000)
    parser.add_argument("--time-series-dtout", type=float)
    parser.add_argument("--profile-dtout", type=float)
    parser.add_argument("--profile-block-time", type=float)
    args = parser.parse_args()

    if args.samples <= 0:
        raise RuntimeError("samples must be positive")
    if args.run_steps <= 0:
        raise RuntimeError("run-steps must be positive")
    if args.thermo_nevery <= 0 or args.time_series_nevery <= 0 or args.profile_nevery <= 0 or args.profile_nblock <= 0:
        raise RuntimeError("nevery and nblock values must be positive")

    run_dir = OUTPUT_ROOT / "runs"
    result_dir = OUTPUT_ROOT / "results"
    run_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)

    for sample in range(args.samples):
        seed = struct.unpack("<I", os.urandom(4))[0] % 2147483646 + 1
        sid = f"{sample:03d}"
        input_path = run_dir / f"input_{sid}.script"
        time_series = result_dir / f"time_series_{sid}.dat"
        profile = result_dir / f"profile_{sid}.dat"
        input_path.write_text(generate_input(args, seed, time_series, profile))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Generate Ohtaka inputs for an order-parameter sine-force response run."""

import argparse
import math
import secrets
from pathlib import Path


MAX_SEED = 2147483646


def random_seed(used_seeds):
    while True:
        seed = secrets.randbelow(MAX_SEED) + 1
        if seed not in used_seeds:
            used_seeds.add(seed)
            return seed


def steps_from_time(time_value, dt, name):
    if time_value is None:
        return None
    steps = round(time_value / dt)
    if steps <= 0 or abs(steps * dt - time_value) > 1.0e-12 * max(1.0, abs(time_value)):
        raise RuntimeError(f"{name} must be a positive integer multiple of dt")
    return steps


def nevery_from_time(time_value, dt, default_steps, name):
    steps = steps_from_time(time_value, dt, name)
    return steps if steps is not None else default_steps


def parse_axis(axis):
    if axis in ("x", "0"):
        return 0
    if axis in ("y", "1"):
        return 1
    raise RuntimeError("force-axis must be x|y|0|1")


def derived_eta(args):
    if args.eta is not None:
        return args.eta
    return f"{float(args.d0) * args.schmidt_number:.16g}"


def derived_force_amplitude(args):
    if args.force_amplitude is not None:
        return float(args.force_amplitude)
    axis_index = parse_axis(args.force_axis)
    length = float(args.length[axis_index])
    wave_number = 2.0 * math.pi * args.force_nk / length
    return float(args.target_amplitude) * float(args.d0) * float(args.free_energy_a) * wave_number * wave_number


def resolve_run_steps(args, dt):
    run_steps = steps_from_time(args.run_time, dt, "run-time")
    relax_steps = steps_from_time(args.relax_time, dt, "relax-time")
    measure_steps = steps_from_time(args.measure_time, dt, "measure-time")

    if run_steps is not None:
        if relax_steps is not None or measure_steps is not None:
            raise RuntimeError("--run-time cannot be combined with --relax-time or --measure-time")
        return 0, run_steps

    if measure_steps is None:
        measure_steps = args.run_steps
    if relax_steps is None:
        relax_steps = 0
    return relax_steps, measure_steps


def generate_input(args, seed, force_amplitude, relax_steps, measure_steps, time_series_nevery, profile_nevery, profile_nblock, time_series_path, profile_path):
    eta = derived_eta(args)
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
        f"fix                 3 order_parameter force/sine on component 0 axis {args.force_axis} nk {args.force_nk} amplitude {force_amplitude:.16e}",
        "",
        f"set                 density uniform value {args.density}",
        "set                 momentum all uniform value 0.0",
        f"set                 order_parameter all uniform value {args.order_parameter}",
        "",
        f"thermo              observe on progress off nevery {args.thermo_nevery}",
        "restart             off",
        "",
    ]

    if relax_steps > 0:
        lines += [
            f"run                 {relax_steps}",
            "",
        ]

    lines += [
        f"measure             1 time_series on nevery {time_series_nevery} file {time_series_path.as_posix()} target E_T E_K E_psi |psi[0]|^2 |d_psi[0]|^2 Jpsi[0]_y",
        f"measure             2 ave/profile on axis y nevery {profile_nevery} nblock {profile_nblock} file {profile_path.as_posix()} average block target psi[0] Jpsi[0]_y vx pi_xy",
        "",
        f"run                 {measure_steps}",
        "",
    ]
    return "\n".join(lines)


def write_config(path, args, eta, force_amplitude, relax_steps, measure_steps, time_series_nevery, profile_nevery, profile_nblock):
    with path.open("w", encoding="utf-8") as config:
        config.write("# Order-parameter Kolmogorov-flow replica configuration\n")
        for key, value in sorted(vars(args).items()):
            config.write(f"{key} {value}\n")
        config.write(f"derived_eta {eta}\n")
        config.write(f"derived_force_amplitude {force_amplitude:.16e}\n")
        config.write(f"relax_steps {relax_steps}\n")
        config.write(f"measure_steps {measure_steps}\n")
        config.write(f"time_series_nevery_effective {time_series_nevery}\n")
        config.write(f"profile_nevery_effective {profile_nevery}\n")
        config.write(f"profile_nblock_effective {profile_nblock}\n")
        config.write("seed_source os_entropy\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--samples", type=int, default=16)
    parser.add_argument("--sample-offset", type=int, default=0)

    parser.add_argument("--grid", nargs=2, type=int, default=[128, 128])
    parser.add_argument("--length", nargs=2, default=["4096", "4096"])
    parser.add_argument("--dealias", default="three_halves")
    parser.add_argument("--dt", default="4.0")
    parser.add_argument("--time-evolution", default="srk3/incompressible")

    parser.add_argument("--D0", dest="d0", default="0.12")
    parser.add_argument("--eta", default=None)
    parser.add_argument("--schmidt-number", type=float, default=1.0)
    parser.add_argument("--free-energy-a", default="1.0")
    parser.add_argument("--kBT", default="1.0")
    parser.add_argument("--density", default="1.0")
    parser.add_argument("--order-parameter", default="0.0")

    parser.add_argument("--force-axis", choices=("x", "y", "0", "1"), default="y")
    parser.add_argument("--force-nk", type=int, default=1)
    parser.add_argument("--force-amplitude", default=None)
    parser.add_argument("--target-amplitude", default="0.1")

    parser.add_argument("--run-steps", type=int, default=10000)
    parser.add_argument("--run-time", type=float)
    parser.add_argument("--relax-time", type=float)
    parser.add_argument("--measure-time", type=float)
    parser.add_argument("--thermo-nevery", type=int, default=100)
    parser.add_argument("--time-series-nevery", type=int, default=250)
    parser.add_argument("--profile-nevery", type=int, default=100)
    parser.add_argument("--profile-nblock", type=int, default=10000)
    parser.add_argument("--time-series-dtout", type=float)
    parser.add_argument("--profile-dtout", type=float)
    parser.add_argument("--profile-block-time", type=float)
    args = parser.parse_args()

    if args.samples <= 0:
        raise RuntimeError("samples must be positive")
    if args.sample_offset < 0:
        raise RuntimeError("sample-offset must be nonnegative")
    if args.run_steps <= 0:
        raise RuntimeError("run-steps must be positive")
    if args.thermo_nevery <= 0 or args.time_series_nevery <= 0 or args.profile_nevery <= 0 or args.profile_nblock <= 0:
        raise RuntimeError("nevery and nblock values must be positive")
    if args.force_nk <= 0:
        raise RuntimeError("force-nk must be positive")

    axis_index = parse_axis(args.force_axis)
    active_size = args.grid[axis_index]
    if args.force_nk >= active_size / 2:
        raise RuntimeError("force-nk must satisfy nk < active_N_axis/2")

    dt = float(args.dt)
    relax_steps, measure_steps = resolve_run_steps(args, dt)
    time_series_nevery = nevery_from_time(args.time_series_dtout, dt, args.time_series_nevery, "time-series-dtout")
    profile_nevery = nevery_from_time(args.profile_dtout, dt, args.profile_nevery, "profile-dtout")
    profile_nblock = nevery_from_time(args.profile_block_time, dt, args.profile_nblock, "profile-block-time")
    if profile_nblock % profile_nevery != 0:
        raise RuntimeError("profile-block-time must be an integer multiple of profile-dtout")

    output_root = Path(args.output_root) / f"nk_{args.force_nk:03d}"
    run_dir = output_root / "runs"
    result_dir = output_root / "results"
    log_dir = output_root / "logs"
    metadata_dir = output_root / "metadata"
    for path in (run_dir, result_dir, log_dir, metadata_dir):
        path.mkdir(parents=True, exist_ok=True)

    eta = derived_eta(args)
    force_amplitude = derived_force_amplitude(args)
    write_config(output_root / "config.dat", args, eta, force_amplitude, relax_steps, measure_steps, time_series_nevery, profile_nevery, profile_nblock)

    used_seeds = set()
    seeds_path = output_root / "seeds.dat"
    with seeds_path.open("w", encoding="utf-8") as seeds:
        seeds.write("# sample noise_seed input time_series profile stdout stderr\n")
        for local_sample in range(args.samples):
            sample = args.sample_offset + local_sample
            sid = f"{sample:03d}"
            noise_seed = random_seed(used_seeds)

            input_path = run_dir / f"input_{sid}.script"
            time_series_path = result_dir / f"time_series_{sid}.dat"
            profile_path = result_dir / f"profile_{sid}.dat"
            stdout_path = log_dir / f"stdout_{sid}.log"
            stderr_path = log_dir / f"stderr_{sid}.log"

            input_path.write_text(
                generate_input(
                    args,
                    noise_seed,
                    force_amplitude,
                    relax_steps,
                    measure_steps,
                    time_series_nevery,
                    profile_nevery,
                    profile_nblock,
                    time_series_path,
                    profile_path,
                ),
                encoding="utf-8",
            )

            seed_record = (
                f"{sid} {noise_seed} {input_path.as_posix()} "
                f"{time_series_path.as_posix()} {profile_path.as_posix()} "
                f"{stdout_path.as_posix()} {stderr_path.as_posix()}\n"
            )
            seeds.write(seed_record)
            (metadata_dir / f"seed_{sid}.dat").write_text(
                "# sample noise_seed input time_series profile stdout stderr\n" + seed_record,
                encoding="utf-8",
            )


if __name__ == "__main__":
    main()

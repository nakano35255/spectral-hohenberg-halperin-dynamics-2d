#!/usr/bin/env python3
"""Generate Kugui input scripts for order-parameter relaxation replicas."""

import argparse
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


def resolve_run_steps(args, dt):
    run_steps = steps_from_time(args.run_time, dt, "run-time")
    return run_steps if run_steps is not None else args.run_steps


def nevery_from_time(time_value, dt, default_steps, name):
    steps = steps_from_time(time_value, dt, name)
    return steps if steps is not None else default_steps


def generate_input(args, noise_seed, snapshot_prefix):
    dt = float(args.dt)
    run_steps = resolve_run_steps(args, dt)
    thermo_nevery = nevery_from_time(args.thermo_dtout, dt, args.thermo_nevery, "thermo-dtout")
    snapshot_nevery = nevery_from_time(args.snapshot_dtout, dt, args.snapshot_nevery, "snapshot-dtout")

    return "\n".join(
        [
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
            f"model transport     constant eta {args.eta} zeta {args.zeta} M[0,0] {args.mobility}",
            "",
            "fix                 1 momentum nonlinear off",
            "fix                 2 order_parameter nonlinear on",
            f"fix                 3 momentum noise on seed {noise_seed} kBT {args.kBT}",
            "",
            f"set                 density uniform value {args.density}",
            "set                 momentum all uniform value 0.0",
            (
                "set                 order_parameter all sine2d "
                f"base {args.order_parameter_base} amplitude {args.order_parameter_amplitude} "
                f"nkx {args.nkx} nky {args.nky}"
            ),
            "",
            f"thermo              observe on progress {args.progress} nevery {thermo_nevery}",
            f"measure             1 snapshot on nevery {snapshot_nevery} file {snapshot_prefix.as_posix()} space physical",
            "restart             off",
            "",
            f"run                 {run_steps}",
            "",
        ]
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--samples", type=int, default=16)
    parser.add_argument("--sample-offset", type=int, default=0)

    parser.add_argument("--grid", nargs=2, type=int, default=[256, 256])
    parser.add_argument("--length", nargs=2, default=["8192", "8192"])
    parser.add_argument("--dealias", default="three_halves")
    parser.add_argument("--dt", default="4.0")
    parser.add_argument("--time-evolution", default="srk3/incompressible")

    parser.add_argument("--free-energy-a", default="1.0")
    parser.add_argument("--eta", default="0.12")
    parser.add_argument("--zeta", default="0.0")
    parser.add_argument("--mobility", default="0.12")
    parser.add_argument("--kBT", default="1.0")
    parser.add_argument("--density", default="1.0")

    parser.add_argument("--order-parameter-base", default="0.0")
    parser.add_argument("--order-parameter-amplitude", default="0.1")
    parser.add_argument("--nkx", type=int, default=1)
    parser.add_argument("--nky", type=int, default=1)

    parser.add_argument("--run-steps", type=int, default=50000)
    parser.add_argument("--run-time", type=float)
    parser.add_argument("--thermo-nevery", type=int, default=5000)
    parser.add_argument("--thermo-dtout", type=float)
    parser.add_argument("--snapshot-nevery", type=int, default=1000)
    parser.add_argument("--snapshot-dtout", type=float)
    parser.add_argument("--progress", choices=("on", "off"), default="on")
    args = parser.parse_args()

    if args.samples <= 0:
        raise RuntimeError("samples must be positive")
    if args.sample_offset < 0:
        raise RuntimeError("sample-offset must be nonnegative")
    if args.run_steps <= 0 or args.thermo_nevery <= 0 or args.snapshot_nevery <= 0:
        raise RuntimeError("step counts must be positive")

    output_root = Path(args.output_root)
    run_dir = output_root / "runs"
    log_dir = output_root / "logs"
    sample_root = output_root / "samples"
    metadata_dir = output_root / "metadata"
    for path in (run_dir, log_dir, sample_root, metadata_dir):
        path.mkdir(parents=True, exist_ok=True)

    config_path = output_root / "config.dat"
    with config_path.open("w", encoding="utf-8") as config:
        config.write("# Order-parameter relaxation replica configuration\n")
        for key, value in sorted(vars(args).items()):
            config.write(f"{key} {value}\n")
        config.write("seed_source os_entropy\n")

    seeds_path = output_root / "seeds.dat"
    used_seeds = set()
    with seeds_path.open("w", encoding="utf-8") as seeds:
        seeds.write("# sample noise_seed input snapshot_prefix stdout stderr\n")
        for local_sample in range(args.samples):
            sample = args.sample_offset + local_sample
            sid = f"{sample:03d}"
            noise_seed = random_seed(used_seeds)

            sample_dir = sample_root / f"sample_{sid}"
            snapshot_dir = sample_dir / "snapshots"
            result_dir = sample_dir / "results"
            for path in (sample_dir, snapshot_dir, result_dir):
                path.mkdir(parents=True, exist_ok=True)

            input_path = run_dir / f"input_{sid}.script"
            stdout_path = log_dir / f"stdout_{sid}.log"
            stderr_path = log_dir / f"stderr_{sid}.log"
            snapshot_prefix = snapshot_dir / "physical"
            input_text = generate_input(args, noise_seed, snapshot_prefix)
            input_path.write_text(input_text, encoding="utf-8")

            seed_record = (
                f"{sid} {noise_seed} {input_path.as_posix()} "
                f"{snapshot_prefix.as_posix()} "
                f"{stdout_path.as_posix()} {stderr_path.as_posix()}\n"
            )
            seeds.write(seed_record)
            (metadata_dir / f"seed_{sid}.dat").write_text(
                "# sample noise_seed input snapshot_prefix stdout stderr\n" + seed_record,
                encoding="utf-8",
            )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Generate Ohtaka input scripts for incompressible YKGK replicas."""

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
    steps = round(time_value / dt)
    if steps <= 0 or abs(steps * dt - time_value) > 1.0e-12 * max(1.0, abs(time_value)):
        raise RuntimeError(f"{name} must be a positive integer multiple of dt")
    return steps


def nevery_from_time(time_value, dt, default_steps, name):
    if time_value is None:
        return default_steps
    return steps_from_time(time_value, dt, name)


def generate_input(args, sample, paths, noise_seed, init_seed):
    global_sample = args.sample_offset + sample
    dt = float(args.dt)
    run_steps = steps_from_time(args.run_time, dt, "run-time")
    ykgk_nevery = nevery_from_time(args.ykgk_dtout, dt, args.ykgk_nevery, "ykgk-dtout")
    ykgk_nblock = nevery_from_time(args.ykgk_block_time, dt, args.ykgk_nblock, "ykgk-block-time")

    if ykgk_nblock % ykgk_nevery != 0:
        raise RuntimeError("ykgk-block-time must be an integer multiple of ykgk-dtout")
    if run_steps % ykgk_nblock != 0:
        raise RuntimeError("run-time must be an integer multiple of ykgk-block-time")

    sid = f"{global_sample:03d}"
    ykgk_file = paths["result_dir"] / f"yokota_green_kubo_{sid}.dat"

    lines = [
        "dimension           2",
        "boundary            p p",
        "",
        f"grid                {args.grid[0]} {args.grid[1]}",
        f"length              {args.length[0]} {args.length[1]}",
        f"dealias             {args.dealias}",
        "",
        "order_parameters    0",
        "",
        f"timestep            {args.dt}",
        f"time_evolution      {args.time_evolution}",
        "",
        f"model transport     constant eta {args.eta} zeta 0.0",
        "",
        "fix                 1 momentum nonlinear on",
        f"fix                 2 momentum noise on seed {noise_seed} kBT {args.kBT}",
        "",
        f"set                 density uniform value {args.density}",
        f"set                 momentum all equilibrium/gaussian/incompressible kBT {args.kBT} seed {init_seed}",
        "",
        f"thermo              observe on progress off nevery {args.thermo_nevery}",
        "restart             off",
        "",
        f"measure             gk yokota_green_kubo on nevery {ykgk_nevery} nblock {ykgk_nblock} file {ykgk_file.as_posix()} mode diagonal",
        "",
        f"run                 {run_steps}",
        "",
    ]
    return "\n".join(lines), noise_seed, init_seed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--samples", type=int, default=32)
    parser.add_argument("--sample-offset", type=int, default=0)
    parser.add_argument("--eta", default="0.1")
    parser.add_argument("--dt", default="0.01")
    parser.add_argument("--run-time", type=float, default=25000.0)
    parser.add_argument("--grid", nargs=2, type=int, default=[256, 256])
    parser.add_argument("--length", nargs=2, default=["256.0", "256.0"])
    parser.add_argument("--dealias", default="three_halves")
    parser.add_argument("--time-evolution", default="srk3/incompressible")
    parser.add_argument("--kBT", default="1.0")
    parser.add_argument("--density", default="1.0")
    parser.add_argument("--thermo-nevery", type=int, default=10000)
    parser.add_argument("--ykgk-nevery", type=int, default=10000)
    parser.add_argument("--ykgk-dtout", type=float, default=100.0)
    parser.add_argument("--ykgk-nblock", type=int, default=500000)
    parser.add_argument("--ykgk-block-time", type=float, default=5000.0)
    parser.add_argument("--seed", type=int, default=12345)
    args = parser.parse_args()

    if args.samples <= 0:
        raise RuntimeError("samples must be positive")
    if args.sample_offset < 0:
        raise RuntimeError("sample-offset must be nonnegative")
    if args.thermo_nevery <= 0 or args.ykgk_nevery <= 0 or args.ykgk_nblock <= 0:
        raise RuntimeError("nevery and nblock values must be positive")

    output_root = Path(args.output_root)
    paths = {
        "run_dir": output_root / "runs",
        "result_dir": output_root / "results",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)

    config_path = output_root / "config.dat"
    with config_path.open("w", encoding="utf-8") as config:
        config.write("# Yokota Green-Kubo incompressible nonlinear run configuration\n")
        for key, value in sorted(vars(args).items()):
            config.write(f"{key} {value}\n")
        config.write("seed_source os_entropy\n")

    used_seeds = set()
    for sample in range(args.samples):
        global_sample = args.sample_offset + sample
        sid = f"{global_sample:03d}"
        input_path = paths["run_dir"] / f"input_{sid}.script"
        noise_seed = random_seed(used_seeds)
        init_seed = random_seed(used_seeds)
        text, noise_seed, init_seed = generate_input(args, sample, paths, noise_seed, init_seed)
        input_path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()

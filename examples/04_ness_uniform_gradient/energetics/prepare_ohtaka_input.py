#!/usr/bin/env python3
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


def label(prefix, value):
    return f"{prefix}_{value}"


def steps_from_time(time_value, dt, name):
    steps = round(time_value / dt)
    if steps <= 0 or abs(steps * dt - time_value) > 1.0e-12 * max(1.0, abs(time_value)):
        raise RuntimeError(f"{name} must be a positive integer multiple of dt")
    return steps


def common_header(args, seed):
    eta = args.eta if args.eta is not None else f"{float(args.d0) * args.schmidt_number:.16g}"
    return [
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
    ]


def initial_condition_block(args):
    return [
        f"set                 density uniform value {args.density}",
        "set                 momentum all uniform value 0.0",
        f"set                 order_parameter all uniform value {args.order_parameter}",
        "",
    ]


def generate_input(args, sample, segment, run_steps, time_series_nevery, paths, seed):
    lines = common_header(args, seed)

    sid = f"{sample:03d}"
    segment_id = f"{segment:03d}"
    prev_segment_id = f"{segment - 1:03d}"
    if args.legacy_single_segment:
        time_series_path = paths["result_dir"] / f"time_series_{sid}.dat"
        restart_write_path = paths["restart_dir"] / f"restart_{sid}.restart"
    else:
        time_series_path = paths["segment_dir"] / f"time_series_{sid}_seg{segment_id}.dat"
        restart_write_path = paths["restart_dir"] / f"restart_{sid}_seg{segment_id}.restart.tmp"

    if segment == 1:
        lines += initial_condition_block(args)
    else:
        lines += [
            f"restart             read file {(paths['restart_dir'] / f'restart_{sid}_seg{prev_segment_id}.restart').as_posix()}",
            "",
        ]

    lines += [
        f"thermo              observe on progress off nevery {args.thermo_nevery}",
        f"measure             1 time_series on nevery {time_series_nevery} file {time_series_path.as_posix()} target E_K |psi[0]|^2 |d_psi[0]|^2 Jpsi[0]_x",
        f"restart             write file {restart_write_path.as_posix()}",
        "",
        f"run                 {run_steps}",
        "",
    ]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--samples", type=int, default=288)
    parser.add_argument("--segments", type=int, default=2)
    parser.add_argument("--legacy-single-segment", action="store_true")
    parser.add_argument("--D0", dest="d0", required=True)
    parser.add_argument("--eta", default=None)
    parser.add_argument("--schmidt-number", type=float, default=1.0)
    parser.add_argument("--dt", required=True)
    parser.add_argument("--run-time-per-segment", type=float, default=200000000.0)
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
    parser.add_argument("--seed", type=int, default=12345)
    args = parser.parse_args()

    if args.samples <= 0:
        raise RuntimeError("samples must be positive")
    if args.legacy_single_segment:
        args.segments = 1
    if args.segments <= 0:
        raise RuntimeError("segments must be positive")
    if args.thermo_nevery <= 0 or args.time_series_nevery <= 0:
        raise RuntimeError("nevery values must be positive")

    dt = float(args.dt)
    run_steps = steps_from_time(args.run_time_per_segment, dt, "run-time-per-segment")
    time_series_nevery = args.time_series_nevery
    if args.time_series_dtout is not None:
        time_series_nevery = steps_from_time(args.time_series_dtout, dt, "time-series-dtout")

    case_dir = Path(label("D0", args.d0)) / label("dt", args.dt)
    output_root = Path(args.output_root)
    paths = {
        "run_dir": output_root / "runs" / case_dir,
        "result_dir": output_root / "results" / case_dir,
        "segment_dir": output_root / "segments" / case_dir,
        "restart_dir": output_root / "restarts" / case_dir,
        "seed_dir": output_root / "seeds" / case_dir,
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)

    used_seeds = set()
    seeds_path = paths["seed_dir"] / "seeds.dat"
    with seeds_path.open("w", encoding="utf-8") as seeds:
        seeds.write("# sample segment noise_seed input\n")
        for sample in range(args.samples):
            sid = f"{sample:03d}"
            for segment in range(1, args.segments + 1):
                segment_id = f"{segment:03d}"
                if args.legacy_single_segment:
                    input_path = paths["run_dir"] / f"input_{sid}.script"
                else:
                    input_path = paths["run_dir"] / f"input_{sid}_seg{segment_id}.script"
                seed = random_seed(used_seeds)
                input_path.write_text(
                    generate_input(args, sample, segment, run_steps, time_series_nevery, paths, seed),
                    encoding="utf-8",
                )
                seeds.write(f"{sid} {segment_id} {seed} {input_path.as_posix()}\n")


if __name__ == "__main__":
    main()

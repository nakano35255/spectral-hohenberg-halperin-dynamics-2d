#!/usr/bin/env python3
import argparse
from pathlib import Path


def steps_from_time(time_value, dt, name):
    steps = round(time_value / dt)
    if steps <= 0 or abs(steps * dt - time_value) > 1.0e-9 * max(1.0, abs(time_value)):
        raise RuntimeError(f"{name} must be a positive integer multiple of dt")
    return steps


def common_header(args, seed):
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
        f"model transport     constant eta {args.eta} zeta 0.0 M[0,0] {args.d0}",
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


def restart_index_name(index):
    return f"relax_{index:03d}.restart"


def generate_relax_input(args, segment, run_steps, nevery):
    output_root = Path(args.output_root)
    seed = args.seed + segment - 1
    lines = common_header(args, seed)

    restart_dir = output_root / "restarts"
    result_dir = output_root / "results"
    prev_restart = restart_dir / restart_index_name(segment - 1)
    next_restart_tmp = restart_dir / f"{restart_index_name(segment)}.tmp"
    time_series = result_dir / f"time_series_relax_{segment:03d}.dat"

    if segment == 1:
        lines += initial_condition_block(args)
    else:
        lines += [
            f"restart             read file {prev_restart.as_posix()}",
            "",
        ]

    lines += [
        f"thermo              observe on progress off nevery {args.thermo_nevery}",
        f"measure             ts time_series on nevery {nevery} file {time_series.as_posix()} target E_K |psi[0]|^2 |d_psi[0]|^2 Jpsi[0]_x",
        f"restart             write file {next_restart_tmp.as_posix()}",
        "",
        f"run                 {run_steps}",
        "",
    ]
    return "\n".join(lines)


def generate_budget_input(args, run_steps, budget_nevery, time_series_nevery):
    output_root = Path(args.output_root)
    seed = args.seed + 100000 + args.restart_index
    lines = common_header(args, seed)

    restart_dir = output_root / "restarts"
    result_dir = output_root / "results"
    source_restart = restart_dir / restart_index_name(args.restart_index)
    final_restart_tmp = restart_dir / f"budget_from_{args.restart_index:03d}.restart.tmp"
    time_series = result_dir / f"time_series_budget_from_{args.restart_index:03d}.dat"
    budget_shell = result_dir / f"budget_shell_from_{args.restart_index:03d}.dat"

    lines += [
        f"restart             read file {source_restart.as_posix()}",
        "",
        f"thermo              observe on progress off nevery {args.thermo_nevery}",
        f"measure             ts time_series on nevery {time_series_nevery} file {time_series.as_posix()} target E_K |psi[0]|^2 |d_psi[0]|^2 Jpsi[0]_x",
        f"measure             bs_shell budget/spectrum on component 0 nevery {budget_nevery} nblock {args.budget_nblock} file {budget_shell.as_posix()} mode shell average running",
        f"restart             write file {final_restart_tmp.as_posix()}",
        "",
        f"run                 {run_steps}",
        "",
    ]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["relax", "budget", "all"], default="all")
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--relax-segments", type=int, default=8)
    parser.add_argument("--restart-index", type=int, default=None)

    parser.add_argument("--D0", dest="d0", default="0.004")
    parser.add_argument("--eta", default="0.004")
    parser.add_argument("--dt", default="16.0")
    parser.add_argument("--grid", nargs=2, type=int, default=[1024, 1024])
    parser.add_argument("--length", nargs=2, default=["32768", "32768"])
    parser.add_argument("--dealias", default="three_halves")
    parser.add_argument("--time-evolution", default="srk3/incompressible")
    parser.add_argument("--free-energy-a", default="1.0")
    parser.add_argument("--gradient-amplitude", default="0.00006103515625")
    parser.add_argument("--kBT", default="1.0")
    parser.add_argument("--density", default="1.0")
    parser.add_argument("--order-parameter", default="0.0")
    parser.add_argument("--seed", type=int, default=12345)

    parser.add_argument("--relax-time-per-segment", type=float, default=100000000.0)
    parser.add_argument("--budget-time", type=float, default=100000000.0)
    parser.add_argument("--time-series-dtout", type=float, default=16384.0)
    parser.add_argument("--budget-nevery", type=int, default=20)
    parser.add_argument("--budget-nblock", type=int, default=200)
    parser.add_argument("--thermo-nevery", type=int, default=1000)
    args = parser.parse_args()

    if args.relax_segments <= 0:
        raise RuntimeError("relax-segments must be positive")
    if args.restart_index is None:
        args.restart_index = args.relax_segments
    if args.restart_index <= 0:
        raise RuntimeError("restart-index must be positive")

    dt = float(args.dt)
    relax_steps = steps_from_time(args.relax_time_per_segment, dt, "relax-time-per-segment")
    budget_steps = steps_from_time(args.budget_time, dt, "budget-time")
    time_series_nevery = steps_from_time(args.time_series_dtout, dt, "time-series-dtout")

    output_root = Path(args.output_root)
    run_dir = output_root / "runs"
    (output_root / "results").mkdir(parents=True, exist_ok=True)
    (output_root / "restarts").mkdir(parents=True, exist_ok=True)
    (output_root / "logs").mkdir(parents=True, exist_ok=True)
    run_dir.mkdir(parents=True, exist_ok=True)

    if args.mode in ("relax", "all"):
        for segment in range(1, args.relax_segments + 1):
            input_path = run_dir / f"input_relax_{segment:03d}.script"
            input_path.write_text(
                generate_relax_input(args, segment, relax_steps, time_series_nevery),
                encoding="utf-8",
            )

    if args.mode in ("budget", "all"):
        input_path = run_dir / f"input_budget_from_{args.restart_index:03d}.script"
        input_path.write_text(
            generate_budget_input(args, budget_steps, args.budget_nevery, time_series_nevery),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()

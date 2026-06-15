#!/usr/bin/env python3
import argparse
import os
import re
import struct
import sys
from pathlib import Path


WORK_EXAMPLE_BASE = Path("/work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/01_best_timestep_equilibrium")


def case_dir_from_context():
    if "SHHD_EXAMPLE_ROOT" in os.environ:
        return Path(os.environ["SHHD_EXAMPLE_ROOT"])
    return Path(sys.argv[0]).resolve().parent


def id_text(value):
    return value.replace("-", "m")


def steps_from_time(time_value, dt, name):
    if time_value is None:
        return None
    steps = round(time_value / dt)
    if steps <= 0 or abs(steps * dt - time_value) > 1.0e-12:
        raise RuntimeError(f"{name} must be a positive integer multiple of dt")
    return steps


def case_defaults(case_name):
    match = re.fullmatch(r"(compressible|incompressible)_transport(.+)", case_name)
    if not match:
        raise RuntimeError(f"cannot infer example/01 settings from case name: {case_name}")
    flow_type, transport = match.groups()
    compressible = flow_type == "compressible"
    return {
        "compressible": compressible,
        "time_evolution": f"srk3/{flow_type}",
        "transport": transport,
    }


def transport_values(args, defaults):
    if args.transport is None and args.eta is None and args.zeta is None and args.mobility is None:
        transport = defaults["transport"]
        eta = transport
        zeta = transport if defaults["compressible"] else None
        mobility = transport
    else:
        transport = args.transport or args.eta or args.zeta or args.mobility
        eta = args.eta if args.eta is not None else transport
        zeta = args.zeta if args.zeta is not None else (transport if defaults["compressible"] else None)
        mobility = args.mobility if args.mobility is not None else transport

    if eta != mobility or (zeta is not None and eta != zeta):
        raise RuntimeError("eta, optional zeta, and M[0,0] must be the same for this example")
    return eta, zeta, mobility


def generate_input(args, defaults, seed, time_series_path):
    dt_text = args.dt if args.dt is not None else "0.005"
    dt = float(dt_text)
    run_steps = steps_from_time(args.tmax, dt, "tmax") if args.tmax is not None else args.run_steps
    nevery = steps_from_time(args.dtout, dt, "dtout") if args.dtout is not None else args.time_series_nevery
    eta, zeta, mobility = transport_values(args, defaults)

    targets = "E_T E_K E_psi E_C" if defaults["compressible"] else "E_T E_K E_psi"
    transport_parts = ["model transport     constant", "eta", eta]
    if zeta is not None:
        transport_parts.extend(["zeta", zeta])
    transport_parts.extend(["M[0,0]", mobility])

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
        f"timestep            {dt_text}",
        f"time_evolution      {defaults['time_evolution']}",
        "",
    ]
    if defaults["compressible"]:
        lines.append(f"model thermo        linear_eos cT {args.cT}")
    lines.extend(
        [
            f"model free_energy   quadratic a[0] {args.free_energy_a}",
            " ".join(transport_parts),
            "",
            f"fix                 1 all noise on seed {seed} kBT {args.kBT}",
            "fix                 2 all nonlinear on",
            "",
            f"set                 density uniform value {args.density}",
            "set                 momentum all uniform value 0.0",
            f"set                 order_parameter all uniform value {args.order_parameter}",
            "",
            f"thermo              observe on progress off nevery {args.thermo_nevery}",
            f"measure             1 time_series on nevery {nevery} file {time_series_path.as_posix()} target {targets}",
            "restart             off",
            "",
            f"run                 {run_steps}",
            "",
        ]
    )
    return "\n".join(lines)


def main():
    case_dir = case_dir_from_context()
    case_name = case_dir.name
    defaults = case_defaults(case_name)
    output_root = Path(os.environ.get("SHHD_OUTPUT_ROOT", WORK_EXAMPLE_BASE / case_name))

    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=16)
    parser.add_argument("--grid", nargs=2, type=int, default=[64, 64])
    parser.add_argument("--length", nargs=2, default=["64.0", "64.0"])
    parser.add_argument("--dealias", default="three_halves")
    parser.add_argument("--dt")
    parser.add_argument("--tmax", type=float)
    parser.add_argument("--dtout", type=float)
    parser.add_argument("--transport")
    parser.add_argument("--eta")
    parser.add_argument("--zeta")
    parser.add_argument("--mobility")
    parser.add_argument("--cT", default="10.0")
    parser.add_argument("--free-energy-a", default="1.0")
    parser.add_argument("--kBT", default="1.0")
    parser.add_argument("--density", default="1.0")
    parser.add_argument("--order-parameter", default="0.0")
    parser.add_argument("--run-steps", type=int, default=10000)
    parser.add_argument("--thermo-nevery", type=int, default=5000)
    parser.add_argument("--time-series-nevery", type=int, default=50)
    args = parser.parse_args()

    if args.samples <= 0:
        raise RuntimeError("samples must be positive")
    if args.run_steps <= 0 or args.thermo_nevery <= 0 or args.time_series_nevery <= 0:
        raise RuntimeError("step counts must be positive")

    run_dir = output_root / "runs"
    result_dir = output_root / "results"
    if args.dt is not None:
        group = "dt" + id_text(args.dt)
        run_dir = run_dir / group
        result_dir = result_dir / group

    run_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)

    for sample in range(args.samples):
        seed = struct.unpack("<I", os.urandom(4))[0] % 2147483646 + 1
        sid = f"{sample:03d}"
        input_path = run_dir / f"input_{sid}.script"
        time_series = result_dir / f"time_series_{sid}.dat"
        input_path.write_text(generate_input(args, defaults, seed, time_series))


if __name__ == "__main__":
    main()

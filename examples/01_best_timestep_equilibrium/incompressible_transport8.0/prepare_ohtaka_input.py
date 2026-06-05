#!/usr/bin/env python3
import argparse
import os
import re
import struct
from pathlib import Path


CASE_NAME = Path(__file__).resolve().parent.name
HOME_EXAMPLE_BASE = Path("/home/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/01_best_timestep_equilibrium")
WORK_EXAMPLE_BASE = Path("/work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/01_best_timestep_equilibrium")
ROOT = Path(os.environ.get("SHHD_EXAMPLE_ROOT", HOME_EXAMPLE_BASE / CASE_NAME))
OUTPUT_ROOT = Path(os.environ.get("SHHD_OUTPUT_ROOT", WORK_EXAMPLE_BASE / CASE_NAME))


def id_text(value: str) -> str:
    return value.replace("-", "m")


def transport_values(text: str):
    match = re.search(r"^model\s+transport\s+constant\b(.*)$", text, flags=re.MULTILINE)
    if not match:
        raise RuntimeError("cannot find constant transport model in input.script")

    tokens = match.group(1).split()
    values = {}
    for index in range(0, len(tokens), 2):
        if index + 1 >= len(tokens):
            raise RuntimeError("transport coefficients must be key-value pairs")
        values[tokens[index]] = tokens[index + 1]

    if "eta" not in values or "M[0,0]" not in values:
        raise RuntimeError("constant transport requires eta and M[0,0] for this example")
    return values["eta"], values.get("zeta"), values["M[0,0]"]


def replace_transport(text: str, eta: str, zeta, mobility: str) -> str:
    if zeta is None:
        replacement = f"model transport     constant eta {eta} M[0,0] {mobility}"
    else:
        replacement = f"model transport     constant eta {eta} zeta {zeta} M[0,0] {mobility}"
    return re.sub(r"^model\s+transport\s+constant\b.*$", replacement, text, flags=re.MULTILINE)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=16)
    parser.add_argument("--dt")
    parser.add_argument("--tmax", type=float)
    parser.add_argument("--dtout", type=float)
    parser.add_argument("--transport")
    parser.add_argument("--eta")
    parser.add_argument("--zeta")
    parser.add_argument("--mobility")
    args = parser.parse_args()

    text = (ROOT / "input.script").read_text()
    eta0, zeta0, mobility0 = transport_values(text)
    if args.transport is not None or args.eta is not None or args.zeta is not None or args.mobility is not None:
        transport = args.transport or args.eta or args.zeta or args.mobility
        eta = args.eta if args.eta is not None else transport
        zeta = args.zeta if args.zeta is not None else (transport if zeta0 is not None else None)
        mobility = args.mobility if args.mobility is not None else transport
    else:
        eta, zeta, mobility = eta0, zeta0, mobility0

    if eta != mobility or (zeta is not None and eta != zeta):
        raise RuntimeError("eta, optional zeta, and M[0,0] must be the same for this example")

    run_dir = OUTPUT_ROOT / "runs"
    result_dir = OUTPUT_ROOT / "results"

    if args.dt:
        dt = float(args.dt)
        group = "dt" + id_text(args.dt)
        run_dir = run_dir / group
        result_dir = result_dir / group
        text = re.sub(r"^timestep\s+\S+", f"timestep            {args.dt}", text, flags=re.MULTILINE)

        if args.tmax is not None:
            steps = round(args.tmax / dt)
            if abs(steps * dt - args.tmax) > 1.0e-12:
                raise RuntimeError("tmax must be an integer multiple of dt")
            text = re.sub(r"^run\s+\S+", f"run                 {steps}", text, flags=re.MULTILINE)

        if args.dtout is not None:
            nevery = round(args.dtout / dt)
            if nevery <= 0 or abs(nevery * dt - args.dtout) > 1.0e-12:
                raise RuntimeError("dtout must be an integer multiple of dt")
            text = re.sub(r"(measure\s+\S+\s+time_series\s+on\s+nevery\s+)\S+", r"\g<1>" + str(nevery), text)

    run_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)

    if args.transport is not None or args.eta is not None or args.zeta is not None or args.mobility is not None:
        text = replace_transport(text, eta, zeta, mobility)

    for sample in range(args.samples):
        seed = struct.unpack("<I", os.urandom(4))[0] % 2147483646 + 1
        sid = f"{sample:03d}"
        input_path = run_dir / f"input_{sid}.script"
        time_series = result_dir / f"time_series_{sid}.dat"

        sample_text = re.sub(r"(fix\s+\S+\s+all\s+noise\s+on\s+.*?seed\s+)\S+", r"\g<1>" + str(seed), text)
        sample_text = re.sub(r"(measure\s+\S+\s+time_series\s+on\s+.*?file\s+)\S+", r"\g<1>" + time_series.as_posix(), sample_text)
        input_path.write_text(sample_text)


if __name__ == "__main__":
    main()

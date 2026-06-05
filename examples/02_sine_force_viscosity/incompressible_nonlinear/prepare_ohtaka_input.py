#!/usr/bin/env python3
import argparse
import os
import re
import struct
from pathlib import Path
from typing import Optional


CASE_NAME = Path(__file__).resolve().parent.name
HOME_EXAMPLE_BASE = Path("/home/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/02_sine_force_viscosity")
WORK_EXAMPLE_BASE = Path("/work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/02_sine_force_viscosity")
ROOT = Path(os.environ.get("SHHD_EXAMPLE_ROOT", HOME_EXAMPLE_BASE / CASE_NAME))
OUTPUT_ROOT = Path(os.environ.get("SHHD_OUTPUT_ROOT", WORK_EXAMPLE_BASE / CASE_NAME))


def replace_required(pattern: str, replacement: str, text: str, name: str) -> str:
    updated, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise RuntimeError(f"cannot replace {name} in input.script")
    return updated


def replace_run_steps(
    text: str,
    relax_steps: Optional[int],
    measure_steps: Optional[int],
    run_steps: Optional[int],
) -> str:
    matches = list(re.finditer(r"^run\s+\S+", text, flags=re.MULTILINE))
    if not matches:
        raise RuntimeError("input.script must contain at least one run command")

    if len(matches) == 1:
        if run_steps is None:
            if relax_steps is not None and measure_steps is not None:
                run_steps = relax_steps + measure_steps
            elif relax_steps is not None:
                run_steps = relax_steps
            else:
                run_steps = measure_steps

        if run_steps is None:
            return text
        match = matches[0]
        return text[:match.start()] + f"run                 {run_steps}" + text[match.end():]

    replacements = []
    if relax_steps is not None:
        replacements.append((matches[0], f"run                 {relax_steps}"))
    if measure_steps is not None:
        replacements.append((matches[1], f"run                 {measure_steps}"))

    for match, replacement in reversed(replacements):
        text = text[:match.start()] + replacement + text[match.end():]
    return text


def steps_from_time(time_value: Optional[float], dt: Optional[float], name: str) -> Optional[int]:
    if time_value is None:
        return None
    if dt is None:
        raise RuntimeError(f"{name} requires --dt")
    steps = round(time_value / dt)
    if steps <= 0 or abs(steps * dt - time_value) > 1.0e-12:
        raise RuntimeError(f"{name} must be a positive integer multiple of dt")
    return steps


def nevery_from_time(time_value: Optional[float], dt: Optional[float], name: str) -> Optional[int]:
    return steps_from_time(time_value, dt, name)


def update_template(text: str, args: argparse.Namespace) -> str:
    dt = float(args.dt) if args.dt is not None else None
    if args.dt is not None:
        text = replace_required(
            r"^timestep\s+\S+",
            f"timestep            {args.dt}",
            text,
            "timestep",
        )

    relax_steps = steps_from_time(args.relax_time, dt, "relax-time")
    measure_steps = steps_from_time(args.measure_time, dt, "measure-time")
    run_steps = steps_from_time(args.run_time, dt, "run-time")
    text = replace_run_steps(text, relax_steps, measure_steps, run_steps)

    time_series_nevery = nevery_from_time(args.time_series_dtout, dt, "time-series-dtout")
    if time_series_nevery is not None:
        text = replace_required(
            r"(measure\s+\S+\s+time_series\s+on\s+nevery\s+)\S+",
            r"\g<1>" + str(time_series_nevery),
            text,
            "time_series nevery",
        )

    profile_nevery = nevery_from_time(args.profile_dtout, dt, "profile-dtout")
    if profile_nevery is not None:
        text = replace_required(
            r"(measure\s+\S+\s+ave/profile\s+on\s+axis\s+\S+\s+nevery\s+)\S+",
            r"\g<1>" + str(profile_nevery),
            text,
            "ave/profile nevery",
        )

    profile_nblock = steps_from_time(args.profile_block_time, dt, "profile-block-time")
    if profile_nblock is not None:
        text = replace_required(
            r"(measure\s+\S+\s+ave/profile\s+on\s+.*?\bnblock\s+)\S+",
            r"\g<1>" + str(profile_nblock),
            text,
            "ave/profile nblock",
        )
        if profile_nevery is not None and profile_nblock % profile_nevery != 0:
            raise RuntimeError("profile-block-time must be an integer multiple of profile-dtout")

    return text


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=16)
    parser.add_argument("--dt")
    parser.add_argument("--run-time", type=float)
    parser.add_argument("--relax-time", type=float)
    parser.add_argument("--measure-time", type=float)
    parser.add_argument("--time-series-dtout", type=float)
    parser.add_argument("--profile-dtout", type=float)
    parser.add_argument("--profile-block-time", type=float)
    args = parser.parse_args()

    if args.samples <= 0:
        raise RuntimeError("samples must be positive")

    text = (ROOT / "input.script").read_text()
    text = update_template(text, args)

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

        sample_text = replace_required(
            r"(fix\s+\S+\s+momentum\s+noise\s+on\s+.*?\bseed\s+)\S+",
            r"\g<1>" + str(seed),
            text,
            "momentum noise seed",
        )
        sample_text = replace_required(
            r"(measure\s+\S+\s+time_series\s+on\s+.*?\bfile\s+)\S+",
            r"\g<1>" + time_series.as_posix(),
            sample_text,
            "time_series file",
        )
        sample_text = replace_required(
            r"(measure\s+\S+\s+ave/profile\s+on\s+.*?\bfile\s+)\S+",
            r"\g<1>" + profile.as_posix(),
            sample_text,
            "ave/profile file",
        )
        input_path.write_text(sample_text)


if __name__ == "__main__":
    main()

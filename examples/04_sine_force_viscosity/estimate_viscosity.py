#!/usr/bin/env python3
"""Estimate viscosity from sine-force ave/profile outputs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math


ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = ROOT / "examples" / "04_sine_force_viscosity"
RESULTS = EXAMPLE / "results"
SUMMARY = RESULTS / "viscosity_estimate.txt"


@dataclass
class CaseConfig:
    name: str
    input_script: Path
    profile_file: Path
    force_axis: str
    nk: int
    amplitude: float
    eta: float
    length_x: float
    length_y: float


@dataclass
class ProfileBlock:
    block: int
    step: int
    time: float
    samples: int
    columns: list[str]
    rows: list[list[float]]


def tokenize(line: str) -> list[str]:
    return line.split("#", 1)[0].split()


def parse_key_values(tokens: list[str]) -> dict[str, str]:
    values: dict[str, str] = {}
    for i in range(0, len(tokens), 2):
        if i + 1 < len(tokens):
            values[tokens[i]] = tokens[i + 1]
    return values


def parse_input(path: Path, name: str) -> CaseConfig:
    length_x = 0.0
    length_y = 0.0
    eta = math.nan
    profile_file: Path | None = None
    force_axis = ""
    nk = 0
    amplitude = math.nan

    with path.open() as handle:
        for raw_line in handle:
            tokens = tokenize(raw_line)
            if not tokens:
                continue

            if tokens[0] == "length":
                length_x = float(tokens[1])
                length_y = float(tokens[2])
                continue

            if tokens[:3] == ["model", "transport", "constant"]:
                values = parse_key_values(tokens[3:])
                eta = float(values["eta"])
                continue

            if (
                len(tokens) >= 6
                and tokens[0] == "fix"
                and tokens[2] in {"momentum", "j"}
                and tokens[3] == "force/sine"
                and tokens[4] == "on"
            ):
                values = parse_key_values(tokens[5:])
                if values.get("component") in {"x", "0"}:
                    force_axis = values["axis"]
                    nk = int(values["nk"])
                    amplitude = float(values["amplitude"])
                continue

            if len(tokens) >= 5 and tokens[0] == "measure" and tokens[2] == "ave/profile" and tokens[3] == "on":
                values = parse_key_values(tokens[4:tokens.index("target")] if "target" in tokens else tokens[4:])
                if "file" in values:
                    profile_file = ROOT / values["file"]
                continue

    if profile_file is None:
        raise ValueError(f"{path} does not define a measure ave/profile file")
    if force_axis not in {"x", "y", "0", "1"}:
        raise ValueError(f"{path} does not define a momentum x force/sine")
    if nk <= 0 or not math.isfinite(amplitude):
        raise ValueError(f"{path} has invalid force/sine parameters")
    if not math.isfinite(eta):
        raise ValueError(f"{path} does not define model transport constant eta")
    if length_x <= 0.0 or length_y <= 0.0:
        raise ValueError(f"{path} has invalid domain length")

    return CaseConfig(
        name=name,
        input_script=path,
        profile_file=profile_file,
        force_axis=force_axis,
        nk=nk,
        amplitude=amplitude,
        eta=eta,
        length_x=length_x,
        length_y=length_y,
    )


def parse_block_metadata(line: str) -> dict[str, str]:
    fields = line[1:].split()
    return {fields[i]: fields[i + 1] for i in range(0, len(fields) - 1, 2)}


def read_profile(path: Path) -> ProfileBlock:
    if not path.exists():
        raise FileNotFoundError(f"missing profile file: {path}")

    blocks: list[ProfileBlock] = []
    current_meta: dict[str, str] | None = None
    current_columns: list[str] | None = None
    current_rows: list[list[float]] = []

    def finish_current() -> None:
        nonlocal current_meta, current_columns, current_rows
        if current_meta is None:
            return
        if current_columns is None:
            raise ValueError(f"profile block is missing columns in {path}")
        if not current_rows:
            raise ValueError(f"profile block is empty in {path}")
        blocks.append(
            ProfileBlock(
                block=int(current_meta["block"]),
                step=int(current_meta["step"]),
                time=float(current_meta["time"]),
                samples=int(current_meta["samples"]),
                columns=current_columns,
                rows=current_rows,
            )
        )
        current_meta = None
        current_columns = None
        current_rows = []

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                finish_current()
                continue
            if line.startswith("# block"):
                finish_current()
                current_meta = parse_block_metadata(line)
                continue
            if line.startswith("# columns"):
                if current_meta is None:
                    raise ValueError(f"columns appeared before block metadata in {path}")
                current_columns = line[1:].split()[1:]
                continue
            if line.startswith("#"):
                continue

            if current_meta is None or current_columns is None:
                raise ValueError(f"profile data appeared before block header in {path}")

            values = [float(value) for value in line.split()]
            if len(values) != len(current_columns):
                raise ValueError(f"unexpected profile row width in {path}: {line}")
            current_rows.append(values)

    finish_current()

    if not blocks:
        raise ValueError(f"no profile blocks found in {path}")
    return blocks[-1]


def projection(values: list[float], basis: list[float]) -> float:
    denominator = sum(b * b for b in basis)
    if denominator == 0.0:
        raise ValueError("zero basis norm")
    return sum(v * b for v, b in zip(values, basis)) / denominator


def safe_ratio(numerator: float, denominator: float, eps: float = 1.0e-12) -> float:
    if abs(denominator) <= eps:
        return math.nan
    return numerator / denominator


def write_profile_diagnostics(
    cfg: CaseConfig,
    block: ProfileBlock,
    coords: list[float],
    vx: list[float],
    pi_xy: list[float],
    k: float,
) -> Path:
    path = RESULTS / f"{cfg.name}_profile_diagnostics.dat"

    with path.open("w") as handle:
        handle.write("# sine-force viscosity profile diagnostics\n")
        handle.write(f"# input {cfg.input_script.relative_to(ROOT)}\n")
        handle.write(f"# source_profile {cfg.profile_file.relative_to(ROOT)}\n")
        handle.write(
            f"# block {block.block} step {block.step} time {block.time:.16e} "
            f"samples {block.samples}\n"
        )
        handle.write(f"# force f_x = {cfg.amplitude:.16e} sin({k:.16e} coord)\n")
        handle.write("# nan means the local denominator is too close to zero\n")
        handle.write(
            "# columns coord f_x vx pi_xy sin_ky cos_ky "
            "f_over_sin vx_over_sin pi_xy_over_cos f_over_v "
            "eta_force_velocity eta_stress_gradient eta_stress_amplitude\n"
        )

        for coord, v_value, stress_value in zip(coords, vx, pi_xy):
            sin_value = math.sin(k * coord)
            cos_value = math.cos(k * coord)
            force_value = cfg.amplitude * sin_value

            f_over_sin = safe_ratio(force_value, sin_value)
            vx_over_sin = safe_ratio(v_value, sin_value)
            pi_over_cos = safe_ratio(stress_value, cos_value)
            f_over_v = safe_ratio(force_value, v_value)

            eta_force_velocity = safe_ratio(force_value, k * k * v_value)

            gradient_from_sine = k * vx_over_sin * cos_value
            eta_stress_gradient = (
                safe_ratio(-stress_value, gradient_from_sine)
                if math.isfinite(gradient_from_sine)
                else math.nan
            )
            eta_stress_amplitude = (
                safe_ratio(-pi_over_cos, k * vx_over_sin)
                if math.isfinite(pi_over_cos) and math.isfinite(vx_over_sin)
                else math.nan
            )

            handle.write(
                f"{coord:.16e} "
                f"{force_value:.16e} "
                f"{v_value:.16e} "
                f"{stress_value:.16e} "
                f"{sin_value:.16e} "
                f"{cos_value:.16e} "
                f"{f_over_sin:.16e} "
                f"{vx_over_sin:.16e} "
                f"{pi_over_cos:.16e} "
                f"{f_over_v:.16e} "
                f"{eta_force_velocity:.16e} "
                f"{eta_stress_gradient:.16e} "
                f"{eta_stress_amplitude:.16e}\n"
            )

    return path


def write_profile_figure(
    cfg: CaseConfig,
    block: ProfileBlock,
    coords: list[float],
    vx: list[float],
    pi_xy: list[float],
    k: float,
    vx_sin: float,
) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    force = [cfg.amplitude * math.sin(k * coord) for coord in coords]
    eta_from_force_velocity = [
        safe_ratio(f_value, k * k * v_value) for f_value, v_value in zip(force, vx)
    ]
    eta_from_stress_velocity = [
        safe_ratio(-stress_value, k * vx_sin * math.cos(k * coord))
        for coord, stress_value in zip(coords, pi_xy)
    ]
    coord_label = "x" if cfg.force_axis in {"x", "0"} else "y"
    path = RESULTS / f"{cfg.name}_profile.png"

    fig, axes = plt.subplots(5, 1, figsize=(7.0, 10.5), sharex=True)
    fig.suptitle(
        f"{cfg.name}: block {block.block}, step {block.step}, samples {block.samples}",
        fontsize=12,
    )

    axes[0].plot(coords, force, "o-", markersize=3.0, linewidth=1.2)
    axes[0].set_ylabel(r"$f_x$")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(coords, vx, "o-", markersize=3.0, linewidth=1.2)
    axes[1].set_ylabel(r"$v_x$")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(coords, pi_xy, "o-", markersize=3.0, linewidth=1.2)
    axes[2].set_ylabel(r"$\Pi_{xy}$")
    axes[2].grid(True, alpha=0.3)

    axes[3].plot(
        coords,
        eta_from_force_velocity,
        "o-",
        markersize=3.0,
        linewidth=1.2,
        label=r"$f_x/(k^2 v_x)$",
    )
    axes[3].axhline(cfg.eta, color="black", linestyle="--", linewidth=1.0, label=r"$\eta$")
    axes[3].set_ylabel(r"$\eta_f$")
    axes[3].grid(True, alpha=0.3)
    axes[3].legend(loc="best", fontsize=9)

    axes[4].plot(
        coords,
        eta_from_stress_velocity,
        "o-",
        markersize=3.0,
        linewidth=1.2,
        label=r"$-\Pi_{xy}/(k A_v \cos ky)$",
    )
    axes[4].axhline(cfg.eta, color="black", linestyle="--", linewidth=1.0, label=r"$\eta$")
    axes[4].set_ylabel(r"$\eta_\Pi$")
    axes[4].set_xlabel(coord_label)
    axes[4].grid(True, alpha=0.3)
    axes[4].legend(loc="best", fontsize=9)

    for axis in axes:
        axis.ticklabel_format(axis="y", style="plain", useOffset=False)

    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)

    return path


def estimate_case(cfg: CaseConfig) -> str:
    block = read_profile(cfg.profile_file)
    if "coord" not in block.columns or "vx" not in block.columns or "pi_xy" not in block.columns:
        raise ValueError(f"{cfg.profile_file} must contain coord, vx, and pi_xy columns")

    coord_index = block.columns.index("coord")
    vx_index = block.columns.index("vx")
    stress_index = block.columns.index("pi_xy")

    coords = [row[coord_index] for row in block.rows]
    vx = [row[vx_index] for row in block.rows]
    pi_xy = [row[stress_index] for row in block.rows]

    length = cfg.length_x if cfg.force_axis in {"x", "0"} else cfg.length_y
    k = 2.0 * math.pi * float(cfg.nk) / length

    sin_basis = [math.sin(k * x) for x in coords]
    cos_basis = [math.cos(k * x) for x in coords]

    vx_sin = projection(vx, sin_basis)
    pi_cos = projection(pi_xy, cos_basis)

    eta_from_velocity = cfg.amplitude / (k * k * vx_sin)
    eta_from_stress = -pi_cos / (k * vx_sin)
    force_from_stress = -k * pi_cos
    diagnostics = write_profile_diagnostics(cfg, block, coords, vx, pi_xy, k)
    figure = write_profile_figure(cfg, block, coords, vx, pi_xy, k, vx_sin)

    return "\n".join(
        [
            f"[{cfg.name}]",
            f"input             : {cfg.input_script.relative_to(ROOT)}",
            f"profile           : {cfg.profile_file.relative_to(ROOT)}",
            f"diagnostics       : {diagnostics.relative_to(ROOT)}",
            f"figure            : {figure.relative_to(ROOT)}",
            f"block             : {block.block}",
            f"step              : {block.step}",
            f"time              : {block.time:.16e}",
            f"samples           : {block.samples}",
            f"eta input         : {cfg.eta:.16e}",
            f"force amplitude   : {cfg.amplitude:.16e}",
            f"k                 : {k:.16e}",
            f"vx sin amplitude  : {vx_sin:.16e}",
            f"pi_xy cos amp     : {pi_cos:.16e}",
            f"force from stress : {force_from_stress:.16e}",
            f"eta from vx       : {eta_from_velocity:.16e}",
            f"eta from pi_xy    : {eta_from_stress:.16e}",
            f"relerr eta(vx)    : {abs(eta_from_velocity - cfg.eta) / abs(cfg.eta):.16e}",
            f"relerr eta(pi_xy) : {abs(eta_from_stress - cfg.eta) / abs(cfg.eta):.16e}",
        ]
    )


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    cases = [
        parse_input(EXAMPLE / "input_incompressible.script", "incompressible"),
        parse_input(EXAMPLE / "input_compressible.script", "compressible"),
    ]

    report = "\n\n".join(estimate_case(case) for case in cases) + "\n"
    SUMMARY.write_text(report)
    print(report)
    print(f"wrote {SUMMARY.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

import argparse
import glob
import math
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


Mode = Tuple[int, int]


@dataclass
class SineInitialCondition:
    base: float = 0.0
    amplitude: float = 0.0
    nkx: int = 0
    nky: int = 0


@dataclass
class Config:
    nx: int = 0
    ny: int = 0
    lx: float = 2.0 * math.pi
    ly: float = 2.0 * math.pi
    dt: float = 0.0
    integrator: str = ""
    a: Dict[int, float] = field(default_factory=dict)
    b: Dict[int, float] = field(default_factory=dict)
    u: Dict[int, float] = field(default_factory=dict)
    mobility: Dict[int, float] = field(default_factory=dict)
    sine: Dict[int, SineInitialCondition] = field(default_factory=dict)
    spectral_prefixes: List[str] = field(default_factory=list)


@dataclass
class Snapshot:
    path: Path
    step: int
    time: float
    nkx: int
    nky: int
    columns: Dict[str, int]
    data: np.ndarray


def strip_comment(line: str) -> str:
    return line.split("#", 1)[0].strip()


def key_values(tokens: List[str], start: int) -> Dict[str, str]:
    if (len(tokens) - start) % 2 != 0:
        raise ValueError("key-value arguments must appear in pairs: " + " ".join(tokens))
    return {tokens[i]: tokens[i + 1] for i in range(start, len(tokens), 2)}


def parse_indexed_key(key: str, name: str) -> int:
    match = re.fullmatch(rf"{re.escape(name)}\[(\d+)\]", key)
    return int(match.group(1)) if match else -1


def parse_mobility_key(key: str) -> Tuple[int, int]:
    match = re.fullmatch(r"M\[(\d+),(\d+)\]", key)
    if not match:
        return (-1, -1)
    return (int(match.group(1)), int(match.group(2)))


def parse_input_script(path: Path) -> Config:
    cfg = Config()

    with path.open() as handle:
        for raw_line in handle:
            line = strip_comment(raw_line)
            if not line:
                continue

            tokens = line.split()
            command = tokens[0]

            if command == "grid":
                cfg.nx = int(tokens[1])
                cfg.ny = int(tokens[2])
            elif command == "length":
                cfg.lx = float(tokens[1])
                cfg.ly = float(tokens[2])
            elif command == "timestep":
                cfg.dt = float(tokens[1])
            elif command == "time_evolution":
                cfg.integrator = tokens[1]
            elif command == "model":
                category = tokens[1]
                model_type = tokens[2]
                values = key_values(tokens, 3)

                if category == "free_energy" and model_type in ("phi4", "coeff"):
                    for key, raw_value in values.items():
                        value = float(raw_value)
                        index = parse_indexed_key(key, "a")
                        if index >= 0:
                            cfg.a[index] = value
                            continue
                        index = parse_indexed_key(key, "b")
                        if index >= 0:
                            cfg.b[index] = value
                            continue
                        index = parse_indexed_key(key, "u")
                        if index >= 0:
                            cfg.u[index] = value
                            continue
                elif category == "transport" and model_type in ("constant", "coeff"):
                    for key, raw_value in values.items():
                        row, col = parse_mobility_key(key)
                        if row == col and row >= 0:
                            cfg.mobility[row] = float(raw_value)
            elif command == "set" and tokens[1] == "order_parameter":
                component = int(tokens[2])
                initial_type = tokens[3]
                if initial_type == "sine":
                    values = key_values(tokens, 4)
                    cfg.sine[component] = SineInitialCondition(
                        base=float(values.get("base", "0.0")),
                        amplitude=float(values["amplitude"]),
                        nkx=int(values["nkx"]),
                        nky=int(values.get("nky", "0")),
                    )
            elif command == "measure" and tokens[2] == "snapshot" and tokens[3] == "on":
                values = key_values(tokens, 4)
                if values.get("space", "physical") == "spectral":
                    cfg.spectral_prefixes.append(values["file"])

    validate_config(cfg)
    return cfg


def validate_config(cfg: Config) -> None:
    if cfg.nx <= 0 or cfg.ny <= 0:
        raise ValueError("grid must be specified")
    if cfg.dt <= 0.0:
        raise ValueError("timestep must be specified")
    if cfg.integrator != "euler/quiescent":
        raise ValueError("this one-step comparison requires time_evolution euler/quiescent")
    if 0 not in cfg.sine:
        raise ValueError("set order_parameter 0 sine is required")
    if cfg.u.get(0, 0.0) == 0.0:
        raise ValueError("this nonlinear comparison requires nonzero u[0]")


def parse_snapshot(path: Path) -> Snapshot:
    step = -1
    time = math.nan
    nkx = -1
    nky = -1
    columns: Dict[str, int] = {}

    with path.open() as handle:
        for raw_line in handle:
            if not raw_line.startswith("#"):
                break

            parts = raw_line[1:].strip().split()
            if not parts:
                continue
            if parts[0] == "step":
                step = int(parts[1])
                time = float(parts[3])
            elif parts[0] == "nkx":
                nkx = int(parts[1])
                nky = int(parts[3])
            elif parts[0] == "kx_index":
                columns = {name: i for i, name in enumerate(parts)}

    if step < 0 or nkx <= 0 or nky <= 0 or not columns:
        raise ValueError("invalid spectral snapshot: " + str(path))

    data = np.loadtxt(path, comments="#", ndmin=2)
    return Snapshot(path, step, time, nkx, nky, columns, data)


def resolve_snapshot_files(cfg: Config, input_path: Path, prefix_arg: str) -> List[Path]:
    repo_root = input_path.resolve().parents[2]
    input_dir = input_path.resolve().parent

    prefixes = [prefix_arg] if prefix_arg else cfg.spectral_prefixes
    if not prefixes:
        raise ValueError("no spectral snapshot prefix was found; pass --prefix explicitly")

    files: List[Path] = []
    for prefix in prefixes:
        raw = Path(prefix)
        candidates = [raw] if raw.is_absolute() else [Path.cwd() / raw, repo_root / raw, input_dir / raw]

        for candidate in candidates:
            matches = [Path(p) for p in glob.glob(str(candidate) + "_step*.snapshot")]
            if matches:
                files.extend(matches)
                break

    files = sorted(set(files), key=lambda p: str(p))
    files.sort(key=lambda p: int(re.search(r"_step(\d+)\.snapshot$", p.name).group(1)))
    return files


def signed_ky(ky_index: np.ndarray, ny: int) -> np.ndarray:
    return np.where(ky_index <= ny // 2, ky_index, ky_index - ny)


def active_mask(cfg: Config, nkx: int, nky: int) -> np.ndarray:
    kx_index = np.arange(nkx)[None, :]
    ky_index = np.arange(nky)[:, None]
    return (kx_index < cfg.nx // 2) & (np.abs(signed_ky(ky_index, cfg.ny)) < cfg.ny // 2)


def wave_numbers(cfg: Config, nkx: int, nky: int) -> np.ndarray:
    kx_index = np.arange(nkx)[None, :]
    ky_index = np.arange(nky)[:, None]
    kx = 2.0 * math.pi * kx_index / cfg.lx
    ky = 2.0 * math.pi * signed_ky(ky_index, cfg.ny) / cfg.ly
    return kx * kx + ky * ky


def initial_sine_hat(cfg: Config, nkx: int, nky: int) -> np.ndarray:
    sine = cfg.sine[0]
    grid_size = float(cfg.nx) * float(cfg.ny)
    out = np.zeros((nky, nkx), dtype=np.complex128)

    out[0, 0] = sine.base * grid_size

    ky = sine.nky if sine.nky >= 0 else cfg.ny + sine.nky
    out[ky, sine.nkx] = -0.5j * sine.amplitude * grid_size

    if sine.nkx == 0:
        conjugate_ky = (-sine.nky) % cfg.ny
        if conjugate_ky != ky:
            out[conjugate_ky, 0] = np.conj(out[ky, 0])

    return out


def expected_one_step(cfg: Config, snapshot: Snapshot) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    init_hat = initial_sine_hat(cfg, snapshot.nkx, snapshot.nky)
    psi = np.fft.irfft2(init_hat, s=(cfg.ny, cfg.nx)).real
    mu_nonlinear_hat = np.fft.rfft2(cfg.u.get(0, 0.0) * psi**3)

    k2 = wave_numbers(cfg, snapshot.nkx, snapshot.nky)
    linear_mu_hat = (cfg.a.get(0, 0.0) + cfg.b.get(0, 0.0) * k2) * init_hat
    rhs_hat = -cfg.mobility.get(0, 0.0) * k2 * (linear_mu_hat + mu_nonlinear_hat)
    rhs_hat = np.where(active_mask(cfg, snapshot.nkx, snapshot.nky), rhs_hat, 0.0)

    return init_hat + cfg.dt * rhs_hat, init_hat, rhs_hat, psi


def snapshot_hat(snapshot: Snapshot) -> np.ndarray:
    values = (
        snapshot.data[:, snapshot.columns["psi_com0_real"]]
        + 1j * snapshot.data[:, snapshot.columns["psi_com0_imag"]]
    )
    return values.reshape((snapshot.nky, snapshot.nkx))


def mode_value(field: np.ndarray, mode: Mode) -> complex:
    return field[mode[1] % field.shape[0], mode[0]]


def write_validation_plot(
    cfg: Config,
    snapshot: Snapshot,
    actual: np.ndarray,
    expected: np.ndarray,
    initial: np.ndarray,
    rhs: np.ndarray,
    output: Path,
) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sine = cfg.sine[0]
    fundamental = (sine.nkx, sine.nky)
    third = (3 * sine.nkx, 3 * sine.nky)
    modes = [fundamental, third]

    actual_rhs = (actual - initial) / cfg.dt
    expected_rhs = rhs
    y_index = 0
    x = np.arange(cfg.nx) * cfg.lx / cfg.nx
    rhs_physical = np.fft.irfft2(expected_rhs, s=(cfg.ny, cfg.nx)).real
    actual_rhs_physical = np.fft.irfft2(actual_rhs, s=(cfg.ny, cfg.nx)).real

    fig = plt.figure(figsize=(12.0, 8.0))
    gs = fig.add_gridspec(2, 2, width_ratios=(1.05, 1.0), height_ratios=(1.0, 1.0))
    ax_rhs = fig.add_subplot(gs[0, 0])
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_spectrum = fig.add_subplot(gs[1, 0])
    ax_error = fig.add_subplot(gs[1, 1])

    ax_rhs.plot(x, actual_rhs_physical[y_index, :], ".", markersize=3, label="snapshot")
    ax_rhs.plot(x, rhs_physical[y_index, :], "-", linewidth=1.0, label="expected")
    ax_rhs.set_xlabel("x")
    ax_rhs.set_ylabel(r"$\partial_t\psi(x,0)$")
    ax_rhs.grid(True, linewidth=0.4, alpha=0.5)
    ax_rhs.legend(fontsize=8)

    scale = max(abs(mode_value(initial, fundamental)), 1.0e-300)
    labels = [f"{mode}" for mode in modes]
    actual_values = [abs(mode_value(actual - initial, mode)) / scale for mode in modes]
    expected_values = [abs(mode_value(expected - initial, mode)) / scale for mode in modes]
    xloc = np.arange(len(modes))
    width = 0.36
    ax_bar.bar(xloc - width / 2, actual_values, width, label="snapshot")
    ax_bar.bar(xloc + width / 2, expected_values, width, label="expected")
    ax_bar.set_xticks(xloc, labels)
    ax_bar.set_ylabel(r"$|\Delta\hat{\psi}| / |\hat{\psi}_{k}(0)|$")
    ax_bar.grid(True, axis="y", linewidth=0.4, alpha=0.5)
    ax_bar.legend(fontsize=8)

    amplitude = np.abs(actual - initial) / scale
    amplitude = np.fft.fftshift(amplitude, axes=0)
    image_spectrum = ax_spectrum.imshow(
        np.log10(amplitude + 1.0e-30),
        origin="lower",
        aspect="auto",
        extent=[0, snapshot.nkx - 1, -snapshot.nky // 2, snapshot.nky // 2 - 1],
        vmin=-16,
        vmax=-1,
    )
    ax_spectrum.set_xlabel("kx_index")
    ax_spectrum.set_ylabel("signed ky_index")
    ax_spectrum.set_title(r"$\log_{10} |\Delta\hat{\psi}| / |\hat{\psi}_{k}(0)|$")
    fig.colorbar(image_spectrum, ax=ax_spectrum, shrink=0.85)

    error = np.abs(actual - expected) / scale
    error = np.fft.fftshift(error, axes=0)
    image_error = ax_error.imshow(
        np.log10(error + 1.0e-30),
        origin="lower",
        aspect="auto",
        extent=[0, snapshot.nkx - 1, -snapshot.nky // 2, snapshot.nky // 2 - 1],
        vmin=-16,
        vmax=-10,
    )
    ax_error.set_xlabel("kx_index")
    ax_error.set_ylabel("signed ky_index")
    ax_error.set_title(r"$\log_{10}$ normalized one-step error")
    fig.colorbar(image_error, ax=ax_error, shrink=0.85)

    fig.suptitle("Example 03: nonlinear phi4 sine harmonic generation")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def analyze(cfg: Config, snapshot_paths: List[Path], args: argparse.Namespace) -> bool:
    snapshot = parse_snapshot(snapshot_paths[-1])
    actual = snapshot_hat(snapshot)
    expected, initial, rhs, _ = expected_one_step(cfg, snapshot)

    error = np.abs(actual - expected)
    abs_error = float(np.max(error))
    fundamental_scale = max(abs(mode_value(initial, (cfg.sine[0].nkx, cfg.sine[0].nky))), 1.0e-300)
    norm_error = float(np.max(error / fundamental_scale))

    sine = cfg.sine[0]
    fundamental = (sine.nkx, sine.nky)
    third = (3 * sine.nkx, 3 * sine.nky)
    generated = np.abs((actual - initial) / fundamental_scale)
    expected_generated = np.abs((expected - initial) / fundamental_scale)

    allowed = np.zeros_like(generated, dtype=bool)
    for mode in (fundamental, third):
        allowed[mode[1] % snapshot.nky, mode[0]] = True
    unexpected = generated.copy()
    unexpected[allowed] = 0.0
    max_unexpected = float(np.max(unexpected))

    print("Input")
    print(f"  grid                 : {cfg.nx} x {cfg.ny}")
    print(f"  length               : {cfg.lx:g} x {cfg.ly:g}")
    print(f"  dt                   : {cfg.dt:g}")
    print(f"  time_evolution       : {cfg.integrator}")
    print(f"  phi4                 : a={cfg.a.get(0, 0.0):g} b={cfg.b.get(0, 0.0):g} u={cfg.u.get(0, 0.0):g}")
    print(f"  mobility             : M={cfg.mobility.get(0, 0.0):g}")
    print(f"  sine mode            : {fundamental}")
    print()
    print("Generated harmonics")
    for mode in (fundamental, third):
        actual_value = abs(mode_value(actual - initial, mode)) / fundamental_scale
        expected_value = abs(mode_value(expected - initial, mode)) / fundamental_scale
        print(f"  {mode!s:12s} snapshot={actual_value:.16e} expected={expected_value:.16e}")
    print(f"  max unexpected       : {max_unexpected:.16e}")
    print()
    print("Worst errors")
    print(f"  abs error            : {abs_error:.16e} (file={snapshot.path})")
    print(f"  normalized error     : {norm_error:.16e}")
    print()

    if args.plot:
        output = Path(args.plot).resolve()
        write_validation_plot(cfg, snapshot, actual, expected, initial, rhs, output)
        print(f"plot: {output}")

    passed = abs_error <= args.tol_abs and norm_error <= args.tol_norm and max_unexpected <= args.tol_unexpected
    print("PASS" if passed else "FAIL")
    return passed


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Validate one-step nonlinear phi4 harmonic generation from a sine field.")
    parser.add_argument("--input", default=str(script_dir / "input.script"))
    parser.add_argument("--prefix", default="")
    parser.add_argument("--plot", default=str(script_dir / "nonlinear_sine_validation.png"))
    parser.add_argument("--tol-abs", type=float, default=1.0e-8)
    parser.add_argument("--tol-norm", type=float, default=1.0e-12)
    parser.add_argument("--tol-unexpected", type=float, default=1.0e-12)
    args = parser.parse_args()

    input_path = Path(args.input).resolve()
    cfg = parse_input_script(input_path)
    snapshots = resolve_snapshot_files(cfg, input_path, args.prefix)
    if not snapshots:
        print("No spectral snapshots found.", file=sys.stderr)
        return 2

    return 0 if analyze(cfg, snapshots, args) else 1


if __name__ == "__main__":
    raise SystemExit(main())

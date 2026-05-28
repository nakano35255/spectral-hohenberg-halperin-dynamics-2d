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
class GaussianInitialCondition:
    base: float = 0.0
    amplitude: float = 0.0
    x0: float = 0.0
    y0: float = 0.0
    sigma_x: float = 1.0
    sigma_y: float = 1.0


@dataclass
class Config:
    nx: int = 0
    ny: int = 0
    lx: float = 2.0 * math.pi
    ly: float = 2.0 * math.pi
    dt: float = 0.0
    integrator: str = ""
    order_parameters: int = 0
    a: Dict[int, float] = field(default_factory=dict)
    b: Dict[int, float] = field(default_factory=dict)
    u: Dict[int, float] = field(default_factory=dict)
    mobility: Dict[int, float] = field(default_factory=dict)
    gaussian: Dict[int, GaussianInitialCondition] = field(default_factory=dict)
    spectral_prefixes: List[str] = field(default_factory=list)


@dataclass
class Snapshot:
    path: Path
    step: int
    time: float
    nkx: int
    nky: int
    order_parameters: int
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

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = strip_comment(raw_line)
            if not line:
                continue

            tokens = line.split()
            head = tokens[0]

            if head == "grid":
                cfg.nx = int(tokens[1])
                cfg.ny = int(tokens[2])
                continue

            if head == "length":
                cfg.lx = float(tokens[1])
                cfg.ly = float(tokens[2])
                continue

            if head == "timestep":
                cfg.dt = float(tokens[1])
                continue

            if head in ("order_parameters", "order_parameter"):
                cfg.order_parameters = int(tokens[1])
                continue

            if head == "time_evolution":
                cfg.integrator = tokens[1]
                continue

            if head == "model" and len(tokens) >= 3:
                category = tokens[1]
                model_type = tokens[2]
                args = key_values(tokens, 3)

                if category == "free_energy" and model_type in ("phi4", "coeff"):
                    for key, value in args.items():
                        for name, target in (("a", cfg.a), ("b", cfg.b), ("u", cfg.u)):
                            index = parse_indexed_key(key, name)
                            if index >= 0:
                                target[index] = float(value)

                if category == "transport" and model_type in ("constant", "coeff"):
                    for key, value in args.items():
                        i, j = parse_mobility_key(key)
                        if i >= 0 and i == j:
                            cfg.mobility[i] = float(value)
                continue

            if head == "set" and len(tokens) >= 5 and tokens[1] in ("order_parameter", "psi"):
                target = tokens[2]
                ic_type = tokens[3]
                if ic_type != "gaussian":
                    continue
                args = key_values(tokens, 4)
                if "sigma" in args:
                    sigma_x = sigma_y = float(args["sigma"])
                else:
                    sigma_x = float(args["sigma_x"])
                    sigma_y = float(args["sigma_y"])
                gaussian = GaussianInitialCondition(
                    base=float(args["base"]),
                    amplitude=float(args["amplitude"]),
                    x0=float(args["x0"]),
                    y0=float(args["y0"]),
                    sigma_x=sigma_x,
                    sigma_y=sigma_y,
                )
                if target == "all":
                    for component in range(cfg.order_parameters):
                        cfg.gaussian[component] = gaussian
                else:
                    cfg.gaussian[int(target)] = gaussian
                continue

            if head == "measure" and len(tokens) >= 5:
                measure_type = tokens[2]
                state = tokens[3]
                if measure_type != "snapshot" or state != "on":
                    continue
                args = key_values(tokens, 4)
                space = args.get("space", "physical")
                if space in ("spectral", "fourier", "spec") and "file" in args:
                    cfg.spectral_prefixes.append(args["file"])

    if cfg.order_parameters != 1:
        raise ValueError("this example analysis expects order_parameters 1")
    if cfg.u.get(0, 0.0) != 0.0:
        raise ValueError("this exact comparison requires phi4 u[0] 0.0")
    if 0 not in cfg.gaussian:
        raise ValueError("this analysis expects a Gaussian initial condition for order_parameter 0")

    return cfg


def signed_ky(ky: np.ndarray, ny: int) -> np.ndarray:
    return np.where(ky <= ny // 2, ky, ky - ny)


def amplification_factor(lam: np.ndarray, cfg: Config, step: int) -> np.ndarray:
    z = lam * cfg.dt
    if cfg.integrator.startswith("srk3/"):
        one_step = 1.0 + z + 0.5 * z * z + (z * z * z) / 6.0
    elif cfg.integrator.startswith("euler/"):
        one_step = 1.0 + z
    else:
        raise ValueError("unsupported integrator: " + cfg.integrator)
    return np.power(one_step, step)


def parse_snapshot(path: Path) -> Snapshot:
    step = -1
    time = math.nan
    nkx = -1
    nky = -1
    order_parameters = -1
    columns: Dict[str, int] = {}

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("# step"):
                parts = line.split()
                step = int(parts[2])
                time = float(parts[4])
                continue
            if line.startswith("# nkx"):
                parts = line.split()
                nkx = int(parts[2])
                nky = int(parts[4])
                order_parameters = int(parts[6])
                continue
            if line.startswith("# kx_index"):
                names = line[2:].split()
                columns = {name: i for i, name in enumerate(names)}
                continue

    if step < 0 or nkx <= 0 or nky <= 0 or not columns:
        raise ValueError("invalid spectral snapshot: " + str(path))

    data = np.loadtxt(path, comments="#", ndmin=2)
    return Snapshot(path, step, time, nkx, nky, order_parameters, columns, data)


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


def spectral_arrays(snapshot: Snapshot, cfg: Config) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    kx_index = snapshot.data[:, snapshot.columns["kx_index"]].astype(int)
    ky_index = snapshot.data[:, snapshot.columns["ky_index"]].astype(int)
    values = (
        snapshot.data[:, snapshot.columns["psi_com0_real"]]
        + 1j * snapshot.data[:, snapshot.columns["psi_com0_imag"]]
    )

    kx = 2.0 * math.pi * kx_index / cfg.lx
    ky = 2.0 * math.pi * signed_ky(ky_index, cfg.ny) / cfg.ly
    return kx_index, ky_index, values, kx, ky


def gaussian_initial_coefficients(snapshot: Snapshot, cfg: Config) -> Tuple[np.ndarray, np.ndarray]:
    gaussian = cfg.gaussian[0]
    kx_index, ky_index, _, kx, ky = spectral_arrays(snapshot, cfg)
    k2 = kx * kx + ky * ky

    grid_size = float(cfg.nx) * float(cfg.ny)
    prefactor = gaussian.amplitude * grid_size * 2.0 * math.pi * gaussian.sigma_x * gaussian.sigma_y / (cfg.lx * cfg.ly)
    envelope = np.exp(-0.5 * (gaussian.sigma_x**2 * kx**2 + gaussian.sigma_y**2 * ky**2))
    phase = -(kx * gaussian.x0 + ky * gaussian.y0)
    coeff0 = prefactor * envelope * np.exp(1j * phase)

    zero = (kx_index == 0) & (ky_index == 0)
    coeff0[zero] += gaussian.base * grid_size
    return coeff0, k2


def expected_coefficients(snapshot: Snapshot, cfg: Config) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    coeff0, k2 = gaussian_initial_coefficients(snapshot, cfg)
    lam = -cfg.mobility.get(0, 0.0) * k2 * (cfg.a.get(0, 0.0) + cfg.b.get(0, 0.0) * k2)
    return coeff0 * amplification_factor(lam, cfg, snapshot.step), coeff0 * np.exp(lam * snapshot.time), coeff0


def mode_index(snapshot: Snapshot, mode: Mode) -> int:
    kx = snapshot.data[:, snapshot.columns["kx_index"]].astype(int)
    ky = snapshot.data[:, snapshot.columns["ky_index"]].astype(int)
    matches = np.flatnonzero((kx == mode[0]) & (ky == mode[1]))
    if matches.size != 1:
        raise ValueError(f"mode {mode} was not found exactly once in {snapshot.path}")
    return int(matches[0])


def selected_modes(snapshot: Snapshot) -> List[Mode]:
    modes = [(0, 0), (1, 0), (4, 0), (8, 0), (12, 0), (16, 16)]
    return [mode for mode in modes if mode[0] < snapshot.nkx and mode[1] < snapshot.nky]


def write_validation_plot(
    snapshots: List[Snapshot],
    cfg: Config,
    output: Path,
    min_initial_amplitude: float,
    min_relative_initial_amplitude: float,
) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    first = snapshots[0]
    final = snapshots[-1]
    modes = selected_modes(final)

    fig = plt.figure(figsize=(12.0, 8.0))
    gs = fig.add_gridspec(2, 2, width_ratios=(1.2, 1.0), height_ratios=(1.0, 1.0))
    ax_decay = fig.add_subplot(gs[:, 0])
    ax_slice = fig.add_subplot(gs[0, 1])
    ax_error = fig.add_subplot(gs[1, 1])

    times = np.asarray([snapshot.time for snapshot in snapshots])

    for mode in modes:
        actual = []
        exact = []
        scale = None
        for snapshot in snapshots:
            _, _, values, _, _ = spectral_arrays(snapshot, cfg)
            _, continuous, coeff0 = expected_coefficients(snapshot, cfg)
            index = mode_index(snapshot, mode)
            if scale is None:
                scale = max(abs(coeff0[index]), 1.0e-300)
            actual.append(abs(values[index]) / scale)
            exact.append(abs(continuous[index]) / scale)

        ax_decay.semilogy(times, actual, ".", markersize=3, label=f"{mode} snapshot")
        ax_decay.semilogy(times, exact, "-", linewidth=1.0, label=f"{mode} exact")

    ax_decay.set_xlabel("time")
    ax_decay.set_ylabel(r"$|\hat{\psi}(t)| / |\hat{\psi}(0)|$")
    ax_decay.grid(True, which="both", linewidth=0.4, alpha=0.5)
    ax_decay.legend(fontsize=7, ncol=2)

    _, _, final_values, _, _ = spectral_arrays(final, cfg)
    final_expected_srk3, final_expected_continuous, final_coeff0 = expected_coefficients(final, cfg)
    spectral_shape = (final.nky, final.nkx)
    sim_hat = final_values.reshape(spectral_shape)
    exact_hat = final_expected_continuous.reshape(spectral_shape)

    sim_phys = np.fft.irfft2(sim_hat, s=(cfg.ny, cfg.nx)).real
    exact_phys = np.fft.irfft2(exact_hat, s=(cfg.ny, cfg.nx)).real

    gaussian = cfg.gaussian[0]
    y_index = int(round((gaussian.y0 / cfg.ly) * cfg.ny)) % cfg.ny
    x = np.arange(cfg.nx) * cfg.lx / cfg.nx
    ax_slice.plot(x, sim_phys[y_index, :], ".", markersize=3, label="snapshot")
    ax_slice.plot(x, exact_phys[y_index, :], "-", linewidth=1.0, label="exact")
    ax_slice.set_xlabel("x")
    ax_slice.set_ylabel(r"$\psi(x, y_0)$")
    ax_slice.grid(True, linewidth=0.4, alpha=0.5)
    ax_slice.legend(fontsize=8)

    amplitude_floor = max(min_initial_amplitude, min_relative_initial_amplitude * float(np.max(np.abs(final_coeff0))))
    denom = np.maximum(np.abs(final_coeff0), amplitude_floor)
    normalized_error = np.abs(final_values - final_expected_srk3) / denom
    error_grid = np.log10(normalized_error.reshape(spectral_shape) + 1.0e-30)
    error_grid = np.fft.fftshift(error_grid, axes=0)

    image = ax_error.imshow(
        error_grid,
        origin="lower",
        aspect="auto",
        extent=[0, final.nkx - 1, -final.nky // 2, final.nky // 2 - 1],
        vmin=-16,
        vmax=-10,
    )
    ax_error.set_xlabel("kx_index")
    ax_error.set_ylabel("signed ky_index")
    ax_error.set_title(r"$\log_{10}$ normalized SRK3 error")
    fig.colorbar(image, ax=ax_error, shrink=0.85)

    fig.suptitle("Example 02: linear phi4 Gaussian validation")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def analyze(cfg: Config, snapshot_paths: List[Path], args: argparse.Namespace) -> bool:
    snapshots = [parse_snapshot(path) for path in snapshot_paths]
    max_abs_error = 0.0
    max_norm_error = 0.0
    max_continuous_error = 0.0
    worst_file = ""
    worst_step = -1

    for snapshot in snapshots:
        _, _, values, _, _ = spectral_arrays(snapshot, cfg)
        expected_srk3, expected_continuous, coeff0 = expected_coefficients(snapshot, cfg)
        abs_error = np.abs(values - expected_srk3)
        continuous_error = np.abs(values - expected_continuous)
        amplitude = np.abs(coeff0)
        amplitude_floor = max(args.min_initial_amplitude, args.min_relative_initial_amplitude * float(np.max(amplitude)))
        significant = amplitude >= amplitude_floor
        if np.any(significant):
            norm_error = abs_error[significant] / amplitude[significant]
            local_norm = float(np.max(norm_error))
        else:
            local_norm = 0.0

        local_abs = float(np.max(abs_error))
        local_cont = float(np.max(continuous_error))

        if local_abs > max_abs_error:
            max_abs_error = local_abs
            worst_file = str(snapshot.path)
            worst_step = snapshot.step
        max_norm_error = max(max_norm_error, local_norm)
        max_continuous_error = max(max_continuous_error, local_cont)

    print("Input")
    print(f"  grid              : {cfg.nx} x {cfg.ny}")
    print(f"  length            : {cfg.lx:g} x {cfg.ly:g}")
    print(f"  dt                : {cfg.dt:g}")
    print(f"  time_evolution    : {cfg.integrator}")
    print(f"  phi4              : a={cfg.a.get(0, 0.0):g} b={cfg.b.get(0, 0.0):g} u={cfg.u.get(0, 0.0):g}")
    print(f"  mobility          : M={cfg.mobility.get(0, 0.0):g}")
    print(f"  snapshots         : {len(snapshots)}")
    print()
    print("Worst errors")
    print(f"  SRK3 abs error       : {max_abs_error:.16e} (step={worst_step}, file={worst_file})")
    print(f"  SRK3 normalized error: {max_norm_error:.16e}")
    print(f"  exact abs error      : {max_continuous_error:.16e}")
    print()

    if args.plot:
        output = Path(args.plot).resolve()
        write_validation_plot(snapshots, cfg, output, args.min_initial_amplitude, args.min_relative_initial_amplitude)
        print(f"plot: {output}")

    passed = max_abs_error <= args.tol_abs or max_norm_error <= args.tol_norm
    print("PASS" if passed else "FAIL")
    return passed


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Validate linear phi4 dynamics from a Gaussian initial condition.")
    parser.add_argument("--input", default=str(script_dir / "input.script"))
    parser.add_argument("--prefix", default="")
    parser.add_argument("--plot", default=str(script_dir / "phi4_gaussian_validation.png"))
    parser.add_argument("--tol-abs", type=float, default=1.0e-7)
    parser.add_argument("--tol-norm", type=float, default=1.0e-10)
    parser.add_argument("--min-initial-amplitude", type=float, default=1.0e-14)
    parser.add_argument("--min-relative-initial-amplitude", type=float, default=1.0e-8)
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

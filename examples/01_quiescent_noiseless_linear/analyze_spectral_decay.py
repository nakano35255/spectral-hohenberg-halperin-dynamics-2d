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
    order_parameters: int = 0
    integrator: str = ""
    a: Dict[int, float] = field(default_factory=dict)
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
    if not match:
        return -1
    return int(match.group(1))


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

                if category == "free_energy" and model_type in ("quadratic", "coeff"):
                    for key, value in args.items():
                        index = parse_indexed_key(key, "a")
                        if index >= 0:
                            cfg.a[index] = float(value)

                if category == "transport" and model_type in ("constant", "coeff"):
                    for key, value in args.items():
                        i, j = parse_mobility_key(key)
                        if i >= 0 and i == j:
                            cfg.mobility[i] = float(value)
                continue

            if head == "set" and len(tokens) >= 5 and tokens[1] in ("order_parameter", "psi"):
                target = tokens[2]
                ic_type = tokens[3]
                if ic_type != "sine":
                    continue
                args = key_values(tokens, 4)
                sine = SineInitialCondition(
                    base=float(args.get("base", "0.0")),
                    amplitude=float(args["amplitude"]),
                    nkx=int(args["nkx"]),
                    nky=int(args.get("nky", "0")),
                )
                if target == "all":
                    for component in range(cfg.order_parameters):
                        cfg.sine[component] = sine
                else:
                    cfg.sine[int(target)] = sine
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

    if cfg.nx <= 0 or cfg.ny <= 0:
        raise ValueError("input script must define a positive grid")
    if cfg.dt <= 0.0:
        raise ValueError("input script must define a positive timestep")
    if cfg.order_parameters <= 0:
        raise ValueError("this analysis expects at least one order parameter")

    return cfg


def wrap_ky(nky: int, ny: int) -> int:
    return nky if nky >= 0 else ny + nky


def signed_ky(ky: int, ny: int) -> int:
    return ky if ky <= ny // 2 else ky - ny


def wave_number_squared(mode: Mode, cfg: Config) -> float:
    kx_index, ky_index = mode
    kx = 2.0 * math.pi * kx_index / cfg.lx
    ky = 2.0 * math.pi * signed_ky(ky_index, cfg.ny) / cfg.ly
    return kx * kx + ky * ky


def initial_coefficients(component: int, cfg: Config) -> Dict[Mode, complex]:
    if component not in cfg.sine:
        return {}

    sine = cfg.sine[component]
    grid_size = float(cfg.nx) * float(cfg.ny)
    coeffs: Dict[Mode, complex] = {}

    if sine.base != 0.0:
        coeffs[(0, 0)] = complex(sine.base * grid_size, 0.0)

    ky = wrap_ky(sine.nky, cfg.ny)
    conjugate_ky = (cfg.ny - ky) % cfg.ny
    positive = complex(0.0, -0.5 * sine.amplitude * grid_size)
    coeffs[(sine.nkx, ky)] = positive
    if sine.nkx == 0 and conjugate_ky != ky:
        coeffs[(sine.nkx, conjugate_ky)] = positive.conjugate()

    return coeffs


def amplification_factor(lam: float, cfg: Config, step: int) -> float:
    z = lam * cfg.dt
    if cfg.integrator.startswith("srk3/"):
        one_step = 1.0 + z + 0.5 * z * z + (z * z * z) / 6.0
    elif cfg.integrator.startswith("euler/"):
        one_step = 1.0 + z
    else:
        raise ValueError("unsupported integrator for this analysis: " + cfg.integrator)
    return one_step ** step


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
        candidates = []
        raw = Path(prefix)
        if raw.is_absolute():
            candidates.append(raw)
        else:
            candidates.append(Path.cwd() / raw)
            candidates.append(repo_root / raw)
            candidates.append(input_dir / raw)

        for candidate in candidates:
            matches = [Path(p) for p in glob.glob(str(candidate) + "_step*.snapshot")]
            if matches:
                files.extend(matches)
                break

    files = sorted(set(files), key=lambda p: str(p))
    files.sort(key=lambda p: int(re.search(r"_step(\d+)\.snapshot$", p.name).group(1)))
    return files


def psi_values(snapshot: Snapshot, component: int) -> np.ndarray:
    real_key = f"psi_com{component}_real"
    imag_key = f"psi_com{component}_imag"
    return snapshot.data[:, snapshot.columns[real_key]] + 1j * snapshot.data[:, snapshot.columns[imag_key]]


def write_decay_plot(
    series: Dict[Tuple[int, Mode], Dict[str, List[float]]],
    output: Path,
    show_srk3: bool,
) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.0, 5.0))

    for (component, mode), values in sorted(series.items()):
        time = np.asarray(values["time"])
        actual = np.asarray(values["actual"])
        srk3 = np.asarray(values["srk3"])
        continuous = np.asarray(values["continuous"])
        label = f"psi{component} {mode}"

        scale = max(float(values["scale"][0]), 1.0)

        ax.semilogy(time, actual / scale, ".", markersize=3, label=f"{label} snapshot")
        ax.semilogy(time, continuous / scale, "-", linewidth=1.0, label=f"{label} exact")
        if show_srk3:
            ax.semilogy(time, srk3 / scale, "--", linewidth=0.8, label=f"{label} SRK3")

    ax.set_xlabel("time")
    ax.set_ylabel(r"$|\hat{\psi}(t)| / |\hat{\psi}(0)|$")
    ax.grid(True, which="both", linewidth=0.4, alpha=0.5)
    ax.legend(fontsize=6, ncol=2)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def analyze(cfg: Config, snapshots: List[Path], args: argparse.Namespace) -> bool:
    expected_by_component = {
        component: initial_coefficients(component, cfg)
        for component in range(cfg.order_parameters)
    }

    print("Input")
    print(f"  grid              : {cfg.nx} x {cfg.ny}")
    print(f"  length            : {cfg.lx:g} x {cfg.ly:g}")
    print(f"  dt                : {cfg.dt:g}")
    print(f"  time_evolution    : {cfg.integrator}")
    print(f"  snapshots         : {len(snapshots)}")
    print()

    print("Expected spectral modes")
    for component, coeffs in expected_by_component.items():
        if not coeffs:
            print(f"  psi_com{component}: none")
            continue
        pieces = []
        for mode, coeff0 in sorted(coeffs.items()):
            lam = -cfg.mobility.get(component, 0.0) * cfg.a.get(component, 0.0) * wave_number_squared(mode, cfg)
            pieces.append(f"{mode} lambda={lam:.16e} |c0|={abs(coeff0):.16e}")
        print(f"  psi_com{component}: " + "; ".join(pieces))
    print()

    worst_unexpected = (0.0, "", -1, -1, (-1, -1))
    worst_srk3_abs = (0.0, "", -1, -1, (-1, -1))
    worst_srk3_norm = (0.0, "", -1, -1, (-1, -1))
    worst_cont_abs = (0.0, "", -1, -1, (-1, -1))
    series: Dict[Tuple[int, Mode], Dict[str, List[float]]] = {}

    selected = snapshots
    if args.max_files > 0:
        selected = selected[: args.max_files]

    target_modes_by_component = {
        component: set(coeffs.keys())
        for component, coeffs in expected_by_component.items()
    }

    for path in selected:
        snapshot = parse_snapshot(path)
        if snapshot.order_parameters != cfg.order_parameters:
            raise ValueError(f"{path}: order_parameter count mismatch")

        kx = snapshot.data[:, snapshot.columns["kx_index"]].astype(int)
        ky = snapshot.data[:, snapshot.columns["ky_index"]].astype(int)

        for component in range(cfg.order_parameters):
            values = psi_values(snapshot, component)
            target_mask = np.zeros(values.shape, dtype=bool)

            for mode in target_modes_by_component[component]:
                target_mask |= (kx == mode[0]) & (ky == mode[1])

            unexpected_indices = np.flatnonzero(~target_mask)
            if unexpected_indices.size > 0:
                unexpected_magnitude = np.abs(values[unexpected_indices])
                local_argmax = int(np.argmax(unexpected_magnitude))
                global_index = int(unexpected_indices[local_argmax])
                magnitude = float(unexpected_magnitude[local_argmax])
                if magnitude > worst_unexpected[0]:
                    mode = (int(kx[global_index]), int(ky[global_index]))
                    worst_unexpected = (magnitude, str(path), snapshot.step, component, mode)

            for mode in target_modes_by_component[component]:
                matches = np.flatnonzero((kx == mode[0]) & (ky == mode[1]))
                if matches.size != 1:
                    raise ValueError(f"{path}: expected mode {mode} for psi_com{component} was not found exactly once")
                value = complex(values[int(matches[0])])
                coeff0 = expected_by_component[component][mode]
                lam = -cfg.mobility.get(component, 0.0) * cfg.a.get(component, 0.0) * wave_number_squared(mode, cfg)
                expected_srk3 = coeff0 * amplification_factor(lam, cfg, snapshot.step)
                expected_cont = coeff0 * math.exp(lam * snapshot.time)

                srk3_abs = abs(value - expected_srk3)
                srk3_norm = srk3_abs / max(abs(coeff0), 1.0)
                cont_abs = abs(value - expected_cont)

                if srk3_abs > worst_srk3_abs[0]:
                    worst_srk3_abs = (srk3_abs, str(path), snapshot.step, component, mode)
                if srk3_norm > worst_srk3_norm[0]:
                    worst_srk3_norm = (srk3_norm, str(path), snapshot.step, component, mode)
                if cont_abs > worst_cont_abs[0]:
                    worst_cont_abs = (cont_abs, str(path), snapshot.step, component, mode)

                key = (component, mode)
                series.setdefault(key, {"time": [], "actual": [], "srk3": [], "continuous": [], "scale": []})
                series[key]["time"].append(snapshot.time)
                series[key]["actual"].append(abs(value))
                series[key]["srk3"].append(abs(expected_srk3))
                series[key]["continuous"].append(abs(expected_cont))
                series[key]["scale"].append(abs(coeff0))

    print("Worst errors")
    print(
        "  unexpected mode amplitude : "
        f"{worst_unexpected[0]:.16e} "
        f"(step={worst_unexpected[2]}, psi_com{worst_unexpected[3]}, mode={worst_unexpected[4]}, file={worst_unexpected[1]})"
    )
    print(
        "  SRK3 abs error            : "
        f"{worst_srk3_abs[0]:.16e} "
        f"(step={worst_srk3_abs[2]}, psi_com{worst_srk3_abs[3]}, mode={worst_srk3_abs[4]}, file={worst_srk3_abs[1]})"
    )
    print(
        "  SRK3 normalized error     : "
        f"{worst_srk3_norm[0]:.16e} "
        f"(step={worst_srk3_norm[2]}, psi_com{worst_srk3_norm[3]}, mode={worst_srk3_norm[4]}, file={worst_srk3_norm[1]})"
    )
    print(
        "  exact abs error           : "
        f"{worst_cont_abs[0]:.16e} "
        f"(step={worst_cont_abs[2]}, psi_com{worst_cont_abs[3]}, mode={worst_cont_abs[4]}, file={worst_cont_abs[1]})"
    )
    print()

    passed = True
    if worst_unexpected[0] > args.tol_unexpected:
        print(f"FAIL: unexpected mode amplitude exceeds {args.tol_unexpected:.3e}")
        passed = False
    if worst_srk3_abs[0] > args.tol_srk3_abs and worst_srk3_norm[0] > args.tol_srk3_norm:
        print(
            "FAIL: target modes deviate from the SRK3 amplification factor "
            f"(abs>{args.tol_srk3_abs:.3e} and normalized>{args.tol_srk3_norm:.3e})"
        )
        passed = False

    if passed:
        print("PASS")

    if args.plot:
        plot_path = Path(args.plot).resolve()
        write_decay_plot(series, plot_path, args.plot_srk3)
        print(f"plot: {plot_path}")

    return passed


def main() -> int:
    script_dir = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(
        description="Check spectral snapshots for noiseless linear quiescent dynamics."
    )
    parser.add_argument(
        "--input",
        default=str(script_dir / "input.script"),
        help="input script used to generate the snapshots",
    )
    parser.add_argument(
        "--prefix",
        default="",
        help="spectral snapshot prefix; defaults to the spectral measure prefix in the input script",
    )
    parser.add_argument(
        "--tol-unexpected",
        type=float,
        default=1.0e-10,
        help="absolute tolerance for non-target spectral modes",
    )
    parser.add_argument(
        "--tol-srk3-abs",
        type=float,
        default=1.0e-9,
        help="absolute tolerance for target-mode agreement with the SRK3 amplification factor",
    )
    parser.add_argument(
        "--tol-srk3-norm",
        type=float,
        default=1.0e-12,
        help="initial-amplitude-normalized tolerance for target-mode agreement with the SRK3 amplification factor",
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=0,
        help="analyze only the first N snapshots; 0 means all snapshots",
    )
    parser.add_argument(
        "--plot",
        default="",
        help="write a semilog decay plot to this file",
    )
    parser.add_argument(
        "--plot-srk3",
        action="store_true",
        help="also draw the SRK3 discrete amplification curve in the plot",
    )
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

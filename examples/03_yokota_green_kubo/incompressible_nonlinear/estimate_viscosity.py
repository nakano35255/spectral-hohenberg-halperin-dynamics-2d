#!/usr/bin/env python3
"""Analyze Yokota Green-Kubo viscosity estimates."""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import re
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


CASE_DIR = Path(__file__).resolve().parent
ROOT = CASE_DIR.parents[2]
DEFAULT_CASE = "eta0_0p1_grid256_L256_dt0p01_T25000_diag_n1152"
RAW_DATA_DIR = CASE_DIR / "raw_data"
PROCESSED_DATA_DIR = CASE_DIR / "processed_data"

DEFAULT_YKGK = RAW_DATA_DIR / DEFAULT_CASE / "results" / "yokota_green_kubo_000.dat"
DEFAULT_TIME_SERIES = RAW_DATA_DIR / DEFAULT_CASE / "results" / "time_series_000.dat"
DEFAULT_INPUT_SCRIPT = RAW_DATA_DIR / DEFAULT_CASE / "runs" / "input_000.script"
DEFAULT_OUTPUT_DIR = PROCESSED_DATA_DIR / DEFAULT_CASE


@dataclass(frozen=True)
class YokotaRow:
    source: Path
    nsample: int
    tau: float
    mode_index: int
    kx: float
    ky: float
    s_mean: float


@dataclass(frozen=True)
class YokotaPoint:
    tau: float
    mode_index: int
    kx: float
    ky: float
    s_mean: float
    s_sem: float
    s_std: float
    count: int


@dataclass(frozen=True)
class ModeSummary:
    mode_index: int
    kx: float
    ky: float
    k_abs: float
    tau_min: float
    tau_max: float
    s_last: float
    s_tail_mean: float
    s_tail_sem: float
    s_tail_std: float
    tail_points: int
    sample_count_min: int


def repo_path(value: str | Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return ROOT / path


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def tokenize(line: str) -> list[str]:
    return line.split("#", 1)[0].split()


def parse_key_values(tokens: list[str]) -> dict[str, str]:
    values: dict[str, str] = {}
    for i in range(0, len(tokens), 2):
        if i + 1 < len(tokens):
            values[tokens[i]] = tokens[i + 1]
    return values


def parse_input_config(input_script: Path) -> dict[str, object]:
    config: dict[str, object] = {
        "eta": None,
        "rho0": 1.0,
        "kBT": 1.0,
        "length": None,
        "grid": None,
    }
    if not input_script.exists():
        return config

    with input_script.open() as handle:
        for raw_line in handle:
            tokens = tokenize(raw_line)
            if not tokens:
                continue
            if tokens[0] == "grid" and len(tokens) >= 3:
                config["grid"] = (int(tokens[1]), int(tokens[2]))
            elif tokens[0] == "length" and len(tokens) >= 3:
                config["length"] = (float(tokens[1]), float(tokens[2]))
            elif len(tokens) >= 4 and tokens[:3] == ["model", "transport", "constant"]:
                values = parse_key_values(tokens[3:])
                if "eta" in values:
                    config["eta"] = float(values["eta"])
            elif len(tokens) >= 6 and tokens[0] == "fix" and tokens[2:5] == ["momentum", "noise", "on"]:
                values = parse_key_values(tokens[5:])
                if "kBT" in values:
                    config["kBT"] = float(values["kBT"])
            elif tokens[:4] == ["set", "density", "uniform", "value"] and len(tokens) >= 5:
                config["rho0"] = float(tokens[4])

    return config


def stderr(values: Iterable[float]) -> float:
    array = np.asarray(list(values), dtype=float)
    if array.size <= 1:
        return 0.0
    return float(array.std(ddof=1) / math.sqrt(array.size))


def fns_eta(q_values, eta0, rho0, kBT, a_uv=1.0):
    q_values = np.asarray(q_values, dtype=float)
    return np.sqrt(eta0 * eta0 + (rho0 * kBT) / (8.0 * np.pi) * np.log(2.0 * np.pi / (a_uv * q_values)))


def discrete_square_eta(q_values, eta0, rho0, kBT, length, a_uv):
    lx, ly = length
    uv_nx = int(np.floor(lx / a_uv))
    uv_ny = int(np.floor(ly / a_uv))
    if uv_nx <= 0 or uv_ny <= 0:
        raise RuntimeError("a_uv is too large for the discrete square cutoff")

    nx = np.arange(-uv_nx, uv_nx + 1, dtype=float)
    ny = np.arange(-uv_ny, uv_ny + 1, dtype=float)
    kx = 2.0 * np.pi * nx / lx
    ky = 2.0 * np.pi * ny / ly
    kx_grid, ky_grid = np.meshgrid(kx, ky, indexing="ij")
    k2 = kx_grid * kx_grid + ky_grid * ky_grid
    nonzero = k2 > 0.0
    integrand = np.zeros_like(k2)
    integrand[nonzero] = (kx_grid[nonzero] ** 2) * (ky_grid[nonzero] ** 2) / (k2[nonzero] ** 3)

    values = []
    for q in np.asarray(q_values, dtype=float):
        mask = nonzero & (np.sqrt(k2) >= q)
        correction = 2.0 * rho0 * kBT * float(np.sum(integrand[mask])) / (lx * ly)
        values.append(np.sqrt(eta0 * eta0 + correction))
    return np.asarray(values)


def expand_input_paths(inputs: list[str] | None, patterns: list[str] | None) -> list[Path]:
    paths: list[Path] = []
    if inputs:
        paths.extend(repo_path(value) for value in inputs)
    if patterns:
        for pattern in patterns:
            paths.extend(Path(value) for value in glob.glob(str(repo_path(pattern))))
    if not paths:
        paths.append(DEFAULT_YKGK)

    unique: dict[str, Path] = {}
    for path in paths:
        unique[str(path)] = path
    return sorted(unique.values())


def read_yokota_table(path: Path) -> list[YokotaRow]:
    if not path.exists():
        raise FileNotFoundError(f"missing Yokota Green-Kubo output: {path}")

    rows: list[YokotaRow] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split()
            if len(fields) != 6:
                raise ValueError(f"unexpected row in {path}: {line}")

            rows.append(
                YokotaRow(
                    source=path,
                    nsample=int(fields[0]),
                    tau=float(fields[1]),
                    mode_index=int(fields[2]),
                    kx=float(fields[3]),
                    ky=float(fields[4]),
                    s_mean=float(fields[5]),
                )
            )

    if not rows:
        raise ValueError(f"no data rows found in {path}")
    return rows


def aggregate_rows(rows: Iterable[YokotaRow]) -> list[YokotaPoint]:
    buckets: dict[tuple[int, float], list[float]] = {}
    metadata: dict[tuple[int, float], tuple[float, float]] = {}

    for row in rows:
        key = (row.mode_index, row.tau)
        if key in metadata:
            kx, ky = metadata[key]
            if not (math.isclose(row.kx, kx, rel_tol=1.0e-11, abs_tol=1.0e-14)
                    and math.isclose(row.ky, ky, rel_tol=1.0e-11, abs_tol=1.0e-14)):
                raise RuntimeError(f"inconsistent wave vector for mode={row.mode_index} tau={row.tau:g}")
        else:
            metadata[key] = (row.kx, row.ky)
        buckets.setdefault(key, []).append(row.s_mean)

    points: list[YokotaPoint] = []
    for mode_index, tau in sorted(buckets, key=lambda item: (item[0], item[1])):
        values = np.asarray(buckets[(mode_index, tau)], dtype=float)
        kx, ky = metadata[(mode_index, tau)]
        s_std = float(values.std(ddof=1)) if values.size > 1 else 0.0
        s_sem = s_std / math.sqrt(values.size) if values.size > 1 else 0.0
        points.append(
            YokotaPoint(
                tau=tau,
                mode_index=mode_index,
                kx=kx,
                ky=ky,
                s_mean=float(values.mean()),
                s_sem=s_sem,
                s_std=s_std,
                count=int(values.size),
            )
        )
    return points


def group_by_mode(points: Iterable[YokotaPoint]) -> dict[int, list[YokotaPoint]]:
    grouped: dict[int, list[YokotaPoint]] = {}
    for point in points:
        grouped.setdefault(point.mode_index, []).append(point)

    for mode_points in grouped.values():
        mode_points.sort(key=lambda point: point.tau)

    return dict(sorted(grouped.items()))


def choose_modes(mode_indices: list[int], requested: list[int] | None, max_modes: int) -> list[int]:
    if requested:
        missing = sorted(set(requested) - set(mode_indices))
        if missing:
            raise ValueError(f"requested mode_index values are absent: {missing}")
        return sorted(dict.fromkeys(requested))

    if max_modes <= 0 or len(mode_indices) <= max_modes:
        return mode_indices

    if max_modes == 1:
        return [mode_indices[0]]

    chosen: set[int] = set()
    n = len(mode_indices)
    for i in range(max_modes):
        index = round(i * (n - 1) / (max_modes - 1))
        chosen.add(mode_indices[index])

    return [mode for mode in mode_indices if mode in chosen]


def summarize_modes(grouped: dict[int, list[YokotaPoint]], tail_fraction: float) -> list[ModeSummary]:
    if not (0.0 < tail_fraction <= 1.0):
        raise ValueError("tail_fraction must satisfy 0 < tail_fraction <= 1")

    summaries: list[ModeSummary] = []
    for mode_index, mode_points in grouped.items():
        if not mode_points:
            continue

        tail_count = max(1, math.ceil(len(mode_points) * tail_fraction))
        tail_points = mode_points[-tail_count:]
        tail_values = [point.s_mean for point in tail_points]
        tail_sem = math.sqrt(sum(point.s_sem * point.s_sem for point in tail_points)) / len(tail_points)
        kx = mode_points[-1].kx
        ky = mode_points[-1].ky
        k_abs = math.hypot(kx, ky)
        tail_std = statistics.stdev(tail_values) if len(tail_values) > 1 else 0.0

        summaries.append(
            ModeSummary(
                mode_index=mode_index,
                kx=kx,
                ky=ky,
                k_abs=k_abs,
                tau_min=mode_points[0].tau,
                tau_max=mode_points[-1].tau,
                s_last=mode_points[-1].s_mean,
                s_tail_mean=statistics.fmean(tail_values),
                s_tail_sem=tail_sem,
                s_tail_std=tail_std,
                tail_points=len(tail_points),
                sample_count_min=min(point.count for point in tail_points),
            )
        )

    return summaries


def write_summary(
    path: Path,
    summaries: list[ModeSummary],
    eta: float | None,
    sources: list[Path],
    tail_fraction: float,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "mode_index",
        "kx",
        "ky",
        "k_abs",
        "tau_min",
        "tau_max",
        "S_last",
        "S_tail_mean",
        "S_tail_sem",
        "S_tail_std",
        "tail_points",
        "sample_count_min",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for summary in summaries:
            writer.writerow(
                {
                    "mode_index": summary.mode_index,
                    "kx": f"{summary.kx:.16e}",
                    "ky": f"{summary.ky:.16e}",
                    "k_abs": f"{summary.k_abs:.16e}",
                    "tau_min": f"{summary.tau_min:.16e}",
                    "tau_max": f"{summary.tau_max:.16e}",
                    "S_last": f"{summary.s_last:.16e}",
                    "S_tail_mean": f"{summary.s_tail_mean:.16e}",
                    "S_tail_sem": f"{summary.s_tail_sem:.16e}",
                    "S_tail_std": f"{summary.s_tail_std:.16e}",
                    "tail_points": summary.tail_points,
                    "sample_count_min": summary.sample_count_min,
                }
            )


def write_metadata(
    path: Path,
    eta: float | None,
    sources: list[Path],
    tail_fraction: float,
    input_script: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    entries = [
        ("case", "yokota_green_kubo"),
        ("source_input_script", display_path(input_script)),
        ("source_files", str(len(sources))),
        ("tail_fraction", f"{tail_fraction:.16e}"),
    ]
    if eta is not None:
        entries.append(("eta_reference", f"{eta:.16e}"))
    for index, source in enumerate(sources[:8]):
        entries.append((f"source_{index}", display_path(source)))
    if len(sources) > 8:
        entries.append(("source_extra", str(len(sources) - 8)))

    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["key", "value"])
        writer.writerows(entries)


def import_pyplot():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required for plotting") from exc
    return plt


def configure_plot_cache(output_dir: Path) -> None:
    cache_root = output_dir / ".plot_cache"
    mpl_config = cache_root / "matplotlib"
    xdg_cache = cache_root / "xdg"
    mpl_config.mkdir(parents=True, exist_ok=True)
    xdg_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config))
    os.environ.setdefault("XDG_CACHE_HOME", str(xdg_cache))


def plot_tau_curves(
    output_path: Path,
    grouped: dict[int, list[YokotaPoint]],
    selected_modes: list[int],
    eta: float | None,
) -> None:
    plt = import_pyplot()
    fig, ax = plt.subplots(figsize=(8.0, 5.0))

    for mode_index in selected_modes:
        points = grouped[mode_index]
        taus = np.asarray([point.tau for point in points], dtype=float)
        means = np.asarray([point.s_mean for point in points], dtype=float)
        sems = np.asarray([point.s_sem for point in points], dtype=float)
        label = f"m={mode_index}, k={points[-1].kx:.3g}"
        ax.plot(taus, means, marker=".", linewidth=1.2, label=label)
        if np.any(sems > 0.0):
            ax.fill_between(taus, means - sems, means + sems, alpha=0.14, linewidth=0)

    if eta is not None:
        ax.axhline(eta, color="black", linestyle="--", linewidth=1.0, label=f"eta0={eta:g}")

    ax.set_xlabel("tau")
    ax.set_ylabel("S(k,tau)")
    ax.set_title("Yokota Green-Kubo: stress integral")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_k_spectrum(
    output_path: Path,
    summaries: list[ModeSummary],
    eta: float | None,
    config: dict[str, object],
    a_uv: float,
    draw_theory: bool,
) -> None:
    plt = import_pyplot()
    q_values = np.asarray([summary.kx for summary in summaries], dtype=float)
    y_values = np.asarray([summary.s_tail_mean for summary in summaries], dtype=float)
    y_errors = np.asarray(
        [summary.s_tail_sem if summary.s_tail_sem > 0.0 else summary.s_tail_std for summary in summaries],
        dtype=float,
    )

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.errorbar(q_values, y_values, yerr=y_errors, marker="o", linestyle="-", capsize=3, label="YKGK diagonal")

    if eta is not None:
        ax.axhline(eta, color="#64748b", linestyle="--", linewidth=1.0, label=rf"$\eta_0={eta:g}$")

    length = config.get("length")
    rho0 = config.get("rho0")
    kBT = config.get("kBT")
    if draw_theory and eta is not None and length is not None and rho0 is not None and kBT is not None:
        q_line = np.geomspace(q_values.min(), q_values.max(), 600)
        eta_fns = fns_eta(q_line, eta, float(rho0), float(kBT), a_uv=a_uv)
        eta_square = discrete_square_eta(q_values, eta, float(rho0), float(kBT), length, a_uv)
        ax.plot(q_line, eta_fns, color="#111827", linestyle=":", linewidth=1.7, label=rf"FNS $a_{{uv}}={a_uv:g}$")
        ax.plot(q_values, eta_square, color="#dc2626", linestyle=":", linewidth=1.5, label="discrete square sum")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$k$")
    ax.set_ylabel(r"$\eta_{\rm eff}(k)$")
    ax.set_title("Yokota Green-Kubo: diagonal mode viscosity")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def read_time_series(path: Path) -> tuple[list[str], list[list[float]]]:
    if not path.exists():
        return [], []

    columns: list[str] = []
    rows: list[list[float]] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("# step"):
                columns = line[1:].split()
                continue
            if line.startswith("#"):
                continue
            rows.append([float(value) for value in line.split()])

    return columns, rows


def plot_time_series(output_path: Path, columns: list[str], rows: list[list[float]]) -> None:
    if not columns or not rows:
        return

    column_index = {name: index for index, name in enumerate(columns)}
    if "time" not in column_index:
        return

    plot_columns = [name for name in ("E_T", "E_K", "pi_xy") if name in column_index]
    if not plot_columns:
        return

    plt = import_pyplot()
    fig, axes = plt.subplots(len(plot_columns), 1, figsize=(8.0, 2.6 * len(plot_columns)), sharex=True)
    if len(plot_columns) == 1:
        axes = [axes]

    time = [row[column_index["time"]] for row in rows]
    for axis, name in zip(axes, plot_columns):
        axis.plot(time, [row[column_index[name]] for row in rows], linewidth=1.0)
        axis.set_ylabel(name)
        axis.grid(True, alpha=0.3)

    axes[-1].set_xlabel("time")
    fig.suptitle("Time series diagnostics")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", default=None, help="Yokota Green-Kubo output file. Repeatable.")
    parser.add_argument("--input-glob", action="append", default=None, help="Glob for many Yokota Green-Kubo files.")
    parser.add_argument("--input-script", default=str(DEFAULT_INPUT_SCRIPT), help="Input script used to read eta and box settings.")
    parser.add_argument("--time-series", default=str(DEFAULT_TIME_SERIES), help="Optional single time_series file for diagnostics.")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR), help="Directory for processed summary CSV output.")
    parser.add_argument("--figure-dir", default=None, help="Directory for PNG figures. Defaults to --output-dir.")
    parser.add_argument("--eta", type=float, default=None, help="Reference eta. Defaults to model transport eta from --input-script.")
    parser.add_argument("--tail-fraction", type=float, default=0.2, help="Fraction of largest tau points used for eta(q).")
    parser.add_argument("--max-modes", type=int, default=6, help="Maximum number of modes to plot in S(k,tau).")
    parser.add_argument("--modes", type=int, nargs="*", default=None, help="Specific mode_index values for S(k,tau).")
    parser.add_argument("--a-uv", type=float, default=1.0, help="Microscopic cutoff used for the theory curves.")
    parser.add_argument("--no-theory", action="store_true", help="Do not draw FNS/discrete-square theory curves.")
    parser.add_argument("--no-plots", action="store_true", help="Only write the processed CSV summary.")
    return parser


def main() -> None:
    args = build_argument_parser().parse_args()

    input_script = repo_path(args.input_script)
    output_dir = repo_path(args.output_dir)
    figure_dir = repo_path(args.figure_dir) if args.figure_dir is not None else output_dir
    time_series_path = repo_path(args.time_series)
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    configure_plot_cache(figure_dir)

    sources = expand_input_paths(args.input, args.input_glob)
    rows: list[YokotaRow] = []
    for source in sources:
        rows.extend(read_yokota_table(source))

    config = parse_input_config(input_script)
    eta = args.eta if args.eta is not None else config.get("eta")
    eta = float(eta) if eta is not None else None

    points = aggregate_rows(rows)
    grouped = group_by_mode(points)
    mode_indices = sorted(grouped)
    selected_modes = choose_modes(mode_indices, args.modes, args.max_modes)
    summaries = summarize_modes(grouped, args.tail_fraction)

    summary_path = output_dir / "yokota_green_kubo_summary.csv"
    metadata_path = output_dir / "metadata.csv"
    write_summary(summary_path, summaries, eta, sources, args.tail_fraction)
    write_metadata(metadata_path, eta, sources, args.tail_fraction, input_script)

    if not args.no_plots:
        tau_plot = figure_dir / "yokota_green_kubo_tau.png"
        spectrum_plot = figure_dir / "yokota_green_kubo_spectrum.png"
        time_series_plot = figure_dir / "yokota_green_kubo_time_series.png"

        plot_tau_curves(tau_plot, grouped, selected_modes, eta)
        plot_k_spectrum(
            spectrum_plot,
            summaries,
            eta,
            config,
            args.a_uv,
            not args.no_theory,
        )

        columns, time_rows = read_time_series(time_series_path)
        if columns and time_rows:
            plot_time_series(time_series_plot, columns, time_rows)

    sample_counts = sorted({point.count for point in points})
    print("Yokota Green-Kubo analysis")
    print(f"  sources      : {len(sources)}")
    print(f"  modes        : {len(mode_indices)}")
    print(f"  tau points   : {len({point.tau for point in points})}")
    print(f"  samples/mode : {sample_counts[0]}..{sample_counts[-1]}" if sample_counts else "  samples/mode : unknown")
    if eta is not None:
        print(f"  eta reference: {eta:g}")
    print(f"  metadata     : {display_path(metadata_path)}")
    print(f"  summary      : {display_path(summary_path)}")
    if not args.no_plots:
        print(f"  tau plot     : {display_path(tau_plot)}")
        print(f"  k plot       : {display_path(spectrum_plot)}")
        if time_series_path.exists():
            print(f"  time series  : {display_path(time_series_plot)}")


if __name__ == "__main__":
    main()

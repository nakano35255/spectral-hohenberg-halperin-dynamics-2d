#!/usr/bin/env python3
"""Analyze Yokota Green-Kubo measure output."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import argparse
import math
import os
import statistics
from typing import Iterable


ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = ROOT / "examples" / "05_yokota_green_kubo"
RESULTS = EXAMPLE / "results"

DEFAULT_YKGK = RESULTS / "yokota_green_kubo.dat"
DEFAULT_TIME_SERIES = RESULTS / "time_series_gaussian.dat"
DEFAULT_INPUT_SCRIPT = EXAMPLE / "input_incompressible.script"


@dataclass(frozen=True)
class YokotaRow:
    nsample: int
    tau: float
    mode_index: int
    kx: float
    ky: float
    s_mean: float


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
    s_tail_std: float
    tail_samples: int


def repo_path(value: str | Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return ROOT / path


def display_path(path: Path) -> str:
    return str(path.relative_to(ROOT) if path.is_relative_to(ROOT) else path)


def tokenize(line: str) -> list[str]:
    return line.split("#", 1)[0].split()


def parse_key_values(tokens: list[str]) -> dict[str, str]:
    values: dict[str, str] = {}
    for i in range(0, len(tokens), 2):
        if i + 1 < len(tokens):
            values[tokens[i]] = tokens[i + 1]
    return values


def parse_eta(input_script: Path) -> float | None:
    if not input_script.exists():
        return None

    with input_script.open() as handle:
        for raw_line in handle:
            tokens = tokenize(raw_line)
            if len(tokens) >= 4 and tokens[:3] == ["model", "transport", "constant"]:
                values = parse_key_values(tokens[3:])
                if "eta" in values:
                    return float(values["eta"])

    return None


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


def group_by_mode(rows: Iterable[YokotaRow]) -> dict[int, list[YokotaRow]]:
    grouped: dict[int, list[YokotaRow]] = {}
    for row in rows:
        grouped.setdefault(row.mode_index, []).append(row)

    for mode_rows in grouped.values():
        mode_rows.sort(key=lambda row: row.tau)

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


def summarize_modes(grouped: dict[int, list[YokotaRow]], tail_fraction: float) -> list[ModeSummary]:
    if not (0.0 < tail_fraction <= 1.0):
        raise ValueError("tail_fraction must satisfy 0 < tail_fraction <= 1")

    summaries: list[ModeSummary] = []
    for mode_index, mode_rows in grouped.items():
        if not mode_rows:
            continue

        tail_count = max(1, math.ceil(len(mode_rows) * tail_fraction))
        tail_values = [row.s_mean for row in mode_rows[-tail_count:]]
        kx = mode_rows[-1].kx
        ky = mode_rows[-1].ky
        k_abs = math.hypot(kx, ky)
        s_tail_std = statistics.stdev(tail_values) if len(tail_values) > 1 else 0.0

        summaries.append(
            ModeSummary(
                mode_index=mode_index,
                kx=kx,
                ky=ky,
                k_abs=k_abs,
                tau_min=mode_rows[0].tau,
                tau_max=mode_rows[-1].tau,
                s_last=mode_rows[-1].s_mean,
                s_tail_mean=statistics.fmean(tail_values),
                s_tail_std=s_tail_std,
                tail_samples=len(tail_values),
            )
        )

    return summaries


def is_diagonal_mode_set(summaries: list[ModeSummary], tolerance: float = 1.0e-12) -> bool:
    return bool(summaries) and all(
        summary.kx > 0.0 and abs(summary.kx - summary.ky) <= tolerance * max(1.0, abs(summary.kx), abs(summary.ky))
        for summary in summaries
    )


def write_summary(path: Path, summaries: list[ModeSummary], eta: float | None, source: Path, tail_fraction: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("# Yokota Green-Kubo analysis summary\n")
        handle.write(f"# source {display_path(source)}\n")
        handle.write(f"# tail_fraction {tail_fraction:.16e}\n")
        if eta is not None:
            handle.write(f"# eta_reference {eta:.16e}\n")
        handle.write("# columns mode_index kx ky |k| tau_min tau_max S_last S_tail_mean S_tail_std tail_samples\n")
        for summary in summaries:
            handle.write(
                f"{summary.mode_index:d} "
                f"{summary.kx:.16e} "
                f"{summary.ky:.16e} "
                f"{summary.k_abs:.16e} "
                f"{summary.tau_min:.16e} "
                f"{summary.tau_max:.16e} "
                f"{summary.s_last:.16e} "
                f"{summary.s_tail_mean:.16e} "
                f"{summary.s_tail_std:.16e} "
                f"{summary.tail_samples:d}\n"
            )


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
    grouped: dict[int, list[YokotaRow]],
    selected_modes: list[int],
    eta: float | None,
) -> None:
    plt = import_pyplot()
    fig, ax = plt.subplots(figsize=(8.0, 5.0))

    for mode_index in selected_modes:
        rows = grouped[mode_index]
        label = f"m={mode_index}, k=({rows[-1].kx:.3g},{rows[-1].ky:.3g})"
        ax.plot([row.tau for row in rows], [row.s_mean for row in rows], marker=".", linewidth=1.2, label=label)

    if eta is not None:
        ax.axhline(eta, color="black", linestyle="--", linewidth=1.0, label=f"eta={eta:g}")

    ax.set_xlabel("tau")
    ax.set_ylabel("S_mean")
    ax.set_title("Yokota Green-Kubo: S_mean(tau)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_k_spectrum(output_path: Path, summaries: list[ModeSummary], eta: float | None) -> None:
    plt = import_pyplot()
    diagonal = is_diagonal_mode_set(summaries)

    if diagonal:
        x_values = [summary.kx for summary in summaries]
        x_label = "kx (= ky)"
    else:
        x_values = [summary.k_abs for summary in summaries]
        x_label = "|k|"

    y_values = [summary.s_tail_mean for summary in summaries]
    y_errors = [summary.s_tail_std for summary in summaries]

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.errorbar(x_values, y_values, yerr=y_errors, marker="o", linestyle="-", capsize=3)

    if eta is not None:
        ax.axhline(eta, color="black", linestyle="--", linewidth=1.0, label=f"eta={eta:g}")
        ax.legend()

    ax.set_xlabel(x_label)
    ax.set_ylabel("tail mean of S_mean")
    ax.set_title("Yokota Green-Kubo: long-tau mode dependence")
    ax.grid(True, alpha=0.3)
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

    plot_columns = [name for name in ("E_K", "pi_xy") if name in column_index]
    if not plot_columns:
        return

    plt = import_pyplot()
    fig, axes = plt.subplots(len(plot_columns), 1, figsize=(8.0, 2.8 * len(plot_columns)), sharex=True)
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
    parser.add_argument("--input", default=str(DEFAULT_YKGK), help="Yokota Green-Kubo output file.")
    parser.add_argument("--input-script", default=str(DEFAULT_INPUT_SCRIPT), help="Input script used to read eta when --eta is omitted.")
    parser.add_argument("--time-series", default=str(DEFAULT_TIME_SERIES), help="Optional time_series file for E_K/pi_xy plot.")
    parser.add_argument("--output-dir", default=str(RESULTS), help="Directory for figures and summary.")
    parser.add_argument("--eta", type=float, default=None, help="Reference eta line. Defaults to model transport eta from --input-script.")
    parser.add_argument("--tail-fraction", type=float, default=0.2, help="Fraction of largest tau points used for the tail summary.")
    parser.add_argument("--max-modes", type=int, default=6, help="Maximum number of modes to plot in S_mean(tau).")
    parser.add_argument("--modes", type=int, nargs="*", default=None, help="Specific mode_index values for the S_mean(tau) plot.")
    parser.add_argument("--no-plots", action="store_true", help="Only write the text summary.")
    return parser


def main() -> None:
    args = build_argument_parser().parse_args()

    input_path = repo_path(args.input)
    input_script = repo_path(args.input_script)
    time_series_path = repo_path(args.time_series)
    output_dir = repo_path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    configure_plot_cache(output_dir)

    eta = args.eta if args.eta is not None else parse_eta(input_script)

    rows = read_yokota_table(input_path)
    grouped = group_by_mode(rows)
    mode_indices = sorted(grouped)
    selected_modes = choose_modes(mode_indices, args.modes, args.max_modes)
    summaries = summarize_modes(grouped, args.tail_fraction)

    summary_path = output_dir / "yokota_green_kubo_summary.dat"
    write_summary(summary_path, summaries, eta, input_path, args.tail_fraction)

    if not args.no_plots:
        plot_tau_curves(output_dir / "yokota_green_kubo_tau.png", grouped, selected_modes, eta)
        plot_k_spectrum(output_dir / "yokota_green_kubo_spectrum.png", summaries, eta)

        columns, time_rows = read_time_series(time_series_path)
        if columns and time_rows:
            plot_time_series(output_dir / "yokota_green_kubo_time_series.png", columns, time_rows)

    nsamples = sorted({row.nsample for row in rows})
    print("Yokota Green-Kubo analysis")
    print(f"  input        : {display_path(input_path)}")
    print(f"  modes        : {len(mode_indices)}")
    print(f"  tau points   : {len({row.tau for row in rows})}")
    print(f"  nsample      : {nsamples[-1] if nsamples else 'unknown'}")
    if eta is not None:
        print(f"  eta reference: {eta:g}")
    print(f"  summary      : {display_path(summary_path)}")
    if not args.no_plots:
        print(f"  tau plot     : {display_path(output_dir / 'yokota_green_kubo_tau.png')}")
        print(f"  k plot       : {display_path(output_dir / 'yokota_green_kubo_spectrum.png')}")
        if time_series_path.exists():
            print(f"  time series  : {display_path(output_dir / 'yokota_green_kubo_time_series.png')}")


if __name__ == "__main__":
    main()

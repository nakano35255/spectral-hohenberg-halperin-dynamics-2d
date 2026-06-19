#!/usr/bin/env python3
"""Compare Yokota Green-Kubo viscosity with Example 02 Kolmogorov-flow data."""

from __future__ import annotations

import argparse
import csv
import glob
import os
from pathlib import Path

import numpy as np


CASE_DIR = Path(__file__).resolve().parent
ROOT = CASE_DIR.parents[2]
DEFAULT_CASE = "eta0_0p1_grid256_L256_dt0p01_T25000_diag_n576"
DEFAULT_YKGK_SUMMARY = CASE_DIR / "processed_data" / DEFAULT_CASE / "yokota_green_kubo_summary.csv"
DEFAULT_KOLOMOGOROV_PROCESSED_GLOB = (
    ROOT
    / "examples"
    / "02_kolomogorov_flow"
    / "incompressible"
    / "viscosity"
    / "processed_data"
    / "*"
    / "steady_response_t30000.csv"
)
DEFAULT_OUTPUT = CASE_DIR / "figures" / "yokota_green_kubo_kolmogorov_compare.png"


def repo_path(value: str | Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return ROOT / path


def configure_plot_cache(output: Path) -> None:
    cache_root = output.parent / ".plot_cache"
    mpl_config = cache_root / "matplotlib"
    xdg_cache = cache_root / "xdg"
    mpl_config.mkdir(parents=True, exist_ok=True)
    xdg_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config))
    os.environ.setdefault("XDG_CACHE_HOME", str(xdg_cache))


def read_ykgk_summary_csv(path: Path) -> dict[str, np.ndarray]:
    rows = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(
                (
                    int(row["mode_index"]),
                    float(row["kx"]),
                    float(row["ky"]),
                    float(row["k_abs"]),
                    float(row["S_tail_mean"]),
                    float(row["S_tail_sem"]),
                    int(row["sample_count_min"]),
                )
            )

    if not rows:
        raise RuntimeError(f"no data rows found in {path}")

    rows.sort(key=lambda row: row[1])
    return {
        "mode": np.asarray([row[0] for row in rows], dtype=int),
        "kx": np.asarray([row[1] for row in rows], dtype=float),
        "ky": np.asarray([row[2] for row in rows], dtype=float),
        "k_abs": np.asarray([row[3] for row in rows], dtype=float),
        "eta": np.asarray([row[4] for row in rows], dtype=float),
        "eta_sem": np.asarray([row[5] for row in rows], dtype=float),
        "samples": np.asarray([row[6] for row in rows], dtype=int),
    }


def read_ykgk_summary_dat(path: Path) -> dict[str, np.ndarray]:
    rows = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 12:
                raise RuntimeError(f"unexpected Yokota Green-Kubo summary row: {line}")
            rows.append(
                (
                    int(fields[0]),
                    float(fields[1]),
                    float(fields[2]),
                    float(fields[3]),
                    float(fields[7]),
                    float(fields[8]),
                    int(fields[11]),
                )
            )

    if not rows:
        raise RuntimeError(f"no data rows found in {path}")

    rows.sort(key=lambda row: row[1])
    return {
        "mode": np.asarray([row[0] for row in rows], dtype=int),
        "kx": np.asarray([row[1] for row in rows], dtype=float),
        "ky": np.asarray([row[2] for row in rows], dtype=float),
        "k_abs": np.asarray([row[3] for row in rows], dtype=float),
        "eta": np.asarray([row[4] for row in rows], dtype=float),
        "eta_sem": np.asarray([row[5] for row in rows], dtype=float),
        "samples": np.asarray([row[6] for row in rows], dtype=int),
    }


def read_ykgk_summary(path: Path) -> dict[str, np.ndarray]:
    if path.suffix == ".csv":
        return read_ykgk_summary_csv(path)
    return read_ykgk_summary_dat(path)


def load_kolmogorov_processed(patterns: list[str]) -> list[dict[str, object]]:
    paths: list[Path] = []
    for pattern in patterns:
        paths.extend(Path(value) for value in glob.glob(str(repo_path(pattern))))
    paths = sorted({path.resolve() for path in paths})
    if not paths:
        joined = ", ".join(patterns)
        raise FileNotFoundError(f"no Example 02 processed data matched: {joined}")

    datasets = []
    for path in paths:
        rows = []
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                rows.append(
                    (
                        float(row["k"]),
                        float(row["eta_eff_mean"]),
                        float(row["eta_eff_sem"]),
                        int(row["n_samples"]),
                    )
                )
        rows.sort(key=lambda row: row[0])
        datasets.append(
            {
                "label": path.parent.name,
                "q": np.asarray([row[0] for row in rows], dtype=float),
                "eta": np.asarray([row[1] for row in rows], dtype=float),
                "eta_sem": np.asarray([row[2] for row in rows], dtype=float),
                "samples": np.asarray([row[3] for row in rows], dtype=int),
            }
        )
    return datasets


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ykgk-summary", default=str(DEFAULT_YKGK_SUMMARY))
    parser.add_argument(
        "--kolmogorov-processed-glob",
        action="append",
        default=None,
        help="Glob for Example 02 processed steady_response CSV files. Repeatable.",
    )
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--title", default="Example 03 YKGK vs Example 02 Kolmogorov flow")
    parser.add_argument(
        "--ykgk-wave-number",
        choices=("component", "abs"),
        default="component",
        help="Use kx=ky for the diagonal YKGK wave number, or |k| for the vector magnitude.",
    )
    return parser


def main() -> None:
    args = build_argument_parser().parse_args()

    ykgk_summary = repo_path(args.ykgk_summary)
    kolmogorov_patterns = args.kolmogorov_processed_glob or [str(DEFAULT_KOLOMOGOROV_PROCESSED_GLOB)]
    output = repo_path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    configure_plot_cache(output)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ykgk = read_ykgk_summary(ykgk_summary)
    kolmogorov = load_kolmogorov_processed(kolmogorov_patterns)
    if args.ykgk_wave_number == "component":
        ykgk_q = ykgk["kx"]
        ykgk_label = "Example 03 YKGK diagonal"
        x_label = r"$k$"
    else:
        ykgk_q = ykgk["k_abs"]
        ykgk_label = r"Example 03 YKGK diagonal, $|\mathbf{k}|$"
        x_label = r"$|\mathbf{k}|$"

    fig, ax = plt.subplots(figsize=(7.5, 5.3))

    ax.errorbar(
        ykgk_q,
        ykgk["eta"],
        yerr=ykgk["eta_sem"],
        fmt="o-",
        ms=4.5,
        lw=1.5,
        capsize=2.5,
        color="#2563eb",
        label=ykgk_label,
        alpha=0.95,
    )

    colors = ["#dc2626", "#16a34a", "#7c3aed"]
    markers = ["s", "^", "D"]
    for index, dataset in enumerate(kolmogorov):
        ax.errorbar(
            dataset["q"],
            dataset["eta"],
            yerr=dataset["eta_sem"],
            fmt=markers[index % len(markers)] + "-",
            ms=4.5,
            lw=1.35,
            capsize=2.5,
            color=colors[index % len(colors)],
            label=f"Example 02 Kolmogorov {dataset['label']}",
            alpha=0.92,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(x_label)
    ax.set_ylabel(r"$\eta_{\rm eff}(k)$")
    ax.set_title(args.title)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=True, fontsize=8.7)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)

    print(f"saved {output}")
    print(
        f"YKGK modes: {len(ykgk_q)}, "
        f"{args.ykgk_wave_number} wave number=[{ykgk_q.min():.6g}, {ykgk_q.max():.6g}]"
    )
    for dataset in kolmogorov:
        print(
            f"Kolmogorov {dataset['label']}: "
            f"points={len(dataset['q'])}, q=[{dataset['q'].min():.6g}, {dataset['q'].max():.6g}]"
        )


if __name__ == "__main__":
    main()

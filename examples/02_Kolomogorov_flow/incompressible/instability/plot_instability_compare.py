#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path

default_mpl_cache = Path("/tmp/mplcache")
default_mpl_cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(default_mpl_cache))

import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
INCOMPRESSIBLE_DIR = CASE_DIR.parent
VISCOSITY_DIR = INCOMPRESSIBLE_DIR / "viscosity"
RAW_DATA_DIR = CASE_DIR / "raw_data"
PROCESSED_DATA_DIR = CASE_DIR / "processed_data"
sys.path.insert(0, str(VISCOSITY_DIR))

import plot_theory_compare as visc


CASES = (
    ("dt=0.01", RAW_DATA_DIR / "eta0_0.1_U0.25_dt0.01", "eta0_0.1_U0.25_dt0.01", "#2563eb", "o"),
    ("dt=0.005", RAW_DATA_DIR / "eta0_0.1_U0.25_dt0.005", "eta0_0.1_U0.25_dt0.005", "#dc2626", "s"),
)


def time_label(value):
    return f"{value:g}".replace(".", "p")


def default_processed_data_root():
    return PROCESSED_DATA_DIR


def default_output_path():
    return CASE_DIR / "figures" / "eta0_0.1_U0.25_instability_compare.png"


def reference_processed_data_dir(steady_start):
    return VISCOSITY_DIR / "processed_data" / "eta0_0.1_U0.025_dt0.01"


def load_reference_rows(steady_start):
    processed_data_dir = reference_processed_data_dir(steady_start)
    if not (processed_data_dir / "metadata.csv").exists() or not (processed_data_dir / "mode_time_series.csv").exists():
        raise RuntimeError(f"missing viscosity reference processed data: {processed_data_dir}")
    return visc.load_processed_case(processed_data_dir, steady_start)


def reference_arrays(reference_rows):
    nks = np.asarray([row["nk"] for row in reference_rows], dtype=float)
    ks = np.asarray([row["k"] for row in reference_rows], dtype=float)
    eta = np.asarray([row["eta_mean"] for row in reference_rows], dtype=float)
    return nks, ks, eta


def force_amplitude(row):
    if "amplitude" in row:
        return row["amplitude"]
    return row["eta_mean"] * row["k"] * row["k"] * row["u_mean"]


def eta_ref_at_nk(nk_values, reference_rows):
    ref_nks, _, ref_eta = reference_arrays(reference_rows)
    return np.interp(np.asarray(nk_values, dtype=float), ref_nks, ref_eta)


def kolmogorov_re_critical(nk):
    nk = np.asarray(nk, dtype=float)
    alpha = 1.0 / nk
    return np.sqrt(2.0) * (1.0 + alpha * alpha) / np.sqrt(1.0 - alpha * alpha)


def unstable_nk_ranges(reference_rows, target_u, max_nk):
    ref_nks, _, ref_eta = reference_arrays(reference_rows)
    length = reference_rows[0]["length"][1]
    rho0 = reference_rows[0]["rho0"]
    grid = np.linspace(1.0001, max_nk, 6000)
    eta_grid = np.interp(grid, ref_nks, ref_eta)
    q_grid = 2.0 * np.pi * grid / length
    critical_u = eta_grid * q_grid * kolmogorov_re_critical(grid) / rho0
    unstable = target_u > critical_u

    ranges = []
    start = None
    previous = grid[0]
    for nk, is_unstable in zip(grid, unstable):
        if is_unstable and start is None:
            start = nk
        if (not is_unstable) and start is not None:
            ranges.append((start, previous))
            start = None
        previous = nk
    if start is not None:
        ranges.append((start, grid[-1]))
    return ranges


def unstable_nk_ranges_constant_eta(reference_rows, eta, target_u, max_nk):
    length = reference_rows[0]["length"][1]
    rho0 = reference_rows[0]["rho0"]
    grid = np.linspace(1.0001, max_nk, 6000)
    q_grid = 2.0 * np.pi * grid / length
    critical_u = eta * q_grid * kolmogorov_re_critical(grid) / rho0
    unstable = target_u > critical_u

    ranges = []
    start = None
    previous = grid[0]
    for nk, is_unstable in zip(grid, unstable):
        if is_unstable and start is None:
            start = nk
        if (not is_unstable) and start is not None:
            ranges.append((start, previous))
            start = None
        previous = nk
    if start is not None:
        ranges.append((start, grid[-1]))
    return ranges


def subtract_ranges(ranges, excluded_ranges):
    remaining = []
    for start, end in ranges:
        pieces = [(start, end)]
        for excluded_start, excluded_end in excluded_ranges:
            next_pieces = []
            for piece_start, piece_end in pieces:
                if excluded_end <= piece_start or piece_end <= excluded_start:
                    next_pieces.append((piece_start, piece_end))
                    continue
                if piece_start < excluded_start:
                    next_pieces.append((piece_start, excluded_start))
                if excluded_end < piece_end:
                    next_pieces.append((excluded_end, piece_end))
            pieces = next_pieces
        remaining.extend(pieces)
    return [(start, end) for start, end in remaining if end > start]


def add_unstable_spans(ax, ranges, length, color, alpha, label=None):
    for index, (left_nk, right_nk) in enumerate(ranges):
        left_q = 2.0 * np.pi * left_nk / length
        right_q = 2.0 * np.pi * right_nk / length
        ax.axvspan(
            left_q,
            right_q,
            color=color,
            alpha=alpha,
            lw=0,
            label=label if label and index == 0 else None,
        )


def draw_relaxation(ax, rows, reference_rows, selected_nks, steady_start, title):
    colors = plt.cm.viridis(np.linspace(0.10, 0.88, len(selected_nks)))
    by_nk = {row["nk"]: row for row in rows}
    max_time = max(float(row["times"][-1]) for row in rows)
    ax.axhline(1.0, color="#111827", lw=1.0, ls="--", alpha=0.75)
    ax.axvspan(steady_start / 1000.0, max_time / 1000.0, color="#dbeafe", alpha=0.45, lw=0)
    for index, (nk, color) in enumerate(zip(selected_nks, colors)):
        if nk not in by_nk:
            continue
        row = by_nk[nk]
        eta_ref = float(eta_ref_at_nk([nk], reference_rows)[0])
        u_ref = force_amplitude(row) / (row["k"] * row["k"] * eta_ref)
        times = row["times"] / 1000.0
        mean = row["u_mean_t"] / u_ref
        sem = row["u_sem_t"] / abs(u_ref)
        zorder = 30 - index
        ax.plot(times, mean, color=color, lw=1.45, label=rf"$n_k={nk}$", zorder=zorder)
        ax.fill_between(times, mean - sem, mean + sem, color=color, alpha=0.12, lw=0, zorder=zorder - 0.5)
    ax.set_title(title)
    ax.set_xlabel(r"$t/10^3$")
    ax.set_ylabel(r"$U_k(t)/U_k^{\rm ref}$")
    ax.set_ylim(0.0, 1.35)
    ax.grid(True, alpha=0.22)
    ax.legend(frameon=False, fontsize=8.2, ncol=3)


def draw_eta_panel(ax, datasets, reference_rows, reference_unstable_ranges, bare_only_unstable_ranges, target_u):
    ref_nks, ref_ks, ref_eta = reference_arrays(reference_rows)
    length = reference_rows[0]["length"][1]
    add_unstable_spans(
        ax,
        bare_only_unstable_ranges,
        length,
        color="#84cc16",
        alpha=0.16,
        label=r"unstable with bare $\eta_0$ only",
    )
    add_unstable_spans(
        ax,
        reference_unstable_ranges,
        length,
        color="#f97316",
        alpha=0.13,
        label="unstable from U=0.025 reference",
    )
    ax.plot(ref_ks, ref_eta, color="#111827", lw=1.7, label="U=0.025 reference")
    for label, rows, color, marker in datasets:
        ks = np.asarray([row["k"] for row in rows])
        eta = np.asarray([row["eta_mean"] for row in rows])
        eta_sem = np.asarray([row["eta_sem"] for row in rows])
        ax.errorbar(ks, eta, yerr=eta_sem, fmt=marker, ms=4.2, capsize=2.2, color=color, label=label, alpha=0.92)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$q = 2\pi n_k/L$")
    ax.set_ylabel(r"$\eta_{\rm eff}(q)$")
    ax.set_title(rf"$U={target_u:g}$ response vs low-$U$ reference")
    ax.grid(True, which="both", alpha=0.22)
    ax.legend(frameon=False, fontsize=8.4)


def draw_ratio_panel(ax, datasets, reference_rows, reference_unstable_ranges, bare_only_unstable_ranges):
    ref_nks, _, _ = reference_arrays(reference_rows)
    length = reference_rows[0]["length"][1]
    add_unstable_spans(ax, bare_only_unstable_ranges, length, color="#84cc16", alpha=0.16)
    add_unstable_spans(ax, reference_unstable_ranges, length, color="#f97316", alpha=0.13)
    ax.axhline(1.0, color="#111827", lw=1.1, ls="--", alpha=0.75)
    for label, rows, color, marker in datasets:
        nks = np.asarray([row["nk"] for row in rows], dtype=float)
        ks = np.asarray([row["k"] for row in rows])
        eta = np.asarray([row["eta_mean"] for row in rows])
        eta_sem = np.asarray([row["eta_sem"] for row in rows])
        eta_ref = eta_ref_at_nk(nks, reference_rows)
        ax.errorbar(
            ks,
            eta / eta_ref,
            yerr=eta_sem / eta_ref,
            fmt=marker,
            ms=4.2,
            capsize=2.2,
            color=color,
            label=label,
            alpha=0.92,
        )
    ax.set_xscale("log")
    ax.set_xlabel(r"$q = 2\pi n_k/L$")
    ax.set_ylabel(r"$\eta_{\rm eff}/\eta_{\rm ref}$")
    ax.set_title("Deviation from U=0.025 reference")
    ax.grid(True, which="both", alpha=0.22)
    ax.legend(frameon=False, fontsize=8.4)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--steady-start", type=float, default=30000.0)
    parser.add_argument("--target-u", type=float, default=0.25)
    parser.add_argument("--relaxation-nks", nargs="*", type=int, default=[1, 2, 4, 8, 16, 20])
    parser.add_argument("--processed-data-root", type=Path)
    parser.add_argument("--rebuild-processed-data", action="store_true")
    parser.add_argument("--output", type=Path, default=default_output_path())
    args = parser.parse_args()
    if args.processed_data_root is None:
        args.processed_data_root = default_processed_data_root()

    datasets = visc.load_or_analyze(args.processed_data_root, CASES, args.steady_start, args.rebuild_processed_data)
    reference_rows = load_reference_rows(args.steady_start)
    max_nk = max(max(row["nk"] for row in rows) for _, rows, _, _ in datasets)
    reference_unstable_ranges = unstable_nk_ranges(reference_rows, args.target_u, max_nk)
    bare_unstable_ranges = unstable_nk_ranges_constant_eta(reference_rows, eta=0.1, target_u=args.target_u, max_nk=max_nk)
    bare_only_unstable_ranges = subtract_ranges(bare_unstable_ranges, reference_unstable_ranges)

    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.1), constrained_layout=True)
    relax_left_ax, relax_right_ax = axes[0]
    eta_ax, ratio_ax = axes[1]

    draw_relaxation(relax_left_ax, datasets[0][1], reference_rows, args.relaxation_nks, args.steady_start, datasets[0][0])
    draw_relaxation(relax_right_ax, datasets[1][1], reference_rows, args.relaxation_nks, args.steady_start, datasets[1][0])
    draw_eta_panel(eta_ax, datasets, reference_rows, reference_unstable_ranges, bare_only_unstable_ranges, args.target_u)
    draw_ratio_panel(ratio_ax, datasets, reference_rows, reference_unstable_ranges, bare_only_unstable_ranges)

    fig.suptitle(rf"$\eta_0=0.1$, $U={args.target_u:g}$, reference from $U=0.025$, steady average $t\geq {args.steady_start:g}$", y=1.02)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight", dpi=180)
    print(f"saved {args.output}")
    print("unstable nk ranges from U=0.025 reference:", ", ".join(f"{a:.3f}-{b:.3f}" for a, b in reference_unstable_ranges))
    print("unstable nk ranges from bare eta0=0.1:", ", ".join(f"{a:.3f}-{b:.3f}" for a, b in bare_unstable_ranges))
    for label, rows, _, _ in datasets:
        print(f"{label}: samples={rows[0]['samples']}")


if __name__ == "__main__":
    main()

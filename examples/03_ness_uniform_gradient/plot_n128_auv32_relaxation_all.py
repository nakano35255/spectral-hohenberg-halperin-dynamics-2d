#!/usr/bin/env python3
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analysis_common import EXAMPLE, numeric_label, parse_input, read_time_series, transform, col


MAIN_ROOT = EXAMPLE / "main"
FIGURE_DIR = EXAMPLE / "figures"
RELAXATION_OUTPUT = FIGURE_DIR / "02_main_relaxation_by_D0.png"
OVERLAY_OUTPUT = FIGURE_DIR / "03_main_observables_vs_time_by_D0.png"


QUANTITIES = [
    (0, r"$E_K$"),
    (1, r"$\langle|\psi|^2\rangle$"),
    (2, r"$D_0\langle|\nabla\delta\psi|^2\rangle/V$"),
    (3, r"$-G\langle J_{\psi,x}\rangle$"),
]


def load_main_case(d0_dir):
    _, d0 = numeric_label(d0_dir, "D0")
    result_dir = d0_dir / "results"
    run_dir = d0_dir / "runs"
    input_path = run_dir / "input_000.script"
    params = parse_input(input_path)

    files = sorted(result_dir.glob("time_series_*.dat"))
    if not files:
        raise FileNotFoundError(f"no time_series files in {result_dir}")

    times = None
    sums = None
    sums2 = None
    for path in files:
        columns, data = read_time_series(path)
        sample_times = data[:, col(columns, "time")]
        if times is None:
            times = sample_times
        elif data.shape[0] != times.shape[0] or not np.allclose(sample_times, times):
            raise ValueError(f"time grid differs in {path}")

        y = transform(data, columns, d0, params)
        if sums is None:
            sums = np.zeros_like(y)
            sums2 = np.zeros_like(y)
        sums += y
        sums2 += y * y

    samples = len(files)
    mean = sums / samples
    if samples > 1:
        variance = np.maximum(sums2 / samples - mean * mean, 0.0)
        sem = np.sqrt(variance / (samples - 1))
    else:
        sem = np.zeros_like(mean)

    return {
        "d0": d0,
        "d0_text": d0_dir.name.split("_", 1)[1],
        "params": params,
        "times": times,
        "mean": mean,
        "sem": sem,
        "samples": samples,
    }


def load_main_cases():
    d0_dirs = sorted(MAIN_ROOT.glob("D0_*"), key=lambda path: numeric_label(path, "D0")[1])
    cases = [load_main_case(d0_dir) for d0_dir in d0_dirs]
    if not cases:
        raise FileNotFoundError(f"no D0 cases found in {MAIN_ROOT}")
    return cases


def plot_indices(size, max_points=2500):
    if size <= max_points:
        return np.arange(size)
    return np.unique(np.linspace(0, size - 1, max_points).astype(int))


def case_color_map(cases):
    cmap = plt.get_cmap("viridis")
    values = np.linspace(0.08, 0.92, len(cases))
    return {id(case): cmap(value) for case, value in zip(cases, values)}


def plot_relaxation_grid(cases, colors):
    fig, axes = plt.subplots(
        len(cases),
        len(QUANTITIES),
        figsize=(16.0, 2.15 * len(cases)),
        sharex=False,
        constrained_layout=True,
    )
    if len(cases) == 1:
        axes = axes[np.newaxis, :]

    for row_index, case in enumerate(cases):
        t = case["times"]
        keep = plot_indices(t.size)
        shade_start = 0.5 * t.max()
        color = colors[id(case)]

        for col_index, (q_index, label) in enumerate(QUANTITIES):
            ax = axes[row_index, col_index]
            y = case["mean"][:, q_index]
            yerr = case["sem"][:, q_index]
            ax.plot(t[keep], y[keep], color=color, lw=1.2)
            ax.fill_between(
                t[keep],
                y[keep] - yerr[keep],
                y[keep] + yerr[keep],
                color=color,
                alpha=0.14,
                linewidth=0,
            )
            ax.axvspan(shade_start, t.max(), color="#dbeafe", alpha=0.18, linewidth=0)
            ax.grid(True, alpha=0.25, linewidth=0.7)
            ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
            if col_index in (2, 3):
                ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
            if row_index == 0:
                ax.set_title(label)
            if col_index == 0:
                ax.set_ylabel(rf"$D_0={case['d0_text']}$" + "\n" + rf"$n={case['samples']}$")
            if row_index == len(cases) - 1:
                ax.set_xlabel(r"$t$")

    params = cases[0]["params"]
    fig.suptitle(
        rf"$N={params['grid'][0]},\ a_{{uv}}={params['a_uv']:.0f},\ L={params['length'][0]:.0f},\ G={params['gradient']:.8g}$: relaxation for all main $D_0$ values",
        fontsize=13,
    )
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(RELAXATION_OUTPUT, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_overlay(cases, colors):
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.4), constrained_layout=True)
    axes = axes.ravel()

    for ax, (q_index, label) in zip(axes, QUANTITIES):
        for case in cases:
            t = case["times"]
            keep = plot_indices(t.size)
            y = case["mean"][:, q_index]
            ax.plot(
                t[keep],
                y[keep],
                color=colors[id(case)],
                lw=1.45,
                label=rf"$D_0={case['d0_text']}$, n={case['samples']}",
            )

        ax.set_title(label)
        ax.set_xlabel(r"$t$")
        ax.set_xlim(left=0.0)
        ax.grid(True, alpha=0.25, linewidth=0.7)
        ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
        if q_index in (2, 3):
            ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncols=5, frameon=False, bbox_to_anchor=(0.5, 1.05))
    params = cases[0]["params"]
    fig.suptitle(
        rf"$N={params['grid'][0]},\ a_{{uv}}={params['a_uv']:.0f},\ L={params['length'][0]:.0f},\ G={params['gradient']:.8g}$: all main $D_0$ time series",
        fontsize=13,
        y=1.11,
    )
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OVERLAY_OUTPUT, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    cases = load_main_cases()
    colors = case_color_map(cases)
    plot_relaxation_grid(cases, colors)
    plot_overlay(cases, colors)
    print(f"saved {RELAXATION_OUTPUT}")
    print(f"saved {OVERLAY_OUTPUT}")


if __name__ == "__main__":
    main()

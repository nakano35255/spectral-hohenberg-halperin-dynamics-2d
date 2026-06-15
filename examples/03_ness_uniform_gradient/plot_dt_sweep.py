#!/usr/bin/env python3
from pathlib import Path
import re

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


EXAMPLE = Path(__file__).resolve().parent
DEFAULT_ROOT = EXAMPLE / "legacy" / "dt_sweep_T10000"
DEFAULT_OUTPUT = EXAMPLE / "figures" / "01_dt_sweep_relaxation_check.png"
GRADIENT = 0.000048828125


def parse_volume(input_path):
    with input_path.open() as handle:
        for raw_line in handle:
            tokens = raw_line.split("#", 1)[0].split()
            if len(tokens) == 3 and tokens[0] == "length":
                return float(tokens[1]) * float(tokens[2])
    raise ValueError(f"missing length command in {input_path}")


def numeric_label(path, prefix):
    match = re.search(rf"{re.escape(prefix)}_([^/]+)$", path.name)
    if not match:
        raise ValueError(f"cannot parse {prefix} value from {path}")
    text = match.group(1).rstrip(",")
    return text, float(text)


def read_time_series(path):
    columns = None
    with path.open() as handle:
        for line in handle:
            if line.startswith("# step time"):
                columns = line[2:].split()
                break
    if columns is None:
        raise ValueError(f"missing time_series header: {path}")

    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] != len(columns):
        raise ValueError(f"column count mismatch in {path}")
    return columns, data


def load_case(case_dir):
    files = sorted(case_dir.glob("time_series_*.dat"))
    if not files:
        raise FileNotFoundError(f"no time_series files in {case_dir}")

    arrays = []
    columns0 = None
    times0 = None
    for path in files:
        columns, data = read_time_series(path)
        if columns0 is None:
            columns0 = columns
            times0 = data[:, 1]
        elif columns != columns0:
            raise ValueError(f"columns differ in {path}")
        elif data.shape[0] != times0.shape[0] or not np.allclose(data[:, 1], times0):
            raise ValueError(f"time grid differs in {path}")
        arrays.append(data)

    stack = np.stack(arrays, axis=0)
    mean = np.mean(stack, axis=0)
    sem = np.std(stack, axis=0, ddof=1) / np.sqrt(stack.shape[0]) if stack.shape[0] > 1 else np.zeros_like(mean)
    return columns0, mean, sem, stack.shape[0]


def discover_cases(root):
    result_root = root / "results"
    if not result_root.exists():
        raise FileNotFoundError(f"missing results directory: {result_root}")

    cases = {}
    for d0_dir in sorted(result_root.glob("D0_*"), key=lambda p: numeric_label(p, "D0")[1]):
        d0_text, d0 = numeric_label(d0_dir, "D0")
        dt_cases = []
        for dt_dir in sorted(d0_dir.glob("dt_*"), key=lambda p: numeric_label(p, "dt")[1]):
            dt_text, dt = numeric_label(dt_dir, "dt")
            columns, mean, sem, samples = load_case(dt_dir)
            dt_cases.append(
                {
                    "dt_text": dt_text,
                    "dt": dt,
                    "columns": columns,
                    "mean": mean,
                    "sem": sem,
                    "samples": samples,
                }
            )
        if dt_cases:
            cases[d0] = {"d0_text": d0_text, "runs": dt_cases}
    if not cases:
        raise FileNotFoundError(f"no D0 cases found under {result_root}")
    return cases


def column_index(columns, name):
    try:
        return columns.index(name)
    except ValueError as exc:
        raise ValueError(f"missing column {name}") from exc


def transformed_series(run, d0, volume, quantity):
    columns = run["columns"]
    mean = run["mean"]
    sem = run["sem"]

    if quantity == "E_K":
        index = column_index(columns, "E_K")
        return mean[:, index], sem[:, index]
    if quantity == "psi2":
        index = column_index(columns, "|psi[0]|^2")
        return mean[:, index], sem[:, index]
    if quantity == "diss":
        index = column_index(columns, "|d_psi[0]|^2")
        return d0 * mean[:, index] / volume, d0 * sem[:, index] / volume
    if quantity == "work":
        index = column_index(columns, "Jpsi[0]_x")
        return -GRADIENT * mean[:, index], GRADIENT * sem[:, index]
    raise ValueError(f"unknown quantity: {quantity}")


def plot(root, output):
    cases = discover_cases(root)
    volume = parse_volume(EXAMPLE / "input.script")
    quantities = [
        ("E_K", r"$E_K$"),
        ("psi2", r"$\langle|\psi|^2\rangle$"),
        ("diss", r"$D_0\langle|\nabla\psi|^2\rangle$"),
        ("work", r"$-G\langle J_{\psi,x}\rangle$"),
    ]
    colors = {
        0.5: "#2563eb",
        1.0: "#16a34a",
        2.0: "#f97316",
        4.0: "#dc2626",
    }

    nrows = len(cases)
    ncols = len(quantities)
    top_pad = 1.18 if nrows == 1 else 1.08
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.0 * ncols, 2.85 * nrows),
        sharex=True,
        squeeze=False,
        constrained_layout=True,
    )

    for row, (d0, case) in enumerate(cases.items()):
        for col, (quantity, ylabel) in enumerate(quantities):
            ax = axes[row][col]
            for run in case["runs"]:
                t = run["mean"][:, 1]
                y, y_sem = transformed_series(run, d0, volume, quantity)
                color = colors.get(run["dt"], None)
                label = rf"$\Delta t={run['dt_text']}$, n={run['samples']}"
                ax.plot(t, y, lw=1.7, color=color, label=label)
                ax.fill_between(t, y - y_sem, y + y_sem, color=color, alpha=0.13, linewidth=0)

            ax.grid(True, alpha=0.25, linewidth=0.7)
            if row == 0:
                ax.set_title(ylabel)
            if col == 0:
                ax.set_ylabel(rf"$D_0=\nu_0={case['d0_text']}$")
            if row == nrows - 1:
                ax.set_xlabel(r"$t$")
            if quantity in ("diss", "work"):
                ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncols=4, frameon=False, bbox_to_anchor=(0.5, top_pad))
    fig.suptitle("Example 03 dt sweep: sample-averaged time evolution", y=top_pad + 0.08, fontsize=13)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    plot(DEFAULT_ROOT, DEFAULT_OUTPUT)


if __name__ == "__main__":
    main()

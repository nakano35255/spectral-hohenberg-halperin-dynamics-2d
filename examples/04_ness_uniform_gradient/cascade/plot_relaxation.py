#!/usr/bin/env python3
import csv
from pathlib import Path
import re

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CASE_DIR = Path(__file__).resolve().parent
DEFAULT_RUN_NAME = "relax_D0_0p004_grid512_L16384_dt16"
DEFAULT_INPUT_ROOT = CASE_DIR / "tmp" / DEFAULT_RUN_NAME
PROCESSED_ROOT = CASE_DIR / "processed_data"
FIGURE_ROOT = CASE_DIR / "figures"


def parse_input(path):
    text = path.read_text()
    length = re.search(r"^length\s+(\S+)\s+(\S+)", text, re.M)
    grid = re.search(r"^grid\s+(\S+)\s+(\S+)", text, re.M)
    transport = re.search(r"^model\s+transport\s+constant\b(.*)$", text, re.M)
    force = re.search(r"^fix\s+\S+\s+order_parameter\s+force/gradient\s+on\s+component\s+\S+\s+direction\s+\S+\s+amplitude\s+(\S+)", text, re.M)
    timestep = re.search(r"^timestep\s+(\S+)", text, re.M)
    if length is None or grid is None or transport is None or force is None or timestep is None:
        raise RuntimeError(f"missing required parameters in {path}")
    tokens = transport.group(1).split()
    transport_values = {tokens[index]: tokens[index + 1] for index in range(0, len(tokens), 2)}
    return {
        "grid": (int(grid.group(1)), int(grid.group(2))),
        "length": (float(length.group(1)), float(length.group(2))),
        "dt": float(timestep.group(1)),
        "d0": float(transport_values["M[0,0]"]),
        "eta": float(transport_values["eta"]),
        "gradient": float(force.group(1)),
    }


def read_time_series(path):
    columns = None
    with path.open() as handle:
        for line in handle:
            if line.startswith("# step time"):
                columns = line[2:].split()
                break
    if columns is None:
        raise RuntimeError(f"missing header in {path}")
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] != len(columns):
        raise RuntimeError(f"column count mismatch in {path}")
    return columns, data


def col(columns, name):
    try:
        return columns.index(name)
    except ValueError as exc:
        raise RuntimeError(f"missing column {name}") from exc


def load_replica(replica_dir, params):
    result_paths = sorted((replica_dir / "results").glob("time_series_relax_*.dat"))
    if not result_paths:
        raise RuntimeError(f"missing time_series files in {replica_dir}")

    rows = []
    for path in result_paths:
        match = re.search(r"time_series_relax_(\d+)\.dat$", path.name)
        if not match:
            raise RuntimeError(f"cannot parse segment from {path}")
        segment = int(match.group(1))
        columns, data = read_time_series(path)
        rows.append(
            np.column_stack(
                [
                    np.full(data.shape[0], segment, dtype=float),
                    data[:, col(columns, "step")],
                    data[:, col(columns, "time")],
                    data[:, col(columns, "E_K")],
                    data[:, col(columns, "|psi[0]|^2")],
                    data[:, col(columns, "|d_psi[0]|^2")],
                    data[:, col(columns, "Jpsi[0]_x")],
                ]
            )
        )

    data = np.concatenate(rows, axis=0)
    order = np.argsort(data[:, 2])
    data = data[order]
    volume = params["length"][0] * params["length"][1]
    d0 = params["d0"]
    gradient = params["gradient"]
    grad2 = data[:, 5]
    return {
        "replica": replica_dir.name,
        "segment": data[:, 0].astype(int),
        "step": data[:, 1],
        "time": data[:, 2],
        "E_K": data[:, 3],
        "psi2": data[:, 4],
        "grad2": grad2,
        "dissipation": d0 * grad2 / volume,
        "production": -gradient * data[:, 6],
        "length": 2.0 * np.pi * np.sqrt(np.maximum(data[:, 4] / grad2, 0.0)),
    }


def mean_sem(values):
    values = np.asarray(values, dtype=float)
    mean = float(np.mean(values))
    if values.size <= 1:
        return mean, 0.0
    return mean, float(np.std(values, ddof=1) / np.sqrt(values.size))


def replica_window_average(replica, mask, quantity):
    if not np.any(mask):
        raise RuntimeError(f"empty window for {replica['replica']} {quantity}")
    return float(np.mean(replica[quantity][mask]))


def write_csv(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_window_summary(replicas, windows):
    quantities = ("E_K", "psi2", "dissipation", "production", "length")
    rows = []
    by_window = {}
    for label, start, end in windows:
        row = {"window": label, "t_start": start, "t_end": end, "n_replica": len(replicas)}
        by_window[label] = {}
        for quantity in quantities:
            values = []
            for replica in replicas:
                mask = (replica["time"] >= start) & (replica["time"] < end)
                values.append(replica_window_average(replica, mask, quantity))
            mean, sem = mean_sem(values)
            row[f"{quantity}_mean"] = f"{mean:.16e}"
            row[f"{quantity}_sem"] = f"{sem:.16e}"
            by_window[label][quantity] = mean
        rows.append(row)
    return rows, by_window


def build_drift_summary(by_window):
    pairs = [
        ("segment2_vs_segment1", "segment_001_all", "segment_002_all"),
        ("segment3_vs_segment2", "segment_002_all", "segment_003_all"),
        ("segment3_last_half_vs_first_half", "segment_003_first_half", "segment_003_last_half"),
        ("segment3_last_quarter_vs_previous", "segment_003_q3", "segment_003_q4"),
    ]
    quantities = ("E_K", "psi2", "dissipation", "production", "length")
    rows = []
    for label, baseline, comparison in pairs:
        for quantity in quantities:
            before = by_window[baseline][quantity]
            after = by_window[comparison][quantity]
            rows.append(
                {
                    "comparison": label,
                    "quantity": quantity,
                    "baseline_window": baseline,
                    "comparison_window": comparison,
                    "baseline_mean": f"{before:.16e}",
                    "comparison_mean": f"{after:.16e}",
                    "relative_change": f"{(after - before) / before:.16e}",
                }
            )
    return rows


def downsample_indices(size, max_points=2500):
    if size <= max_points:
        return np.arange(size)
    return np.unique(np.linspace(0, size - 1, max_points).astype(int))


def plot_time_series(path, replicas, params):
    quantities = [
        ("E_K", r"$E_K$", "#1f2937"),
        ("psi2", r"$\langle|\psi|^2\rangle$", "#b45309"),
        ("dissipation", r"$D_0\langle|\nabla\psi|^2\rangle/V$", "#047857"),
        ("production", r"$-G\langle J_{\psi,x}\rangle$", "#be123c"),
        ("length", r"$2\pi\sqrt{\langle|\psi|^2\rangle/\langle|\nabla\psi|^2\rangle}$", "#2563eb"),
    ]

    fig, axes = plt.subplots(3, 2, figsize=(13.5, 10.0), constrained_layout=True)
    axes = axes.ravel()
    time = replicas[0]["time"]
    time_scale = 1.0e8
    stack_by_quantity = {quantity: np.stack([replica[quantity] for replica in replicas], axis=0) for quantity, _, _ in quantities}
    keep = downsample_indices(time.size)

    for ax, (quantity, label, color) in zip(axes, quantities):
        for replica in replicas:
            ax.plot(replica["time"][keep] / time_scale, replica[quantity][keep], color=color, alpha=0.10, lw=0.8)
        stack = stack_by_quantity[quantity]
        mean = np.mean(stack, axis=0)
        sem = np.std(stack, axis=0, ddof=1) / np.sqrt(stack.shape[0])
        ax.plot(time[keep] / time_scale, mean[keep], color=color, lw=1.8, label="replica mean")
        ax.fill_between(time[keep] / time_scale, mean[keep] - sem[keep], mean[keep] + sem[keep], color=color, alpha=0.18, linewidth=0)
        for boundary in (4.0e8, 8.0e8):
            ax.axvline(boundary / time_scale, color="#64748b", lw=0.9, ls=":")
        ax.axvspan(1.0e9 / time_scale, 1.2e9 / time_scale, color="#dbeafe", alpha=0.22, linewidth=0)
        ax.set_title(label)
        ax.set_xlabel(r"$t/10^8$")
        ax.grid(True, alpha=0.25)
        if quantity in ("dissipation", "production"):
            ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))

    axes[-1].axis("off")
    axes[0].legend(frameon=False, loc="best")
    fig.suptitle(
        rf"Kugui relaxation, $D_0=\eta={params['d0']:.3g}$, "
        rf"$N={params['grid'][0]}$, $L={params['length'][0]:.0f}$, "
        rf"$G={params['gradient']:.8g}$; segments 1-3, $n={len(replicas)}$",
        fontsize=13,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_window_drift(path, window_rows):
    quantities = [
        ("psi2", r"$\langle|\psi|^2\rangle$", "#b45309"),
        ("dissipation", r"$D_0\langle|\nabla\psi|^2\rangle/V$", "#047857"),
        ("production", r"$-G\langle J_{\psi,x}\rangle$", "#be123c"),
        ("length", r"$\ell$", "#2563eb"),
    ]
    labels = [row["window"] for row in window_rows]
    x = np.arange(len(labels))
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 7.8), constrained_layout=True)
    axes = axes.ravel()
    for ax, (quantity, ylabel, color) in zip(axes, quantities):
        means = np.asarray([float(row[f"{quantity}_mean"]) for row in window_rows])
        sems = np.asarray([float(row[f"{quantity}_sem"]) for row in window_rows])
        ax.errorbar(x, means, yerr=sems, fmt="o-", color=color, capsize=3, lw=1.4)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=8)
        ax.set_ylabel(ylabel)
        ax.grid(True, axis="y", alpha=0.25)
        if quantity in ("dissipation", "production"):
            ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    fig.suptitle("Relaxation drift by window")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    input_root = DEFAULT_INPUT_ROOT
    if not input_root.exists():
        raise RuntimeError(f"missing input root: {input_root}")

    first_input = sorted(input_root.glob("replica_*/runs/input_relax_001.script"))[0]
    params = parse_input(first_input)
    replicas = [load_replica(replica_dir, params) for replica_dir in sorted(input_root.glob("replica_*"))]
    if not replicas:
        raise RuntimeError(f"no replicas in {input_root}")

    reference_time = replicas[0]["time"]
    for replica in replicas[1:]:
        if replica["time"].shape != reference_time.shape or not np.allclose(replica["time"], reference_time):
            raise RuntimeError(f"time grid differs for {replica['replica']}")

    windows = [
        ("segment_001_all", 0.0, 4.0e8),
        ("segment_002_all", 4.0e8, 8.0e8),
        ("segment_003_all", 8.0e8, 1.2e9),
        ("segment_003_first_half", 8.0e8, 1.0e9),
        ("segment_003_last_half", 1.0e9, 1.2e9 + 1.0),
        ("segment_003_q3", 1.0e9, 1.1e9),
        ("segment_003_q4", 1.1e9, 1.2e9 + 1.0),
    ]

    processed_dir = PROCESSED_ROOT / DEFAULT_RUN_NAME
    window_rows, by_window = build_window_summary(replicas, windows)
    drift_rows = build_drift_summary(by_window)
    write_csv(
        processed_dir / "relaxation_window_summary.csv",
        list(window_rows[0].keys()),
        window_rows,
    )
    write_csv(
        processed_dir / "relaxation_drift_summary.csv",
        list(drift_rows[0].keys()),
        drift_rows,
    )

    plot_time_series(FIGURE_ROOT / "02_kugui_relaxation_D0_0p004_grid512_segments1_to3.png", replicas, params)
    plot_window_drift(FIGURE_ROOT / "03_kugui_relaxation_D0_0p004_grid512_window_drift.png", window_rows)

    print(f"replicas={len(replicas)} points_per_replica={reference_time.size}")
    for row in drift_rows:
        if row["comparison"] in ("segment3_vs_segment2", "segment3_last_quarter_vs_previous"):
            print(
                row["comparison"],
                row["quantity"],
                f"{100.0 * float(row['relative_change']):+.3f}%",
            )
    print(f"wrote {processed_dir / 'relaxation_window_summary.csv'}")
    print(f"wrote {processed_dir / 'relaxation_drift_summary.csv'}")


if __name__ == "__main__":
    main()

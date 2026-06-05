#!/usr/bin/env python3
"""Create results.png from a case-local Ohtaka dt sweep."""

import argparse
import math
import os
import re
import tempfile
from pathlib import Path

default_mpl_cache = Path(tempfile.gettempdir()) / "shhd-matplotlib-cache"
default_mpl_cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(default_mpl_cache))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ENERGY_COLUMNS = ("E_T", "E_K", "E_psi", "E_C")
CASE_DIR = Path(__file__).resolve().parent
CASE_NAME = CASE_DIR.name
WORK_EXAMPLE_BASE = Path("/work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/01_best_timestep_equilibrium")


def parse_float_id(text):
    return float(text.replace("m", "-").replace("p", "."))


def read_required(pattern, text, name):
    match = re.search(pattern, text, flags=re.MULTILINE)
    if not match:
        raise RuntimeError("cannot read {} from input script".format(name))
    return match


def read_config(run_dir):
    text = (run_dir / "input_000.script").read_text()
    grid = read_required(r"^grid\s+(\d+)\s+(\d+)", text, "grid")
    time_evolution = read_required(r"^time_evolution\s+(\S+)", text, "time_evolution")
    timestep = read_required(r"^timestep\s+(\S+)", text, "timestep")
    kbt = read_required(r"\bfix\s+\S+\s+\S+\s+noise\s+on\b.*\bkBT\s+(\S+)", text, "kBT")
    return {
        "active_nx": int(grid.group(1)),
        "active_ny": int(grid.group(2)),
        "time_evolution": time_evolution.group(1),
        "dt": float(timestep.group(1)),
        "kBT": float(kbt.group(1)),
    }


def scalar_nonzero_dof(active_nx, active_ny):
    dof = 0
    ky_cut = active_ny // 2
    for gx in range(active_nx // 2):
        weight = 1 if gx == 0 else 2
        for ky in range(-ky_cut + 1, ky_cut):
            if gx == 0 and ky == 0:
                continue
            dof += weight
    return dof


def theory_values(config):
    scalar_energy = 0.5 * scalar_nonzero_dof(config["active_nx"], config["active_ny"]) * config["kBT"]
    if "incompressible" in config["time_evolution"]:
        return {"E_T": 2.0 * scalar_energy, "E_K": scalar_energy, "E_psi": scalar_energy}
    return {"E_T": 4.0 * scalar_energy, "E_K": 2.0 * scalar_energy, "E_psi": scalar_energy, "E_C": scalar_energy}


def read_time_series(path):
    header = None
    rows = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("# step"):
                header = line[1:].split()
            elif not line.startswith("#"):
                values = line.split()
                if header is not None and len(values) != len(header):
                    continue
                rows.append([float(value) for value in values])
    if header is None:
        raise RuntimeError("missing time_series header: {}".format(path))
    if not rows:
        raise RuntimeError("no complete time_series rows: {}".format(path))
    return header, np.asarray(rows, dtype=float)


def group_dt_from_name(name):
    if name.startswith("dt"):
        return parse_float_id(name[2:])
    match = re.search(r"_dt(.+)$", name)
    if match:
        return parse_float_id(match.group(1))
    return math.inf


def group_sort_key(path, runs_dir):
    input_path = runs_dir / path.name / "input_000.script"
    if input_path.exists():
        try:
            return (read_config(input_path.parent)["dt"], path.name)
        except Exception:
            pass
    return (group_dt_from_name(path.name), path.name)


def load_group(result_dir, runs_dir):
    run_dir = runs_dir / result_dir.name
    if not run_dir.exists():
        raise RuntimeError("missing run directory: {}".format(run_dir))
    config = read_config(run_dir)
    files = sorted(result_dir.glob("time_series_*.dat"))
    if not files:
        raise RuntimeError("no time_series files in {}".format(result_dir))

    headers = []
    arrays = []
    for path in files:
        header, data = read_time_series(path)
        headers.append(header)
        arrays.append(data)

    header0 = headers[0]
    if any(header != header0 for header in headers):
        raise RuntimeError("header mismatch in {}".format(result_dir))
    if any(data.shape[1] != arrays[0].shape[1] for data in arrays):
        raise RuntimeError("column-count mismatch in {}".format(result_dir))
    min_rows = min(data.shape[0] for data in arrays)
    arrays = [data[:min_rows, :] for data in arrays]
    if any(not np.allclose(data[:, :2], arrays[0][:, :2]) for data in arrays):
        raise RuntimeError("time-grid mismatch in {}".format(result_dir))

    return {
        "name": result_dir.name,
        "config": config,
        "header": header0,
        "time": arrays[0][:, 1],
        "stack": np.stack(arrays, axis=0),
    }


def finite_mean_and_stderr(values):
    finite = np.isfinite(values)
    counts = finite.sum(axis=0)
    sums = np.where(finite, values, 0.0).sum(axis=0)
    mean = np.full(values.shape[1], np.nan)
    np.divide(sums, counts, out=mean, where=(counts > 0))

    centered = np.where(finite, values - mean[None, :], 0.0)
    variance = np.full(values.shape[1], np.nan)
    denom = counts - 1
    np.divide((centered * centered).sum(axis=0), denom, out=variance, where=(denom > 0))
    stderr = np.sqrt(variance) / np.sqrt(np.maximum(counts, 1))
    return mean, stderr


def colors_for_groups(groups):
    cmap = plt.get_cmap("viridis")
    if len(groups) == 1:
        return [cmap(0.35)]
    return [cmap(index / float(len(groups) - 1)) for index in range(len(groups))]


def discover_groups(results_dir, runs_dir, requested):
    if requested:
        result_dirs = [results_dir / name for name in requested]
    else:
        result_dirs = [
            path for path in results_dir.iterdir()
            if path.is_dir() and list(path.glob("time_series_*.dat"))
        ]
        result_dirs = sorted(result_dirs, key=lambda path: group_sort_key(path, runs_dir))
    return [load_group(path, runs_dir) for path in result_dirs]


def plot_results(path, groups, columns, theory, with_stderr):
    colors = colors_for_groups(groups)
    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.2), sharex=True)
    axes_flat = list(axes.flat)

    for ax, column in zip(axes_flat, columns):
        for group, color in zip(groups, colors):
            idx = group["header"].index(column)
            values = group["stack"][:, :, idx]
            mean, stderr = finite_mean_and_stderr(values)
            label = "dt={:g}".format(group["config"]["dt"])
            ax.plot(group["time"], mean, color=color, linewidth=1.9, label=label)
            if with_stderr:
                ax.fill_between(group["time"], mean - stderr, mean + stderr, color=color, alpha=0.14, linewidth=0)
        if column in theory:
            ax.axhline(theory[column], color="black", linestyle="--", linewidth=1.0, alpha=0.65)
        ax.set_title(column)
        ax.set_ylabel("sample mean energy")
        ax.grid(True, alpha=0.25)

    for ax in axes_flat[len(columns):]:
        ax.set_visible(False)
    for ax in axes[-1, :]:
        ax.set_xlabel("time")
    axes_flat[min(1, len(columns) - 1)].legend(frameon=False, fontsize=8, ncol=1)
    fig.suptitle("{} sample-averaged time series by dt".format(CASE_NAME))
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", type=Path, default=Path(os.environ.get("SHHD_OUTPUT_ROOT", WORK_EXAMPLE_BASE / CASE_NAME)))
    parser.add_argument("--runs-dir", type=Path)
    parser.add_argument("--results-dir", type=Path)
    parser.add_argument("--output", type=Path, default=CASE_DIR / "results.png")
    parser.add_argument("--group", action="append", default=None)
    parser.add_argument("--no-stderr", action="store_true")
    args = parser.parse_args()

    runs_dir = args.runs_dir or args.output_root / "runs"
    results_dir = args.results_dir or args.output_root / "results"
    groups = discover_groups(results_dir, runs_dir, args.group)
    if not groups:
        raise RuntimeError("no dt-sweep result groups found")

    columns = [column for column in ENERGY_COLUMNS if all(column in group["header"] for group in groups)]
    if not columns:
        raise RuntimeError("no plottable energy columns found")

    theory = theory_values(groups[0]["config"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    plot_results(args.output, groups, columns, theory, not args.no_stderr)

    print("output {}".format(args.output))
    print("groups " + " ".join("{}:dt={:g},samples={}".format(group["name"], group["config"]["dt"], group["stack"].shape[0]) for group in groups))


if __name__ == "__main__":
    main()

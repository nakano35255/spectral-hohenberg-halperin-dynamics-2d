#!/usr/bin/env python3
import argparse
import csv
import os
import re
from pathlib import Path

default_mpl_cache = Path("/tmp/mplcache")
default_mpl_cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(default_mpl_cache))

import matplotlib.pyplot as plt
import numpy as np


CASE_DIR = Path(__file__).resolve().parent
RAW_DATA_DIR = CASE_DIR / "raw_data"
PROCESSED_DATA_DIR = CASE_DIR / "processed_data"


def cases_for_eta0(eta0_label):
    if eta0_label == "0.1":
        return (
            ("U=0.01", RAW_DATA_DIR / "eta0_0.1_U0.01", "eta0_0.1_U0.01_dt0.01", "#2563eb", "o"),
            ("U=0.025", RAW_DATA_DIR / "eta0_0.1_U0.025", "eta0_0.1_U0.025_dt0.01", "#dc2626", "s"),
        )
    if eta0_label == "0.5":
        return (
            ("U=0.025", RAW_DATA_DIR / "eta0_0.5_U0.025", "eta0_0.5_U0.025_dt0.01", "#dc2626", "s"),
        )
    raise RuntimeError(f"unknown eta0 label: {eta0_label}")


def eta0_file_label(eta0_label):
    return eta0_label


def time_label(value):
    return f"{value:g}".replace(".", "p")


def extract_name_value(name, key):
    for pattern in (rf"(?:^|_){re.escape(key)}_([^_]+)", rf"(?:^|_){re.escape(key)}([^_]+)"):
        match = re.search(pattern, name)
        if match is not None:
            return match.group(1)
    return ""


def stderr(values):
    values = np.asarray(values, dtype=float)
    if values.size <= 1:
        return np.nan
    return values.std(ddof=1) / np.sqrt(values.size)


def parse_grid(path):
    text = path.read_text()
    match = re.search(r"^grid\s+(\S+)\s+(\S+)", text, re.M)
    if match is None:
        raise RuntimeError(f"cannot parse grid in {path}")
    return tuple(int(value) for value in match.groups())


def parse_input(path):
    text = path.read_text()
    length = tuple(float(value) for value in re.search(r"^length\s+(\S+)\s+(\S+)", text, re.M).groups())
    eta0 = float(re.search(r"^model\s+transport\s+constant\s+eta\s+(\S+)", text, re.M).group(1))
    dt = float(re.search(r"^timestep\s+(\S+)", text, re.M).group(1))
    dealias = re.search(r"^dealias\s+(\S+)", text, re.M)
    time_evolution = re.search(r"^time_evolution\s+(\S+)", text, re.M)
    force = re.search(
        r"^fix\s+\S+\s+momentum\s+force/sine\s+on\s+component\s+(\S+)\s+axis\s+(\S+)\s+nk\s+(\S+)\s+amplitude\s+(\S+)",
        text,
        re.M,
    )
    if force is None:
        raise RuntimeError(f"cannot parse force/sine in {path}")
    component, axis, nk, amplitude = force.groups()
    if component != "x" or axis != "y":
        raise RuntimeError("this plot assumes force/sine component x axis y")
    noise = re.search(r"^fix\s+\S+\s+momentum\s+noise\s+on\s+.*?\bkBT\s+(\S+)", text, re.M)
    density = re.search(r"^set\s+density\s+uniform\s+value\s+(\S+)", text, re.M)
    return {
        "length": length,
        "dt": dt,
        "eta0": eta0,
        "nk": int(nk),
        "amplitude": float(amplitude),
        "dealias": dealias.group(1) if dealias else "",
        "time_evolution": time_evolution.group(1) if time_evolution else "",
        "kBT": float(noise.group(1)) if noise else 1.0,
        "rho0": float(density.group(1)) if density else 1.0,
    }


def fns_eta(k, eta0, rho0, kBT, a_uv=1.0):
    return np.sqrt(eta0 * eta0 + (rho0 * kBT) / (8.0 * np.pi) * np.log(2.0 * np.pi / (a_uv * k)))


def read_profile_observables(path, k, steady_start):
    times = []
    u_modes = []
    coord = None
    sin_ky = None
    steady_vx_sum = None
    steady_count = 0
    current_time = None
    current_coord = []
    current_vx = []
    vx_index = None

    def finalize_block():
        nonlocal coord, sin_ky, steady_vx_sum, steady_count
        nonlocal current_time, current_coord, current_vx
        if current_time is None or not current_vx:
            return
        vx = np.asarray(current_vx, dtype=float)
        if coord is None:
            coord = np.asarray(current_coord, dtype=float)
            sin_ky = np.sin(k * coord)
            steady_vx_sum = np.zeros_like(coord)
        elif vx.size != coord.size:
            raise RuntimeError(f"{path}: block has {vx.size} rows, expected {coord.size}")

        times.append(current_time)
        u_modes.append(2.0 * float(np.dot(vx, sin_ky)) / coord.size)
        if current_time >= steady_start:
            steady_vx_sum += vx
            steady_count += 1

    with path.open() as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                finalize_block()
                current_time = None
                current_coord = []
                current_vx = []
                continue
            if line.startswith("# block"):
                finalize_block()
                parts = line.split()
                current_time = float(parts[parts.index("time") + 1])
                current_coord = []
                current_vx = []
                continue
            if line.startswith("# columns"):
                columns = line.split()[2:]
                vx_index = columns.index("vx")
                continue
            if line.startswith("#"):
                continue
            if current_time is None:
                raise RuntimeError(f"{path}: data row before block header")
            values = line.split()
            if coord is None:
                current_coord.append(float(values[0]))
            current_vx.append(float(values[vx_index]))

    finalize_block()
    if steady_count == 0:
        raise RuntimeError(f"{path}: no blocks in steady window")
    return np.asarray(times), np.asarray(u_modes), coord, steady_vx_sum / steady_count


def analyze_case(root, steady_start):
    rows = []
    for nk_dir in sorted(root.glob("nk_*")):
        input_path = nk_dir / "runs" / "input_000.script"
        config = parse_input(input_path)
        nk = config["nk"]
        k = 2.0 * np.pi * nk / config["length"][1]
        eta_samples = []
        u_samples = []
        mode_arrays = []
        time_grid = None
        for profile in sorted((nk_dir / "results").glob("profile_*.dat")):
            times, u_modes, _, _ = read_profile_observables(profile, k, steady_start)
            if time_grid is None:
                time_grid = times
            elif not np.allclose(times, time_grid):
                raise RuntimeError(f"time grid mismatch in {profile}")
            mode_arrays.append(u_modes)
            mask = times >= steady_start
            if not np.any(mask):
                raise RuntimeError(f"{profile}: no data for t >= {steady_start}")
            u_mean = float(np.mean(u_modes[mask]))
            u_samples.append(u_mean)
            eta_samples.append(config["amplitude"] / (k * k * u_mean))
        if not eta_samples:
            raise RuntimeError(f"no profile files in {nk_dir}")
        modes = np.stack(mode_arrays, axis=0)
        rows.append(
            {
                **config,
                "grid": parse_grid(input_path),
                "k": k,
                "times": time_grid,
                "u_mean_t": modes.mean(axis=0),
                "u_sem_t": modes.std(axis=0, ddof=1) / np.sqrt(modes.shape[0]),
                "eta_mean": float(np.mean(eta_samples)),
                "eta_sem": stderr(eta_samples),
                "u_mean": float(np.mean(u_samples)),
                "u_sem": stderr(u_samples),
                "samples": len(eta_samples),
            }
        )
    return sorted(rows, key=lambda row: row["nk"])


def format_float(value):
    return f"{float(value):.17g}"


def metadata_path(processed_data_dir):
    return processed_data_dir / "metadata.csv"


def mode_time_series_path(processed_data_dir):
    return processed_data_dir / "mode_time_series.csv"


def steady_response_path(processed_data_dir, steady_start):
    return processed_data_dir / f"steady_response_t{time_label(steady_start)}.csv"


def relative_to_case(path):
    try:
        return path.relative_to(CASE_DIR).as_posix()
    except ValueError:
        return path.as_posix()


def write_metadata(path, entries):
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["key", "value"])
        for key, value in entries:
            writer.writerow([key, value])


def read_metadata(path):
    with path.open(newline="") as fh:
        return {record["key"]: record["value"] for record in csv.DictReader(fh)}


def save_processed_case(processed_data_dir, rows, steady_start, source_raw_data_dir, case_name, processed_run_name):
    if not rows:
        raise RuntimeError("cannot save empty processed data")
    first = rows[0]
    sample_counts = [row["samples"] for row in rows]
    processed_data_dir.mkdir(parents=True, exist_ok=True)
    metadata_entries = [
        ("case", case_name),
        ("source_raw_data_dir", relative_to_case(source_raw_data_dir)),
        ("source_run_name", source_raw_data_dir.name),
        ("processed_run_name", processed_run_name),
        ("eta0", format_float(first["eta0"])),
        ("target_U", extract_name_value(processed_run_name, "U")),
        ("dt", format_float(first["dt"])),
        ("grid_x", first["grid"][0]),
        ("grid_y", first["grid"][1]),
        ("length_x", format_float(first["length"][0])),
        ("length_y", format_float(first["length"][1])),
        ("dealias", first["dealias"]),
        ("time_evolution", first["time_evolution"]),
        ("rho0", format_float(first["rho0"])),
        ("kBT", format_float(first["kBT"])),
        ("n_samples_min", min(sample_counts)),
        ("n_samples_max", max(sample_counts)),
        ("steady_start_default", format_float(steady_start)),
    ]
    write_metadata(metadata_path(processed_data_dir), metadata_entries)

    with mode_time_series_path(processed_data_dir).open("w", newline="") as fh:
        fieldnames = ["nk", "k", "amplitude", "time", "u_mean", "u_sem", "n_samples"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            for time, u_mean_t, u_sem_t in zip(row["times"], row["u_mean_t"], row["u_sem_t"]):
                writer.writerow(
                    {
                        "nk": row["nk"],
                        "k": format_float(row["k"]),
                        "amplitude": format_float(row["amplitude"]),
                        "time": format_float(time),
                        "u_mean": format_float(u_mean_t),
                        "u_sem": format_float(u_sem_t),
                        "n_samples": row["samples"],
                    }
                )

    with steady_response_path(processed_data_dir, steady_start).open("w", newline="") as fh:
        fieldnames = [
            "steady_start",
            "nk",
            "k",
            "amplitude",
            "u_steady_mean",
            "u_steady_sem",
            "eta_eff_mean",
            "eta_eff_sem",
            "n_samples",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "steady_start": format_float(steady_start),
                    "nk": row["nk"],
                    "k": format_float(row["k"]),
                    "amplitude": format_float(row["amplitude"]),
                    "u_steady_mean": format_float(row["u_mean"]),
                    "u_steady_sem": format_float(row["u_sem"]),
                    "eta_eff_mean": format_float(row["eta_mean"]),
                    "eta_eff_sem": format_float(row["eta_sem"]),
                    "n_samples": row["samples"],
                }
            )


def load_processed_case(processed_data_dir, steady_start):
    meta_path = metadata_path(processed_data_dir)
    series_path = mode_time_series_path(processed_data_dir)
    response_path = steady_response_path(processed_data_dir, steady_start)
    if not meta_path.exists() or not series_path.exists() or not response_path.exists():
        missing = [str(path) for path in (meta_path, series_path, response_path) if not path.exists()]
        raise RuntimeError(f"missing processed data: {', '.join(missing)}")

    metadata = read_metadata(meta_path)

    time_series = {}
    with series_path.open(newline="") as fh:
        for record in csv.DictReader(fh):
            nk = int(record["nk"])
            if nk not in time_series:
                time_series[nk] = {
                    "k": float(record["k"]),
                    "amplitude": float(record["amplitude"]),
                    "samples": int(record["n_samples"]),
                    "values": [],
                }
            time_series[nk]["values"].append((float(record["time"]), float(record["u_mean"]), float(record["u_sem"])))

    rows = []
    with response_path.open(newline="") as fh:
        for record in csv.DictReader(fh):
            stored_steady_start = float(record["steady_start"])
            if abs(stored_steady_start - steady_start) > 1.0e-12:
                raise RuntimeError(f"{response_path}: steady_start={stored_steady_start:g}, expected {steady_start:g}")
            nk = int(record["nk"])
            if nk not in time_series:
                raise RuntimeError(f"{series_path}: missing mode time series for nk={nk}")
            values = sorted(time_series[nk]["values"], key=lambda item: item[0])
            rows.append(
                {
                    "nk": nk,
                    "k": float(record["k"]),
                    "amplitude": float(record["amplitude"]),
                    "u_mean": float(record["u_steady_mean"]),
                    "u_sem": float(record["u_steady_sem"]),
                    "eta_mean": float(record["eta_eff_mean"]),
                    "eta_sem": float(record["eta_eff_sem"]),
                    "samples": int(record["n_samples"]),
                    "eta0": float(metadata["eta0"]),
                    "rho0": float(metadata["rho0"]),
                    "kBT": float(metadata["kBT"]),
                    "length": (float(metadata["length_x"]), float(metadata["length_y"])),
                    "grid": (int(metadata["grid_x"]), int(metadata["grid_y"])),
                    "dt": float(metadata["dt"]),
                    "dealias": metadata["dealias"],
                    "time_evolution": metadata["time_evolution"],
                    "times": np.asarray([item[0] for item in values], dtype=float),
                    "u_mean_t": np.asarray([item[1] for item in values], dtype=float),
                    "u_sem_t": np.asarray([item[2] for item in values], dtype=float),
                }
            )
    return sorted(rows, key=lambda row: row["nk"])


def load_or_analyze(processed_data_root, cases, steady_start, rebuild_processed_data):
    datasets = []
    for label, root, processed_name, color, marker in cases:
        processed_data_dir = processed_data_root / processed_name
        if not rebuild_processed_data:
            try:
                rows = load_processed_case(processed_data_dir, steady_start)
                datasets.append((label, rows, color, marker))
                continue
            except RuntimeError as error:
                print(f"rebuilding processed data for {processed_name}: {error}")

        rows = analyze_case(root, steady_start)
        save_processed_case(processed_data_dir, rows, steady_start, root, CASE_DIR.name, processed_name)
        datasets.append((label, rows, color, marker))
    return datasets


def default_processed_data_root():
    return PROCESSED_DATA_DIR


def default_output_path(eta0_label):
    return CASE_DIR / "figures" / f"eta0_{eta0_file_label(eta0_label)}_theory_compare.png"


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


def draw_dataset(ax, datasets):
    for label, rows, color, marker in datasets:
        ks = np.asarray([row["k"] for row in rows])
        eta = np.asarray([row["eta_mean"] for row in rows])
        eta_sem = np.asarray([row["eta_sem"] for row in rows])
        ax.errorbar(
            ks,
            eta,
            yerr=eta_sem,
            fmt=marker,
            ms=4.3,
            capsize=2.3,
            color=color,
            label=label,
            alpha=0.92,
        )


def draw_relaxation(ax, rows, selected_nks, steady_start, title):
    colors = plt.cm.viridis(np.linspace(0.12, 0.86, len(selected_nks)))
    by_nk = {row["nk"]: row for row in rows}
    max_time = max(float(row["times"][-1]) for row in rows)
    ax.axhline(1.0, color="#111827", lw=1.0, ls="--", alpha=0.75)
    ax.axvspan(steady_start / 1000.0, max_time / 1000.0, color="#dbeafe", alpha=0.45, lw=0)
    for index, (nk, color) in enumerate(zip(selected_nks, colors)):
        if nk not in by_nk:
            continue
        row = by_nk[nk]
        scale = row["u_mean"]
        times = row["times"] / 1000.0
        mean = row["u_mean_t"] / scale
        sem = row["u_sem_t"] / abs(scale)
        zorder = 20 - index
        ax.plot(times, mean, color=color, lw=1.45, label=rf"$n_k={nk}$", zorder=zorder)
        ax.fill_between(times, mean - sem, mean + sem, color=color, alpha=0.13, lw=0, zorder=zorder - 0.5)
    ax.set_title(title)
    ax.set_xlabel(r"$t/10^3$")
    ax.set_ylabel(r"$U_k(t)/U_k^{\rm steady}$")
    ax.set_ylim(0.0, 1.25)
    ax.grid(True, alpha=0.22)
    ax.legend(frameon=False, fontsize=8.4, ncol=2)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--eta0", choices=["0.1", "0.5"], default="0.1")
    parser.add_argument("--steady-start", type=float, default=30000.0)
    parser.add_argument("--relaxation-nks", nargs="*", type=int, default=[1, 2, 4, 8])
    parser.add_argument("--a-uv", type=float, default=1.0)
    parser.add_argument("--processed-data-root", type=Path)
    parser.add_argument("--rebuild-processed-data", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cases = cases_for_eta0(args.eta0)
    if args.processed_data_root is None:
        args.processed_data_root = default_processed_data_root()
    if args.output is None:
        args.output = default_output_path(args.eta0)

    datasets = load_or_analyze(args.processed_data_root, cases, args.steady_start, args.rebuild_processed_data)

    first_rows = datasets[0][1]
    first = first_rows[0]
    all_ks = np.asarray([row["k"] for _, rows, _, _ in datasets for row in rows])
    k_line = np.linspace(all_ks.min(), all_ks.max(), 600)
    k_data = np.asarray([row["k"] for row in first_rows])
    eta_fns_a1 = fns_eta(k_line, first["eta0"], first["rho0"], first["kBT"], a_uv=args.a_uv)
    eta_square = discrete_square_eta(k_data, first["eta0"], first["rho0"], first["kBT"], first["length"], args.a_uv)

    if len(datasets) == 1:
        fig, axes = plt.subplot_mosaic(
            [["relax", "relax"], ["continuum", "square"]],
            figsize=(11.6, 8.0),
            constrained_layout=True,
        )
        draw_relaxation(axes["relax"], datasets[0][1], args.relaxation_nks, args.steady_start, datasets[0][0])
        left_ax = axes["continuum"]
        right_ax = axes["square"]
    else:
        fig, axes = plt.subplot_mosaic(
            [["relax0", "relax1"], ["continuum", "square"]],
            figsize=(11.6, 8.0),
            constrained_layout=True,
        )
        draw_relaxation(axes["relax0"], datasets[0][1], args.relaxation_nks, args.steady_start, datasets[0][0])
        draw_relaxation(axes["relax1"], datasets[1][1], args.relaxation_nks, args.steady_start, datasets[1][0])
        left_ax = axes["continuum"]
        right_ax = axes["square"]

    draw_dataset(left_ax, datasets)
    left_ax.plot(k_line, eta_fns_a1, color="#111827", lw=1.7, label=rf"FNS $a_{{uv}}={args.a_uv:g}$")
    left_ax.set_title("Continuum circular cutoff")

    draw_dataset(right_ax, datasets)
    right_ax.plot(k_data, eta_square, color="#111827", lw=1.7, label="Discrete square sum")
    right_ax.set_title(r"Square cutoff $\Lambda_i=2\pi/a_{uv}$")

    for ax in (left_ax, right_ax):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$q = 2\pi n_k/L$")
        ax.grid(True, which="both", alpha=0.22)
        ax.legend(frameon=False, fontsize=8.7)

    left_ax.set_ylabel(r"$\eta_{\rm eff}(q)$")
    fig.suptitle(rf"$\eta_0={args.eta0}$, steady average $t\geq {args.steady_start:g}$", y=1.02)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight", dpi=180)
    print(f"saved {args.output}")
    for label, rows, _, _ in datasets:
        print(f"{label}: samples={rows[0]['samples']}")
    print(
        "discrete square: "
        f"a_uv={args.a_uv:g} "
        f"length={first['length'][0]:g}x{first['length'][1]:g} "
        f"processed_data_root={args.processed_data_root}"
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse
import os
import re
from pathlib import Path

default_mpl_cache = Path("/tmp/mplcache")
default_mpl_cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(default_mpl_cache))

import matplotlib.pyplot as plt
import numpy as np


CASE_DIR = Path(__file__).resolve().parent
CACHE_VERSION = 2


def cases_for_eta0(eta0_label):
    if eta0_label == "0.1":
        return (
            ("U=0.01", CASE_DIR / "main" / "eta0_0.1_U0.01", "#2563eb", "o"),
            ("U=0.025", CASE_DIR / "main" / "eta0_0.1_U0.025", "#dc2626", "s"),
        )
    if eta0_label == "0.5":
        return (
            ("U=0.025", CASE_DIR / "main" / "eta0_0.5_U0.025", "#dc2626", "s"),
        )
    raise RuntimeError(f"unknown eta0 label: {eta0_label}")


def eta0_file_label(eta0_label):
    return eta0_label


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


def save_cache(path, datasets, steady_start):
    payload = {
        "version": np.asarray([CACHE_VERSION], dtype=int),
        "steady_start": np.asarray([steady_start], dtype=float),
        "num_cases": np.asarray([len(datasets)], dtype=int),
    }
    for index, (label, rows, color, marker) in enumerate(datasets):
        payload[f"label_{index}"] = np.asarray(label)
        payload[f"color_{index}"] = np.asarray(color)
        payload[f"marker_{index}"] = np.asarray(marker)
        payload[f"nk_{index}"] = np.asarray([row["nk"] for row in rows], dtype=int)
        payload[f"k_{index}"] = np.asarray([row["k"] for row in rows], dtype=float)
        payload[f"eta_mean_{index}"] = np.asarray([row["eta_mean"] for row in rows], dtype=float)
        payload[f"eta_sem_{index}"] = np.asarray([row["eta_sem"] for row in rows], dtype=float)
        payload[f"u_mean_{index}"] = np.asarray([row["u_mean"] for row in rows], dtype=float)
        payload[f"u_sem_{index}"] = np.asarray([row["u_sem"] for row in rows], dtype=float)
        payload[f"samples_{index}"] = np.asarray([row["samples"] for row in rows], dtype=int)
        payload[f"eta0_{index}"] = np.asarray([rows[0]["eta0"]], dtype=float)
        payload[f"rho0_{index}"] = np.asarray([rows[0]["rho0"]], dtype=float)
        payload[f"kBT_{index}"] = np.asarray([rows[0]["kBT"]], dtype=float)
        payload[f"length_{index}"] = np.asarray(rows[0]["length"], dtype=float)
        payload[f"grid_{index}"] = np.asarray(rows[0]["grid"], dtype=int)
        for row_index, row in enumerate(rows):
            payload[f"times_{index}_{row_index}"] = np.asarray(row["times"], dtype=float)
            payload[f"u_mean_t_{index}_{row_index}"] = np.asarray(row["u_mean_t"], dtype=float)
            payload[f"u_sem_t_{index}_{row_index}"] = np.asarray(row["u_sem_t"], dtype=float)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, **payload)


def load_cache(path, steady_start):
    data = np.load(path)
    if "version" not in data or int(data["version"][0]) != CACHE_VERSION:
        raise RuntimeError(f"{path}: incompatible cache version")
    cached_steady_start = float(data["steady_start"][0])
    if abs(cached_steady_start - steady_start) > 1.0e-12:
        raise RuntimeError(f"{path}: steady_start={cached_steady_start:g}, expected {steady_start:g}")

    datasets = []
    for index in range(int(data["num_cases"][0])):
        label = str(data[f"label_{index}"])
        color = str(data[f"color_{index}"])
        marker = str(data[f"marker_{index}"])
        rows = []
        for j, nk in enumerate(data[f"nk_{index}"]):
            rows.append(
                {
                    "nk": int(nk),
                    "k": float(data[f"k_{index}"][j]),
                    "eta_mean": float(data[f"eta_mean_{index}"][j]),
                    "eta_sem": float(data[f"eta_sem_{index}"][j]),
                    "u_mean": float(data[f"u_mean_{index}"][j]),
                    "u_sem": float(data[f"u_sem_{index}"][j]),
                    "samples": int(data[f"samples_{index}"][j]),
                    "eta0": float(data[f"eta0_{index}"][0]),
                    "rho0": float(data[f"rho0_{index}"][0]),
                    "kBT": float(data[f"kBT_{index}"][0]),
                    "length": tuple(float(value) for value in data[f"length_{index}"]),
                    "grid": tuple(int(value) for value in data[f"grid_{index}"]),
                    "times": np.asarray(data[f"times_{index}_{j}"], dtype=float),
                    "u_mean_t": np.asarray(data[f"u_mean_t_{index}_{j}"], dtype=float),
                    "u_sem_t": np.asarray(data[f"u_sem_t_{index}_{j}"], dtype=float),
                }
            )
        datasets.append((label, rows, color, marker))
    return datasets


def load_or_analyze(cache_path, cases, steady_start, rebuild_cache):
    if cache_path.exists() and not rebuild_cache:
        try:
            return load_cache(cache_path, steady_start)
        except RuntimeError as error:
            print(f"rebuilding cache: {error}")

    datasets = []
    for label, root, color, marker in cases:
        rows = analyze_case(root, steady_start)
        datasets.append((label, rows, color, marker))
    save_cache(cache_path, datasets, steady_start)
    return datasets


def default_cache_path(eta0_label, steady_start):
    time_label = f"{steady_start:g}".replace(".", "p")
    return CASE_DIR / "cache" / f"eta0_{eta0_file_label(eta0_label)}_theory_compare_t{time_label}.npz"


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
    parser.add_argument("--cache", type=Path)
    parser.add_argument("--rebuild-cache", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cases = cases_for_eta0(args.eta0)
    if args.cache is None:
        args.cache = default_cache_path(args.eta0, args.steady_start)
    if args.output is None:
        args.output = default_output_path(args.eta0)

    datasets = load_or_analyze(args.cache, cases, args.steady_start, args.rebuild_cache)

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
        f"cache={args.cache}"
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse
import csv
import math
import os
from pathlib import Path
import re

default_mpl_cache = Path("/tmp/mplcache")
default_mpl_cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(default_mpl_cache))

import matplotlib.pyplot as plt
import numpy as np


CASE_DIR = Path(__file__).resolve().parent
RAW_DATA_DIR = CASE_DIR / "raw_data"
PROCESSED_DATA_DIR = CASE_DIR / "processed_data"
FIGURE_DIR = CASE_DIR / "figures"
DEFAULT_RUN_NAME = "D0_0p12_targetpsi0p1_grid128_L4096_dt4_T20000000_n72"


def time_label(value):
    return f"{value:g}".replace(".", "p").replace("+", "")


def format_float(value):
    return f"{float(value):.17g}"


def stderr(values):
    values = np.asarray(values, dtype=float)
    if values.size <= 1:
        return np.nan
    return float(values.std(ddof=1) / np.sqrt(values.size))


def parse_value(text):
    text = text.strip()
    if text == "None":
        return None
    if text.startswith("[") and text.endswith("]"):
        inner = text[1:-1].strip()
        if not inner:
            return []
        return [parse_value(item.strip().strip("'\"")) for item in inner.split(",")]
    try:
        if re.fullmatch(r"[+-]?\d+", text):
            return int(text)
        return float(text)
    except ValueError:
        return text


def parse_config(path):
    values = {}
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            key, value = line.split(maxsplit=1)
            values[key] = parse_value(value)
    return values


def read_sample_modes(path):
    by_nk = {}
    with path.open(newline="") as handle:
        for record in csv.DictReader(handle):
            nk = int(record["nk"])
            by_nk.setdefault(nk, []).append(record)
    return by_nk


def read_mode_time_series(path):
    by_nk = {}
    with path.open(newline="") as handle:
        for record in csv.DictReader(handle):
            nk = int(record["nk"])
            by_nk.setdefault(nk, []).append(record)
    for rows in by_nk.values():
        rows.sort(key=lambda record: float(record["time"]))
    return by_nk


def eq18_delta_2d(k, rho0, kBT, a_uv):
    cutoff = 2.0 * math.pi / a_uv
    if np.any(np.asarray(k) >= cutoff):
        raise RuntimeError("Eq. (18) requires k < 2*pi/a_uv")
    return kBT / (4.0 * math.pi * rho0) * np.log(cutoff / k)


def eq18_dr(k, d0, nu0, rho0, kBT, a_uv):
    delta = eq18_delta_2d(k, rho0, kBT, a_uv)
    return 0.5 * (d0 - nu0) + np.sqrt(0.25 * (d0 + nu0) ** 2 + delta)


def solve_scalar_self_consistent(correction, d0):
    low = max(0.0, 0.5 * d0)
    high = max(1.0, 2.0 * d0)

    def residual(value):
        return value - d0 - correction(value)

    low_residual = residual(low)
    high_residual = residual(high)
    while high_residual <= 0.0:
        high *= 2.0
        high_residual = residual(high)
        if high > 1.0e8:
            raise RuntimeError("failed to bracket self-consistent diffusion")

    for _ in range(120):
        mid = 0.5 * (low + high)
        mid_residual = residual(mid)
        if mid_residual <= 0.0:
            low = mid
            low_residual = mid_residual
        else:
            high = mid
            high_residual = mid_residual
    return 0.5 * (low + high)


def square_wave_vectors(length, a_uv):
    lx, ly = length
    uv_nx = int(np.floor(lx / a_uv))
    uv_ny = int(np.floor(ly / a_uv))
    if uv_nx <= 0 or uv_ny <= 0:
        raise RuntimeError("a_uv is too large for the discrete square cutoff")

    nx = np.arange(-uv_nx, uv_nx + 1, dtype=float)
    ny = np.arange(-uv_ny, uv_ny + 1, dtype=float)
    qx = 2.0 * np.pi * nx / lx
    qy = 2.0 * np.pi * ny / ly
    return np.meshgrid(qx, qy, indexing="ij")


def active_wave_vectors(length, grid):
    lx, ly = length
    nx_active, ny_active = grid
    nx = np.arange(-nx_active // 2 + 1, nx_active // 2, dtype=float)
    ny = np.arange(-ny_active // 2 + 1, ny_active // 2, dtype=float)
    nx_grid, ny_grid = np.meshgrid(nx, ny, indexing="ij")
    qx = 2.0 * np.pi * nx_grid / lx
    qy = 2.0 * np.pi * ny_grid / ly
    return nx_grid, ny_grid, qx, qy


def discrete_square_dr(q_values, d0, nu0, rho0, kBT, length, a_uv, force_axis):
    lx, ly = length
    qx_grid, qy_grid = square_wave_vectors(length, a_uv)
    k2 = qx_grid * qx_grid + qy_grid * qy_grid
    nonzero = k2 > 0.0
    if force_axis in ("y", "1"):
        transverse2 = qx_grid * qx_grid
    elif force_axis in ("x", "0"):
        transverse2 = qy_grid * qy_grid
    else:
        raise RuntimeError(f"unknown force_axis: {force_axis}")

    integrand = np.zeros_like(k2)
    integrand[nonzero] = transverse2[nonzero] / (k2[nonzero] * k2[nonzero])

    values = []
    volume = lx * ly
    for q in np.asarray(q_values, dtype=float):
        mask = nonzero & (np.sqrt(k2) >= q)
        delta = kBT * float(np.sum(integrand[mask])) / (rho0 * volume)
        values.append(0.5 * (d0 - nu0) + math.sqrt(0.25 * (d0 + nu0) ** 2 + delta))
    return np.asarray(values)


def discrete_s54_dr(q_values, d0, nu0, rho0, kBT, length, a_uv, force_axis):
    lx, ly = length
    volume = lx * ly
    qx_grid, qy_grid = square_wave_vectors(length, a_uv)
    q2 = qx_grid * qx_grid + qy_grid * qy_grid
    nonzero_q = q2 > 0.0
    values = []

    for q_external in np.asarray(q_values, dtype=float):
        if force_axis in ("y", "1"):
            kx_external = 0.0
            ky_external = q_external
        elif force_axis in ("x", "0"):
            kx_external = q_external
            ky_external = 0.0
        else:
            raise RuntimeError(f"unknown force_axis: {force_axis}")

        k2_external = q_external * q_external
        q_dot_k = qx_grid * kx_external + qy_grid * ky_external
        numerator = q2 * k2_external - q_dot_k * q_dot_k
        px = kx_external - qx_grid
        py = ky_external - qy_grid
        p2 = px * px + py * py
        mask = nonzero_q & (p2 > 1.0e-30) & (numerator > 0.0)

        numerator_values = numerator[mask]
        p2_values = p2[mask]
        q2_values = q2[mask]

        def correction(d_renormalized):
            denominator = p2_values * (nu0 * p2_values + d_renormalized * q2_values)
            integral = float(np.sum(numerator_values / denominator)) / volume
            return kBT * integral / (rho0 * k2_external)

        values.append(solve_scalar_self_consistent(correction, d0))

    return np.asarray(values)


def external_wave_vector(q_external, force_axis):
    if force_axis in ("y", "1"):
        return 0.0, q_external
    if force_axis in ("x", "0"):
        return q_external, 0.0
    raise RuntimeError(f"unknown force_axis: {force_axis}")


def s46_radial_grid(q_values, q_abs, radial_count):
    q_min = float(np.min(q_abs[q_abs > 0.0]))
    q_max = float(np.max(q_abs))
    log_grid = np.geomspace(q_min, q_max, radial_count)
    return np.unique(np.concatenate([log_grid, np.asarray(q_values, dtype=float)]))


def s46_update_grid(k_grid, d_of_q, qx, qy, q2, q_abs, d0, nu0, rho0, kBT, volume, force_axis):
    d_q = np.interp(q_abs, k_grid, d_of_q, left=d_of_q[0], right=d_of_q[-1])
    updated = np.empty_like(k_grid)
    for index, k_external in enumerate(k_grid):
        kx_external, ky_external = external_wave_vector(k_external, force_axis)
        k2_external = k_external * k_external
        q_dot_k = qx * kx_external + qy * ky_external
        numerator = q2 * k2_external - q_dot_k * q_dot_k
        px = kx_external - qx
        py = ky_external - qy
        p2 = px * px + py * py
        mask = (p2 > 1.0e-30) & (numerator > 0.0)
        denominator = p2[mask] * (nu0 * p2[mask] + d_q[mask] * q2[mask])
        integral = float(np.sum(numerator[mask] / denominator)) / volume
        updated[index] = d0 + kBT * integral / (rho0 * k2_external)
    return updated


def s46_update_grid_active_mask(
    k_grid,
    d_of_q,
    nx,
    ny,
    qx,
    qy,
    q2,
    q_abs,
    d0,
    nu0,
    rho0,
    kBT,
    length,
    grid,
    force_axis,
):
    lx, ly = length
    nx_active, ny_active = grid
    volume = lx * ly
    d_q = np.interp(q_abs, k_grid, d_of_q, left=d_of_q[0], right=d_of_q[-1])
    updated = np.empty_like(k_grid)
    for index, k_external in enumerate(k_grid):
        kx_external, ky_external = external_wave_vector(k_external, force_axis)
        nx_external = kx_external * lx / (2.0 * np.pi)
        ny_external = ky_external * ly / (2.0 * np.pi)
        k2_external = k_external * k_external
        q_dot_k = qx * kx_external + qy * ky_external
        numerator = q2 * k2_external - q_dot_k * q_dot_k
        px = kx_external - qx
        py = ky_external - qy
        p2 = px * px + py * py
        px_index = nx_external - nx
        py_index = ny_external - ny
        p_active = (np.abs(px_index) < nx_active / 2.0) & (np.abs(py_index) < ny_active / 2.0)
        mask = (p2 > 1.0e-30) & (numerator > 0.0) & p_active
        denominator = p2[mask] * (nu0 * p2[mask] + d_q[mask] * q2[mask])
        integral = float(np.sum(numerator[mask] / denominator)) / volume
        updated[index] = d0 + kBT * integral / (rho0 * k2_external)
    return updated


def discrete_s46_active_dr(
    q_values,
    d0,
    nu0,
    rho0,
    kBT,
    length,
    grid,
    force_axis,
    radial_count=384,
    max_iter=300,
    tolerance=1.0e-8,
    damping=0.45,
):
    nx_grid, ny_grid, qx_grid, qy_grid = active_wave_vectors(length, grid)
    q2_grid = qx_grid * qx_grid + qy_grid * qy_grid
    mask = q2_grid > 0.0
    nx = nx_grid[mask]
    ny = ny_grid[mask]
    qx = qx_grid[mask]
    qy = qy_grid[mask]
    q2 = q2_grid[mask]
    q_abs = np.sqrt(q2)
    k_grid = s46_radial_grid(q_values, q_abs, radial_count)
    d_current = np.full_like(k_grid, d0, dtype=float)

    last_error = np.inf
    iterations = 0
    for iteration in range(1, max_iter + 1):
        d_target = s46_update_grid_active_mask(
            k_grid,
            d_current,
            nx,
            ny,
            qx,
            qy,
            q2,
            q_abs,
            d0,
            nu0,
            rho0,
            kBT,
            length,
            grid,
            force_axis,
        )
        d_next = (1.0 - damping) * d_current + damping * d_target
        scale = np.maximum(np.abs(d_next), 1.0e-12)
        last_error = float(np.max(np.abs(d_next - d_current) / scale))
        d_current = d_next
        iterations = iteration
        if last_error < tolerance:
            break

    values = np.interp(np.asarray(q_values, dtype=float), k_grid, d_current)
    return values, iterations, last_error


def read_steady_profiles(run_name):
    path = PROCESSED_DATA_DIR / run_name / "steady_profiles.csv"
    by_nk = {}
    with path.open(newline="") as handle:
        for record in csv.DictReader(handle):
            nk = int(record["nk"])
            by_nk.setdefault(nk, []).append(record)
    for nk, records in by_nk.items():
        records.sort(key=lambda record: float(record["coord"]))
    return by_nk


def profile_modes(record_rows):
    coord = np.asarray([float(record["coord"]) for record in record_rows])
    psi = np.asarray([float(record["psi[0]_mean"]) for record in record_rows])
    psi_sem = np.asarray([float(record["psi[0]_sem"]) for record in record_rows])
    k = float(record_rows[0]["k"])
    sin_ky = np.sin(k * coord)
    cos_ky = np.cos(k * coord)
    psi_sin = 2.0 * float(np.dot(psi, sin_ky)) / coord.size
    psi_cos = 2.0 * float(np.dot(psi, cos_ky)) / coord.size
    fit = psi_sin * sin_ky + psi_cos * cos_ky
    residual = psi - fit
    fit_rms = float(np.sqrt(np.mean(fit * fit)))
    residual_rms = float(np.sqrt(np.mean(residual * residual)))
    return {
        "coord": coord,
        "psi": psi,
        "psi_sem": psi_sem,
        "k": k,
        "psi_sin": psi_sin,
        "psi_cos": psi_cos,
        "fit": fit,
        "residual": residual,
        "fit_rms": fit_rms,
        "residual_rms": residual_rms,
        "mode_amplitude": float(np.hypot(psi_sin, psi_cos)),
    }


def build_rows(run_name, steady_start):
    processed_dir = PROCESSED_DATA_DIR / run_name
    sample_path = processed_dir / "sample_steady_modes.csv"
    series_path = processed_dir / "mode_time_series.csv"
    if not sample_path.exists() or not series_path.exists():
        raise RuntimeError("run the profile processing first; sample_steady_modes.csv and mode_time_series.csv are required")

    sample_modes = read_sample_modes(sample_path)
    time_series = read_mode_time_series(series_path)
    rows = []
    for nk in sorted(sample_modes):
        config = parse_config(RAW_DATA_DIR / run_name / f"nk_{nk:03d}" / "config.dat")
        length = (float(config["length"][0]), float(config["length"][1]))
        grid = (int(config["grid"][0]), int(config["grid"][1]))
        k = 2.0 * math.pi * nk / length[1]
        free_energy_a = float(config["free_energy_a"])
        force_amplitude = float(config["derived_force_amplitude"])
        psi_samples = np.asarray([float(record["psi[0]_sin"]) for record in sample_modes[nk]], dtype=float)
        d_eff_samples = force_amplitude / (free_energy_a * k * k * psi_samples)
        row = {
            "nk": nk,
            "k": k,
            "amplitude": force_amplitude,
            "psi_steady_mean": float(np.mean(psi_samples)),
            "psi_steady_sem": stderr(psi_samples),
            "d_eff_mean": float(np.mean(d_eff_samples)),
            "d_eff_sem": stderr(d_eff_samples),
            "n_samples": len(psi_samples),
            "d0": float(config["d0"]),
            "nu0": float(config["eta"]) / float(config["density"]),
            "rho0": float(config["density"]),
            "kBT": float(config["kBT"]),
            "length": length,
            "grid": grid,
            "force_axis": config["force_axis"],
            "steady_start": steady_start,
        }
        rows.append(row)
    return rows, time_series


def save_response(run_name, rows, steady_start, s46_radial_count, s46_max_iter, s46_tolerance, s46_damping):
    processed_dir = PROCESSED_DATA_DIR / run_name
    output = processed_dir / f"theory_compare_steady_response_t{time_label(steady_start)}_s46_active.csv"
    with output.open("w", newline="") as handle:
        fieldnames = [
            "steady_start",
            "nk",
            "k",
            "amplitude",
            "psi_steady_mean",
            "psi_steady_sem",
            "d_eff_mean",
            "d_eff_sem",
            "d_s46_active",
            "s46_iterations",
            "s46_error",
            "n_samples",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        q_values = np.asarray([row["k"] for row in rows])
        first = rows[0]
        s46_active, s46_iterations, s46_error = discrete_s46_active_dr(
            q_values,
            first["d0"],
            first["nu0"],
            first["rho0"],
            first["kBT"],
            first["length"],
            first["grid"],
            first["force_axis"],
            radial_count=s46_radial_count,
            max_iter=s46_max_iter,
            tolerance=s46_tolerance,
            damping=s46_damping,
        )
        for row, d_s46_active in zip(rows, s46_active):
            writer.writerow(
                {
                    "steady_start": format_float(steady_start),
                    "nk": row["nk"],
                    "k": format_float(row["k"]),
                    "amplitude": format_float(row["amplitude"]),
                    "psi_steady_mean": format_float(row["psi_steady_mean"]),
                    "psi_steady_sem": format_float(row["psi_steady_sem"]),
                    "d_eff_mean": format_float(row["d_eff_mean"]),
                    "d_eff_sem": format_float(row["d_eff_sem"]),
                    "d_s46_active": format_float(d_s46_active),
                    "s46_iterations": s46_iterations,
                    "s46_error": format_float(s46_error),
                    "n_samples": row["n_samples"],
                }
            )
    theory_values = {
        "d_s46_active": s46_active,
        "s46_iterations": s46_iterations,
        "s46_error": s46_error,
    }
    return output, theory_values


def draw_relaxation(ax, rows, time_series, selected_nks, steady_start):
    colors = plt.cm.viridis(np.linspace(0.12, 0.86, len(selected_nks)))
    by_nk = {row["nk"]: row for row in rows}
    max_time = max(float(records[-1]["time"]) for records in time_series.values())
    ax.axhline(1.0, color="#111827", lw=1.0, ls="--", alpha=0.75)
    ax.axvspan(steady_start / 1.0e6, max_time / 1.0e6, color="#dbeafe", alpha=0.45, lw=0)
    for index, (nk, color) in enumerate(zip(selected_nks, colors)):
        if nk not in by_nk or nk not in time_series:
            continue
        records = time_series[nk]
        scale = by_nk[nk]["psi_steady_mean"]
        times = np.asarray([float(record["time"]) for record in records]) / 1.0e6
        mean = np.asarray([float(record["psi[0]_sin_mean"]) for record in records]) / scale
        sem = np.asarray([float(record["psi[0]_sin_sem"]) for record in records]) / abs(scale)
        zorder = 20 - index
        ax.plot(times, mean, color=color, lw=1.45, label=rf"$n_k={nk}$", zorder=zorder)
        ax.fill_between(times, mean - sem, mean + sem, color=color, alpha=0.13, lw=0, zorder=zorder - 0.5)
    ax.set_title("Order-parameter mode relaxation")
    ax.set_xlabel(r"$t/10^6$")
    ax.set_ylabel(r"$\Psi_k(t)/\Psi_k^{\rm steady}$")
    ax.set_ylim(0.0, 1.25)
    ax.grid(True, alpha=0.22)
    ax.legend(frameon=False, fontsize=8.4, ncol=3)


def draw_dataset(ax, rows):
    k = np.asarray([row["k"] for row in rows])
    d_eff = np.asarray([row["d_eff_mean"] for row in rows])
    d_sem = np.asarray([row["d_eff_sem"] for row in rows])
    ax.errorbar(k, d_eff, yerr=d_sem, fmt="o", ms=4.3, capsize=2.3, color="#2563eb", label="simulation", alpha=0.92)


def save_profile_metrics(run_name, rows):
    profiles = read_steady_profiles(run_name)
    output = PROCESSED_DATA_DIR / run_name / "profile_linearity_metrics.csv"
    with output.open("w", newline="") as handle:
        fieldnames = [
            "nk",
            "k",
            "target_amplitude",
            "psi_sin",
            "psi_cos",
            "mode_amplitude",
            "mode_amplitude_over_target",
            "fit_rms",
            "residual_rms",
            "residual_over_fit",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            nk = row["nk"]
            metrics = profile_modes(profiles[nk])
            config = parse_config(RAW_DATA_DIR / run_name / f"nk_{nk:03d}" / "config.dat")
            target = float(config["target_amplitude"])
            residual_over_fit = metrics["residual_rms"] / metrics["fit_rms"] if metrics["fit_rms"] > 0.0 else np.nan
            writer.writerow(
                {
                    "nk": nk,
                    "k": format_float(row["k"]),
                    "target_amplitude": format_float(target),
                    "psi_sin": format_float(metrics["psi_sin"]),
                    "psi_cos": format_float(metrics["psi_cos"]),
                    "mode_amplitude": format_float(metrics["mode_amplitude"]),
                    "mode_amplitude_over_target": format_float(metrics["mode_amplitude"] / target),
                    "fit_rms": format_float(metrics["fit_rms"]),
                    "residual_rms": format_float(metrics["residual_rms"]),
                    "residual_over_fit": format_float(residual_over_fit),
                }
            )
    return output


def make_profile_figure(run_name, rows, selected_nks, output):
    profiles = read_steady_profiles(run_name)
    selected = [nk for nk in selected_nks if nk in profiles]
    if not selected:
        raise RuntimeError("no selected nk values are available in steady_profiles.csv")

    target_config = parse_config(RAW_DATA_DIR / run_name / f"nk_{selected[0]:03d}" / "config.dat")
    target_amplitude = float(target_config["target_amplitude"])
    by_nk = {row["nk"]: row for row in rows}
    colors = plt.cm.viridis(np.linspace(0.12, 0.86, len(selected)))
    fig = plt.figure(figsize=(12.0, max(7.2, 1.18 * len(selected))), constrained_layout=True)
    gridspec = fig.add_gridspec(len(selected), 2, width_ratios=[3.0, 1.28])
    left_axes = []
    first_left = None
    for index, (nk, color) in enumerate(zip(selected, colors)):
        ax = fig.add_subplot(gridspec[index, 0], sharex=first_left)
        if first_left is None:
            first_left = ax
        left_axes.append(ax)
        metrics = profile_modes(profiles[nk])
        coord = metrics["coord"]
        target = target_amplitude * np.sin(metrics["k"] * coord)
        ax.plot(coord, target, color="#6b7280", lw=1.0, ls="--", label="bare target" if index == 0 else None)
        ax.plot(coord, metrics["fit"], color="#111827", lw=1.0, ls=":", label="first harmonic" if index == 0 else None)
        ax.plot(coord, metrics["psi"], color=color, lw=1.55, label="simulation" if index == 0 else None)
        ax.fill_between(
            coord,
            metrics["psi"] - metrics["psi_sem"],
            metrics["psi"] + metrics["psi_sem"],
            color=color,
            alpha=0.16,
            lw=0,
        )
        ax.axhline(0.0, color="#d1d5db", lw=0.8)
        ax.set_ylabel(rf"$n_k={nk}$")
        ax.grid(True, alpha=0.18)
        if index + 1 < len(selected):
            ax.tick_params(labelbottom=False)
    left_axes[-1].set_xlabel(r"$y$")
    left_axes[0].legend(frameon=False, fontsize=8.0, ncol=3, loc="upper right")

    all_k = []
    amplitude_ratio = []
    residual_ratio = []
    for row in rows:
        metrics = profile_modes(profiles[row["nk"]])
        all_k.append(row["k"])
        amplitude_ratio.append(metrics["mode_amplitude"] / target_amplitude)
        residual_ratio.append(metrics["residual_rms"] / metrics["fit_rms"] if metrics["fit_rms"] > 0.0 else np.nan)

    split = max(1, len(selected) // 2)
    ax_amp = fig.add_subplot(gridspec[:split, 1])
    ax_res = fig.add_subplot(gridspec[split:, 1], sharex=ax_amp)
    ax_amp.plot(all_k, amplitude_ratio, "o-", color="#2563eb", ms=4.0, lw=1.25)
    ax_amp.set_xscale("log")
    ax_amp.set_ylabel(r"$|\Psi_1|/\Psi_{\rm target}$")
    ax_amp.set_title("Profile diagnostics")
    ax_amp.grid(True, which="both", alpha=0.22)

    ax_res.plot(all_k, residual_ratio, "o-", color="#dc2626", ms=4.0, lw=1.25)
    ax_res.set_xscale("log")
    ax_res.set_xlabel(r"$q = 2\pi n_k/L$")
    ax_res.set_ylabel("residual RMS / first harmonic RMS")
    ax_res.grid(True, which="both", alpha=0.22)

    fig.suptitle("Steady order-parameter profiles and first-harmonic residuals", y=1.01)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight", dpi=180)
    return output


def make_figure(run_name, rows, time_series, steady_start, a_uv, selected_nks, output, theory_values, response_y_scale):
    k_data = np.asarray([row["k"] for row in rows])
    d_eff = np.asarray([row["d_eff_mean"] for row in rows])
    d_eff_sem = np.asarray([row["d_eff_sem"] for row in rows])
    d_s46_active = theory_values["d_s46_active"]
    relative_difference = (d_eff - d_s46_active) / d_s46_active

    fig, axes = plt.subplot_mosaic(
        [["relax", "relax"], ["continuum", "square"]],
        figsize=(11.6, 8.0),
        constrained_layout=True,
    )
    draw_relaxation(axes["relax"], rows, time_series, selected_nks, steady_start)

    draw_dataset(axes["continuum"], rows)
    axes["continuum"].plot(k_data, d_s46_active, color="#111827", lw=1.9, label="Eq. (S46), active spectral mask")
    axes["continuum"].set_title("S46 with simulation active spectral mask")

    axes["square"].axhline(0.0, color="#111827", lw=1.0, alpha=0.75)
    axes["square"].errorbar(
        k_data,
        relative_difference,
        yerr=d_eff_sem / d_s46_active,
        fmt="o-",
        ms=4.2,
        lw=1.1,
        capsize=2.2,
        color="#2563eb",
        label="simulation - S46 active",
    )
    axes["square"].set_title("Relative deviation from S46 active")

    for ax in (axes["continuum"], axes["square"]):
        ax.set_xscale("log")
        ax.set_xlabel(r"$q = 2\pi n_k/L$")
        ax.grid(True, which="both", alpha=0.22)
        ax.legend(frameon=False, fontsize=8.7)

    axes["continuum"].set_yscale(response_y_scale)
    axes["continuum"].set_ylabel(r"$D_{\rm eff}(q)$")
    axes["square"].set_ylabel(r"$(D_{\rm sim}-D_{\rm S46})/D_{\rm S46}$")

    fig.suptitle(rf"Order-parameter response, steady average $t\geq {steady_start:g}$", y=1.02)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight", dpi=180)
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-name", default=DEFAULT_RUN_NAME)
    parser.add_argument("--steady-start", type=float, default=1.0e7)
    parser.add_argument("--a-uv", type=float)
    parser.add_argument("--relaxation-nks", nargs="*", type=int, default=[1, 2, 4, 8, 16, 32, 48])
    parser.add_argument("--s46-radial-count", type=int, default=384)
    parser.add_argument("--s46-max-iter", type=int, default=300)
    parser.add_argument("--s46-tolerance", type=float, default=1.0e-8)
    parser.add_argument("--s46-damping", type=float, default=0.45)
    parser.add_argument("--response-y-scale", choices=("linear", "log"), default="linear")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--profile-output", type=Path)
    args = parser.parse_args()

    rows, time_series = build_rows(args.run_name, args.steady_start)
    first = rows[0]
    if args.a_uv is None:
        args.a_uv = first["length"][1] / first["grid"][1]
    if args.output is None:
        args.output = FIGURE_DIR / args.run_name / f"order_parameter_s46_active_theory_compare_y{args.response_y_scale}.png"
    if args.profile_output is None:
        args.profile_output = FIGURE_DIR / args.run_name / "order_parameter_profile_linearity_check.png"

    response_path, theory_values = save_response(
        args.run_name,
        rows,
        args.steady_start,
        args.s46_radial_count,
        args.s46_max_iter,
        args.s46_tolerance,
        args.s46_damping,
    )
    output = make_figure(
        args.run_name,
        rows,
        time_series,
        args.steady_start,
        args.a_uv,
        args.relaxation_nks,
        args.output,
        theory_values,
        args.response_y_scale,
    )
    profile_metrics_path = save_profile_metrics(args.run_name, rows)
    profile_output = make_profile_figure(args.run_name, rows, args.relaxation_nks, args.profile_output)
    print(f"saved {response_path}")
    print(f"saved {output}")
    print(f"saved {profile_metrics_path}")
    print(f"saved {profile_output}")
    print(f"S46 iterations={theory_values['s46_iterations']}, error={theory_values['s46_error']:.3e}")
    print(f"a_uv={args.a_uv:g}, samples={rows[0]['n_samples']}, nk_count={len(rows)}")


if __name__ == "__main__":
    main()

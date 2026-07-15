#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import re
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw_data"
PROCESSED = ROOT / "processed_data"
FIGURES = ROOT / "figures"

CASES = [
    ("relax_D0_4p00_grid512_L16384_dt4", "4p00", 1, 4, 7),
    ("relax_D0_0p12_grid512_L16384_dt4", "0p12", 2, 5, 8),
    ("relax_D0_0p004_grid512_L16384_dt16", "0p004", 3, 6, 9),
]

COLORS = {
    "4p00": "#7c3aed",
    "0p12": "#0891b2",
    "0p004": "#dc2626",
    "E_K": "#0f172a",
    "psi2": "#b45309",
    "dissipation": "#047857",
    "production": "#be123c",
    "length": "#2563eb",
    "transfer": "#2563eb",
    "total": "#0f172a",
    "fit": "#64748b",
}

PAPER_FIGSIZE = (9.6, 8.0)
MAIN_LINEWIDTH = 3.0
GUIDE_LINEWIDTH = 2.5

plt.rcParams.update(
    {
        "font.size": 20,
        "axes.titlesize": 20,
        "axes.labelsize": 24,
        "xtick.labelsize": 20,
        "ytick.labelsize": 20,
        "legend.fontsize": 18,
        "figure.titlesize": 20,
        "axes.linewidth": 1.3,
        "xtick.major.size": 7,
        "ytick.major.size": 7,
        "xtick.major.width": 1.3,
        "ytick.major.width": 1.3,
        "xtick.minor.size": 4,
        "ytick.minor.size": 4,
        "xtick.minor.width": 1.0,
        "ytick.minor.width": 1.0,
    }
)


def clean_outputs() -> None:
    FIGURES.mkdir(exist_ok=True)
    for path in FIGURES.glob("*.png"):
        path.unlink()

    PROCESSED.mkdir(exist_ok=True)
    for path in PROCESSED.iterdir():
        if path.name == ".gitkeep":
            continue
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def fmt(value: object) -> object:
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.16e}"
    return value


def label_to_float(label: str) -> float:
    return float(label.replace("p", "."))


def dt_label(dt: str) -> str:
    value = float(dt)
    return str(int(value)) if value.is_integer() else str(value).replace(".", "p")


def shell_delta_k(params: dict[str, str]) -> float:
    return min(2.0 * math.pi / float(params["length_x"]), 2.0 * math.pi / float(params["length_y"]))


def read_input(path: Path) -> dict[str, str]:
    params: dict[str, str] = {}
    with path.open() as handle:
        for raw_line in handle:
            parts = raw_line.strip().split()
            if not parts:
                continue
            if parts[0] == "grid":
                params["grid_x"], params["grid_y"] = parts[1], parts[2]
            elif parts[0] == "length":
                params["length_x"], params["length_y"] = parts[1], parts[2]
            elif parts[0] == "timestep":
                params["dt"] = parts[1]
            elif parts[:3] == ["model", "transport", "constant"]:
                transport = dict(zip(parts[3::2], parts[4::2]))
                params["eta"] = transport["eta"]
                params["d0"] = transport.get("M[0,0]", parts[-1])
            elif parts[:3] == ["fix", "3", "order_parameter"]:
                params["gradient"] = parts[-1]
    return params


def read_table(path: Path, columns: int) -> np.ndarray:
    rows: list[list[float]] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            values = [float(value) for value in line.split()]
            if len(values) != columns:
                raise RuntimeError(f"unexpected column count in {path}: {line}")
            rows.append(values)
    if not rows:
        raise RuntimeError(f"no numeric rows in {path}")
    return np.asarray(rows, dtype=float)


def read_last_block(path: Path, columns: int) -> np.ndarray:
    current: list[list[float]] = []
    seen_block = False
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("# block"):
                current = []
                seen_block = True
                continue
            if line.startswith("#"):
                continue
            values = [float(value) for value in line.split()]
            if len(values) != columns:
                raise RuntimeError(f"unexpected column count in {path}: {line}")
            current.append(values)
    if not seen_block or not current:
        raise RuntimeError(f"no data block in {path}")
    return np.asarray(current, dtype=float)


def mean_sem(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = values.mean(axis=0)
    sem = values.std(axis=0, ddof=1) / math.sqrt(values.shape[0]) if values.shape[0] > 1 else np.zeros_like(mean)
    return mean, sem


def same_shell_grid(blocks: list[np.ndarray]) -> None:
    reference = blocks[0][:, :2]
    for block in blocks[1:]:
        if block.shape != blocks[0].shape or not np.allclose(block[:, :2], reference, rtol=0.0, atol=0.0):
            raise RuntimeError("shell grids do not match across replicas")


def load_relaxation(case_name: str) -> dict[str, object]:
    case_root = RAW / case_name
    replicas = []
    segment_bounds: list[float] = []

    for replica_dir in sorted(case_root.glob("replica_*")):
        chunks = []
        bounds = []
        for path in sorted((replica_dir / "results").glob("time_series_relax_*.dat")):
            data = read_table(path, 6)
            chunks.append(data)
            bounds.append(float(data[-1, 1]))
        if chunks:
            replicas.append(np.vstack(chunks))
            segment_bounds = bounds

    if not replicas:
        raise RuntimeError(f"no relaxation data under {case_root}")

    time = replicas[0][:, 1]
    for data in replicas[1:]:
        if data.shape != replicas[0].shape or not np.allclose(data[:, 1], time):
            raise RuntimeError(f"relaxation time grids differ in {case_name}")

    first_input = sorted((case_root / "replica_000" / "runs").glob("input_relax_*.script"))[0]
    params = read_input(first_input)
    stack = np.stack(replicas, axis=0)
    d0 = float(params["d0"])
    gradient = float(params["gradient"])
    volume = float(params["length_x"]) * float(params["length_y"])
    metrics = {
        "E_K": stack[:, :, 2],
        "psi2": stack[:, :, 3],
        "dissipation": d0 * stack[:, :, 4] / volume,
        "production": -gradient * stack[:, :, 5],
        "length": 2.0 * math.pi * np.sqrt(stack[:, :, 3] / stack[:, :, 4]),
    }
    return {
        "params": params,
        "time": time,
        "metrics": metrics,
        "segment_bounds": segment_bounds,
        "n_replica": len(replicas),
    }


def plot_relaxation(case_name: str, label: str, number: int) -> None:
    data = load_relaxation(case_name)
    params = data["params"]
    time = data["time"] / 1.0e8
    metrics = data["metrics"]

    rows = []
    means: dict[str, np.ndarray] = {}
    sems: dict[str, np.ndarray] = {}
    for key, values in metrics.items():
        means[key], sems[key] = mean_sem(values)
    for index, t_value in enumerate(data["time"]):
        rows.append([fmt(t_value), *[fmt(means[key][index]) for key in metrics], *[fmt(sems[key][index]) for key in metrics]])
    write_csv(
        PROCESSED / case_name / "relaxation_mean.csv",
        ["time", *[f"{key}_mean" for key in metrics], *[f"{key}_sem" for key in metrics]],
        rows,
    )

    panels = [
        ("E_K", r"$E_K$"),
        ("psi2", r"$\langle|\psi|^2\rangle$"),
        ("dissipation", r"$D_0\langle|\nabla\psi|^2\rangle/V$"),
        ("production", r"$-G\langle J_{\psi,x}\rangle$"),
        ("length", r"$\ell_\psi$"),
    ]
    fig, axes = plt.subplots(3, 2, figsize=(12.0, 10.0), constrained_layout=True)
    axes = axes.ravel()
    for ax, (key, ylabel) in zip(axes, panels):
        color = COLORS[key]
        for replica_values in metrics[key]:
            ax.plot(time, replica_values, color=color, alpha=0.12, linewidth=0.8)
        ax.plot(time, means[key], color=color, linewidth=MAIN_LINEWIDTH)
        ax.fill_between(time, means[key] - sems[key], means[key] + sems[key], color=color, alpha=0.16, linewidth=0)
        for bound in data["segment_bounds"][:-1]:
            ax.axvline(bound / 1.0e8, color="#64748b", linestyle=":", linewidth=1.0)
        ax.set_xlabel(r"$t/10^8$")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.25)
    axes[-1].axis("off")
    axes[-1].legend(
        [
            plt.Line2D([0], [0], color=COLORS["E_K"], linewidth=MAIN_LINEWIDTH),
            plt.Line2D([0], [0], color="#64748b", linestyle=":", linewidth=1.0),
        ],
        ["replica mean", "segment boundary"],
        loc="center",
        frameon=False,
    )
    fig.suptitle(
        rf"Relaxation, $D_0=\eta={float(params['d0']):g}$, $n={data['n_replica']}$ replicas",
    )
    path = FIGURES / (
        f"{number:02d}_kugui_relaxation_D0_{label}_grid{params['grid_x']}_L{int(float(params['length_x']))}"
        f"_dt{dt_label(params['dt'])}_relaxation.png"
    )
    fig.savefig(path, dpi=200)
    plt.close(fig)


def load_shell_blocks(case_name: str, pattern: str, columns: int) -> tuple[dict[str, str], np.ndarray, np.ndarray]:
    case_root = RAW / case_name
    paths = sorted(case_root.glob(f"replica_*/results/{pattern}"))
    if not paths:
        raise RuntimeError(f"no shell files for {case_name}: {pattern}")
    blocks = [read_last_block(path, columns) for path in paths]
    same_shell_grid(blocks)
    first_path = paths[0]
    match = re.search(r"from_(\d+)\.dat$", first_path.name)
    if match is None:
        raise RuntimeError(f"cannot parse restart index from {first_path}")
    mode = "budget" if "budget" in first_path.name else "structure"
    input_path = first_path.parents[1] / "runs" / f"input_{mode}_from_{int(match.group(1)):03d}.script"
    params = read_input(input_path)
    return params, blocks[0][:, :2], np.stack(blocks, axis=0)


def symlog_linthresh(values: np.ndarray) -> float:
    nonzero = np.abs(values[np.abs(values) > 0.0])
    return 1.0e-16 if nonzero.size == 0 else 10.0 ** (math.floor(math.log10(float(nonzero.max()))) - 5)


def budget_figure_path(params: dict[str, str], label: str, number: int, suffix: str) -> Path:
    return FIGURES / (
        f"{number:02d}_kugui_budget_D0_{label}_grid{params['grid_x']}_L{int(float(params['length_x']))}"
        f"_dt{dt_label(params['dt'])}_shell_spectrum{suffix}.png"
    )


def draw_budget_figure(
    path: Path,
    k: np.ndarray,
    means: dict[str, np.ndarray],
    sems: dict[str, np.ndarray],
    params: dict[str, str],
    n_replica: int,
    *,
    xscale: str,
    yscale: str,
    ylabel: str = "budget spectrum",
    title_quantity: str = "budget shell spectrum",
    xlim: tuple[float, float] | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=PAPER_FIGSIZE, constrained_layout=True)
    ax.set_xscale(xscale)
    ax.axhline(0, color="#6b7280", linewidth=0.9)
    ax.grid(True, which="both", alpha=0.25)

    if yscale == "symlog":
        ax.set_yscale("symlog", linthresh=symlog_linthresh(np.r_[means["transfer"], means["dissipation"], means["production"]]))

    for key in ("transfer", "dissipation", "production"):
        ax.plot(k, means[key], color=COLORS[key], linewidth=MAIN_LINEWIDTH, label=key)
        ax.fill_between(k, means[key] - sems[key], means[key] + sems[key], color=COLORS[key], alpha=0.13, linewidth=0)
    if xlim is not None:
        ax.set_xlim(*xlim)
        zoom = (k >= xlim[0]) & (k <= xlim[1])
        main_values = np.r_[
            means["transfer"][zoom] - sems["transfer"][zoom],
            means["transfer"][zoom] + sems["transfer"][zoom],
            means["dissipation"][zoom] - sems["dissipation"][zoom],
            means["dissipation"][zoom] + sems["dissipation"][zoom],
            means["production"][zoom] - sems["production"][zoom],
            means["production"][zoom] + sems["production"][zoom],
        ]
        span = float(np.nanmax(main_values) - np.nanmin(main_values))
        pad = 0.08 * span if span > 0.0 else 1.0
        ax.set_ylim(float(np.nanmin(main_values)) - pad, float(np.nanmax(main_values)) + pad)
    ax.set_ylabel(ylabel)
    if yscale == "linear":
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    ax.set_xlabel(r"$k$")
    ax.legend(frameon=False, ncol=3, loc="upper right")
    title_map = {
        "budget shell spectrum": "budget",
        "premultiplied budget shell spectrum": "premultiplied budget",
        "premultiplied budget shell spectrum, low-$k$ zoom": "premultiplied budget, low-$k$ zoom",
    }
    ax.set_title(rf"{title_map.get(title_quantity, title_quantity)}, $D_0={float(params['d0']):g}$, $n={n_replica}$ replicas")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def low_k_zoom_xlim(k: np.ndarray, transfer: np.ndarray) -> tuple[float, float]:
    negative = transfer < 0.0
    positive = transfer > 0.0
    crossings = np.where(negative[:-1] & positive[1:])[0]
    if crossings.size:
        upper = float(k[int(crossings[0]) + 1])
    else:
        log_min = math.log(float(k.min()))
        log_max = math.log(float(k.max()))
        upper = math.exp(log_min + 0.55 * (log_max - log_min))
    return float(k.min()), upper


def plot_budget(case_name: str, label: str, number: int) -> None:
    params, shells, stack = load_shell_blocks(case_name, "budget_shell_from_*.dat", 6)
    k = shells[:, 0]
    terms = {
        "transfer": stack[:, :, 2],
        "dissipation": stack[:, :, 3],
        "production": stack[:, :, 4],
        "total": stack[:, :, 5],
    }
    means = {}
    sems = {}
    for key, values in terms.items():
        means[key], sems[key] = mean_sem(values)
    write_csv(
        PROCESSED / case_name / "budget_shell_mean.csv",
        ["k", "count", *[f"{key}_mean" for key in terms], *[f"{key}_sem" for key in terms]],
        [
            [fmt(shells[i, 0]), fmt(shells[i, 1]), *[fmt(means[key][i]) for key in terms], *[fmt(sems[key][i]) for key in terms]]
            for i in range(k.size)
        ],
    )

    draw_budget_figure(
        budget_figure_path(params, label, number, ""),
        k,
        means,
        sems,
        params,
        stack.shape[0],
        xscale="log",
        yscale="symlog",
    )

    log_density_factor = k / shell_delta_k(params)
    premultiplied_means = {key: log_density_factor * values for key, values in means.items()}
    premultiplied_sems = {key: log_density_factor * values for key, values in sems.items()}
    draw_budget_figure(
        budget_figure_path(params, label, number, "_premultiplied"),
        k,
        premultiplied_means,
        premultiplied_sems,
        params,
        stack.shape[0],
        xscale="log",
        yscale="linear",
        ylabel=r"$k\,B_{\rm shell}(k)/\Delta k$",
        title_quantity="premultiplied budget shell spectrum",
    )
    draw_budget_figure(
        budget_figure_path(params, label, number, "_premultiplied_zoom"),
        k,
        premultiplied_means,
        premultiplied_sems,
        params,
        stack.shape[0],
        xscale="log",
        yscale="linear",
        ylabel=r"$k\,B_{\rm shell}(k)/\Delta k$",
        title_quantity="premultiplied budget shell spectrum, low-$k$ zoom",
        xlim=low_k_zoom_xlim(k, means["transfer"]),
    )


def positive_band(mean: np.ndarray, sem: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return np.maximum(mean - sem, np.nextafter(0.0, 1.0)), np.maximum(mean + sem, np.nextafter(0.0, 1.0))


def guide_window(k: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float, float]:
    valid = np.isfinite(k) & np.isfinite(y) & (k > 0) & (y > 0)
    k_valid = k[valid]
    y_valid = y[valid]
    low = math.exp(math.log(float(k_valid.min())) + 0.10 * math.log(float(k_valid.max() / k_valid.min())))
    high = math.exp(math.log(float(k_valid.min())) + 0.50 * math.log(float(k_valid.max() / k_valid.min())))
    mask = (k_valid >= low) & (k_valid <= high)
    if np.count_nonzero(mask) < 3:
        mask = np.arange(k_valid.size) < min(12, k_valid.size)
        low = float(k_valid[mask].min())
        high = float(k_valid[mask].max())
    return valid, low, high


def add_power_guide(ax: plt.Axes, k: np.ndarray, y: np.ndarray, power: float, label: str, offset: float = 2.2) -> None:
    valid, low, high = guide_window(k, y)
    k_valid = k[valid]
    y_valid = y[valid]
    mask = (k_valid >= low) & (k_valid <= high)
    coefficient = offset * float(np.median(y_valid[mask] / k_valid[mask] ** power))
    guide_k = np.geomspace(low, high, 120)
    ax.loglog(guide_k, coefficient * guide_k**power, color="#111827", linestyle="--", linewidth=GUIDE_LINEWIDTH, label=label)


def add_power_fit(ax: plt.Axes, k: np.ndarray, y: np.ndarray) -> None:
    valid, low, high = guide_window(k, y)
    k_valid = k[valid]
    y_valid = y[valid]
    mask = (k_valid >= low) & (k_valid <= high)
    slope, intercept = np.polyfit(np.log(k_valid[mask]), np.log(y_valid[mask]), 1)
    fit_k = np.geomspace(low, high, 120)
    ax.loglog(
        fit_k,
        np.exp(intercept) * fit_k**slope,
        color=COLORS["fit"],
        linestyle="--",
        linewidth=GUIDE_LINEWIDTH,
        label=rf"fit $\propto k^{{{slope:.2f}}}$",
    )


def structure_length(case_name: str) -> tuple[float, float]:
    lengths = []
    case_root = RAW / case_name
    restart = sorted(case_root.glob("replica_*/results/static_corr_shell_from_*.dat"))[0].stem.rsplit("_", 1)[-1]
    for path in sorted(case_root.glob(f"replica_*/results/time_series_structure_from_{restart}.dat")):
        data = read_table(path, 6)
        lengths.append(math.sqrt(float(data[:, 3].mean()) / float(data[:, 4].mean())))
    values = np.asarray(lengths, dtype=float)
    sem = values.std(ddof=1) / math.sqrt(values.size) if values.size > 1 else 0.0
    return float(values.mean()), float(sem)


def plot_structure(case_name: str, label: str, number: int) -> dict[str, object]:
    params, shells, stack = load_shell_blocks(case_name, "static_corr_shell_from_*.dat", 3)
    k = shells[:, 0]
    count = shells[:, 1]
    mean, sem = mean_sem(stack[:, :, 2])
    length_mean, length_sem = structure_length(case_name)
    write_csv(
        PROCESSED / case_name / "structure_shell_mean.csv",
        ["k", "count", "Spsi_mean", "Spsi_sem"],
        [[fmt(k[i]), fmt(count[i]), fmt(mean[i]), fmt(sem[i])] for i in range(k.size)],
    )

    fig, ax = plt.subplots(figsize=PAPER_FIGSIZE, constrained_layout=True)
    lower, upper = positive_band(mean, sem)
    ax.loglog(k, mean, color=COLORS[label], linewidth=MAIN_LINEWIDTH, marker="o", markersize=7.5, markeredgewidth=0, label=rf"$S_\psi(k)$")
    ax.fill_between(k, lower, upper, color=COLORS[label], alpha=0.16, linewidth=0)
    add_power_guide(ax, k, mean, -4.0, r"$\propto k^{-4}$")
    if label == "0p004":
        add_power_fit(ax, k, mean)
    ax.set_title(
        rf"Structure factor, $D_0={float(params['d0']):g}$, $n={stack.shape[0]}$ replicas",
    )
    ax.set_xlabel(r"$k$")
    ax.set_ylabel(r"$S_\psi(k)$")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    path = FIGURES / (
        f"{number:02d}_kugui_structure_D0_{label}_grid{params['grid_x']}_L{int(float(params['length_x']))}"
        f"_dt{dt_label(params['dt'])}_shell_spectrum.png"
    )
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return {"case": case_name, "label": label, "params": params, "k": k, "mean": mean, "sem": sem, "length_mean": length_mean, "length_sem": length_sem}


def plot_structure_comparison(results: list[dict[str, object]]) -> None:
    fig, ax = plt.subplots(figsize=PAPER_FIGSIZE, constrained_layout=True)
    for result in sorted(results, key=lambda item: label_to_float(str(item["label"])), reverse=True):
        label = str(result["label"])
        k = result["k"]
        mean = result["mean"]
        sem = result["sem"]
        lower, upper = positive_band(mean, sem)
        ax.loglog(k, mean, color=COLORS[label], linewidth=MAIN_LINEWIDTH, marker="o", markersize=7.0, markeredgewidth=0, label=rf"$D_0={float(result['params']['d0']):g}$")
        ax.fill_between(k, lower, upper, color=COLORS[label], alpha=0.12, linewidth=0)
    guide = max(results, key=lambda item: float(np.nanmedian(item["mean"][:8])))
    add_power_guide(ax, guide["k"], guide["mean"], -4.0, r"linear theory $\propto k^{-4}$", offset=2.4)
    ax.set_title(r"Kugui structure factor comparison, $N=512$, $L=16384$")
    ax.set_xlabel(r"$k$")
    ax.set_ylabel(r"$S_\psi(k)$")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(FIGURES / "10_kugui_structure_comparison_grid512_L16384_shell_spectrum.png", dpi=200)
    plt.close(fig)


def plot_length_scale(results: list[dict[str, object]]) -> None:
    ordered = sorted(results, key=lambda item: float(item["params"]["d0"]))
    d0 = np.asarray([float(item["params"]["d0"]) for item in ordered])
    length = np.asarray([float(item["length_mean"]) for item in ordered])
    sem = np.asarray([float(item["length_sem"]) for item in ordered])
    slope, intercept = np.polyfit(np.log(d0), np.log(length), 1)
    write_csv(
        PROCESSED / "structure_length_scale_summary.csv",
        ["case", "D0", "length_mean", "length_sem"],
        [[item["case"], fmt(float(item["params"]["d0"])), fmt(float(item["length_mean"])), fmt(float(item["length_sem"]))] for item in ordered],
    )

    fig, ax = plt.subplots(figsize=PAPER_FIGSIZE, constrained_layout=True)
    ax.errorbar(d0, length, yerr=sem, fmt="o", color="#0f172a", ecolor="#64748b", capsize=4, markersize=8, label=r"time-series ratio")
    fit_x = np.geomspace(float(d0.min()), float(d0.max()), 200)
    ax.loglog(fit_x, np.exp(intercept) * fit_x**slope, color="#dc2626", linewidth=MAIN_LINEWIDTH, label=rf"fit $\alpha={slope:.2f}$")
    ax.set_title(r"Structure length scale")
    ax.set_xlabel(r"$D_0$")
    ax.set_ylabel(r"$\ell_\psi$")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(FIGURES / "11_kugui_structure_length_scale_vs_D0.png", dpi=200)
    plt.close(fig)


def main() -> None:
    clean_outputs()
    structure_results = []
    for case_name, label, relax_number, budget_number, structure_number in CASES:
        plot_relaxation(case_name, label, relax_number)
        plot_budget(case_name, label, budget_number)
        structure_results.append(plot_structure(case_name, label, structure_number))
    plot_structure_comparison(structure_results)
    plot_length_scale(structure_results)
    print(f"wrote {len(list(FIGURES.glob('*.png')))} figures and minimal processed CSV files")


if __name__ == "__main__":
    main()

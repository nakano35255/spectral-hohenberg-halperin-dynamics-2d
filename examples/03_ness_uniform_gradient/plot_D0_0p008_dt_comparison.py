#!/usr/bin/env python3
from pathlib import Path
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analysis_common import EXAMPLE, parse_input, read_time_series, transform, col


FIGURE_DIR = EXAMPLE / "figures"
CACHE_DIR = EXAMPLE / "cache"
CACHE_PATH = CACHE_DIR / "D0_0p008_dt_comparison_timeseries.npz"
OUTPUT = FIGURE_DIR / "05_D0_0p008_dt_comparison.png"
D0 = 0.008
OVERLAP_TMAX = 1.0e8


CASES = [
    {
        "key": "legacy_dt4",
        "label": r"legacy: $\Delta t=4,\ T=10^8,\ n=576$",
        "result_dir": EXAMPLE
        / "legacy"
        / "production_D0_0p008_dt4_T100000000_n576"
        / "results"
        / "D0_0.008"
        / "dt_4.0",
        "input_path": EXAMPLE
        / "legacy"
        / "production_D0_0p008_dt4_T100000000_n576"
        / "runs"
        / "D0_0.008"
        / "dt_4.0"
        / "input_000.script",
        "color": "#f97316",
        "linestyle": "--",
    },
    {
        "key": "main_dt8",
        "label": r"main: $\Delta t=8,\ T=2\times10^8,\ n=288$",
        "result_dir": EXAMPLE / "main" / "D0_0.008" / "results",
        "input_path": EXAMPLE / "main" / "D0_0.008" / "runs" / "input_000.script",
        "color": "#2563eb",
        "linestyle": "-",
    },
]

QUANTITIES = [
    (0, r"$E_K$"),
    (1, r"$\langle|\psi|^2\rangle$"),
    (2, r"$D_0\langle|\nabla\delta\psi|^2\rangle/V$"),
    (3, r"$-G\langle J_{\psi,x}\rangle$"),
]


def source_signature():
    entries = []
    for case in CASES:
        files = sorted(case["result_dir"].glob("time_series_*.dat"))
        stats = [path.stat() for path in files]
        input_stat = case["input_path"].stat()
        entries.append(
            {
                "key": case["key"],
                "count": len(files),
                "total_size": sum(stat.st_size for stat in stats),
                "max_mtime_ns": max([stat.st_mtime_ns for stat in stats] + [input_stat.st_mtime_ns]),
                "input_size": input_stat.st_size,
            }
        )
    return {"version": 1, "overlap_tmax": OVERLAP_TMAX, "entries": entries}


def load_case(case):
    params = parse_input(case["input_path"])
    files = sorted(case["result_dir"].glob("time_series_*.dat"))
    if not files:
        raise FileNotFoundError(f"no time_series files in {case['result_dir']}")

    times = None
    sums = None
    sums2 = None
    for path in files:
        columns, data = read_time_series(path)
        sample_times = data[:, col(columns, "time")]
        mask = sample_times <= OVERLAP_TMAX
        sample_times = sample_times[mask]
        data = data[mask]
        if times is None:
            times = sample_times
        elif data.shape[0] != times.shape[0] or not np.allclose(sample_times, times):
            raise ValueError(f"time grid differs in {path}")

        y = transform(data, columns, D0, params)
        if sums is None:
            sums = np.zeros_like(y)
            sums2 = np.zeros_like(y)
        sums += y
        sums2 += y * y

    samples = len(files)
    mean = sums / samples
    variance = np.maximum(sums2 / samples - mean * mean, 0.0)
    sem = np.sqrt(variance / max(samples - 1, 1))
    return {"times": times, "mean": mean, "sem": sem, "samples": samples}


def load_cache(signature):
    if not CACHE_PATH.exists():
        return None
    with np.load(CACHE_PATH) as data:
        if json.loads(str(data["signature_json"])) != signature:
            return None
        loaded = {}
        for case in CASES:
            key = case["key"]
            loaded[key] = {
                "times": data[f"{key}_times"],
                "mean": data[f"{key}_mean"],
                "sem": data[f"{key}_sem"],
                "samples": int(data[f"{key}_samples"]),
            }
        return loaded


def save_cache(results, signature):
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    payload = {"signature_json": np.asarray(json.dumps(signature, sort_keys=True))}
    for key, result in results.items():
        payload[f"{key}_times"] = result["times"]
        payload[f"{key}_mean"] = result["mean"]
        payload[f"{key}_sem"] = result["sem"]
        payload[f"{key}_samples"] = np.asarray(result["samples"], dtype=int)
    np.savez(CACHE_PATH, **payload)
    print(f"saved cache {CACHE_PATH}")


def load_results():
    signature = source_signature()
    cached = load_cache(signature)
    if cached is not None:
        print(f"loaded cache {CACHE_PATH}")
        return cached
    results = {case["key"]: load_case(case) for case in CASES}
    save_cache(results, signature)
    return results


def plot(results):
    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.2), constrained_layout=True)
    axes = axes.ravel()

    for ax, (q_index, ylabel) in zip(axes, QUANTITIES):
        for case in CASES:
            result = results[case["key"]]
            t = result["times"]
            y = result["mean"][:, q_index]
            yerr = result["sem"][:, q_index]
            label = case["label"]
            ax.plot(t, y, color=case["color"], ls=case["linestyle"], lw=1.5, label=label)
            ax.fill_between(t, y - yerr, y + yerr, color=case["color"], alpha=0.12, linewidth=0)

        ax.set_xlim(0.0, OVERLAP_TMAX)
        ax.set_title(ylabel)
        ax.set_xlabel(r"$t$")
        ax.grid(True, alpha=0.25, linewidth=0.7)
        ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
        if q_index in (2, 3):
            ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncols=2, frameon=False, bbox_to_anchor=(0.5, 1.05))
    fig.suptitle(r"$D_0=0.008$: $\Delta t=4$ legacy vs $\Delta t=8$ main, overlap window", y=1.11)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {OUTPUT}")


def main():
    plot(load_results())


if __name__ == "__main__":
    main()

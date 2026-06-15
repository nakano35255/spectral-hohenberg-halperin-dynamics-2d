#!/usr/bin/env python3
from pathlib import Path
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analysis_common import (
    EXAMPLE,
    col,
    mct_induced_values,
    numeric_label,
    parse_input,
    read_time_series,
    transform,
)


MAIN_ROOT = EXAMPLE / "main"
FIGURE_DIR = EXAMPLE / "figures"
CACHE_DIR = EXAMPLE / "cache"
CACHE_PATH = CACHE_DIR / "main_N128_auv32_steady_last_half.npz"
OUTPUT = FIGURE_DIR / "04_main_steady_D0_dependence.png"


def summarize_case(d0_dir):
    _, d0 = numeric_label(d0_dir, "D0")
    result_dir = d0_dir / "results"
    input_path = d0_dir / "runs" / "input_000.script"
    params = parse_input(input_path)
    files = sorted(result_dir.glob("time_series_*.dat"))
    if not files:
        raise FileNotFoundError(f"no time_series files in {result_dir}")

    sample_avgs = []
    sample_lengths = []
    steady_start = None
    steady_end = None
    for path in files:
        columns, data = read_time_series(path)
        times = data[:, col(columns, "time")]
        if steady_start is None:
            steady_start = 0.5 * times[-1]
            steady_end = times[-1]
        mask = times >= steady_start
        if not np.any(mask):
            raise ValueError(f"steady window has no data in {path}")

        y = transform(data, columns, d0, params)
        avg = y[mask].mean(axis=0)
        sample_avgs.append(avg)

        grad2 = avg[2] * params["volume"] / d0
        sample_lengths.append(2.0 * np.pi * np.sqrt(avg[1] / grad2))

    sample_avgs = np.asarray(sample_avgs)
    sample_lengths = np.asarray(sample_lengths)
    mean = sample_avgs.mean(axis=0)
    sem = sample_avgs.std(axis=0, ddof=1) / np.sqrt(sample_avgs.shape[0])
    return {
        "d0": d0,
        "d0_text": d0_dir.name.split("_", 1)[1],
        "params": params,
        "samples": sample_avgs.shape[0],
        "steady_start": steady_start,
        "steady_end": steady_end,
        "psi2": mean[1],
        "psi2_sem": sem[1],
        "induced": mean[2],
        "induced_sem": sem[2],
        "length": sample_lengths.mean(),
        "length_sem": sample_lengths.std(ddof=1) / np.sqrt(sample_lengths.size),
    }


def load_rows():
    signature = source_signature()
    cached = load_cache(signature)
    if cached is not None:
        print(f"loaded cache {CACHE_PATH}")
        return cached

    rows = [
        summarize_case(d0_dir)
        for d0_dir in sorted(MAIN_ROOT.glob("D0_*"), key=lambda path: numeric_label(path, "D0")[1])
    ]
    if not rows:
        raise FileNotFoundError(f"no D0 directories in {MAIN_ROOT}")
    save_cache(rows, signature)
    return rows


def source_signature():
    entries = []
    for d0_dir in sorted(MAIN_ROOT.glob("D0_*"), key=lambda path: numeric_label(path, "D0")[1]):
        files = sorted((d0_dir / "results").glob("time_series_*.dat"))
        stats = [path.stat() for path in files]
        input_path = d0_dir / "runs" / "input_000.script"
        input_stat = input_path.stat()
        entries.append(
            {
                "d0_dir": d0_dir.name,
                "count": len(files),
                "total_size": sum(stat.st_size for stat in stats),
                "max_mtime_ns": max([stat.st_mtime_ns for stat in stats] + [input_stat.st_mtime_ns]),
                "input_size": input_stat.st_size,
            }
        )
    return {"version": 1, "main_root": str(MAIN_ROOT), "entries": entries}


def serializable_params(params):
    result = {}
    for key, value in params.items():
        if isinstance(value, tuple):
            result[key] = list(value)
        else:
            result[key] = value
    return result


def restore_params(params):
    result = dict(params)
    for key in ("length", "grid"):
        if key in result:
            result[key] = tuple(result[key])
    return result


def load_cache(signature):
    if not CACHE_PATH.exists():
        return None
    with np.load(CACHE_PATH) as data:
        cached_signature = json.loads(str(data["signature_json"]))
        if cached_signature != signature:
            return None
        params = restore_params(json.loads(str(data["params_json"])))
        rows = []
        for i, d0 in enumerate(data["d0"]):
            rows.append(
                {
                    "d0": float(d0),
                    "d0_text": str(data["d0_text"][i]),
                    "params": params,
                    "samples": int(data["samples"][i]),
                    "steady_start": float(data["steady_start"][i]),
                    "steady_end": float(data["steady_end"][i]),
                    "psi2": float(data["psi2"][i]),
                    "psi2_sem": float(data["psi2_sem"][i]),
                    "induced": float(data["induced"][i]),
                    "induced_sem": float(data["induced_sem"][i]),
                    "length": float(data["length"][i]),
                    "length_sem": float(data["length_sem"][i]),
                }
            )
        return rows


def save_cache(rows, signature):
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    np.savez(
        CACHE_PATH,
        signature_json=np.asarray(json.dumps(signature, sort_keys=True)),
        params_json=np.asarray(json.dumps(serializable_params(rows[0]["params"]), sort_keys=True)),
        d0=np.asarray([row["d0"] for row in rows], dtype=float),
        d0_text=np.asarray([row["d0_text"] for row in rows]),
        samples=np.asarray([row["samples"] for row in rows], dtype=int),
        steady_start=np.asarray([row["steady_start"] for row in rows], dtype=float),
        steady_end=np.asarray([row["steady_end"] for row in rows], dtype=float),
        psi2=np.asarray([row["psi2"] for row in rows], dtype=float),
        psi2_sem=np.asarray([row["psi2_sem"] for row in rows], dtype=float),
        induced=np.asarray([row["induced"] for row in rows], dtype=float),
        induced_sem=np.asarray([row["induced_sem"] for row in rows], dtype=float),
        length=np.asarray([row["length"] for row in rows], dtype=float),
        length_sem=np.asarray([row["length_sem"] for row in rows], dtype=float),
    )
    print(f"saved cache {CACHE_PATH}")


def power_law_fit(x, y, mask):
    coeff = np.polyfit(np.log(x[mask]), np.log(y[mask]), 1)
    return coeff[0], np.exp(coeff[1])


def main():
    rows = load_rows()
    d0 = np.asarray([row["d0"] for row in rows])
    psi2 = np.asarray([row["psi2"] for row in rows])
    psi2_sem = np.asarray([row["psi2_sem"] for row in rows])
    induced = np.asarray([row["induced"] for row in rows])
    induced_sem = np.asarray([row["induced_sem"] for row in rows])
    length = np.asarray([row["length"] for row in rows])
    length_sem = np.asarray([row["length_sem"] for row in rows])

    params = rows[0]["params"]
    theory_d0 = np.exp(np.linspace(np.log(1.0e-3), np.log(max(d0.max() * 1.15, 4.5)), 700))
    mct = mct_induced_values(theory_d0, params)
    mct_zero = mct[0]

    fit_mask = d0 <= 0.5
    fit = power_law_fit(d0, length, fit_mask)
    fit_x = np.exp(np.linspace(np.log(d0[fit_mask].min()), np.log(d0[fit_mask].max()), 240))

    fig, axes = plt.subplots(1, 3, figsize=(15.8, 4.6), constrained_layout=True)

    ax = axes[0]
    ax.errorbar(d0, psi2, yerr=psi2_sem, fmt="o", color="#2563eb", capsize=3)
    ax.set_xscale("log")
    ax.set_xlabel(r"$D_0=\nu_0$")
    ax.set_ylabel(r"$\langle|\psi|^2\rangle$")
    ax.set_title(r"(a) Static fluctuation amplitude")
    ax.grid(True, which="both", alpha=0.25)

    ax = axes[1]
    ax.errorbar(d0, induced, yerr=induced_sem, fmt="o", color="#2563eb", capsize=3, label="simulation")
    ax.plot(theory_d0, mct, color="#111827", ls="--", lw=1.6, label="self-consistent MCT")
    ax.axhline(mct_zero, color="#64748b", ls=":", lw=1.3, label=r"$D_0\to0$")
    ax.set_xscale("log")
    ax.set_ylim(0.0, 1.10 * max(float(mct.max()), float((induced + induced_sem).max())))
    ax.set_xlabel(r"$D_0=\nu_0$")
    ax.set_ylabel(r"$D_0\langle|\nabla\delta\psi|^2\rangle/V$")
    ax.set_title(r"(b) Induced part, MCT zoom")
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, fontsize=8.3)

    ax = axes[2]
    ax.errorbar(d0, length, yerr=length_sem, fmt="o", color="#2563eb", capsize=3)
    ax.plot(fit_x, fit[1] * fit_x ** fit[0], color="#111827", ls="--", lw=1.3, label=rf"$D_0\leq0.5$: $\ell\propto D_0^{{{fit[0]:.2f}}}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$D_0=\nu_0$")
    ax.set_ylabel(r"$2\pi\sqrt{\langle|\psi|^2\rangle/\langle|\nabla\psi|^2\rangle}$")
    ax.set_title(r"(c) Characteristic length")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, fontsize=8.1)

    fig.suptitle(
        rf"$N={params['grid'][0]},\ a_{{uv}}={params['a_uv']:.0f},\ L={params['length'][0]:.0f},\ G={params['gradient']:.8g}$; "
        r"steady averages over the last half of each run",
        fontsize=12.5,
    )
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=240, bbox_inches="tight")
    plt.close(fig)

    print("D0 samples steady_start steady_end psi2 psi2_sem induced induced_sem length length_sem")
    for row in rows:
        print(
            f"{row['d0']:.6g} {row['samples']} {row['steady_start']:.8e} {row['steady_end']:.8e} "
            f"{row['psi2']:.8e} {row['psi2_sem']:.2e} "
            f"{row['induced']:.8e} {row['induced_sem']:.2e} "
            f"{row['length']:.8e} {row['length_sem']:.2e}"
        )
    print(f"log-log fit length D0<=0.5: exponent={fit[0]:.6f}, prefactor={fit[1]:.8e}")
    print(f"saved {OUTPUT}")


if __name__ == "__main__":
    main()

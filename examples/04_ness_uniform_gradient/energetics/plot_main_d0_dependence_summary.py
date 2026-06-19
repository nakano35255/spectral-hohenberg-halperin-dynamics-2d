#!/usr/bin/env python3
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analysis_common import (
    FIGURE_DIR,
    load_or_build_processed_energetics,
    mct_induced_values,
    mct_psi2_values,
)


SC_LABELS = ("Sc1", "Sc4")
OUTPUT = FIGURE_DIR / "04_main_steady_D0_dependence.png"
SC4_OUTPUT = FIGURE_DIR / "04_main_steady_D0_dependence_sc4.png"
COMPARISON_OUTPUT = FIGURE_DIR / "05_main_steady_D0_dependence_sc1_sc4.png"


def load_rows(sc_label):
    return load_or_build_processed_energetics(sc_label)


def power_law_fit(x, y, mask):
    coeff = np.polyfit(np.log(x[mask]), np.log(y[mask]), 1)
    return coeff[0], np.exp(coeff[1])


def schmidt_number(params):
    return params["eta"] / (params["density"] * params["mobility"])


def rows_to_arrays(rows):
    d0 = np.asarray([row["d0"] for row in rows])
    return {
        "d0": d0,
        "psi2": np.asarray([row["psi2"] for row in rows]),
        "psi2_sem": np.asarray([row["psi2_sem"] for row in rows]),
        "induced": np.asarray([row["induced"] for row in rows]),
        "induced_sem": np.asarray([row["induced_sem"] for row in rows]),
        "length": np.asarray([row["length"] for row in rows]),
        "length_sem": np.asarray([row["length_sem"] for row in rows]),
    }


def theory_for_rows(rows):
    d0 = np.asarray([row["d0"] for row in rows])
    params = rows[0]["params"]
    theory_d0 = np.exp(np.linspace(np.log(1.0e-3), np.log(max(d0.max() * 1.15, 4.5)), 700))
    return {
        "params": params,
        "d0": theory_d0,
        "induced": mct_induced_values(theory_d0, params),
        "psi2": mct_psi2_values(theory_d0, params),
    }


def length_fit(arrays):
    fit_mask = arrays["d0"] <= 0.5
    fit = power_law_fit(arrays["d0"], arrays["length"], fit_mask)
    fit_x = np.exp(np.linspace(np.log(arrays["d0"][fit_mask].min()), np.log(arrays["d0"][fit_mask].max()), 240))
    return fit, fit_x


def plot_single(rows, output, sc_label):
    arrays = rows_to_arrays(rows)
    theory = theory_for_rows(rows)
    params = theory["params"]
    fit, fit_x = length_fit(arrays)

    fig, axes = plt.subplots(1, 3, figsize=(15.8, 4.6), constrained_layout=True)

    ax = axes[0]
    ax.errorbar(arrays["d0"], arrays["psi2"], yerr=arrays["psi2_sem"], fmt="o", color="#2563eb", capsize=3, label="simulation")
    ax.plot(theory["d0"], theory["psi2"], color="#111827", ls="--", lw=1.6, label="self-consistent MCT")
    ax.set_xscale("log")
    ax.set_xlabel(r"$D_0$")
    ax.set_ylabel(r"$\langle|\psi|^2\rangle$")
    ax.set_title(r"(a) Static fluctuation amplitude")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, fontsize=8.3)

    ax = axes[1]
    ax.errorbar(arrays["d0"], arrays["induced"], yerr=arrays["induced_sem"], fmt="o", color="#2563eb", capsize=3, label="simulation")
    ax.plot(theory["d0"], theory["induced"], color="#111827", ls="--", lw=1.6, label="self-consistent MCT")
    ax.axhline(theory["induced"][0], color="#64748b", ls=":", lw=1.3, label=r"$D_0\to0$")
    ax.set_xscale("log")
    ax.set_ylim(0.0, 1.10 * max(float(theory["induced"].max()), float((arrays["induced"] + arrays["induced_sem"]).max())))
    ax.set_xlabel(r"$D_0$")
    ax.set_ylabel(r"$D_0\langle|\nabla\delta\psi|^2\rangle/V$")
    ax.set_title(r"(b) Induced part, MCT zoom")
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, fontsize=8.3)

    ax = axes[2]
    ax.errorbar(arrays["d0"], arrays["length"], yerr=arrays["length_sem"], fmt="o", color="#2563eb", capsize=3)
    ax.plot(fit_x, fit[1] * fit_x ** fit[0], color="#111827", ls="--", lw=1.3, label=rf"$D_0\leq0.5$: $\ell\propto D_0^{{{fit[0]:.2f}}}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$D_0$")
    ax.set_ylabel(r"$2\pi\sqrt{\langle|\psi|^2\rangle/\langle|\nabla\psi|^2\rangle}$")
    ax.set_title(r"(c) Characteristic length")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, fontsize=8.1)

    fig.suptitle(
        rf"$S_c={schmidt_number(params):.0f},\ N={params['grid'][0]},\ a_{{uv}}={params['a_uv']:.0f},\ L={params['length'][0]:.0f},\ G={params['gradient']:.8g}$; "
        r"steady averages over the last half of each run",
        fontsize=12.5,
    )
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(fig)
    return fit


def plot_comparison(rows_by_sc, output):
    styles = {
        "Sc1": {"color": "#2563eb", "marker": "o"},
        "Sc4": {"color": "#dc2626", "marker": "s"},
    }
    fig, axes = plt.subplots(1, 3, figsize=(16.6, 4.8), constrained_layout=True)
    induced_max = 0.0

    for sc_label, rows in rows_by_sc.items():
        arrays = rows_to_arrays(rows)
        theory = theory_for_rows(rows)
        params = theory["params"]
        fit, fit_x = length_fit(arrays)
        style = styles[sc_label]
        label = rf"$S_c={schmidt_number(params):.0f}$"

        ax = axes[0]
        ax.errorbar(
            arrays["d0"], arrays["psi2"], yerr=arrays["psi2_sem"],
            fmt=style["marker"] + "-", color=style["color"], capsize=3, lw=1.25,
            label=label + " simulation",
        )
        ax.plot(theory["d0"], theory["psi2"], color=style["color"], ls=":", lw=1.8, label=label + " MCT")

        ax = axes[1]
        ax.errorbar(
            arrays["d0"], arrays["induced"], yerr=arrays["induced_sem"],
            fmt=style["marker"] + "-", color=style["color"], capsize=3, lw=1.25,
            label=label + " simulation",
        )
        ax.plot(theory["d0"], theory["induced"], color=style["color"], ls=":", lw=1.8, label=label + " MCT")
        ax.axhline(theory["induced"][0], color=style["color"], ls=":", lw=1.1, alpha=0.75)
        induced_max = max(
            induced_max,
            float(theory["induced"].max()),
            float((arrays["induced"] + arrays["induced_sem"]).max()),
        )

        ax = axes[2]
        ax.errorbar(
            arrays["d0"], arrays["length"], yerr=arrays["length_sem"],
            fmt=style["marker"] + "-", color=style["color"], capsize=3, lw=1.25,
            label=label + " simulation",
        )
        ax.plot(
            fit_x,
            fit[1] * fit_x ** fit[0],
            color=style["color"],
            ls=":",
            lw=1.8,
            label=label + rf": $\ell\propto D_0^{{{fit[0]:.2f}}}$",
        )

    axes[0].set_xscale("log")
    axes[0].set_xlabel(r"$D_0$")
    axes[0].set_ylabel(r"$\langle|\psi|^2\rangle$")
    axes[0].set_title(r"(a) Static fluctuation amplitude")
    axes[0].grid(True, which="both", alpha=0.25)
    axes[0].legend(frameon=False, fontsize=7.7)

    axes[1].set_xscale("log")
    axes[1].set_ylim(0.0, 1.10 * induced_max)
    axes[1].set_xlabel(r"$D_0$")
    axes[1].set_ylabel(r"$D_0\langle|\nabla\delta\psi|^2\rangle/V$")
    axes[1].set_title(r"(b) Induced part")
    axes[1].ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    axes[1].grid(True, which="both", alpha=0.25)
    axes[1].legend(frameon=False, fontsize=7.7)

    axes[2].set_xscale("log")
    axes[2].set_yscale("log")
    axes[2].set_xlabel(r"$D_0$")
    axes[2].set_ylabel(r"$2\pi\sqrt{\langle|\psi|^2\rangle/\langle|\nabla\psi|^2\rangle}$")
    axes[2].set_title(r"(c) Characteristic length")
    axes[2].grid(True, which="both", alpha=0.25)
    axes[2].legend(frameon=False, fontsize=7.4)

    params = next(iter(rows_by_sc.values()))[0]["params"]
    fig.suptitle(
        rf"$N={params['grid'][0]},\ a_{{uv}}={params['a_uv']:.0f},\ L={params['length'][0]:.0f},\ G={params['gradient']:.8g}$; "
        r"Sc comparison, steady averages over the last half of each run",
        fontsize=12.5,
    )
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(fig)


def print_summary(sc_label, rows, fit):
    print(f"{sc_label}: D0 samples steady_start steady_end psi2 psi2_sem induced induced_sem length length_sem")
    for row in rows:
        print(
            f"{row['d0']:.6g} {row['samples']} {row['steady_start']:.8e} {row['steady_end']:.8e} "
            f"{row['psi2']:.8e} {row['psi2_sem']:.2e} "
            f"{row['induced']:.8e} {row['induced_sem']:.2e} "
            f"{row['length']:.8e} {row['length_sem']:.2e}"
        )
    print(f"{sc_label}: log-log fit length D0<=0.5: exponent={fit[0]:.6f}, prefactor={fit[1]:.8e}")


def main():
    rows_by_sc = {sc_label: load_rows(sc_label) for sc_label in SC_LABELS}
    outputs = {"Sc1": OUTPUT, "Sc4": SC4_OUTPUT}
    fits = {}
    for sc_label, rows in rows_by_sc.items():
        fits[sc_label] = plot_single(rows, outputs[sc_label], sc_label)
    plot_comparison(rows_by_sc, COMPARISON_OUTPUT)

    for sc_label in SC_LABELS:
        print_summary(sc_label, rows_by_sc[sc_label], fits[sc_label])
    for output in (OUTPUT, SC4_OUTPUT, COMPARISON_OUTPUT):
        print(f"saved {output}")


if __name__ == "__main__":
    main()

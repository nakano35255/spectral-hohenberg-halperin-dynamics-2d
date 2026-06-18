#!/usr/bin/env python3
from pathlib import Path
import re

import numpy as np


EXAMPLE = Path(__file__).resolve().parent


def numeric_label(path, prefix):
    match = re.fullmatch(rf"{re.escape(prefix)}_(.+)", path.name)
    if not match:
        raise ValueError(f"cannot parse {prefix} value from {path}")
    text = match.group(1)
    return text, float(text)


def parse_input(path):
    values = {"density": 1.0, "kBT": 1.0, "free_energy_a": 1.0}
    text = path.read_text()
    length = re.search(r"^length\s+(\S+)\s+(\S+)", text, re.M)
    grid = re.search(r"^grid\s+(\S+)\s+(\S+)", text, re.M)
    dealias = re.search(r"^dealias\s+(\S+)", text, re.M)
    transport = re.search(r"^model\s+transport\s+constant\b(.*)$", text, re.M)
    force = re.search(
        r"^fix\s+\S+\s+order_parameter\s+force/gradient\s+on\s+component\s+(\S+)\s+direction\s+(\S+)\s+amplitude\s+(\S+)",
        text,
        re.M,
    )
    density = re.search(r"^set\s+density\s+uniform\s+value\s+(\S+)", text, re.M)
    noise = re.search(r"^fix\s+\S+\s+momentum\s+noise\s+on\s+.*?\bkBT\s+(\S+)", text, re.M)
    free_energy = re.search(r"^model\s+free_energy\s+quadratic\s+a\[0\]\s+(\S+)", text, re.M)

    if length is None or grid is None or transport is None or force is None:
        raise ValueError(f"missing required parameters in {path}")

    tokens = transport.group(1).split()
    transport_values = {tokens[i]: tokens[i + 1] for i in range(0, len(tokens), 2)}
    if force.group(1) != "0":
        raise ValueError("this analysis assumes force/gradient component 0")

    values.update(
        {
            "length": (float(length.group(1)), float(length.group(2))),
            "grid": (int(grid.group(1)), int(grid.group(2))),
            "dealias": dealias.group(1) if dealias else "none",
            "eta": float(transport_values["eta"]),
            "mobility": float(transport_values["M[0,0]"]),
            "gradient_direction": force.group(2),
            "gradient": float(force.group(3)),
        }
    )
    if density:
        values["density"] = float(density.group(1))
    if noise:
        values["kBT"] = float(noise.group(1))
    if free_energy:
        values["free_energy_a"] = float(free_energy.group(1))
    values["chi"] = 1.0 / values["free_energy_a"]
    values["volume"] = values["length"][0] * values["length"][1]
    values["a_uv"] = values["length"][0] / values["grid"][0]
    return values


def read_time_series(path):
    columns = None
    with path.open() as handle:
        for line in handle:
            if line.startswith("# step time"):
                columns = line[2:].split()
                break
    if columns is None:
        raise ValueError(f"missing header in {path}")
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] != len(columns):
        raise ValueError(f"column count mismatch in {path}")
    return columns, data


def col(columns, name):
    try:
        return columns.index(name)
    except ValueError as exc:
        raise ValueError(f"missing column {name}") from exc


def transform(data, columns, d0, params):
    return np.column_stack(
        [
            data[:, col(columns, "E_K")],
            data[:, col(columns, "|psi[0]|^2")],
            d0 * data[:, col(columns, "|d_psi[0]|^2")] / params["volume"],
            -params["gradient"] * data[:, col(columns, "Jpsi[0]_x")],
        ]
    )


def compute_grid_size(active_size, dealias):
    if dealias in ("none", "off"):
        return active_size
    if dealias in ("three_halves", "3/2"):
        return 3 * active_size // 2
    if dealias in ("two", "2"):
        return 2 * active_size
    raise ValueError(f"unknown dealias rule: {dealias}")


def signed_index(index, n):
    return index if index <= n // 2 else index - n


def active_mode_arrays(params):
    active_nx, active_ny = params["grid"]
    compute_ny = compute_grid_size(active_ny, params["dealias"])
    lx, ly = params["length"]
    k_values = []
    angles = []
    weights = []
    for gx in range(0, active_nx // 2):
        kx = 2.0 * np.pi * gx / lx
        for gy_index in range(compute_ny):
            gy = signed_index(gy_index, compute_ny)
            if abs(gy) >= active_ny // 2:
                continue
            if gx == 0 and gy == 0:
                continue
            ky = 2.0 * np.pi * gy / ly
            k2 = kx * kx + ky * ky
            if params["gradient_direction"] == "x":
                angle = ky * ky / (k2 * k2)
            elif params["gradient_direction"] == "y":
                angle = kx * kx / (k2 * k2)
            else:
                raise ValueError(f"unknown gradient direction: {params['gradient_direction']}")
            k_values.append(np.sqrt(k2))
            angles.append(angle)
            weights.append(1.0 if gx == 0 else 2.0)
    return np.asarray(k_values), np.asarray(angles), np.asarray(weights)


def kinematic_viscosity_for_d0(d0_values, params):
    mobility = params["mobility"]
    if mobility == 0.0:
        raise ValueError("MCT comparison requires nonzero mobility.")
    schmidt = params["eta"] / (params["density"] * mobility)
    return schmidt * np.asarray(d0_values, dtype=float)


def mct_renormalized_diffusion(d0_values, k, params):
    d0_values = np.asarray(d0_values, dtype=float)
    nu0 = kinematic_viscosity_for_d0(d0_values, params)
    cutoff = 2.0 * np.pi / params["a_uv"]
    delta = params["kBT"] * np.log(cutoff / k) / (4.0 * np.pi * params["density"])
    d0_grid = d0_values[:, None]
    nu0_grid = nu0[:, None]
    return 0.5 * (d0_grid - nu0_grid) + np.sqrt(0.25 * (d0_grid + nu0_grid) ** 2 + delta[None, :])


def mct_induced(d0, params):
    k, angle, weight = active_mode_arrays(params)
    dr = mct_renormalized_diffusion(np.asarray([d0]), k, params)[0]
    nu0 = kinematic_viscosity_for_d0(np.asarray([d0]), params)[0]
    total = np.sum(weight * angle / (nu0 + dr))
    return params["kBT"] * params["gradient"] ** 2 * total / (
        params["density"] * params["chi"] * params["volume"]
    )


def mct_induced_values(d0_values, params):
    d0_values = np.asarray(d0_values, dtype=float)
    k, angle, weight = active_mode_arrays(params)
    dr = mct_renormalized_diffusion(d0_values, k, params)
    nu0 = kinematic_viscosity_for_d0(d0_values, params)
    total = np.sum(weight[None, :] * angle[None, :] / (nu0[:, None] + dr), axis=1)
    return params["kBT"] * params["gradient"] ** 2 * total / (
        params["density"] * params["chi"] * params["volume"]
    )


def mct_psi2_values(d0_values, params):
    d0_values = np.asarray(d0_values, dtype=float)
    k, angle, weight = active_mode_arrays(params)
    dr = mct_renormalized_diffusion(d0_values, k, params)
    nu0 = kinematic_viscosity_for_d0(d0_values, params)
    total = np.sum(
        weight[None, :] * angle[None, :]
        / (k[None, :] ** 2 * dr * (nu0[:, None] + dr)),
        axis=1,
    )
    return params["kBT"] * params["gradient"] ** 2 * total / (
        params["density"] * params["chi"]
    )

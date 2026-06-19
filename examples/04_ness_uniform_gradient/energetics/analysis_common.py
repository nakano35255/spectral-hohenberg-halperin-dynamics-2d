#!/usr/bin/env python3
import csv
from pathlib import Path
import re

import numpy as np


EXAMPLE = Path(__file__).resolve().parent
RAW_DATA_ROOT = EXAMPLE / "raw_data"
PROCESSED_DATA_ROOT = EXAMPLE / "processed_data"
FIGURE_DIR = EXAMPLE / "figures"

QUANTITY_KEYS = ("E_K", "psi2", "induced", "production")
PROCESSED_TIME_SERIES_MAX_POINTS = 2500


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


def format_float(value):
    return f"{float(value):.17g}"


def sampled_indices(size, max_points=PROCESSED_TIME_SERIES_MAX_POINTS):
    if size <= max_points:
        return np.arange(size)
    return np.unique(np.linspace(0, size - 1, max_points).astype(int))


def metadata_path(sc_label, processed_root=PROCESSED_DATA_ROOT):
    return processed_root / sc_label / "metadata.csv"


def time_series_path(sc_label, processed_root=PROCESSED_DATA_ROOT):
    return processed_root / sc_label / "time_series.csv"


def steady_response_path(sc_label, processed_root=PROCESSED_DATA_ROOT):
    return processed_root / sc_label / "steady_response_last_half.csv"


def write_metadata(path, entries):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["key", "value"])
        for key, value in entries:
            writer.writerow([key, value])


def read_metadata(path):
    with path.open(newline="") as handle:
        return {record["key"]: record["value"] for record in csv.DictReader(handle)}


def raw_energetics_case(d0_dir):
    _, d0 = numeric_label(d0_dir, "D0")
    result_dir = d0_dir / "results"
    input_path = d0_dir / "runs" / "input_000.script"
    params = parse_input(input_path)
    files = sorted(result_dir.glob("time_series_*.dat"))
    if not files:
        raise FileNotFoundError(f"no time_series files in {result_dir}")

    times = None
    sums = None
    sums2 = None
    steady_avgs = []
    length_samples = []
    steady_start = None
    steady_end = None

    for path in files:
        columns, data = read_time_series(path)
        sample_times = data[:, col(columns, "time")]
        if times is None:
            times = sample_times
            steady_start = 0.5 * times[-1]
            steady_end = times[-1]
        elif data.shape[0] != times.shape[0] or not np.allclose(sample_times, times):
            raise ValueError(f"time grid differs in {path}")

        transformed = transform(data, columns, d0, params)
        if sums is None:
            sums = np.zeros_like(transformed)
            sums2 = np.zeros_like(transformed)
        sums += transformed
        sums2 += transformed * transformed

        mask = sample_times >= steady_start
        if not np.any(mask):
            raise ValueError(f"steady window has no data in {path}")
        avg = transformed[mask].mean(axis=0)
        steady_avgs.append(avg)
        grad2 = avg[2] * params["volume"] / d0
        length_samples.append(2.0 * np.pi * np.sqrt(avg[1] / grad2))

    samples = len(files)
    mean = sums / samples
    if samples > 1:
        variance = np.maximum(sums2 / samples - mean * mean, 0.0)
        sem = np.sqrt(variance / (samples - 1))
    else:
        sem = np.zeros_like(mean)

    steady_avgs = np.asarray(steady_avgs)
    length_samples = np.asarray(length_samples)
    steady_mean = steady_avgs.mean(axis=0)
    steady_sem = steady_avgs.std(axis=0, ddof=1) / np.sqrt(steady_avgs.shape[0])
    length_sem = length_samples.std(ddof=1) / np.sqrt(length_samples.size) if length_samples.size > 1 else 0.0

    return {
        "d0": d0,
        "d0_text": d0_dir.name.split("_", 1)[1],
        "params": params,
        "times": times,
        "mean": mean,
        "sem": sem,
        "samples": samples,
        "steady_start": steady_start,
        "steady_end": steady_end,
        "psi2": float(steady_mean[1]),
        "psi2_sem": float(steady_sem[1]),
        "induced": float(steady_mean[2]),
        "induced_sem": float(steady_sem[2]),
        "length": float(length_samples.mean()),
        "length_sem": float(length_sem),
    }


def load_raw_energetics_cases(sc_label, raw_root=RAW_DATA_ROOT):
    data_root = raw_root / sc_label
    d0_dirs = sorted(data_root.glob("D0_*"), key=lambda path: numeric_label(path, "D0")[1])
    cases = [raw_energetics_case(d0_dir) for d0_dir in d0_dirs]
    if not cases:
        raise FileNotFoundError(f"no D0 cases found in {data_root}")
    return cases


def schmidt_number_from_params(params):
    return params["eta"] / (params["density"] * params["mobility"])


def params_to_metadata(params):
    return [
        ("grid_x", params["grid"][0]),
        ("grid_y", params["grid"][1]),
        ("length_x", format_float(params["length"][0])),
        ("length_y", format_float(params["length"][1])),
        ("dealias", params["dealias"]),
        ("density", format_float(params["density"])),
        ("kBT", format_float(params["kBT"])),
        ("free_energy_a", format_float(params["free_energy_a"])),
        ("chi", format_float(params["chi"])),
        ("gradient_direction", params["gradient_direction"]),
        ("gradient", format_float(params["gradient"])),
        ("a_uv", format_float(params["a_uv"])),
        ("volume", format_float(params["volume"])),
    ]


def metadata_to_base_params(metadata):
    params = {
        "grid": (int(metadata["grid_x"]), int(metadata["grid_y"])),
        "length": (float(metadata["length_x"]), float(metadata["length_y"])),
        "dealias": metadata["dealias"],
        "density": float(metadata["density"]),
        "kBT": float(metadata["kBT"]),
        "free_energy_a": float(metadata["free_energy_a"]),
        "chi": float(metadata["chi"]),
        "gradient_direction": metadata["gradient_direction"],
        "gradient": float(metadata["gradient"]),
        "a_uv": float(metadata["a_uv"]),
        "volume": float(metadata["volume"]),
    }
    return params


def write_processed_energetics(sc_label, cases, processed_root=PROCESSED_DATA_ROOT):
    if not cases:
        raise RuntimeError("cannot write empty processed energetics data")

    out_dir = processed_root / sc_label
    out_dir.mkdir(parents=True, exist_ok=True)
    first = cases[0]
    sample_counts = [case["samples"] for case in cases]
    steady_starts = [case["steady_start"] for case in cases]
    steady_ends = [case["steady_end"] for case in cases]
    entries = [
        ("case", "energetics"),
        ("sc_label", sc_label),
        ("source_raw_data_dir", f"raw_data/{sc_label}"),
        ("n_cases", len(cases)),
        ("n_samples_min", min(sample_counts)),
        ("n_samples_max", max(sample_counts)),
        ("steady_window", "last_half"),
        ("steady_start_min", format_float(min(steady_starts))),
        ("steady_start_max", format_float(max(steady_starts))),
        ("steady_end_min", format_float(min(steady_ends))),
        ("steady_end_max", format_float(max(steady_ends))),
        ("schmidt_number", format_float(schmidt_number_from_params(first["params"]))),
        ("time_series_sampling", f"uniform_index_max_{PROCESSED_TIME_SERIES_MAX_POINTS}"),
    ]
    entries.extend(params_to_metadata(first["params"]))
    write_metadata(metadata_path(sc_label, processed_root), entries)

    with time_series_path(sc_label, processed_root).open("w", newline="") as handle:
        fieldnames = [
            "d0",
            "d0_label",
            "time",
            "E_K_mean",
            "E_K_sem",
            "psi2_mean",
            "psi2_sem",
            "induced_mean",
            "induced_sem",
            "production_mean",
            "production_sem",
            "n_samples",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for case in cases:
            for row_index in sampled_indices(case["times"].size):
                time = case["times"][row_index]
                values = {
                    "d0": format_float(case["d0"]),
                    "d0_label": case["d0_text"],
                    "time": format_float(time),
                    "n_samples": case["samples"],
                }
                for quantity_index, key in enumerate(QUANTITY_KEYS):
                    values[f"{key}_mean"] = format_float(case["mean"][row_index, quantity_index])
                    values[f"{key}_sem"] = format_float(case["sem"][row_index, quantity_index])
                writer.writerow(values)

    with steady_response_path(sc_label, processed_root).open("w", newline="") as handle:
        fieldnames = [
            "d0",
            "d0_label",
            "steady_start",
            "steady_end",
            "psi2",
            "psi2_sem",
            "induced",
            "induced_sem",
            "length",
            "length_sem",
            "n_samples",
            "eta",
            "mobility",
            "schmidt_number",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for case in cases:
            params = case["params"]
            writer.writerow(
                {
                    "d0": format_float(case["d0"]),
                    "d0_label": case["d0_text"],
                    "steady_start": format_float(case["steady_start"]),
                    "steady_end": format_float(case["steady_end"]),
                    "psi2": format_float(case["psi2"]),
                    "psi2_sem": format_float(case["psi2_sem"]),
                    "induced": format_float(case["induced"]),
                    "induced_sem": format_float(case["induced_sem"]),
                    "length": format_float(case["length"]),
                    "length_sem": format_float(case["length_sem"]),
                    "n_samples": case["samples"],
                    "eta": format_float(params["eta"]),
                    "mobility": format_float(params["mobility"]),
                    "schmidt_number": format_float(schmidt_number_from_params(params)),
                }
            )


def load_processed_energetics(sc_label, processed_root=PROCESSED_DATA_ROOT):
    metadata = read_metadata(metadata_path(sc_label, processed_root))
    base_params = metadata_to_base_params(metadata)

    steady_rows = {}
    with steady_response_path(sc_label, processed_root).open(newline="") as handle:
        for record in csv.DictReader(handle):
            d0 = float(record["d0"])
            params = dict(base_params)
            params["eta"] = float(record["eta"])
            params["mobility"] = float(record["mobility"])
            steady_rows[d0] = {
                "d0": d0,
                "d0_text": record["d0_label"],
                "params": params,
                "samples": int(record["n_samples"]),
                "steady_start": float(record["steady_start"]),
                "steady_end": float(record["steady_end"]),
                "psi2": float(record["psi2"]),
                "psi2_sem": float(record["psi2_sem"]),
                "induced": float(record["induced"]),
                "induced_sem": float(record["induced_sem"]),
                "length": float(record["length"]),
                "length_sem": float(record["length_sem"]),
                "times": [],
                "mean_rows": [],
                "sem_rows": [],
            }

    with time_series_path(sc_label, processed_root).open(newline="") as handle:
        for record in csv.DictReader(handle):
            d0 = float(record["d0"])
            if d0 not in steady_rows:
                raise ValueError(f"time_series has D0={d0:g}, but steady response is missing it")
            steady_rows[d0]["times"].append(float(record["time"]))
            steady_rows[d0]["mean_rows"].append([float(record[f"{key}_mean"]) for key in QUANTITY_KEYS])
            steady_rows[d0]["sem_rows"].append([float(record[f"{key}_sem"]) for key in QUANTITY_KEYS])

    cases = []
    for d0 in sorted(steady_rows):
        case = steady_rows[d0]
        order = np.argsort(np.asarray(case["times"], dtype=float))
        case["times"] = np.asarray(case["times"], dtype=float)[order]
        case["mean"] = np.asarray(case.pop("mean_rows"), dtype=float)[order]
        case["sem"] = np.asarray(case.pop("sem_rows"), dtype=float)[order]
        cases.append(case)
    return cases


def load_or_build_processed_energetics(sc_label, rebuild=False, raw_root=RAW_DATA_ROOT, processed_root=PROCESSED_DATA_ROOT):
    paths = [metadata_path(sc_label, processed_root), time_series_path(sc_label, processed_root), steady_response_path(sc_label, processed_root)]
    if not rebuild and all(path.exists() for path in paths):
        return load_processed_energetics(sc_label, processed_root)

    cases = load_raw_energetics_cases(sc_label, raw_root)
    write_processed_energetics(sc_label, cases, processed_root)
    return cases

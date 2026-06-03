#!/usr/bin/env python3
"""Check the FDR stationary distribution from time_series output.

This checker is intentionally narrow: it assumes the example uses a quadratic
Hamiltonian and tests only whether the measured stationary distribution obeys
equipartition.  For each quadratic sector,

    <E> = (kBT / 2) * N_dof.

The conserved zero modes start from zero in input.script and are excluded.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from statistics import mean, pstdev
import math


ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = ROOT / "examples" / "01_thermal_fluctuations"
INPUT_SCRIPT = EXAMPLE / "input.script"
DEFAULT_RESULTS = EXAMPLE / "results"
SUMMARY_FILE = DEFAULT_RESULTS / "thermal_fluctuation_check.txt"

BURN_IN_FRACTION = 1.0 / 3.0
MEAN_RELATIVE_TOLERANCE = 0.30
STD_RELATIVE_TOLERANCE = 0.60
TOTAL_CLOSURE_RELATIVE_TOLERANCE = 1.0e-10


@dataclass
class InputConfig:
    active_nx: int = 128
    active_ny: int = 128
    order_parameters: int = 0
    time_evolution: str = "srk3/compressible"
    free_energy_type: str = ""
    free_energy_coefficients: dict[str, float] = field(default_factory=dict)
    thermo_type: str = ""
    transport_coefficients: dict[str, float] = field(default_factory=dict)
    kbt: float = 0.0
    noise_enabled: bool = False
    measure_nevery: int | None = None
    time_series_file: Path = DEFAULT_RESULTS / "time_series.dat"
    total_run_steps: int = 0


@dataclass
class EnergyCheck:
    name: str
    measured_mean: float
    measured_std: float
    expected_mean: float
    expected_std: float
    mean_relative_error: float
    std_relative_error: float
    passed: bool


def tokenize(line: str) -> list[str]:
    return line.split("#", 1)[0].split()


def parse_key_values(tokens: list[str]) -> dict[str, float]:
    values: dict[str, float] = {}
    for i in range(0, len(tokens), 2):
        if i + 1 < len(tokens):
            values[tokens[i]] = float(tokens[i + 1])
    return values


def parse_input(path: Path) -> InputConfig:
    cfg = InputConfig()

    with path.open() as handle:
        for raw_line in handle:
            tokens = tokenize(raw_line)
            if not tokens:
                continue

            if tokens[0] == "grid":
                cfg.active_nx = int(tokens[1])
                cfg.active_ny = int(tokens[2])
                continue

            if tokens[0] in {"order_parameters", "order_parameter"}:
                cfg.order_parameters = int(tokens[1])
                continue

            if tokens[0] == "time_evolution":
                cfg.time_evolution = tokens[1]
                continue

            if tokens[:2] == ["model", "thermo"]:
                cfg.thermo_type = tokens[2]
                continue

            if tokens[:2] == ["model", "free_energy"]:
                if tokens[2] == "coeff":
                    entries = tokens[3:]
                else:
                    cfg.free_energy_type = tokens[2]
                    entries = tokens[3:]
                cfg.free_energy_coefficients.update(parse_key_values(entries))
                continue

            if tokens[:2] == ["model", "transport"]:
                entries = tokens[3:]
                cfg.transport_coefficients.update(parse_key_values(entries))
                continue

            if len(tokens) >= 5 and tokens[0] == "fix" and tokens[2] == "all" and tokens[3] == "noise":
                cfg.noise_enabled = tokens[4] == "on"
                if cfg.noise_enabled:
                    values = parse_key_values(tokens[5:])
                    cfg.kbt = values.get("kBT", values.get("temperature", 1.0))
                continue

            if tokens[0] == "measure" and len(tokens) >= 4 and tokens[2] == "time_series":
                if tokens[3] == "on":
                    entries = tokens[4:]
                    i = 0
                    while i < len(entries):
                        key = entries[i]
                        if key == "target":
                            break
                        if i + 1 >= len(entries):
                            break
                        value = entries[i + 1]
                        if key == "nevery":
                            cfg.measure_nevery = int(value)
                        elif key == "file":
                            cfg.time_series_file = ROOT / value
                        i += 2
                continue

            if tokens[0] == "run":
                cfg.total_run_steps += int(tokens[1])

    return cfg


def active_real_dofs_excluding_zero(active_nx: int, active_ny: int) -> int:
    """Count active real Fourier degrees of freedom except the zero mode.

    This mirrors SpectralMask2D:
      gx = 0, ..., active_nx / 2 - 1
      abs(signed_ky) < active_ny / 2
    """

    dofs = 0
    for gx in range(active_nx // 2):
        for ky in range(-(active_ny // 2) + 1, active_ny // 2):
            if gx == 0:
                if ky > 0:
                    dofs += 2
            else:
                dofs += 2
    return dofs


def read_time_series_energy(path: Path) -> list[dict[str, float]]:
    if not path.exists():
        raise FileNotFoundError(f"missing time_series file: {path}")

    rows: list[dict[str, float]] = []
    header: list[str] | None = None
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                fields = line[1:].split()
                if len(fields) >= 6 and fields[0] == "step" and fields[1] == "time":
                    header = fields
                continue

            if header is None:
                raise ValueError(f"missing time_series header in {path}")

            parts = line.split()
            if len(parts) != len(header):
                raise ValueError(f"unexpected time_series row: {line}")

            values = dict(zip(header, parts))
            for required in ("E_T", "E_K", "E_psi", "E_C"):
                if required not in values:
                    raise ValueError(f"time_series output is missing required column: {required}")

            rows.append(
                {
                    "step": int(values["step"]),
                    "time": float(values["time"]),
                    "total": float(values["E_T"]),
                    "kinetic": float(values["E_K"]),
                    "free": float(values["E_psi"]),
                    "compressibility": float(values["E_C"]),
                }
            )

    if not rows:
        raise ValueError(f"no samples found in {path}")
    return rows


def late_samples(rows: list[dict[str, float]], cfg: InputConfig) -> list[dict[str, float]]:
    burn_in_step = int(round(cfg.total_run_steps * BURN_IN_FRACTION))
    samples = [row for row in rows if row["step"] >= burn_in_step]
    if not samples:
        samples = rows[len(rows) // 2 :]
    return samples


def finite_rows(rows: list[dict[str, float]]) -> bool:
    return all(math.isfinite(value) for row in rows for value in row.values())


def total_closure_ok(rows: list[dict[str, float]]) -> tuple[bool, float, float]:
    max_error = max(
        abs(row["total"] - (row["kinetic"] + row["free"] + row["compressibility"]))
        for row in rows
    )
    scale = max(1.0, max(abs(row["total"]) for row in rows))
    tolerance = TOTAL_CLOSURE_RELATIVE_TOLERANCE * scale
    return max_error <= tolerance, max_error, tolerance


def relative_error(measured: float, expected: float) -> float:
    if expected == 0.0:
        return 0.0 if measured == 0.0 else math.inf
    return abs(measured - expected) / abs(expected)


def make_energy_check(name: str, values: list[float], dofs: int, kbt: float) -> EnergyCheck:
    measured_mean = mean(values)
    measured_std = pstdev(values)
    expected_mean = 0.5 * kbt * dofs
    expected_std = kbt * math.sqrt(0.5 * dofs)
    mean_err = relative_error(measured_mean, expected_mean)
    std_err = relative_error(measured_std, expected_std)
    return EnergyCheck(
        name=name,
        measured_mean=measured_mean,
        measured_std=measured_std,
        expected_mean=expected_mean,
        expected_std=expected_std,
        mean_relative_error=mean_err,
        std_relative_error=std_err,
        passed=(
            mean_err <= MEAN_RELATIVE_TOLERANCE
            and std_err <= STD_RELATIVE_TOLERANCE
        ),
    )


def expected_checks(cfg: InputConfig, samples: list[dict[str, float]]) -> list[EnergyCheck]:
    if not cfg.noise_enabled or cfg.kbt <= 0.0:
        raise ValueError("FDR check requires fix noise on with kBT > 0.")

    if cfg.free_energy_type != "quadratic":
        raise ValueError("FDR check requires model free_energy quadratic.")

    scalar_dofs = active_real_dofs_excluding_zero(cfg.active_nx, cfg.active_ny)
    checks: list[EnergyCheck] = []

    active_order_parameters = 0
    for op in range(cfg.order_parameters):
        coefficient = cfg.free_energy_coefficients.get(f"a[{op}]", 0.0)
        mobility = cfg.transport_coefficients.get(f"M[{op},{op}]", 0.0)
        if coefficient > 0.0 and mobility > 0.0:
            active_order_parameters += 1

    if active_order_parameters > 0:
        checks.append(make_energy_check(
            "order-parameter free energy",
            [row["free"] for row in samples],
            scalar_dofs * active_order_parameters,
            cfg.kbt,
        ))

    if "compressible" in cfg.time_evolution:
        if cfg.thermo_type == "linear_eos":
            checks.append(make_energy_check(
                "density compressibility energy",
                [row["compressibility"] for row in samples],
                scalar_dofs,
                cfg.kbt,
            ))

        eta = cfg.transport_coefficients.get("eta", cfg.transport_coefficients.get("shear_viscosity", 0.0))
        zeta = cfg.transport_coefficients.get("zeta", cfg.transport_coefficients.get("bulk_viscosity", 0.0))
        if eta > 0.0 and zeta > 0.0:
            checks.append(make_energy_check(
                "momentum kinetic energy",
                [row["kinetic"] for row in samples],
                2 * scalar_dofs,
                cfg.kbt,
            ))

    if checks:
        total_dofs = sum(int(round(2.0 * check.expected_mean / cfg.kbt)) for check in checks)
        checks.append(make_energy_check(
            "total quadratic energy",
            [row["total"] for row in samples],
            total_dofs,
            cfg.kbt,
        ))

    return checks


def format_summary(
    cfg: InputConfig,
    rows: list[dict[str, float]],
    samples: list[dict[str, float]],
    checks: list[EnergyCheck],
    closure: tuple[bool, float, float],
) -> str:
    closure_ok, closure_error, closure_tolerance = closure
    finite_ok = finite_rows(rows)
    all_checks_ok = all(check.passed for check in checks)
    ok = finite_ok and closure_ok and all_checks_ok

    lines = [
        "Thermal FDR stationary-distribution check",
        f"  input                 : {INPUT_SCRIPT.relative_to(ROOT)}",
        f"  time_series           : {cfg.time_series_file.relative_to(ROOT)}",
        f"  time evolution        : {cfg.time_evolution}",
        f"  free energy           : {cfg.free_energy_type}",
        f"  active grid           : {cfg.active_nx} x {cfg.active_ny}",
        f"  active real dofs      : {active_real_dofs_excluding_zero(cfg.active_nx, cfg.active_ny)}",
        f"  kBT                   : {cfg.kbt:.8e}",
        f"  samples               : {len(rows)}",
        f"  stationary samples    : {len(samples)}",
        "",
        "Preconditions:",
        f"  [{'PASS' if finite_ok else 'FAIL'}] finite time_series energy columns",
        (
            f"  [{'PASS' if closure_ok else 'FAIL'}] total column closure: "
            f"max error {closure_error:.3e}, tolerance {closure_tolerance:.3e}"
        ),
        "",
        "Equipartition checks:",
    ]

    for check in checks:
        status = "PASS" if check.passed else "FAIL"
        lines.extend([
            f"  [{status}] {check.name}",
            f"      mean measured/expected : {check.measured_mean:.8e} / {check.expected_mean:.8e}",
            f"      std  measured/expected : {check.measured_std:.8e} / {check.expected_std:.8e}",
            f"      measured/expected ratio: mean {check.measured_mean / check.expected_mean:.4f}, std {check.measured_std / check.expected_std:.4f}",
            f"      relative errors        : mean {check.mean_relative_error:.4f}, std {check.std_relative_error:.4f}",
        ])

    if not checks:
        lines.append("  [FAIL] no quadratic stochastic sector was found")
        ok = False

    lines.extend([
        "",
        f"tolerances              : mean {MEAN_RELATIVE_TOLERANCE:.2f}, std {STD_RELATIVE_TOLERANCE:.2f}",
        f"result                  : {'PASS' if ok else 'FAIL'}",
    ])

    return "\n".join(lines) + "\n"


def main() -> int:
    cfg = parse_input(INPUT_SCRIPT)
    rows = read_time_series_energy(cfg.time_series_file)
    samples = late_samples(rows, cfg)
    closure = total_closure_ok(rows)
    checks = expected_checks(cfg, samples)

    summary = format_summary(cfg, rows, samples, checks, closure)
    print(summary, end="")
    SUMMARY_FILE.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_FILE.write_text(summary)

    return 0 if "result                  : PASS" in summary else 1


if __name__ == "__main__":
    raise SystemExit(main())

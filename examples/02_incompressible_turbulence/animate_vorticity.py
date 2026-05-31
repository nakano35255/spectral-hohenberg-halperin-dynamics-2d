#!/usr/bin/env python3
"""Create an animation from incompressible turbulence spectral snapshots."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation


ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = ROOT / "examples" / "02_incompressible_turbulence"
DEFAULT_INPUT = EXAMPLE / "input.script"
DEFAULT_SNAPSHOT_GLOB = EXAMPLE / "results" / "spectral_step*.snapshot"
DEFAULT_OUTPUT = EXAMPLE / "results" / "vorticity.mp4"


@dataclass
class SpectralSnapshot:
    path: Path
    step: int
    time: float
    nx: int
    ny: int
    rho_hat: np.ndarray
    jx_hat: np.ndarray
    jy_hat: np.ndarray


def parse_input_lengths(path: Path) -> tuple[float | None, float | None]:
    if not path.exists():
        return None, None

    lx = None
    ly = None
    with path.open() as handle:
        for raw_line in handle:
            tokens = raw_line.split("#", 1)[0].split()
            if len(tokens) == 3 and tokens[0] == "length":
                lx = float(tokens[1])
                ly = float(tokens[2])
                break

    return lx, ly


def parse_header(path: Path) -> tuple[int, float, int, int, list[str]]:
    step = -1
    time = 0.0
    nkx = -1
    nky = -1
    columns: list[str] = []

    with path.open() as handle:
        for line in handle:
            if not line.startswith("#"):
                break

            parts = line[1:].split()
            if len(parts) >= 4 and parts[0] == "step":
                step = int(parts[1])
                time = float(parts[3])
            elif len(parts) >= 4 and parts[0] == "nkx":
                nkx = int(parts[1])
                nky = int(parts[3])
            elif len(parts) >= 2 and parts[0] == "kx_index":
                columns = parts

    if step < 0 or nkx <= 0 or nky <= 0 or not columns:
        raise ValueError(f"invalid spectral snapshot header: {path}")

    return step, time, nkx, nky, columns


def read_spectral_snapshot(path: Path) -> SpectralSnapshot:
    step, time, nkx, nky, columns = parse_header(path)
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)

    column_index = {name: i for i, name in enumerate(columns)}
    required = [
        "kx_index",
        "ky_index",
        "rho_real",
        "rho_imag",
        "jx_real",
        "jx_imag",
        "jy_real",
        "jy_imag",
    ]
    for name in required:
        if name not in column_index:
            raise ValueError(f"missing column {name} in {path}")

    rho_hat = np.zeros((nky, nkx), dtype=np.complex128)
    jx_hat = np.zeros((nky, nkx), dtype=np.complex128)
    jy_hat = np.zeros((nky, nkx), dtype=np.complex128)

    gx_values = data[:, column_index["kx_index"]].astype(np.int64)
    gy_values = data[:, column_index["ky_index"]].astype(np.int64)

    rho_hat[gy_values, gx_values] = (
        data[:, column_index["rho_real"]]
        + 1j * data[:, column_index["rho_imag"]]
    )
    jx_hat[gy_values, gx_values] = (
        data[:, column_index["jx_real"]]
        + 1j * data[:, column_index["jx_imag"]]
    )
    jy_hat[gy_values, gx_values] = (
        data[:, column_index["jy_real"]]
        + 1j * data[:, column_index["jy_imag"]]
    )

    nx = 2 * (nkx - 1)
    return SpectralSnapshot(path, step, time, nx, nky, rho_hat, jx_hat, jy_hat)


def field_from_snapshot(snapshot: SpectralSnapshot, lx: float, ly: float, field: str) -> np.ndarray:
    rho0 = snapshot.rho_hat[0, 0].real / (snapshot.nx * snapshot.ny)
    if rho0 == 0.0:
        raise ValueError(f"zero reference density in {snapshot.path}")

    ux_hat = snapshot.jx_hat / rho0
    uy_hat = snapshot.jy_hat / rho0

    if field == "speed":
        ux = np.fft.irfft2(ux_hat, s=(snapshot.ny, snapshot.nx))
        uy = np.fft.irfft2(uy_hat, s=(snapshot.ny, snapshot.nx))
        return np.sqrt(ux * ux + uy * uy)

    gx = np.arange(snapshot.nx // 2 + 1, dtype=np.float64)
    gy_index = np.arange(snapshot.ny)
    signed_gy = np.where(gy_index <= snapshot.ny // 2, gy_index, gy_index - snapshot.ny)
    kx = (2.0 * np.pi / lx) * gx[np.newaxis, :]
    ky = (2.0 * np.pi / ly) * signed_gy[:, np.newaxis]

    omega_hat = 1j * (kx * uy_hat - ky * ux_hat)
    return np.fft.irfft2(omega_hat, s=(snapshot.ny, snapshot.nx))


def natural_step_key(path: Path) -> tuple[int, str]:
    match = re.search(r"_step(\d+)\.snapshot$", path.name)
    if match:
        return int(match.group(1)), path.name
    return -1, path.name


def load_frames(
    snapshot_glob: Path,
    input_script: Path,
    field: str,
    stride: int,
    max_frames: int | None,
) -> tuple[list[np.ndarray], list[SpectralSnapshot], float, float]:
    paths = sorted(snapshot_glob.parent.glob(snapshot_glob.name), key=natural_step_key)
    paths = paths[::stride]
    if max_frames is not None:
        paths = paths[:max_frames]
    if not paths:
        raise FileNotFoundError(f"no snapshots matched: {snapshot_glob}")

    snapshots = [read_spectral_snapshot(path) for path in paths]
    lx, ly = parse_input_lengths(input_script)
    if lx is None:
        lx = float(snapshots[0].nx)
    if ly is None:
        ly = float(snapshots[0].ny)

    frames = [field_from_snapshot(snapshot, lx, ly, field) for snapshot in snapshots]
    return frames, snapshots, lx, ly


def save_animation(
    frames: list[np.ndarray],
    snapshots: list[SpectralSnapshot],
    lx: float,
    ly: float,
    field: str,
    output: Path,
    fps: int,
    dpi: int,
    percentile: float,
) -> Path:
    output.parent.mkdir(parents=True, exist_ok=True)

    stack = np.asarray(frames)
    if field == "speed":
        vmin = 0.0
        vmax = float(np.percentile(stack, percentile))
        cmap = "viridis"
        colorbar_label = r"$|\mathbf{u}|$"
    else:
        vmax = float(np.percentile(np.abs(stack), percentile))
        vmin = -vmax
        cmap = "RdBu_r"
        colorbar_label = r"$\omega$"

    if vmax == 0.0:
        vmax = 1.0
        if field == "speed":
            vmin = 0.0
        else:
            vmin = -1.0

    fig, ax = plt.subplots(figsize=(7.0, 6.0), constrained_layout=True)
    image = ax.imshow(
        frames[0],
        origin="lower",
        extent=(0.0, lx, 0.0, ly),
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        interpolation="bilinear",
    )
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    title = ax.set_title("")
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label(colorbar_label)

    def update(frame_index: int):
        snapshot = snapshots[frame_index]
        image.set_data(frames[frame_index])
        title.set_text(f"{field}  step={snapshot.step}  t={snapshot.time:.6g}")
        return image, title

    movie = animation.FuncAnimation(
        fig,
        update,
        frames=len(frames),
        interval=1000.0 / fps,
        blit=False,
    )

    suffix = output.suffix.lower()
    writer = "pillow" if suffix == ".gif" else "ffmpeg"

    try:
        movie.save(output, writer=writer, fps=fps, dpi=dpi)
        saved = output
    except Exception:
        if suffix == ".gif":
            raise
        fallback = output.with_suffix(".gif")
        movie.save(fallback, writer="pillow", fps=fps, dpi=dpi)
        saved = fallback

    plt.close(fig)
    return saved


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Animate vorticity or speed from spectral snapshots."
    )
    parser.add_argument("--snapshots", type=Path, default=DEFAULT_SNAPSHOT_GLOB)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--field", choices=("vorticity", "speed"), default="vorticity")
    parser.add_argument("--fps", type=int, default=12)
    parser.add_argument("--dpi", type=int, default=140)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--max-frames", type=int, default=None)
    parser.add_argument("--percentile", type=float, default=99.0)
    args = parser.parse_args()

    if args.stride <= 0:
        raise ValueError("--stride must be positive")
    if not (0.0 < args.percentile <= 100.0):
        raise ValueError("--percentile must be in (0, 100]")

    frames, snapshots, lx, ly = load_frames(
        args.snapshots,
        args.input,
        args.field,
        args.stride,
        args.max_frames,
    )
    saved = save_animation(
        frames,
        snapshots,
        lx,
        ly,
        args.field,
        args.output,
        args.fps,
        args.dpi,
        args.percentile,
    )

    print(f"saved {saved}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

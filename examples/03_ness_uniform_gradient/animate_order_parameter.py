#!/usr/bin/env python3
"""Animate an order-parameter field from physical snapshots."""

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
EXAMPLE = ROOT / "examples" / "03_ness_uniform_gradient"
DEFAULT_INPUT = EXAMPLE / "input.script"
DEFAULT_SNAPSHOT_GLOB = EXAMPLE / "results" / "physical_step*.snapshot"
DEFAULT_OUTPUT = EXAMPLE / "results" / "order_parameter.mp4"


@dataclass
class PhysicalSnapshot:
    path: Path
    step: int
    time: float
    nx: int
    ny: int
    field: np.ndarray


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


def parse_header(path: Path) -> tuple[int, float, int, int, int, list[str]]:
    space = ""
    step = -1
    time = 0.0
    nx = -1
    ny = -1
    order_parameters = -1
    columns: list[str] = []

    with path.open() as handle:
        for line in handle:
            if not line.startswith("#"):
                break

            parts = line[1:].split()
            if len(parts) >= 3 and parts[0] == "snapshot" and parts[1] == "space":
                space = parts[2]
            elif len(parts) >= 4 and parts[0] == "step":
                step = int(parts[1])
                time = float(parts[3])
            elif len(parts) >= 6 and parts[0] == "nx":
                nx = int(parts[1])
                ny = int(parts[3])
                order_parameters = int(parts[5])
            elif len(parts) >= 2 and parts[0] == "x" and parts[1] == "y":
                columns = parts

    if space != "physical":
        raise ValueError(f"{path} is not a physical snapshot")
    if step < 0 or nx <= 0 or ny <= 0 or order_parameters < 0 or not columns:
        raise ValueError(f"invalid physical snapshot header: {path}")

    return step, time, nx, ny, order_parameters, columns


def read_physical_snapshot(path: Path, component: int) -> PhysicalSnapshot:
    step, time, nx, ny, order_parameters, columns = parse_header(path)
    if component < 0 or component >= order_parameters:
        raise ValueError(
            f"component {component} is out of range for {path} "
            f"(order_parameters={order_parameters})"
        )

    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)

    column_index = {name: i for i, name in enumerate(columns)}
    field_name = f"psi_com{component}"
    for name in ("x", "y", field_name):
        if name not in column_index:
            raise ValueError(f"missing column {name} in {path}")

    if data.shape[0] != nx * ny:
        raise ValueError(f"unexpected row count in {path}: got {data.shape[0]}, expected {nx * ny}")

    gx_values = data[:, column_index["x"]].astype(np.int64)
    gy_values = data[:, column_index["y"]].astype(np.int64)

    field = np.zeros((ny, nx), dtype=np.float64)
    field[gy_values, gx_values] = data[:, column_index[field_name]]

    return PhysicalSnapshot(path, step, time, nx, ny, field)


def natural_step_key(path: Path) -> tuple[int, str]:
    match = re.search(r"_step(\d+)\.snapshot$", path.name)
    if match:
        return int(match.group(1)), path.name
    return -1, path.name


def load_frames(
    snapshot_glob: Path,
    input_script: Path,
    component: int,
    mode: str,
    stride: int,
    max_frames: int | None,
) -> tuple[list[np.ndarray], list[PhysicalSnapshot], float, float]:
    paths = sorted(snapshot_glob.parent.glob(snapshot_glob.name), key=natural_step_key)
    paths = paths[::stride]
    if max_frames is not None:
        paths = paths[:max_frames]
    if not paths:
        raise FileNotFoundError(f"no snapshots matched: {snapshot_glob}")

    snapshots = [read_physical_snapshot(path, component) for path in paths]
    frames = [snapshot.field.copy() for snapshot in snapshots]

    if mode == "fluctuation":
        frames = [frame - np.mean(frame) for frame in frames]

    lx, ly = parse_input_lengths(input_script)
    if lx is None:
        lx = float(snapshots[0].nx)
    if ly is None:
        ly = float(snapshots[0].ny)

    return frames, snapshots, lx, ly


def color_limits(frames: list[np.ndarray], mode: str, percentile: float) -> tuple[float, float, str]:
    stack = np.asarray(frames)

    if mode == "fluctuation":
        vmax = float(np.percentile(np.abs(stack), percentile))
        if vmax == 0.0:
            vmax = 1.0
        return -vmax, vmax, "RdBu_r"

    lower = 100.0 - percentile
    upper = percentile
    vmin, vmax = np.percentile(stack, [lower, upper])
    vmin = float(vmin)
    vmax = float(vmax)
    if vmin == vmax:
        vmin -= 0.5
        vmax += 0.5

    return vmin, vmax, "viridis"


def save_animation(
    frames: list[np.ndarray],
    snapshots: list[PhysicalSnapshot],
    lx: float,
    ly: float,
    component: int,
    mode: str,
    output: Path,
    fps: int,
    dpi: int,
    percentile: float,
) -> Path:
    output.parent.mkdir(parents=True, exist_ok=True)

    vmin, vmax, cmap = color_limits(frames, mode, percentile)
    label = rf"$\psi_{{{component}}}$"
    if mode == "fluctuation":
        label = rf"$\psi_{{{component}}} - \langle \psi_{{{component}}} \rangle$"

    fig, ax = plt.subplots(figsize=(7.0, 6.0), constrained_layout=True)
    image = ax.imshow(
        frames[0],
        origin="lower",
        extent=(0.0, lx, 0.0, ly),
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
    )
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    title = ax.set_title("")
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label(label)

    def update(frame_index: int):
        snapshot = snapshots[frame_index]
        image.set_data(frames[frame_index])
        title.set_text(f"{label}  step={snapshot.step}  t={snapshot.time:.6g}")
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
        description="Animate an order-parameter field from physical snapshots."
    )
    parser.add_argument("--snapshots", type=Path, default=DEFAULT_SNAPSHOT_GLOB)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--component", type=int, default=0)
    parser.add_argument("--mode", choices=("raw", "fluctuation"), default="raw")
    parser.add_argument("--fps", type=int, default=12)
    parser.add_argument("--dpi", type=int, default=140)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--max-frames", type=int, default=None)
    parser.add_argument("--percentile", type=float, default=99.5)
    args = parser.parse_args()

    if args.stride <= 0:
        raise ValueError("--stride must be positive")
    if args.fps <= 0:
        raise ValueError("--fps must be positive")
    if not (50.0 <= args.percentile <= 100.0):
        raise ValueError("--percentile must be in [50, 100]")

    frames, snapshots, lx, ly = load_frames(
        args.snapshots,
        args.input,
        args.component,
        args.mode,
        args.stride,
        args.max_frames,
    )
    saved = save_animation(
        frames,
        snapshots,
        lx,
        ly,
        args.component,
        args.mode,
        args.output,
        args.fps,
        args.dpi,
        args.percentile,
    )

    print(f"saved {saved}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

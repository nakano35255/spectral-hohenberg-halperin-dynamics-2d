#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


def parse_header(path: Path) -> tuple[int, int, int, float]:
    nx = ny = step = None
    time = None
    with path.open() as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            parts = line.strip().split()
            if len(parts) >= 5 and parts[1] == "step":
                step = int(parts[2])
                time = float(parts[4])
            if len(parts) >= 5 and parts[1] == "nx":
                nx = int(parts[2])
                ny = int(parts[4])

    if nx is None or ny is None or step is None or time is None:
        raise RuntimeError(f"Could not parse snapshot header: {path}")

    return nx, ny, step, time


def load_psi(path: Path) -> tuple[np.ndarray, int, float]:
    nx, ny, step, time = parse_header(path)
    data = np.loadtxt(path, comments="#")
    expected_columns = 5
    if data.ndim != 2 or data.shape[1] < expected_columns:
        raise RuntimeError(f"Unexpected snapshot columns: {path}")

    psi = data[:, 3].reshape((ny, nx))
    return psi, step, time


def read_domain_length(path: Path) -> tuple[float, float] | None:
    if not path.exists():
        return None

    with path.open() as handle:
        for line in handle:
            tokens = line.split()
            if len(tokens) >= 3 and tokens[0] == "length":
                return float(tokens[1]), float(tokens[2])

    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-glob", default="raw_data/sine2d_grid256_L8192_dt4_T2400000/samples/sample_000/snapshots/physical_step*.snapshot")
    parser.add_argument("--output", default="figures/sample000_order_parameter.gif")
    parser.add_argument("--fps", type=int, default=10)
    parser.add_argument("--dpi", type=int, default=140)
    args = parser.parse_args()

    base_dir = Path(__file__).resolve().parent
    snapshot_paths = sorted(base_dir.glob(args.input_glob))
    if not snapshot_paths:
        raise RuntimeError(f"No snapshots matched: {base_dir / args.input_glob}")

    frames = []
    steps = []
    times = []
    for path in snapshot_paths:
        psi, step, time = load_psi(path)
        frames.append(psi)
        steps.append(step)
        times.append(time)

    vmin = min(float(np.min(frame)) for frame in frames)
    vmax = max(float(np.max(frame)) for frame in frames)
    if vmin == vmax:
        vmin -= 1.0
        vmax += 1.0

    output_path = base_dir / args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)

    domain_length = read_domain_length(base_dir / "input.script")
    imshow_kwargs = {
        "origin": "lower",
        "cmap": "turbo",
        "vmin": vmin,
        "vmax": vmax,
    }
    if domain_length is not None:
        lx, ly = domain_length
        imshow_kwargs["extent"] = (0.0, lx, 0.0, ly)

    fig, ax = plt.subplots(figsize=(5.2, 4.8), constrained_layout=True)
    image = ax.imshow(frames[0], **imshow_kwargs)
    title = ax.set_title(f"step {steps[0]}  time {times[0]:.3g}")
    if domain_length is None:
        ax.set_xlabel("x grid")
        ax.set_ylabel("y grid")
    else:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    fig.colorbar(image, ax=ax, label="psi[0]")

    def update(index: int):
        image.set_data(frames[index])
        title.set_text(f"step {steps[index]}  time {times[index]:.3g}")
        return image, title

    movie = animation.FuncAnimation(
        fig,
        update,
        frames=len(frames),
        interval=1000 / max(args.fps, 1),
        blit=False,
    )
    writer = animation.PillowWriter(fps=args.fps)
    movie.save(output_path, writer=writer, dpi=args.dpi)
    plt.close(fig)

    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

from pathlib import Path
import sys


def read_last_block(path):
    blocks = []
    current = []

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("# block"):
                if current:
                    blocks.append(current)
                    current = []
                continue
            if line.startswith("#"):
                continue

            values = line.split()
            if len(values) != 6:
                raise RuntimeError(f"Unexpected row in {path}: {line}")

            current.append([float(value) for value in values])

    if current:
        blocks.append(current)

    if not blocks:
        raise RuntimeError(f"No data blocks found in {path}")

    return blocks[-1]


def write_summary(path, rows):
    transfer = sum(row[2] for row in rows)
    dissipation = sum(row[3] for row in rows)
    production = sum(row[4] for row in rows)
    total = sum(row[5] for row in rows)

    with path.open("w") as handle:
        handle.write("# shell-integrated budget from the last block\n")
        handle.write(f"transfer {transfer:.16e}\n")
        handle.write(f"dissipation {dissipation:.16e}\n")
        handle.write(f"production {production:.16e}\n")
        handle.write(f"total {total:.16e}\n")

    print(f"Wrote {path}")
    print(f"transfer   {transfer:.8e}")
    print(f"dissipation {dissipation:.8e}")
    print(f"production {production:.8e}")
    print(f"total      {total:.8e}")


def write_plot(path, rows):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not available; skipped figure output.")
        return

    k = [row[0] for row in rows]
    transfer = [row[2] for row in rows]
    dissipation = [row[3] for row in rows]
    production = [row[4] for row in rows]
    total = [row[5] for row in rows]

    fig, ax = plt.subplots(figsize=(7.0, 4.5), constrained_layout=True)
    ax.axhline(0.0, color="0.6", linewidth=0.8)
    ax.plot(k, transfer, label="transfer")
    ax.plot(k, dissipation, label="dissipation")
    ax.plot(k, production, label="production")
    ax.plot(k, total, color="black", linewidth=1.6, label="total")
    ax.set_xlabel("k")
    ax.set_ylabel("budget spectrum")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.25)
    fig.savefig(path, dpi=180)
    print(f"Wrote {path}")


def main():
    root = Path(__file__).resolve().parent
    input_path = Path(sys.argv[1]) if len(sys.argv) > 1 else root / "results" / "budget_shell.dat"
    summary_path = Path(sys.argv[2]) if len(sys.argv) > 2 else input_path.with_name("budget_shell_summary.txt")
    figure_path = Path(sys.argv[3]) if len(sys.argv) > 3 else input_path.with_suffix(".png")

    rows = read_last_block(input_path)
    write_summary(summary_path, rows)
    write_plot(figure_path, rows)


if __name__ == "__main__":
    main()

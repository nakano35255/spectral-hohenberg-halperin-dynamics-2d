#!/usr/bin/env python3
import argparse
import os
import struct
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = ROOT / "examples" / "05_yokota_green_kubo"
TEMPLATE = EXAMPLE / "input_ohtaka_linear_sample.template"


def random_seed() -> int:
    raw = struct.unpack("<I", os.urandom(4))[0]
    return raw % 2147483646 + 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", default="prepared")
    parser.add_argument("--samples", type=int, default=128)
    args = parser.parse_args()

    run_dir = EXAMPLE / "runs" / "ohtaka_linear_{}".format(args.run_id)
    result_dir = EXAMPLE / "results" / "ohtaka_linear_{}".format(args.run_id)
    run_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)

    template = TEMPLATE.read_text()
    seeds_path = run_dir / "seeds.dat"
    with seeds_path.open("w") as seeds:
        seeds.write("# sample noise_seed init_seed input output\n")
        for sample in range(args.samples):
            input_path = run_dir / "input_{:03d}.script".format(sample)
            output_path = result_dir / "yokota_green_kubo_{:03d}.dat".format(sample)
            noise_seed = random_seed()
            init_seed = random_seed()

            text = (
                template
                .replace("@NOISE_SEED@", str(noise_seed))
                .replace("@INIT_SEED@", str(init_seed))
                .replace("@YKGK_FILE@", str(output_path.relative_to(ROOT)))
            )
            input_path.write_text(text)
            seeds.write(
                "{:03d} {:d} {:d} {} {}\n".format(
                    sample,
                    noise_seed,
                    init_seed,
                    input_path.relative_to(ROOT),
                    output_path.relative_to(ROOT),
                )
            )

    print("run id     : {}".format(args.run_id))
    print("inputs     : {}".format(run_dir.relative_to(ROOT)))
    print("results    : {}".format(result_dir.relative_to(ROOT)))
    print("seeds      : {}".format(seeds_path.relative_to(ROOT)))


if __name__ == "__main__":
    main()

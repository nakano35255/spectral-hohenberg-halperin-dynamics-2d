#!/bin/bash
#SBATCH -p i8cpu
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -c 1
#SBATCH -t 00:12:00
#SBATCH -J ykgk-linear
#SBATCH -o ykgk-linear-%j.out
#SBATCH -e ykgk-linear-%j.err

set -euo pipefail

find_repo_root() {
    local dir=$1
    while [ "$dir" != "/" ]; do
        if [ -f "$dir/Makefile.ohtaka" ] && [ -d "$dir/src" ]; then
            printf "%s\n" "$dir"
            return 0
        fi
        dir=$(dirname "$dir")
    done
    return 1
}

REPO_ROOT=$(find_repo_root "${SLURM_SUBMIT_DIR:-$PWD}")
cd "$REPO_ROOT"

module purge
module load oneapi_compiler/2023.0.0 oneapi_mpi/2023.0.0

export HEFFTE_ROOT=${HEFFTE_ROOT:-$HOME/local/heffte-oneapi-fftw}
export FFTW_ROOT=${FFTW_ROOT:-$HOME/local/fftw-oneapi}
export LD_LIBRARY_PATH=$HEFFTE_ROOT/lib:$FFTW_ROOT/lib:${LD_LIBRARY_PATH:-}

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export MKL_DYNAMIC=FALSE

export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=core

EXAMPLE_DIR=examples/05_yokota_green_kubo
TEMPLATE=$EXAMPLE_DIR/input_ohtaka_linear_sample.template
JOB_TAG=${SLURM_JOB_ID:-manual}
RUN_DIR=$EXAMPLE_DIR/runs/ohtaka_linear_$JOB_TAG
RESULT_DIR=$EXAMPLE_DIR/results/ohtaka_linear_$JOB_TAG
NSAMPLES=128

mkdir -p "$RUN_DIR" "$RESULT_DIR"

echo "repo       : $REPO_ROOT"
echo "template   : $TEMPLATE"
echo "run_dir    : $RUN_DIR"
echo "result_dir : $RESULT_DIR"
echo "nsamples   : $NSAMPLES"

PIDS=()

for SAMPLE_ID in $(seq 0 $((NSAMPLES - 1))); do
    INPUT_FILE=$(printf "%s/input_%03d.script" "$RUN_DIR" "$SAMPLE_ID")
    RESULT_FILE=$(printf "%s/yokota_green_kubo_%03d.dat" "$RESULT_DIR" "$SAMPLE_ID")
    STDOUT_FILE=$(printf "%s/run_%03d.out" "$RESULT_DIR" "$SAMPLE_ID")
    STDERR_FILE=$(printf "%s/run_%03d.err" "$RESULT_DIR" "$SAMPLE_ID")

    NOISE_SEED=$((12345000 + 2 * SAMPLE_ID))
    INIT_SEED=$((12345000 + 2 * SAMPLE_ID + 1))

    sed \
        -e "s|@NOISE_SEED@|$NOISE_SEED|g" \
        -e "s|@INIT_SEED@|$INIT_SEED|g" \
        -e "s|@YKGK_FILE@|$RESULT_FILE|g" \
        "$TEMPLATE" > "$INPUT_FILE"

    srun --exclusive --cpu-bind=cores --distribution=block:block \
        -n 1 -c 1 \
        ./src/out.exe "$INPUT_FILE" \
        > "$STDOUT_FILE" 2> "$STDERR_FILE" &
    PIDS+=("$!")
done

STATUS=0
for PID in "${PIDS[@]}"; do
    if ! wait "$PID"; then
        STATUS=1
    fi
done

if [ "$STATUS" -ne 0 ]; then
    echo "at least one sample run failed; see $RESULT_DIR/run_*.err" >&2
    exit "$STATUS"
fi

python3 - "$RESULT_DIR" "$RESULT_DIR/yokota_green_kubo_average.dat" <<'PY'
from __future__ import annotations

import collections
import pathlib
import sys

result_dir = pathlib.Path(sys.argv[1])
output_path = pathlib.Path(sys.argv[2])

acc = collections.defaultdict(lambda: [0.0, 0.0])
meta = {}

for path in sorted(result_dir.glob("yokota_green_kubo_[0-9][0-9][0-9].dat")):
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 6:
                raise RuntimeError(f"unexpected row in {path}: {line.rstrip()}")

            nsample = int(fields[0])
            tau = float(fields[1])
            mode_index = int(fields[2])
            kx = float(fields[3])
            ky = float(fields[4])
            s_mean = float(fields[5])

            key = (tau, mode_index)
            if key in meta and meta[key] != (kx, ky):
                raise RuntimeError(f"inconsistent wave vector for key {key}")
            meta[key] = (kx, ky)
            acc[key][0] += nsample * s_mean
            acc[key][1] += nsample

if not acc:
    raise RuntimeError(f"no yokota_green_kubo sample files found in {result_dir}")

with output_path.open("w") as out:
    out.write("# nsample tau mode_index kx ky S_mean\n")
    for tau, mode_index in sorted(acc):
        weighted_sum, weight = acc[(tau, mode_index)]
        kx, ky = meta[(tau, mode_index)]
        out.write(
            f"{int(weight)} "
            f"{tau:.16e} "
            f"{mode_index:d} "
            f"{kx:.16e} "
            f"{ky:.16e} "
            f"{weighted_sum / weight:.16e}\n"
        )

print(f"wrote {output_path} with total nsample={int(next(iter(acc.values()))[1])}")
PY

python3 "$EXAMPLE_DIR/estimate_viscosity.py" \
    --input "$RESULT_DIR/yokota_green_kubo_average.dat" \
    --input-script "$TEMPLATE" \
    --output-dir "$RESULT_DIR/analysis" \
    --eta 0.1

echo "average    : $RESULT_DIR/yokota_green_kubo_average.dat"
echo "analysis   : $RESULT_DIR/analysis"

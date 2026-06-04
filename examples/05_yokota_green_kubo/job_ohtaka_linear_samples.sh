#!/bin/bash
#SBATCH -p i8cpu
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -c 1
#SBATCH -t 00:30:00
#SBATCH -J ykgk-linear
#SBATCH -o examples/05_yokota_green_kubo/results/ykgk-linear-%j.out
#SBATCH -e examples/05_yokota_green_kubo/results/ykgk-linear-%j.err

set -eu

module purge
module load oneapi_compiler/2023.0.0 oneapi_mpi/2023.0.0

export LD_LIBRARY_PATH=$HOME/local/heffte-oneapi-fftw/lib:$HOME/local/fftw-oneapi/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "$HOME/spectral-hohenberg-halperin-dynamics-2d"

tag=${YKGK_RUN_ID:-prepared}
run_dir=examples/05_yokota_green_kubo/runs/ohtaka_linear_$tag
result_dir=examples/05_yokota_green_kubo/results/ohtaka_linear_$tag

for sample in $(seq 0 127); do
    input=$(printf "%s/input_%03d.script" "$run_dir" "$sample")
    stdout=$(printf "%s/run_%03d.out" "$result_dir" "$sample")
    stderr=$(printf "%s/run_%03d.err" "$result_dir" "$sample")

    srun --exclusive --mem-per-cpu=1840 --cpu-bind=cores -n 1 -c 1 -N 1 \
    bash -c 'echo "start $(date +%s) sample=$1" >&2; ./src/out.exe "$2"; status=$?; echo "end $(date +%s) sample=$1 status=$status" >&2; exit "$status"' \
    _ "$sample" "$input" > "$stdout" 2> "$stderr" &
done

wait

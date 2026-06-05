#!/bin/bash
#SBATCH -p i8cpu
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -c 1
#SBATCH -t 00:30:00
#SBATCH -J shhd02-visc
#SBATCH -o /work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/02_sine_force_viscosity/%x-%j.out
#SBATCH -e /work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/02_sine_force_viscosity/%x-%j.err

set -eu

module purge
module load oneapi_compiler/2023.0.0 oneapi_mpi/2023.0.0

export LD_LIBRARY_PATH=$HOME/local/heffte-oneapi-fftw/lib:$HOME/local/fftw-oneapi/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export MKL_DYNAMIC=FALSE
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=core

repo=/home/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d
work_base=/work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/02_sine_force_viscosity
case_name=incompressible_nonlinear
case_dir=examples/02_sine_force_viscosity/$case_name
output_example=$work_base/$case_name

cd "$repo"
mkdir -p "$output_example"

samples=16
tasks_per_sample=$((128 / samples))

if [ $((128 % samples)) -ne 0 ]; then
    echo "samples must divide 128" >&2
    exit 1
fi

dt=0.01
run_time=10000.0
time_series_dtout=10.0
profile_dtout=1.0
profile_block_time=10.0

SHHD_EXAMPLE_ROOT="$repo/$case_dir" \
SHHD_OUTPUT_ROOT="$output_example" \
python3 "$case_dir/prepare_ohtaka_input.py" \
    --samples "$samples" \
    --dt "$dt" \
    --run-time "$run_time" \
    --time-series-dtout "$time_series_dtout" \
    --profile-dtout "$profile_dtout" \
    --profile-block-time "$profile_block_time"

for sample in $(seq 0 $((samples - 1))); do
    sample_id=$(printf "%03d" "$sample")
    input=$output_example/runs/input_${sample_id}.script
    stdout=$output_example/results/stdout_${sample_id}.log
    stderr=$output_example/results/stderr_${sample_id}.log

    srun --exclusive --mem-per-cpu=1840 --cpu-bind=cores -n "$tasks_per_sample" -c 1 -N 1 ./src/out.exe "$input" > "$stdout" 2> "$stderr" &
done

wait

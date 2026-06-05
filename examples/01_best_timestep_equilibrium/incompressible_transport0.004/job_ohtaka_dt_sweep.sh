#!/bin/bash
#SBATCH -p i8cpu
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -c 1
#SBATCH -t 00:30:00
#SBATCH -J shhd01-dt
#SBATCH -o /work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/01_best_timestep_equilibrium/%x-%j.out
#SBATCH -e /work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/01_best_timestep_equilibrium/%x-%j.err

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
work_base=/work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/01_best_timestep_equilibrium
case_name=incompressible_transport0.004
case_dir=examples/01_best_timestep_equilibrium/$case_name
output_example=$work_base/$case_name

cd "$repo"
mkdir -p "$output_example"

samples=16
tasks_per_sample=$((128 / samples))

if [ $((128 % samples)) -ne 0 ]; then
    echo "samples must divide 128" >&2
    exit 1
fi

tmax=200.0
dtout=0.40
dts="0.001 0.002 0.005 0.01 0.02 0.04 0.08"

for dt in $dts; do
    group=dt${dt}

    SHHD_EXAMPLE_ROOT="$repo/$case_dir" \
    SHHD_OUTPUT_ROOT="$output_example" \
    python3 "$case_dir/prepare_ohtaka_input.py" \
        --samples "$samples" \
        --dt "$dt" \
        --tmax "$tmax" \
        --dtout "$dtout"

    for sample in $(seq 0 $((samples - 1))); do
        sample_id=$(printf "%03d" "$sample")
        input=$output_example/runs/$group/input_${sample_id}.script
        stdout=$output_example/results/$group/stdout_${sample_id}.log
        stderr=$output_example/results/$group/stderr_${sample_id}.log

        srun --exclusive --mem-per-cpu=1840 --cpu-bind=cores -n "$tasks_per_sample" -c 1 -N 1 ./src/out.exe "$input" > "$stdout" 2> "$stderr" &
    done

    wait
done

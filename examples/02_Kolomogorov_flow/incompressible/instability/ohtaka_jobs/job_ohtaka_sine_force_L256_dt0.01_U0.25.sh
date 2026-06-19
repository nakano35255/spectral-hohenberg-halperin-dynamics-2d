#!/bin/bash
#SBATCH -p F72cpu
#SBATCH -N 72
#SBATCH -n 9216
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH -J shhd02-L256-dt001
#SBATCH -o /work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/02_kolomogorov_flow/incompressible/instability/ohtaka_jobs/%x-%j.out
#SBATCH -e /work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/02_kolomogorov_flow/incompressible/instability/ohtaka_jobs/%x-%j.err

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
work_base=/work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/02_kolomogorov_flow/incompressible/instability/raw_data
case_dir=examples/02_kolomogorov_flow/incompressible
job_dir=examples/02_kolomogorov_flow/incompressible/instability/ohtaka_jobs
run_name=eta0_0.1_U0.25_dt0.01
output_root=$work_base/$run_name

cd "$repo"
mkdir -p "$output_root"

samples=16
total_tasks=9216
nks=(1 2 3 4 5 6 8 10 12 16 20 24 28 32 36 40 48 64)

# Force amplitudes chosen from the FNS estimate so that the target steady
# velocity amplitude is approximately U=0.25 for L=256.
amplitudes=(
    0.000072324092
    0.00027144877
    0.00058599503
    0.0010093680
    0.0015367204
    0.0021641624
    0.0037066351
    0.0056150852
    0.0078719295
    0.013372873
    0.020109385
    0.027999134
    0.036971134
    0.046962310
    0.057915420
    0.069777689
    0.096035544
    0.15745869
)

if [ "${#nks[@]}" -ne "${#amplitudes[@]}" ]; then
    echo "nks and amplitudes must have the same length" >&2
    exit 1
fi

cases=${#nks[@]}
tasks_per_sample=$((total_tasks / (cases * samples)))

if [ $((total_tasks % (cases * samples))) -ne 0 ]; then
    echo "cases * samples must divide total_tasks" >&2
    exit 1
fi

if [ "$tasks_per_sample" -ne 32 ]; then
    echo "expected 32 MPI tasks per sample, got $tasks_per_sample" >&2
    exit 1
fi

dt=0.01
run_time=50000.0
time_series_dtout=10.0
profile_dtout=1.0
profile_block_time=100.0

for index in "${!nks[@]}"; do
    nk=${nks[$index]}
    amplitude=${amplitudes[$index]}
    nk_label=$(printf "nk_%03d" "$nk")
    case_output=$output_root/$nk_label

    SHHD_EXAMPLE_ROOT="$repo/$case_dir" \
    SHHD_OUTPUT_ROOT="$case_output" \
    python3 "$job_dir/prepare_ohtaka_input.py" \
        --samples "$samples" \
        --grid 256 256 \
        --length 256.0 256.0 \
        --dt "$dt" \
        --force-nk "$nk" \
        --force-amplitude "$amplitude" \
        --run-time "$run_time" \
        --time-series-dtout "$time_series_dtout" \
        --profile-dtout "$profile_dtout" \
        --profile-block-time "$profile_block_time"
done

for nk in "${nks[@]}"; do
    nk_label=$(printf "nk_%03d" "$nk")
    case_output=$output_root/$nk_label

    for sample in $(seq 0 $((samples - 1))); do
        sample_id=$(printf "%03d" "$sample")
        input=$case_output/runs/input_${sample_id}.script
        stdout=$case_output/results/stdout_${sample_id}.log
        stderr=$case_output/results/stderr_${sample_id}.log

        srun --exclusive --mem-per-cpu=1840 --cpu-bind=cores -n "$tasks_per_sample" -c 1 -N 1 ./src/out.exe "$input" > "$stdout" 2> "$stderr" &
    done
done

wait

#!/bin/bash
#SBATCH -p F72cpu
#SBATCH -N 72
#SBATCH -n 9216
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH -J shhd05-op-sine-D012
#SBATCH -o /work/k0565/k056500/spectral-hohenberg-halperin-dynamics-2d/examples/05_order_parameter_kolmogorov_flow/ohtaka_jobs/%x-%j.out
#SBATCH -e /work/k0565/k056500/spectral-hohenberg-halperin-dynamics-2d/examples/05_order_parameter_kolmogorov_flow/ohtaka_jobs/%x-%j.err

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

repo=/home/k0565/k056500/spectral-hohenberg-halperin-dynamics-2d
case_dir=examples/05_order_parameter_kolmogorov_flow
job_dir=$case_dir/ohtaka_jobs
work_base=/work/k0565/k056500/spectral-hohenberg-halperin-dynamics-2d/examples/05_order_parameter_kolmogorov_flow/raw_data
run_name=D0_0p12_targetpsi0p1_grid128_L4096_dt4_T20000000_n72
output_root=$work_base/$run_name

cd "$repo"
mkdir -p "$output_root"

samples=72
total_tasks=9216
nks=(1 2 3 4 5 6 8 10 12 16 20 24 28 32 40 48)

cases=${#nks[@]}
tasks_per_sample=8

active_tasks=$((cases * samples * tasks_per_sample))

if [ "$active_tasks" -ne "$total_tasks" ]; then
    echo "active MPI tasks must equal allocation: $active_tasks != $total_tasks" >&2
    exit 1
fi

dt=4.0
D0=0.12
eta0=0.12
target_amplitude=0.1
run_time=20000000.0
time_series_dtout=1000.0
profile_nevery=100
profile_nblock=10000

for nk in "${nks[@]}"; do
    python3 "$job_dir/prepare_ohtaka_input.py" \
        --output-root "$output_root" \
        --samples "$samples" \
        --grid 128 128 \
        --length 4096 4096 \
        --dt "$dt" \
        --D0 "$D0" \
        --eta "$eta0" \
        --force-nk "$nk" \
        --target-amplitude "$target_amplitude" \
        --run-time "$run_time" \
        --time-series-dtout "$time_series_dtout" \
        --profile-nevery "$profile_nevery" \
        --profile-nblock "$profile_nblock"
done

pids=()

for nk in "${nks[@]}"; do
    nk_label=$(printf "nk_%03d" "$nk")
    case_output=$output_root/$nk_label

    for sample in $(seq 0 $((samples - 1))); do
        sample_id=$(printf "%03d" "$sample")
        input=$case_output/runs/input_${sample_id}.script
        stdout=$case_output/logs/stdout_${sample_id}.log
        stderr=$case_output/logs/stderr_${sample_id}.log

        srun --exclusive --mem-per-cpu=1840 --cpu-bind=cores -n "$tasks_per_sample" -c 1 -N 1 ./src/out.exe "$input" > "$stdout" 2> "$stderr" &
        pids+=($!)
    done
done

status=0
for pid in "${pids[@]}"; do
    wait "$pid" || status=1
done
exit "$status"

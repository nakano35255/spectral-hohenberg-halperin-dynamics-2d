#!/bin/bash
#SBATCH -p F72cpu
#SBATCH -N 72
#SBATCH -n 9216
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH -J shhd04-Dlow-dt8
#SBATCH -o /work/k0565/k056500/spectral-hohenberg-halperin-dynamics-2d/examples/04_ness_uniform_gradient/energetics/ohtaka_jobs/%x-%j.out
#SBATCH -e /work/k0565/k056500/spectral-hohenberg-halperin-dynamics-2d/examples/04_ness_uniform_gradient/energetics/ohtaka_jobs/%x-%j.err

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
case_dir=examples/04_ness_uniform_gradient/energetics
work_base=/work/k0565/k056500/spectral-hohenberg-halperin-dynamics-2d/examples/04_ness_uniform_gradient/energetics/raw_data
run_name=production_D0_0p008_0p015_dt8_T200000000_n288
output_root=$work_base/$run_name

cd "$repo"
mkdir -p "$output_root"

D0s=(0.008 0.015)
dt=8.0
samples=288
total_tasks=9216
tasks_per_sample=16
run_time=200000000.0
time_series_dtout=1000.0

active_tasks=$((${#D0s[@]} * samples * tasks_per_sample))
if [ "$active_tasks" -gt "$total_tasks" ]; then
    echo "active MPI tasks exceed allocation: $active_tasks > $total_tasks" >&2
    exit 1
fi

for D0 in "${D0s[@]}"; do
    python3 "$case_dir/prepare_ohtaka_input.py" \
        --legacy-single-segment \
        --output-root "$output_root" \
        --samples "$samples" \
        --D0 "$D0" \
        --dt "$dt" \
        --run-time-per-segment "$run_time" \
        --time-series-dtout "$time_series_dtout"
done

pids=()

for D0 in "${D0s[@]}"; do
    case_path=D0_${D0}/dt_${dt}

    for sample in $(seq 0 $((samples - 1))); do
        sample_id=$(printf "%03d" "$sample")
        input=$output_root/runs/$case_path/input_${sample_id}.script
        stdout=$output_root/results/$case_path/stdout_${sample_id}.log
        stderr=$output_root/results/$case_path/stderr_${sample_id}.log

        srun --exclusive --mem-per-cpu=1840 --cpu-bind=cores -n "$tasks_per_sample" -c 1 -N 1 ./src/out.exe "$input" > "$stdout" 2> "$stderr" &
        pids+=($!)
    done
done

status=0
for pid in "${pids[@]}"; do
    wait "$pid" || status=1
done

exit "$status"

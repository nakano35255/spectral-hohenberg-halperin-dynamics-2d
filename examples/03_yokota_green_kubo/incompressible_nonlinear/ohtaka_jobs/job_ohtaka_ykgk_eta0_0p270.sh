#!/bin/bash
#SBATCH -p F72cpu
#SBATCH -N 72
#SBATCH -n 9216
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH -J shhd03-ykgk-e0270
#SBATCH -o /work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/03_yokota_green_kubo/incompressible_nonlinear/ohtaka_jobs/%x-%j.out
#SBATCH -e /work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/03_yokota_green_kubo/incompressible_nonlinear/ohtaka_jobs/%x-%j.err

set -eu

module purge
module load oneapi_compiler/2023.0.0 oneapi_mpi/2023.0.0

export HEFFTE_ROOT=$HOME/local/heffte-oneapi-fftw
export FFTW_ROOT=$HOME/local/fftw-oneapi
export LD_LIBRARY_PATH=$HEFFTE_ROOT/lib:$FFTW_ROOT/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export MKL_DYNAMIC=FALSE
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=core

repo=${REPO_ROOT:-/home/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d}
job_dir=examples/03_yokota_green_kubo/incompressible_nonlinear/ohtaka_jobs
work_base=${WORK_BASE:-/work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/03_yokota_green_kubo/incompressible_nonlinear/raw_data}
run_name=${RUN_NAME:-eta0_0p270_grid256_L256_dt0p01_T25000_diag_n576}
output_root=$work_base/$run_name

cd "$repo"
mkdir -p "$output_root" "$work_base"

eta0=${ETA0:-0.270}
dt=${DT:-0.01}
run_time=${RUN_TIME:-25000.0}
ykgk_dtout=${YKGK_DTOUT:-100.0}
ykgk_block_time=${YKGK_BLOCK_TIME:-5000.0}
samples=${SAMPLES:-576}
sample_offset=${SAMPLE_OFFSET:-0}
total_tasks=9216
tasks_per_sample=${TASKS_PER_SAMPLE:-16}
base_seed=${BASE_SEED:-270123}

active_tasks=$((samples * tasks_per_sample))
if [ "$active_tasks" -gt "$total_tasks" ]; then
    echo "active MPI tasks exceed allocation: $active_tasks > $total_tasks" >&2
    exit 1
fi

python3 "$job_dir/prepare_ohtaka_inputs.py" \
    --output-root "$output_root" \
    --samples "$samples" \
    --sample-offset "$sample_offset" \
    --eta "$eta0" \
    --dt "$dt" \
    --run-time "$run_time" \
    --ykgk-dtout "$ykgk_dtout" \
    --ykgk-block-time "$ykgk_block_time" \
    --seed "$base_seed"

pids=()
for sample in $(seq 0 $((samples - 1))); do
    global_sample=$((sample_offset + sample))
    sample_id=$(printf "%03d" "$global_sample")
    input=$output_root/runs/input_${sample_id}.script
    stdout=$output_root/results/stdout_${sample_id}.log
    stderr=$output_root/results/stderr_${sample_id}.log

    srun --exclusive --mem-per-cpu=1840 --cpu-bind=cores -n "$tasks_per_sample" -c 1 -N 1 ./src/out.exe "$input" > "$stdout" 2> "$stderr" &
    pids+=($!)
done

status=0
for pid in "${pids[@]}"; do
    wait "$pid" || status=1
done

exit "$status"

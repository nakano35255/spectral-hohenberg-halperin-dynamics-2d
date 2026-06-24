#!/bin/bash
#SBATCH -p i8cpu
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH -J shhd05-relax
#SBATCH -o examples/05_order_parameter_relaxation/ohtaka_jobs/%x-%j.out
#SBATCH -e examples/05_order_parameter_relaxation/ohtaka_jobs/%x-%j.err

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
case_dir=examples/05_order_parameter_relaxation
job_dir=$case_dir/ohtaka_jobs
work_base=${WORK_BASE:-/work/i0019/i001900/spectral-hohenberg-halperin-dynamics-2d/examples/05_order_parameter_relaxation/raw_data}
run_name=${RUN_NAME:-sine2d_grid256_L8192_dt4_T200000}
output_root=$work_base/$run_name

samples=${SAMPLES:-16}
sample_offset=${SAMPLE_OFFSET:-0}
tasks_per_sample=${TASKS_PER_SAMPLE:-8}
total_tasks=${SLURM_NTASKS:-128}
build_before_run=${BUILD_BEFORE_RUN:-0}

grid_x=${GRID_X:-256}
grid_y=${GRID_Y:-256}
length_x=${LENGTH_X:-8192}
length_y=${LENGTH_Y:-8192}
dt=${DT:-4.0}
run_time=${RUN_TIME:-200000.0}
snapshot_dtout=${SNAPSHOT_DTOUT:-4000.0}
thermo_dtout=${THERMO_DTOUT:-20000.0}
eta=${ETA:-0.12}
mobility=${MOBILITY:-0.12}
zeta=${ZETA:-0.0}
kBT=${KBT:-1.0}
amplitude=${AMPLITUDE:-0.1}
nkx=${NKX:-1}
nky=${NKY:-1}

active_tasks=$((samples * tasks_per_sample))
if [ "$active_tasks" -gt "$total_tasks" ]; then
    echo "active MPI tasks exceed allocation: $active_tasks > $total_tasks" >&2
    exit 1
fi

cd "$repo"
mkdir -p "$output_root"

if [ "$build_before_run" -ne 0 ]; then
    make -f Makefile.ohtaka -j 8
fi

python3 "$job_dir/prepare_ohtaka_inputs.py" \
    --output-root "$output_root" \
    --samples "$samples" \
    --sample-offset "$sample_offset" \
    --grid "$grid_x" "$grid_y" \
    --length "$length_x" "$length_y" \
    --dt "$dt" \
    --run-time "$run_time" \
    --snapshot-dtout "$snapshot_dtout" \
    --thermo-dtout "$thermo_dtout" \
    --eta "$eta" \
    --mobility "$mobility" \
    --zeta "$zeta" \
    --kBT "$kBT" \
    --order-parameter-amplitude "$amplitude" \
    --nkx "$nkx" \
    --nky "$nky"

pids=()
for local_sample in $(seq 0 $((samples - 1))); do
    sample=$((sample_offset + local_sample))
    sample_id=$(printf "%03d" "$sample")
    input=$output_root/runs/input_${sample_id}.script
    stdout=$output_root/logs/stdout_${sample_id}.log
    stderr=$output_root/logs/stderr_${sample_id}.log

    srun --exclusive --mem-per-cpu=1840 --cpu-bind=cores -n "$tasks_per_sample" -c 1 -N 1 ./src/out.exe "$input" > "$stdout" 2> "$stderr" &
    pids+=($!)
done

status=0
for pid in "${pids[@]}"; do
    wait "$pid" || status=1
done

exit "$status"

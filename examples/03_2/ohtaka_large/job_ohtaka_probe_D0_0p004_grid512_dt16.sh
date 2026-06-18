#!/bin/bash
#SBATCH -p F72cpu
#SBATCH -N 72
#SBATCH -n 9216
#SBATCH -c 1
#SBATCH -t 01:00:00
#SBATCH -J shhd03-2-probe512-D0004
#SBATCH -o /home/k0565/k056500/spectral-hohenberg-halperin-dynamics-2d/examples/03_2/ohtaka_large/%x-%j.out
#SBATCH -e /home/k0565/k056500/spectral-hohenberg-halperin-dynamics-2d/examples/03_2/ohtaka_large/%x-%j.err

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

repo=/home/k0565/k056500/spectral-hohenberg-halperin-dynamics-2d
case_dir=examples/03_2/ohtaka_large
work_base=/work/k0565/k056500/spectral-hohenberg-halperin-dynamics-2d/examples/03_2/ohtaka_large
run_name=${RUN_NAME:-f72_probe_D0_0p004_grid512_L16384_dt16}
output_root=$work_base/$run_name

cd "$repo"
mkdir -p "$output_root" "$work_base"

if [ "${BUILD_BEFORE_RUN:-0}" = "1" ]; then
    make -f Makefile.ohtaka yes-PASSIVE-SCALAR
    make -f Makefile.ohtaka clean
    make -f Makefile.ohtaka -j 8
fi

relax_time=${RELAX_TIME:-51200.0}
time_series_dtout=${TIME_SERIES_DTOUT:-1024.0}
replicas=${REPLICAS:-36}
tasks_per_run=${TASKS_PER_RUN:-256}
nodes_per_run=${NODES_PER_RUN:-2}
total_tasks=9216
total_nodes=72
base_seed=${BASE_SEED:-12345}
replica_ids=${REPLICA_IDS:-}

if [ "$tasks_per_run" -gt 385 ]; then
    echo "TASKS_PER_RUN=$tasks_per_run exceeds grid512 spectral x size 385." >&2
    exit 1
fi

if [ -z "$replica_ids" ]; then
    replica_ids=$(seq 0 $((replicas - 1)))
fi

active_replicas=0
for replica in $replica_ids; do
    active_replicas=$((active_replicas + 1))
done

active_tasks=$((active_replicas * tasks_per_run))
active_nodes=$((active_replicas * nodes_per_run))
if [ "$active_tasks" -gt "$total_tasks" ]; then
    echo "active MPI tasks exceed allocation: $active_tasks > $total_tasks" >&2
    exit 1
fi
if [ "$active_nodes" -gt "$total_nodes" ]; then
    echo "active nodes exceed allocation: $active_nodes > $total_nodes" >&2
    exit 1
fi

for replica in $replica_ids; do
    replica_num=$((10#$replica))
    rid=$(printf "%03d" "$replica_num")
    replica_root=$output_root/replica_${rid}
    seed=$((base_seed + 1000000 * replica_num))

    python3 "$case_dir/prepare_ohtaka_large_inputs.py" \
        --mode relax \
        --output-root "$replica_root" \
        --relax-segments 1 \
        --seed "$seed" \
        --D0 0.004 \
        --eta 0.004 \
        --dt 16.0 \
        --grid 512 512 \
        --length 16384 16384 \
        --gradient-amplitude 0.0001220703125 \
        --relax-time-per-segment "$relax_time" \
        --time-series-dtout "$time_series_dtout"
done

run_replica() {
    replica=$1
    replica_num=$((10#$replica))
    rid=$(printf "%03d" "$replica_num")
    replica_root=$output_root/replica_${rid}
    input=$replica_root/runs/input_relax_001.script
    stdout=$replica_root/logs/stdout_relax_001.log
    stderr=$replica_root/logs/stderr_relax_001.log
    final_restart=$replica_root/restarts/relax_001.restart
    tmp_restart=${final_restart}.tmp

    if [ -s "$final_restart" ]; then
        echo "skip completed replica ${rid} probe run: $final_restart"
        return 0
    fi

    rm -f "$tmp_restart"
    echo "run replica ${rid} grid512 probe"
    srun --exclusive --mem-per-cpu=1840 --cpu-bind=cores --distribution=block:block \
        -n "$tasks_per_run" -c 1 -N "$nodes_per_run" \
        ./src/out.exe "$input" > "$stdout" 2> "$stderr"
    mv "$tmp_restart" "$final_restart"
}

pids=()
for replica in $replica_ids; do
    run_replica "$replica" &
    pids+=($!)
done

status=0
for pid in "${pids[@]}"; do
    wait "$pid" || status=1
done

exit "$status"

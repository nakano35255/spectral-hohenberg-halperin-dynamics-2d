#!/bin/bash
#SBATCH -p F72cpu
#SBATCH -N 72
#SBATCH -n 9216
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH -J shhd03-2-budget-D0004
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
run_name=relax_D0_0p004_grid1024_L32768_dt16
output_root=$work_base/$run_name

cd "$repo"
mkdir -p "$output_root" "$work_base"

if [ "${BUILD_BEFORE_RUN:-0}" = "1" ]; then
    make -f Makefile.ohtaka yes-PASSIVE-SCALAR
    make -f Makefile.ohtaka clean
    make -f Makefile.ohtaka -j 8
fi

restart_index=${RESTART_INDEX:-1}
budget_time=${BUDGET_TIME:-100000000.0}
replicas=${REPLICAS:-18}
tasks_per_run=${TASKS_PER_RUN:-512}
nodes_per_run=${NODES_PER_RUN:-4}
total_tasks=9216
total_nodes=72
base_seed=${BASE_SEED:-12345}
replica_ids=${REPLICA_IDS:-}

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
        --mode budget \
        --output-root "$replica_root" \
        --restart-index "$restart_index" \
        --seed "$seed" \
        --D0 0.004 \
        --eta 0.004 \
        --dt 16.0 \
        --grid 1024 1024 \
        --length 32768 32768 \
        --gradient-amplitude 0.00006103515625 \
        --budget-time "$budget_time" \
        --time-series-dtout 16384.0 \
        --budget-nevery 20 \
        --budget-nblock 200
done

run_replica() {
    replica=$1
    replica_num=$((10#$replica))
    rid=$(printf "%03d" "$replica_num")
    restart_id=$(printf "%03d" "$restart_index")
    replica_root=$output_root/replica_${rid}
    source_restart=$replica_root/restarts/relax_${restart_id}.restart
    final_restart=$replica_root/restarts/budget_from_${restart_id}.restart
    tmp_restart=${final_restart}.tmp
    input=$replica_root/runs/input_budget_from_${restart_id}.script
    stdout=$replica_root/logs/stdout_budget_from_${restart_id}.log
    stderr=$replica_root/logs/stderr_budget_from_${restart_id}.log

    if [ ! -s "$source_restart" ]; then
        echo "missing source restart: $source_restart" >&2
        return 1
    fi

    if [ -s "$final_restart" ]; then
        echo "skip completed replica ${rid} budget run: $final_restart"
        return 0
    fi

    rm -f "$tmp_restart"
    echo "run replica ${rid} budget from relax restart ${restart_id}"
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

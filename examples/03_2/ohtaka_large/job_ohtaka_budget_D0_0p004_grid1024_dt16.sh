#!/bin/bash
#SBATCH -p F72cpu
#SBATCH -N 4
#SBATCH -n 512
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

restart_index=${RESTART_INDEX:-8}
budget_time=100000000.0
tasks_per_run=512
nodes_per_run=4

rid=$(printf "%03d" "$restart_index")
source_restart=$output_root/restarts/relax_${rid}.restart
final_restart=$output_root/restarts/budget_from_${rid}.restart
tmp_restart=${final_restart}.tmp
input=$output_root/runs/input_budget_from_${rid}.script
stdout=$output_root/logs/stdout_budget_from_${rid}.log
stderr=$output_root/logs/stderr_budget_from_${rid}.log

if [ ! -s "$source_restart" ]; then
    echo "missing source restart: $source_restart" >&2
    exit 1
fi

if [ -s "$final_restart" ]; then
    echo "skip completed budget run: $final_restart"
    exit 0
fi

python3 "$case_dir/prepare_ohtaka_large_inputs.py" \
    --mode budget \
    --output-root "$output_root" \
    --restart-index "$restart_index" \
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

rm -f "$tmp_restart"
echo "run budget from relax restart ${rid}"
srun --exclusive --mem-per-cpu=1840 --cpu-bind=cores --distribution=block:block \
    -n "$tasks_per_run" -c 1 -N "$nodes_per_run" \
    ./src/out.exe "$input" > "$stdout" 2> "$stderr"
mv "$tmp_restart" "$final_restart"

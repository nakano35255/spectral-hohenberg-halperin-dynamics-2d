#!/bin/bash
#SBATCH -p F72cpu
#SBATCH -N 4
#SBATCH -n 512
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH -J shhd03-2-relax-D0004
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

relax_segments=8
relax_time_per_segment=100000000.0
tasks_per_run=512
nodes_per_run=4

python3 "$case_dir/prepare_ohtaka_large_inputs.py" \
    --mode relax \
    --output-root "$output_root" \
    --relax-segments "$relax_segments" \
    --D0 0.004 \
    --eta 0.004 \
    --dt 16.0 \
    --grid 1024 1024 \
    --length 32768 32768 \
    --gradient-amplitude 0.00006103515625 \
    --relax-time-per-segment "$relax_time_per_segment" \
    --time-series-dtout 16384.0

for segment in $(seq 1 "$relax_segments"); do
    sid=$(printf "%03d" "$segment")
    prev_id=$(printf "%03d" $((segment - 1)))
    input=$output_root/runs/input_relax_${sid}.script
    stdout=$output_root/logs/stdout_relax_${sid}.log
    stderr=$output_root/logs/stderr_relax_${sid}.log
    final_restart=$output_root/restarts/relax_${sid}.restart
    tmp_restart=${final_restart}.tmp

    if [ -s "$final_restart" ]; then
        echo "skip completed relax segment ${sid}: $final_restart"
        continue
    fi

    if [ "$segment" -gt 1 ]; then
        prev_restart=$output_root/restarts/relax_${prev_id}.restart
        if [ ! -s "$prev_restart" ]; then
            echo "missing previous restart: $prev_restart" >&2
            exit 1
        fi
    fi

    rm -f "$tmp_restart"
    echo "run relax segment ${sid}"
    srun --exclusive --mem-per-cpu=1840 --cpu-bind=cores --distribution=block:block \
        -n "$tasks_per_run" -c 1 -N "$nodes_per_run" \
        ./src/out.exe "$input" > "$stdout" 2> "$stderr"
    mv "$tmp_restart" "$final_restart"
done

#!/bin/bash
#SBATCH -p F72cpu
#SBATCH -N 72
#SBATCH -n 9216
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH -J shhd04-rst-Sc14-D0004
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

cd "$repo"

D0=0.004
dt=16.0
samples=288
max_segments=2
total_tasks=9216
tasks_per_sample=16
run_time_per_segment=200000000.0
time_series_dtout=1024.0
sc_labels=(Sc1 Sc4)
sc_values=(1.0 4.0)

active_tasks=$((${#sc_values[@]} * samples * tasks_per_sample))
if [ "$active_tasks" -gt "$total_tasks" ]; then
    echo "active MPI tasks exceed allocation: $active_tasks > $total_tasks" >&2
    exit 1
fi

output_root_for() {
    sc_label=$1
    echo "$work_base/production_restart_${sc_label}_D0_0p004_dt16_T200000000x2_n288"
}

prepare_inputs() {
    for i in "${!sc_values[@]}"; do
        sc_label=${sc_labels[$i]}
        sc_value=${sc_values[$i]}
        output_root=$(output_root_for "$sc_label")
        mkdir -p "$output_root"

        python3 "$case_dir/prepare_ohtaka_input.py" \
            --output-root "$output_root" \
            --samples "$samples" \
            --segments "$max_segments" \
            --D0 "$D0" \
            --schmidt-number "$sc_value" \
            --dt "$dt" \
            --run-time-per-segment "$run_time_per_segment" \
            --time-series-dtout "$time_series_dtout"
    done
}

check_previous_restarts() {
    segment=$1
    if [ "$segment" -le 1 ]; then
        return 0
    fi

    prev_segment_id=$(printf "%03d" $((segment - 1)))
    case_path=D0_${D0}/dt_${dt}
    for sc_label in "${sc_labels[@]}"; do
        restart_dir=$(output_root_for "$sc_label")/restarts/$case_path
        for sample in $(seq 0 $((samples - 1))); do
            sample_id=$(printf "%03d" "$sample")
            prev_restart=$restart_dir/restart_${sample_id}_seg${prev_segment_id}.restart
            if [ ! -s "$prev_restart" ]; then
                echo "missing previous restart: $prev_restart" >&2
                return 1
            fi
        done
    done
}

run_segment() {
    segment=$1
    segment_id=$(printf "%03d" "$segment")
    case_path=D0_${D0}/dt_${dt}

    check_previous_restarts "$segment"

    pids=()
    for sc_label in "${sc_labels[@]}"; do
        output_root=$(output_root_for "$sc_label")
        run_dir=$output_root/runs/$case_path
        segment_dir=$output_root/segments/$case_path
        restart_dir=$output_root/restarts/$case_path
        result_dir=$output_root/results/$case_path

        for sample in $(seq 0 $((samples - 1))); do
            sample_id=$(printf "%03d" "$sample")
            input=$run_dir/input_${sample_id}_seg${segment_id}.script
            stdout=$result_dir/stdout_${sample_id}_seg${segment_id}.log
            stderr=$result_dir/stderr_${sample_id}_seg${segment_id}.log
            segment_series=$segment_dir/time_series_${sample_id}_seg${segment_id}.dat
            final_restart=$restart_dir/restart_${sample_id}_seg${segment_id}.restart
            tmp_restart=${final_restart}.tmp

            if [ -s "$final_restart" ] && [ -s "$segment_series" ]; then
                echo "skip completed ${sc_label} sample ${sample_id} segment ${segment_id}: $final_restart"
                continue
            fi

            rm -f "$tmp_restart"
            (
                srun --exclusive --mem-per-cpu=1840 --cpu-bind=cores -n "$tasks_per_sample" -c 1 -N 1 \
                    ./src/out.exe "$input" > "$stdout" 2> "$stderr"
                test -s "$tmp_restart"
                mv "$tmp_restart" "$final_restart"
            ) &
            pids+=($!)
        done
    done

    status=0
    for pid in "${pids[@]}"; do
        wait "$pid" || status=1
    done
    return "$status"
}

merge_time_series() {
    case_path=D0_${D0}/dt_${dt}
    for sc_label in "${sc_labels[@]}"; do
        output_root=$(output_root_for "$sc_label")
        segment_dir=$output_root/segments/$case_path
        result_dir=$output_root/results/$case_path
        mkdir -p "$result_dir"

        for sample in $(seq 0 $((samples - 1))); do
            sample_id=$(printf "%03d" "$sample")
            output=$result_dir/time_series_${sample_id}.dat
            tmp_output=${output}.tmp
            first=1
            : > "$tmp_output"
            for current_segment in $(seq 1 "$max_segments"); do
                current_segment_id=$(printf "%03d" "$current_segment")
                source=$segment_dir/time_series_${sample_id}_seg${current_segment_id}.dat
                if [ ! -s "$source" ]; then
                    echo "missing time series segment: $source" >&2
                    return 1
                fi
                if [ "$first" -eq 1 ]; then
                    cat "$source" >> "$tmp_output"
                    first=0
                else
                    awk 'substr($0,1,1) != "#"' "$source" >> "$tmp_output"
                fi
            done
            mv "$tmp_output" "$output"
        done
    done
}

prepare_inputs

for segment in $(seq 1 "$max_segments"); do
    echo "run segment ${segment}/${max_segments}: D0=$D0 Sc=${sc_values[*]} samples=$samples tasks_per_sample=$tasks_per_sample"
    run_segment "$segment"
done

merge_time_series

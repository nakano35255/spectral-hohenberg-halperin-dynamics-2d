#!/bin/bash
set -euo pipefail

segment=1
job_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cases=("0p008 08" "0p015 08" "0p03 04" "0p06 04" "0p25 04" "0p50 04" "1p00 04" "2p00 04")

for spec in "${cases[@]}"; do
    read -r d0 dt <<< "$spec"
    for replica in {0..7}; do
        qsub -v REPLICA_ID="$replica",RELAX_SEGMENT="$segment" \
            "$job_dir/job_kugui_relax_D0_${d0}_grid512_dt${dt}.pbs"
    done
done

#!/bin/bash
#SBATCH --job-name=NAARC-ZIP
#SBATCH --time=24:00:00
#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial
#SBATCH --mem=124G

. /work/n01/n01/jdha/cfgs/NAARC/scripts/env/archer2/gnu-mpich-ucx

set -euo pipefail

cd OUTPUTS || exit 1

roots=($(ls *_0000.nc | sed 's/_0000\.nc$//' | sort -u))

for fname in "${roots[@]}"; do
    echo "Processing ${fname}.nc"
    date

    log_file="rebuild_${fname}.log"
    rebuild_cmd="/work/n01/n01/jdha/working_mes/nemo_main/tools/REBUILD_NEMO/rebuild_nemo"

    if $rebuild_cmd -d 1 -x 256 -y 256 -z 75 -t 1 -p 24 "$fname" 96 > "$log_file" 2>&1; then
        echo "✅ Rebuild completed successfully for ${fname}.nc"
        mv ${fname}.nc ../OUTPUTS_ZIP/
        mv ${fname}_????.nc ../OUTPUTS_PROCESSED/
        echo "Moved files for ${fname}.nc"
    else
        echo "❌ Rebuild failed for ${fname}.nc — check ${log_file}"
    fi

    date
    echo "----------------------------------------"
done


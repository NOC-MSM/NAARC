#!/bin/bash
#SBATCH --job-name=NAARC-ZIP
#SBATCH --time=24:00:00
#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial
#SBATCH --mem=24G

. XXX_ENV_XXX

cd OUTPUTS

for fname in NAARC_1h_*_0000.nc NAARC_1d_*_0000.nc NAARC_5d_*_0000.nc NAARC_1m_*_0000.nc
do

/work/n01/n01/jdha/working_mes/nemo_main/tools/REBUILD_NEMO/rebuild_nemo -d 1 -x 256 -y 256 -z 75 -t 1 -p 12 ${fname/_0000.nc/} 16

mv ${fname/_0000.nc/}.nc ../OUTPUTS_ZIP/
mv ${fname/_0000.nc/}_????.nc ../OUTPUTS_PROCESSED/

echo "processed ${fname/_0000.nc/}.nc"

done

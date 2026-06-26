#!/bin/bash
#SBATCH --job-name=NAARC-ZIP
#SBATCH --time=24:00:00
#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial
#SBATCH --mem=124G

. XXX_ENV_XXX

if [ ! -d OUTPUTS_ZIP ]
then
	mkdir OUTPUTS_ZIP
fi

if [ ! -d OUTPUTS_PROCESSED ]
then
	mkdir OUTPUTS_PROCESSED
fi

cd OUTPUTS

for fname in NAARC_1d_* NAARC_1m_*
do
        if [[ $fname == *"grid_T"* ]]; then
	   nccopy -4 -d 1 -c x/256,y/256,deptht/75 $fname ../OUTPUTS_ZIP/$fname
	elif [[ $fname == *"grid_U"* ]]; then
           nccopy -4 -d 1 -c x/256,y/256,depthu/75 $fname ../OUTPUTS_ZIP/$fname
	elif [[ $fname == *"grid_V"* ]]; then
           nccopy -4 -d 1 -c x/256,y/256,depthv/75 $fname ../OUTPUTS_ZIP/$fname
	elif [[ $fname == *"grid_W"* ]]; then
           nccopy -4 -d 1 -c x/256,y/256,depthw/75 $fname ../OUTPUTS_ZIP/$fname
	fi

	mv $fname ../OUTPUTS_PROCESSED/

done


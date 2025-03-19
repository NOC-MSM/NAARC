# North Atlantic and ARCtic NEMO configuration

The setup script has been tested and will checkout, compile and run the NAARC (NEMO 4.2.2) code on: ARCHER2 for Cray-MPICH and GNU-MPICH, and Anemone for iFort. 

[Documentation](https://noc-msm.github.io/NAARC/)

## Quick Start:
On ARCHER2
```
git clone git@github.com:NOC-MSM/NAARC.git
./NAARC/scripts/setup/NAARC_setup -p $PWD/NAARC_RUNS  -r $PWD/NAARC -n 4.2.2 -x 2 -m archer2 -a mpich -c gnu
cd NAARC_RUNS/nemo/cfgs/NAARC//
cp -rP EXPREF EXP_MYRUN
cd EXP_MYRUN
ln -s ../INPUTS/domain_cfg_mes.nc domain_cfg.nc
```
or if using ANEMONE, replace use options:
```
-m anemone -a impi -c ifort
```
Edit the project code and options in  `runscript_continuous.slurm` then:
```
sbatch runscript.slurm -y 1979 -s 1
```
This will produce a 5 day mean output from the beginning of 1979. The run should take 15 minutes to complete once in the machine.

### Forcing data:

[NAARC](https://gws-access.jasmin.ac.uk/public/jmmp/NAARC/)

_this is automatically transferred when the setup script is executed_

For ARCHER2 users these data are held under `/work/n01/shared/NAARC` and `/work/n01/shared/nemo/FORCING` and are linked during the setup.

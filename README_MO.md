# North Atlantic and ARCtic NEMO configuration - MO project

The setup script has been tested and will checkout, compile and run the NAARC (NEMO 4.2.2) code on Anemone using ifort.


<img width="541" alt="Screenshot 2025-06-19 at 16 40 30" src="https://github.com/user-attachments/assets/1c681919-c59d-4750-a92e-9e0b9a0d5411" />

[Documentation](https://noc-msm.github.io/NAARC/)

## Quick Start:
On Anemone
```
git clone git@github.com:NOC-MSM/NAARC.git
./NAARC/scripts/setup/NAARC_setup -p $PWD/NAARC_MO  -r $PWD/NAARC -n 4.2.2 -x 2 -m aemone -a impi -c ifort
cd NAARC_MO/nemo/cfgs/NAARC/
cp -rP EXPREF EXP_MES
cd EXP_MES
ln -s ../INPUTS/domain_cfg_mes.nc domain_cfg.nc
```
Before submitting a job, ammend the namelist:
```
sed -i "s/sn_uoatm   =  'uos',/sn_voatm   =  'NOT USED',/" namelist_cfg_template
sed -i "s/sn_voatm   =  'vos',/sn_voatm   =  'NOT USED',/" namelist_cfg_template
```
Submit a test job:
```
sbatch runscript_MES_MO.slurm -y 1979 -s 1
```
This will produce a 1 day mean output from the beginning of 1979. The run should take 15 minutes to complete once in the machine.

You can create `EXP_MES_TIDE` and `EXP_ZPS` for the addtional experiments. Remember to link either the `domain_cfg_mes.nc` or `domain_cfg_zps.nc`
files before running. You will then have to `cp runscript_MES_MO.slurm runscript_MES_TIDE_MO.slurm` etc and edit according to the experiment you
are running:

_MES_
```
ln_tide='.false.'
ln_cdmin2d='.true.'
ln_loglayer='.true.'
sed -i "s/rn_boost =  2./rn_boost =  1./" namelist_cfg_template
```

_MES TIDE_
```
ln_tide='.true.'
ln_cdmin2d='.true.'
ln_loglayer='.true.'
```

_ZPS_
```
ln_tide='.false.'
ln_cdmin2d='.false.'
ln_loglayer='.false.'
sed -i "s/rn_boost =  2./rn_boost =  50./" namelist_cfg_template
```
I am still yet to optimise the mpp layout. You can have a play by changing the following in the runscript:

`SRUN_CMD='mpiexec.hydra -print-rank-map -ppn 1 -np 140 ./xios_server.exe : -np 1326 ./nemo'`

NB: at the moment the MES out-of-the-box case will run for 1d with no output. You will have to go into the `file_def*`
and switch the 1m/1mo to 1d and also `type="one_file"` to `type="mulitple_file"`.

### Forcing data:

[NAARC](https://gws-access.jasmin.ac.uk/public/jmmp/NAARC/)

_this is automatically transferred when the setup script is executed_

Or in the case of Aneome, linked to the following directories: `/dssgfs01/working/acc/FORCING/JRA` and `/dssgfs01/working/jdha/NAARC/4.2.2`

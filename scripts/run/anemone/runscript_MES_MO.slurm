#!/bin/bash
#SBATCH --job-name=NAARC
#SBATCH --partition=compute
#SBATCH --time=03:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=26	

echo $PWD

export OMP_NUM_THREADS=1
export I_MPI_SHM=off
unset I_MPI_PMI_LIBRARY

display_usage() {
   echo
   echo "  Runscript: NAARC on ARCHER2/Anemone"
   echo "  ***************************************"
   echo
   echo "  usage: ${0##*/} -y year -s spin"
   echo
   echo "  flags:  -y first year to be run"
   echo "          -s spin cycle (only valid for first year of the run)"
   echo
   exit 1
        }
# if less than two arguments supplied, display usage
        if [  $# -le 3 ]
        then
                display_usage
                exit 1
        fi
# if more than two arguments supplied, display usage
        if [  $# -ge 5 ]
        then
                display_usage
                exit 1
        fi
# check whether user had supplied -h or --help . If yes display usage
        if [[ ( $# == "--help") ||  $# == "-h" ]]
        then
                display_usage
                exit 0
        fi

while getopts y:s: option
  do
  case "${option}"
  in
  y) export YEAR=${OPTARG};;
  s) export SPIN=${OPTARG};;
  esac
done

echo $YEAR $SPIN 

#################### nemo runscript options ############################
# For info on the parameters see namelist_ref                          #
########################################################################
rn_rdt=300          ; ln_zps='.false.'       ; ln_tmx_itf='.false.'
ln_bt_auto='.true.' ; rn_bt_cmax=0.8         ; nn_baro=30
nn_mxlice=2         ; nn_z0_ice=1            ; ln_rnf_new='.false.'
ln_rstdate='.true.' ; ln_shlat2d='.true.'   ; nn_diaharm=1991
rn_Cd0=2.5e-3       ; ln_loglayer='.true.'   ; ln_tide='.false.'
ln_boost='.true.'  ; ln_gls='.true.'        ; ln_int_wave_drag='.false.' 
ln_hpg_djc='.true.' ; ln_calc_tdiss='.false.'; ln_crt_fbk='.false.'
ln_frm_drg='.false.'; rn_frm_ht0=0.0         ; ln_cdmin2d='.true.'
########################################################################


#################### runscript options #################################
# Controls for the simulation                                          #
########################################################################
ln_spin='.false.'   # Do we require a spin up phase for the ice
nn_spin=-1          # How many years spin up
year_st=1979        # Start year of the simulation (ramains constant)
year_en=1979        # End year of this job submission
offset=0            # If rdt has changed, an offset to nn_it000 is needed
nn_spin_cycle=$SPIN # Change to nn_spin+1 if year is .gt. year_st
year=$YEAR          # Change if starting job part way through simulation
########################################################################


########################################################################
#                No need to edit below this line                       #
########################################################################

if [ "$ln_zps" = ".false." ]
then

   #################### nemo runscript options #########################
   # The following options are for the MES domain_cfg.nc               #
   #####################################################################
   
   SRUN_CMD='mpiexec.hydra -print-rank-map -ppn 1 -np 140 ./xios_server.exe : -np 1326 ./nemo'

else

   #################### nemo runscript options #########################
   # The following options are for the ZPS domain_cfg.nc               #
   #####################################################################

   SRUN_CMD='mpiexec.hydra -print-rank-map -ppn 1 -np 16 ./xios_server.exe : -np 1326 ./nemo'

fi

if [ "$ln_frm_drg" = ".false." ]
then
   ln_cst_drg='.true.'
else
   ln_cst_drg='.false.'
fi
 
if [ "$ln_zps" = ".true." ]
then
   ln_trabbl='.true.'
   ln_traldf_triads='.false.'
   ln_traldf_hor='.false.'
   ln_traldf_iso='.true.'
   ln_hpg_prj='.false.'
   ln_hpg_sco='.true.'
   ln_hpg_djc='.false.' 
else
   ln_trabbl='.false.'
   ln_traldf_triads='.true.'
   ln_traldf_hor='.false.'
   ln_traldf_iso='.false.'
   if [ "$ln_hpg_djc" = ".true." ]
   then
      ln_hpg_prj='.false.'
   else
      ln_hpg_prj='.true.'
   fi
   ln_hpg_sco='.false.'
fi

if [ "$ln_tide" = ".true." ]
then
   nn_dyn2d_dta=3 
else
   nn_dyn2d_dta=1
fi

if [ "$ln_loglayer" = ".true." ]
then
   ln_non_lin='.false.'
else
   ln_non_lin='.true.'
fi

if [ "$ln_gls" = ".true." ]
then
   ln_tke='.false.'
   ln_zdfevd='.false.'
else
   ln_tke='.true.'
   ln_zdfevd='.true.'
fi

if [ "$ln_rnf_new" = ".true." ]
then
	sn_rnf="'eORCA12_rivers',24,'rorunoff',.false.,.false.,'yearly'"
	sn_cnf="'JRA55_eORCA12coef'   ,       -12.        , 'socoefr',   .false.   , .true. , 'yearly'"
else
	sn_rnf="'runoff_1m_nomask',-1,'sornficb',.true.,.true.,'yearly'"
	sn_cnf="'runoff_1m_nomask'   ,       -12.        , 'socoefr',   .false.   , .true. , 'yearly'"
fi

echo $SLURM_NNODES

. XXX_ENV_XXX

export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1

chmod u+x ./myscript_wrapper.sh

nn_jpni=-1
nn_jpnj=-1

nn_date0=$year_st\0101
nextyear=true

while [ $year -le $year_en ] && $nextyear
do

   echo $year

   # Array pretending to be a Pythonic dictionary {EXP_NAME:NZ}
   ARRAY=( "NAARC:75"
         )

   for ens in "${ARRAY[@]}"
   do

      nam="${ens%%:*}"
      jpk="${ens##*:}"

      if [ $year -eq $nn_diaharm ] && [ "$ln_tide" = ".true." ]
      then
	 ln_diaharm='.true.'
      else
	 ln_diaharm='.false.'
      fi

      if [ $year -eq $year_st ]
      then
         nn_it000=1
	 if [ $nn_spin_cycle -eq 1 ] || [ "$ln_spin" = ".false." ]
         then
            ln_rstart='.false.'	  
	    ln_tide_ramp='.true.'
	 else
            ln_rstart='.true.'	  
	    ln_tide_ramp='.false.'
         fi
	 nn_rstctl=0
         ln_reset_ts='.true.'
      else
         nn_it000=`./end_time_step $(( $year - 1 )) $rn_rdt $year_st` ; rs0=`printf "%08d" $nn_it000`; nn_it000=$(( $nn_it000 + 1 + $offset ))
         ln_rstart=".true."
	 nn_rstctl=2
	 ln_tide_ramp='.false.'
         ln_reset_ts='.false.'
      fi
       
      if [ "$ln_rstdate" = ".true." ]
      then
         if [ $year -eq $year_st ]
         then
            rs0=$(($year+1))\0101
         else
            rs0=$year\0101
         fi
      fi 

      nn_itend=`./end_time_step $year $rn_rdt $year_st`; nn_itend=$(( $nn_itend + $offset ))
      if [ $nn_spin_cycle -le $nn_spin ] && [ "$ln_spin" = ".true." ] 
      then
         sed    "s/XXX_EXP_XXX/$nam\_S$nn_spin_cycle/g" namelist_cfg_template > namelist_cfg
         sed -i "s/XXX_RS0_XXX/$nam\_S$(($nn_spin_cycle-1))_$rs0\_restart/g" namelist_cfg
         sed "s/XXX_RS0_XXX/$nam\_S$(($nn_spin_cycle-1))_$rs0\_restart_ice/g" namelist_ice_cfg_template > namelist_ice_cfg
      else
         sed    "s/XXX_EXP_XXX/$nam/g" namelist_cfg_template > namelist_cfg
	 if [ $year -eq $year_st ] && [ "$ln_spin" = ".true." ]
         then
            sed -i "s/XXX_RS0_XXX/$nam\_S$(($nn_spin_cycle-1))_$rs0\_restart/g" namelist_cfg
            sed "s/XXX_RS0_XXX/$nam\_S$(($nn_spin_cycle-1))_$rs0\_restart_ice/g" namelist_ice_cfg_template > namelist_ice_cfg
         else
            sed -i "s/XXX_RS0_XXX/$nam\_$rs0\_restart/g" namelist_cfg
            sed "s/XXX_RS0_XXX/$nam\_$rs0\_restart_ice/g" namelist_ice_cfg_template > namelist_ice_cfg
         fi
      fi
      sed -i "s/XXX_TST_XXX/$nn_it000/g" namelist_cfg
      sed -i "s/XXX_YEAR_XXX/$year/g" namelist_cfg
      sed -i "s/XXX_TEN_XXX/$nn_itend/g" namelist_cfg
      sed -i "s/XXX_RDT_XXX/$rn_rdt/g" namelist_cfg
      sed -i "s/XXX_RST_XXX/$ln_rstart/g" namelist_cfg
      sed -i "s/XXX_RCT_XXX/$nn_rstctl/g" namelist_cfg
      sed -i "s/XXX_RSD_XXX/$ln_rstdate/g" namelist_cfg
      sed -i "s/XXX_RTS_XXX/$ln_reset_ts/g" namelist_cfg
      sed -i "s/XXX_RMP_XXX/$ln_tide_ramp/g" namelist_cfg
      sed -i "s/XXX_BTA_XXX/$ln_bt_auto/g" namelist_cfg
      sed -i "s/XXX_BTC_XXX/$rn_bt_cmax/g" namelist_cfg
      sed -i "s/XXX_RNF_XXX/$sn_rnf/g" namelist_cfg
      sed -i "s/XXX_CNF_XXX/$sn_cnf/g" namelist_cfg
      sed -i "s/XXX_TID_XXX/$ln_tide/g" namelist_cfg
      sed -i "s/XXX_BTD_XXX/$nn_dyn2d_dta/g" namelist_cfg
      sed -i "s/XXX_ITF_XXX/$ln_tmx_itf/g" namelist_cfg
      sed -i "s/XXX_DHM_XXX/$ln_diaharm/g" namelist_cfg
      sed -i "s/XXX_BOO_XXX/$ln_boost/g" namelist_cfg
      sed -i "s/XXX_SHA_XXX/$ln_shlat2d/g" namelist_cfg
      sed -i "s/XXX_CD0_XXX/$rn_Cd0/g" namelist_cfg
      sed -i "s/XXX_GLS_XXX/$ln_gls/g" namelist_cfg
      sed -i "s/XXX_TKE_XXX/$ln_tke/g" namelist_cfg
      sed -i "s/XXX_EVD_XXX/$ln_zdfevd/g" namelist_cfg
      sed -i "s/XXX_BBL_XXX/$ln_trabbl/g" namelist_cfg
      sed -i "s/XXX_HOR_XXX/$ln_traldf_hor/g" namelist_cfg
      sed -i "s/XXX_ISO_XXX/$ln_traldf_iso/g" namelist_cfg
      sed -i "s/XXX_TRI_XXX/$ln_traldf_triads/g" namelist_cfg
      sed -i "s/XXX_PRJ_XXX/$ln_hpg_prj/g" namelist_cfg
      sed -i "s/XXX_DJC_XXX/$ln_hpg_djc/g" namelist_cfg
      sed -i "s/XXX_SCO_XXX/$ln_hpg_sco/g" namelist_cfg
      sed -i "s/XXX_LOG_XXX/$ln_loglayer/g" namelist_cfg
      sed -i "s/XXX_NLN_XXX/$ln_non_lin/g" namelist_cfg
      sed -i "s/XXX_Z0I_XXX/$nn_z0_ice/g" namelist_cfg
      sed -i "s/XXX_IMX_XXX/$nn_mxlice/g" namelist_cfg
      sed -i "s/XXX_PNI_XXX/$nn_jpni/g" namelist_cfg
      sed -i "s/XXX_PNJ_XXX/$nn_jpnj/g" namelist_cfg
      sed -i "s/XXX_TDG_XXX/$ln_int_wave_drag/g" namelist_cfg
      sed -i "s/XXX_CTD_XXX/$ln_calc_tdiss/g" namelist_cfg
      sed -i "s/XXX_FBK_XXX/$ln_crt_fbk/g" namelist_cfg
      sed -i "s/XXX_C2D_XXX/$ln_cdmin2d/g" namelist_cfg
      sed -i "s/XXX_FRM_XXX/$ln_frm_drg/g" namelist_cfg
      sed -i "s/XXX_CST_XXX/$ln_cst_drg/g" namelist_cfg
      sed -i "s/XXX_HT0_XXX/$rn_frm_ht0/g" namelist_cfg

   done

   for ens in "${ARRAY[@]}"
   do
      date
      eval $SRUN_CMD &
      echo $SRUN_CMD &
   #  ./time_step_chk $SLURM_JOB_ID $nn_itend &
   done

   wait
   date
   for ens in "${ARRAY[@]}"
   do
      if [ $year -eq $year_st ] && [ $nn_spin_cycle -le $nn_spin ] && [ "$ln_spin" = ".true." ]
      then
         suf=$year\_S$nn_spin_cycle
      else
         suf=$year
      fi
      #xp="${ens%%:*}"
      #cd $PBS_O_WORKDIR/../ENSEMBLE_MEMBERS/ENS_$xp
      mv *_??_*grid*.nc OUTPUTS &
      mv *_??_*icemod*.nc OUTPUTS &
      cp ocean.output meta_out/ocean.output.$suf
      cp namelist_cfg meta_out/namelist_cfg.$suf
      cp namelist_ice_cfg meta_out/namelist_ice_cfg.$suf
      cp run.stat meta_out/run.stat.$suf
      cp time.step meta_out/time.step.$suf
      mv volume_transport TRANSPORTS/volume_transport.$suf
      mv heat_transport TRANSPORTS/heat_transport.$suf
      mv salt_transport TRANSPORTS/salt_transport.$suf
   done
   wait

   echo "Moved output"

   current_stp=`sed -n 1,1p time.step`
   if [ ! $current_stp -eq $nn_itend ]
   then
      exit
   fi

   echo "current_stp"
   echo $current_stp
   echo $nn_itend

   if [ $nn_spin_cycle -gt $nn_spin ] || [ "$ln_spin" = ".false." ]
   then
      year=$(($year+1))
   else
      nn_spin_cycle=$(($nn_spin_cycle+1))
   fi

   nextyear=`./check_runtime $SLURM_JOB_ID $YEAR $year`

done

# chunk/zip data
echo "Zipping"
sbatch runzip.slurm

# if year -gt year_end no need to go again
if [ $year -gt $year_en ]
then
   exit
else
   echo 'Submitting runscript ' $year $nn_spin_cycle
   sbatch runscript.slurm -y $year -s $nn_spin_cycle
fi

exit

Added an option for NEMO version 5.0 in scripts/setup/NAARC_setup using the case
for NEMO_VER.
This was added to the catch case, externals and XIOS.
An error messgae for XIOS was added in the case that NEMO_VER is less 4.0.

Copied the NEMO version 4.2.2 files to 5.0 for the folders /arch, /CPP, /MY_SRC, EXPREF.

*Tired to compile but error:
E R R O R : key key_isf is not found in /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs//WORK routines...*

key-isf is in the NEMO v5 namelist_ref (in /cfgs/EXP00).
Removed key_isf from CPP/5.0/cpp_NAARC.fcm and CPP/5.0.1/cpp_NAARC.fcm.

*in_out_manager.F90:192:2: fatal error: do_loop_substitute.h90: No such file or directory
  192 |    !!----------------------------------------------------------------------
      |  ^ ~~~~~~~~~~~~~~~~~~~~~~
compilation terminated.
'cpp -Dkey_nosignedzero -P -traditional' -P -traditional -I/work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs//WORK -Dkey_xios -Dkey_qco -Dkey_si3 -I/work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs/BLD/inc in_out_manager.F90 failed (1) at /mnt/lustre/a2fs-work1/work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/ext/FCM/bin/../lib/Fcm/BuildSrc.pm line 741.*

Comment out !#  include "do_loop_substitute.h90" in MY_SRC/5.0/in_out_manager.F90 and MY_SRC/5.0.1/in_out_manager.F90

Commenting out the line wasn't the right approach so undid it.

Tried renaming WORK to WORK_old as Stefanie said it looked out of place and could be conflicting with v5. Same error persists.

Tried renaming MY_SRC/5.0 to MY_SRC/5.0_old.

"do_loop_substitute.h90" error not coming up. Instead:

*E R R O R : key key_xios is not found in /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs//WORK routines...*

James says remove key_xios from CPP/5.0/cpp_NAARC.fcm and CPP/5.0.1/cpp_NAARC.fcm.

*E R R O R : key key_qco is not found in /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs//WORK routines...*

Removed key_qco from CPP/5.0/cpp_NAARC.fcm and CPP/5.0.1/cpp_NAARC.fcm.

*E R R O R : key key_si3 is not found in /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs//WORK routines...*

Removed key_si3 from CPP/5.0/cpp_NAARC.fcm and CPP/5.0.1/cpp_NAARC.fcm.

*ERROR: /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs//BLD/bld.cfg: LINE 44:
       bld::pp::nemo: invalid sub-package in declaration.
ERROR: /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs//BLD/bld.cfg: LINE 47:
       bld::tool::fppflags::nemo: invalid sub-package in declaration.*

The problem was the compile command needs "-n NAARC" so "./makenemo -m archer2-gnu-mpich -r NAARC -n NAARC -j 16". Edited scripts/setup/NAARC_setup


Restored previous keys in CPP/5.0/cpp_NAARC.fcm and CPP/5.0.1/cpp_NAARC.fcm.
Added key_vco_3d in CPP/5.0/cpp_NAARC.fcm and CPP/5.0.1/cpp_NAARC.fcm.

NAARC compiled with v5 but does not run.

*Error with not finding domain_cfg.nc*

INPUTS was empty because a folder for 5.0 did not exist in /work/n01/shared/NAARC/ to link to. The folder /work/n01/shared/NAARC/4.2.2 is now linked to /work/n01/shared/NAARC/5.0 and /work/n01/shared/NAARC/5.0.1.

*Namelist related error:
misspelled variable in namelist namcfg (ref) iostat =  5010 
misspelled variable in namelist namtsd (ref) iostat =  5010 
misspelled variable in namelist namlbc (ref) iostat =  5010
misspelled variable in namelist namrun (ref) iostat =  5010
misspelled variable in namelist namdom (ref) iostat =  5010*

Changes have been made to translate the namelist_cfg_template and namelist_ice_cfg_template to v5 namelists by removing variables not present in the v5 namelist or changing the name of variables that changed. Values of variables have been left unchanged with the exeption of rn_alb_dpnd which was set to 0.30 but a new comment suggest that is outside the obs range 0.12 -- 0.25 so the new default 0.18 has be used for rn_alb_dpnd. 
rn_frm_ht0 was left out of the v5 namelist_cfg_template because ice form drag is included in nemo v5. This means the ice\*.F90 files probably won't need to be carried over in MY_SRC

Copied namelist_ref from nemo_v5 SHARED to EXPREF/5.0/ and EXPREF/5.0.1/. 

*Namelist related errors persits except namdom:
misspelled variable in namelist namsbc_blk (cfg) iostat =  5010
misspelled variable in namelist namrun (cfg) iostat =  5010
misspelled variable in namelist namlbc (cfg) iostat =  5010
misspelled variable in namelist namtsd (cfg) iostat =  5010
 File ./data_1m_potential_temperature_nomask.nc* not found
 File ./data_1m_salinity_nomask.nc* not found
 File ./data_1m_potential_temperature_nomask.nc* not found
 File ./data_1m_salinity_nomask.nc* not found
rn_tide_ramp_dt must be lower than run duration*

For namsbc_blk, ln_NCAR variable, work on MY_SRC files sbcblk.F90. Done.
sbcblk.F90 also needs fldread.F90, sbc_oce.F90, iom.F90. Done iom.F90 didn't need updating.
 
For namrun, ln_rstdate and ln_reset_ts variables, work on MY_SRC files domain.F90, icerst.F90, in_out_manager.F90, restart.F90, dtatsd.F90, istate.F90. Done, dtatsd.F90 needed the gdept() comment editing.
No additions needed.

For namlbc, ln_shlat2d, cn_shlat2d_file, cn_shlat2d_var variables, work on MY_SRC files dommsk.F90. Done.
No additons needed.

For namtsd, ln_tsd_interp, sn_dep, sn_msk variables, work on MY_SRC files dtatsd.F90, sbcrnf.F90. Done.
No additions needed.

*MY_SRC errors reduced to:
misspelled variable in namelist namsbc_blk (cfg) iostat =  5010
misspelled variable in namelist namrun (cfg) iostat =  5010
misspelled variable in namelist namrun (ref) iostat =  5010
rn_tide_ramp_dt must be lower than run duration*

Started converting other files but lraving out momentum bit for now:
diahth.F90. Done
diapea.F90 not in v5 src so copied from v4.2.2. Done
diawri.F90. Done
nemogcm.F90. Done
sbcrnf.F90. Done
sbcssm.F90. Done
tide.h90. Done
traldf.F90. Didn't need updating, no changes.
traldf_triad.F90. Done, added "ldfull=.TRUE." to end of added function in line with v5 changes that occured in the previous line.
traqsr.F90. Didn't need updating, no changes because the variable is defined in v5 as: nksr = nkV       ! name of max level of light extinction used in traatf(\_qco).F90
trdini.F90. Done
trdmxl.F90. Done
trdmxl_rst.F90. Done
trd_oce.F90. Didn't need updating, no changes.
trdtra.F90. Done
zdfdrg.F90 Didn't need updating, no changes.
zdfgls.F90 Didn't need updating, no changes.
zdftke.F90. Done.

Momentum files:
bdydyn.F90 the key statement "IF ( l_trddyn )" is not present in v5 src. Added in the 4.2.2 MY_SRC version. Done.
dynatf_qco.F90 not in v5 src so copied from v4.2.2. Done
dynhpg.F90 Didn't need updating, no changes.
dynspg.F90 Didn't need updating, no changes.
dynvor.F90 Didn't need updating, no changes.
dynzdf.F90 Didn't need updating, no changes.
sbcmod.F90 Didn't need updating, no changes.
trddyn.F90 Didn't need updating, no changes.
trdini.F90. Done.
dynspg_ts.F90 updated changes from 4.2.2 MY_SRC.
iceupdate.F90 zrhoco was not present in v5, no changes.



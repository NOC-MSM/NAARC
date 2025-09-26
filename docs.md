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
rn_frm_ht0 was left out of the v5 namelist_cfg_template because ice form drag is included in nemo v5. Changes in sbc_oce.F90 suggest the ice\*.F90 files probably won't need to be carried over in MY_SRC

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

Tried compiling without Momentum bits first. Error:
/work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs/NAARC/BLD/ppsrc/nemo/trdtra.f90:164:97:

  164 |          CASE( jptra_xad  )   ;   CALL trd_tra_adv( ptrd , pu  , ptra, 'X'  , ztrds, Kmm, Krhs )
      |                                                                                                 1
Error: More actual than formal arguments in procedure call at (1)
/work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs/NAARC/BLD/ppsrc/nemo/trdtra.f90:165:83:

  165 |                                   CALL trd_tra_mng( trdtx, ztrds, ktrd, kt, Kmm   )
      |                                                                                   1
Error: Missing actual argument for argument 'krhs' at (1)

Changed trdtra.F90 in MY_SRC (I'd added the Krhs arg to the wrong function).

Compilation successful.
Run Errors:
misspelled variable in namelist namrun (ref) iostat =  5010
misspelled variable in namelist namrun (cfg) iostat =  5010
misspelled variable in namelist namsbc_blk (cfg) iostat =  5010
rn_tide_ramp_dt must be lower than run duration

Tried again at updating the namelist and namelist_ice from 4.2.2 to 5.0, I think I did it wrong before.
I've set nn_fct_imp =  2 (default is 1) because nn_fct_h and nn_ft_v = 4 which is non-default
nn_fct_imp left out of namelist_cfg_template

In namelist_ref namrun, added ln_rstdate=XXX_RSD_XXX and ln_reset_ts=.false. to namelist_ref
In namelist_cfg_template namsbc_blk, changed csw, csa, cra, crw, cfa, cfw to rn_Cs_io, rn_Cs_ia, rn_Cr_ia, rn_Cr_io, rn_Cf_ia, rn_Cf_io

Run Errors:
misspelled variable in namelist namrun (ref) iostat =  5010
misspelled variable in namelist namrun (ref) iostat =  5010
misspelled variable in namelist namsbc_blk (cfg) iostat =  5010
rn_tide_ramp_dt must be lower than run duration
nemo_gcm: a total of            1  errors have been found

The nemo_gcm error is because of rn_tide_ramp_dt. Tide has been turned off in the runscript.slurm for now.

Run Errors:
misspelled variable in namelist namrun (ref) iostat =  5010
misspelled variable in namelist namrun (ref) iostat =  5010
misspelled variable in namelist namsbc_blk (cfg) iostat =  5010

rn_frm_ht0 was left out of the v5 namelist_cfg_template because ice form drag is included in nemo v5. This means the ice\*.F90 files probably won't need to be carried over in MY_SRC
nn_hls is another one to come back to because the default ahs changes from 1 to 2.

Run Errors:
misspelled variable in namelist namrun (ref) iostat =  5010
misspelled variable in namelist namrun (ref) iostat =  5010

Tired adding ln_1st_euler back in to both namelists but it didn't change the errors. ln_1st_euler in domain.F90 and istate.F90
ln_rst_eos in domain.F90, in_out_manager.F90, restart.F90
ln_1st_euler is in v5 namelist_ref on gitlab but not in the nemo_v5 folder locally. 

I've been comparing with nemo-main instead of nemo v5 so I need to check everything again. When making comparison get two tkdiff windows up. One with src v4.2.2 to src v5.0 and one with MY_SRC v4.2.2 to MY_SRC v5.0.

Done:
namelist_cfg_template, namelist_ref, namelist_ice_cfg_template, namelist_ice_ref
Done:
bdydyn.F90, diahth.F90, diawri.F90, domain.F90, dommsk.F90, dtatsd.F90, dynatf_qco.F90, fldread.F90, icerst.F90, in_out_manager.F90, istate.F90, nemogcm.F90, restart.F90, sbc_oce.F90, sbcblk.F90, sbcrnf.nc, sbcssm.F90, tide.h90, tide_mod.F90, trdini.F90, trdmxl_rst.F90, trdtra.F90, zdftke.F90
dynspg_ts.F90 - not sure "IF( l_trddyn ) THEN" bit line 963 4.2.2 MY_SRC, should be in the if loop next to pvv_b
traldf_triad.F90 - line 134 to 144 in 4.2.2 MY_SRC, not sure if this should be included with v5 change or instead? The 5.0 My_SRC is a fair combination of v5 and MY_SRC.
trdmxl.F90 - the subroutine trd_mxl_zint is mssing in 4.2.2 MY_SRC, this has been continued in 5.0 MY_SRC
MY_SRC only:
diapea.F90

I think I'm getting a segmentation fault.
Adding debug flag (-O0 -g) in arch
Trying gdb4hpc.
gdb nemo core
bt
#10 0x00000000008e8a7c in iom::iom_use (cdname=..., _cdname=9) at /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs/NAARC/BLD/ppsrc/nemo/iom.f90:2893
#11 0x0000000000cf96ec in icedyn_rdgrft::ice_dyn_rdgrft_init () at /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs/NAARC/BLD/ppsrc/nemo/icedyn_rdgrft.f90:1276
#12 0x0000000000c8ebdd in icedyn::ice_dyn_init () at /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs/NAARC/BLD/ppsrc/nemo/icedyn.f90:382
#13 0x0000000000880051 in icestp::ice_init (kbb=1, kmm=2, kaa=3) at /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs/NAARC/BLD/ppsrc/nemo/icestp.f90:336
#14 0x000000000048ee5f in sbcmod::sbc_init (kbb=1, kmm=2, kaa=3) at /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs/NAARC/BLD/ppsrc/nemo/sbcmod.f90:382
#15 0x0000000000466fd5 in nemogcm::nemo_init () at /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs/NAARC/BLD/ppsrc/nemo/nemogcm.f90:409
#16 0x0000000000467210 in nemogcm::nemo_gcm () at /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs/NAARC/BLD/ppsrc/nemo/nemogcm.f90:161
#17 0x00000000004632df in nemo () at /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs/NAARC/WORK/nemo.f90:17
#18 0x0000000000463321 in main (argc=1, argv=0x7ffc1dbd6de3) at /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs/NAARC/WORK/nemo.f90:11
#19 0x000015376b1e829d in __libc_start_main () from /lib64/libc.so.6
#20 0x000000000046321a in _start () at ../sysdeps/x86_64/start.S:120

It maybe something related to the MY_SRC icestp.f90 and ice files I haven't updated to v5.0.

Isolate Ice and Momentum files.

#10 0x00000000008e8a7c in iom::iom_use (cdname=..., _cdname=9)
    at /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs/NAARC/BLD/ppsrc/nemo/iom.f90:2893
#11 0x0000000000cf96ec in icedyn_rdgrft::ice_dyn_rdgrft_init ()
    at /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs/NAARC/BLD/ppsrc/nemo/icedyn_rdgrft.f90:1276
#12 0x0000000000c8ebdd in icedyn::ice_dyn_init ()
    at /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs/NAARC/BLD/ppsrc/nemo/icedyn.f90:382
#13 0x0000000000880051 in icestp::ice_init (kbb=1, kmm=2, kaa=3)
--Type <RET> for more, q to quit, c to continue without paging--
tp.f90:336
#14 0x000000000048ee5f in sbcmod::sbc_init (kbb=1, kmm=2, kaa=3) at /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs/NAARC/BLD/ppsrc/nemo/sbcmod.f90:382

Ice files that weren't obviouse:
sbcmod.F90, sbcblk.F90


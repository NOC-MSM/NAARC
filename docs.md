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
rn_frm_ht0 was left out of the v5 namelist_cfg_template because ice form drag is included in nemo v5. Changes in sbcblk.F90 suggest the ice\*.F90 files probably won't need to be carried over in MY_SRC

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
sbcmod.F90 Done.
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
Trying: 
module load gdb4hpc.
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
    at /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs/NAARC/BLD/ppsrc/nemo/icestp.f90:336
#14 0x000000000048ee5f in sbcmod::sbc_init (kbb=1, kmm=2, kaa=3) at /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs/NAARC/BLD/ppsrc/nemo/sbcmod.f90:382

Ice files that weren't obvious and have been moved to the Ice folder:
sbcmod.F90, sbcblk.F90

Still same error.

Time to compare the MY_SRC ice files to v5.0 ice. The ice.F90 is substantially different in nemo v5.0, I'm not sure where to start with making changes. It is more different between v4.2.2 than it is the same.

namelist_ice_ref in EXPREF is the same as nemo v5. namelist_ice_cfg_template has one additional variable rn_alb_lpnd and rn_cio has changed to run_Cd_io.

Updating EXPREF/field_def_nemo-ice.xml to v5. EXPREF/field_def_nemo-ice.xml only adds lines relative to v4.2.2 SHARED/field_def_nemo-ice.xml, only adding sbcssm variables and not ice related ones.
EXPREF/file_def_nemo-ice.xml doesn't have any of the custom variables added so I left as it was in v4.2.2 EXPREF.

Updating EXPREF/field_def_nemo-oce.xml to v5. The tidal harmonics have been changed to v5 and not copied from EXPREF v4.2.2 because nn_tide_var = 1 in the namelist_cfg_template.

#6  0x0000000000423c13 in xios::CField::solveGridReference() [clone .cold] ()
#7  0x0000000001049692 in xios::CField::solveOnlyReferenceEnabledField(bool) ()
#8  0x000000000105c167 in xios::CFile::solveOnlyRefOfEnabledFields(bool) ()
#9  0x0000000000fd8f57 in xios::CContext::solveOnlyRefOfEnabledFields(bool) ()
#10 0x0000000000fe082b in xios::CContext::postProcessing() ()
#11 0x0000000000fe6d64 in xios::CContext::postProcessingGlobalAttributes() ()
#12 0x0000000000fe7258 in xios::CContext::closeDefinition() ()
#13 0x000000000117dd59 in cxios_context_close_definition ()
#14 0x00000000009120af in iom::iom_init_closedef (cdname=..., _cdname=0)
    at /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs/NAARC/BLD/ppsrc/nemo/iom.f90:322
#15 0x00000000004cf896 in stpmlf::stp_mlf (kstp=1)
    at /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs/NAARC/BLD/ppsrc/nemo/stpmlf.f90:124
#16 0x0000000000467340 in nemogcm::nemo_gcm ()
    at /work/n01/n01/benbar/NAARC/NAARC_RUNS/nemo/cfgs/NAARC/BLD/ppsrc/nemo/nemogcm.f90:186

Try using v5 default .xml files instead of custom ones.
Error didn't change so I guess keep the custom ones.

According to James the nn_ice bits in sbcmod should be included in sbcmod.F90.

I've re-checked the files I previously marked as "Didn't need updating" and some of them do need updating now I'm comparing the right nemo versions.
I've made changes to main MY_SRC:
diadct.F90, dynspg_ts.F90, iom.F90, ldftra.F90, traldf.F90, traqsr.F90

Momentum MY_SRC:
dynhpg.F90, dynspg.F90, dynvor.F90, dynzdf.F90, sbcmod.F90, trd_oce.F90, trddyn.F90

Still having the same segmentation fault in backtrace.
Error in slurm-\*.out:

In file "field.cpp", function "void xios::CField::solveGridReference()",  line 1
604 -> A grid must be defined for field 'fsitherm' .

Commented out fsitherm in file_def_nemo-oce.xml.
Deleted line fsitherm in file_def_nemo-oce.xml.

Slurm-\*.out error:
STOP from timing: try to stop stp_MLF but we point toward dia_wri, MPI rank: 10

Moved dia_wri.F90 out of MY_SRC to see if it makes a difference.
It did make a difference, it still doesn't run but it made some output.abort.nc files.
slurm error:
-> report :  Memory report : Context <nemo> : client side : total memory used fo
r buffer 43198880 bytes

ocean.output error:
  ===>>> : E R R O R

          ===========

   stp_ctl: |ssh| > 20 m  or  |U| > 10 m/s  or  S <= 0  or  S >= 100  or  NaN encounter in the tests
 
 kt 2 |ssh| max    Infinity at i j   3564 1967    found in   24 MPI tasks, spread out among ranks  156 to 1110
 kt 2 |U|   max   9.561     at i j k 2732 2671  1 MPI rank  586
 kt 2 |V|   max   10.33     at i j k 2732 2670  1 MPI rank  586
 kt 2 Sal   min   1.020     at i j k  674 3527  1 MPI rank 1389
 kt 2 Sal   max   39.21     at i j k 3848 2492 30 MPI rank  454
 
        ===> output of last computed fields in output.abort* files


Perhaps the problem is the diapea.F90 file called by diawri.F90. Commented out diapea.

Same Slurm-\*.out error as earlier:
STOP from timing: try to stop stp_MLF but we point toward dia_wri, MPI rank: 10

diapea was false in the namelist anyway. ln_zdftke is also false. ln_ldfeiv is false. Could be teos10. Edit namelist_cfg_template to turn teos10 off.

Didn't help.

Starting with fresh v5.0 diawri.F90 and slowly adding bit of MY_SRC.
It got past the diawri issue.

Adding zdftke htau to diawri.F90.
It got past the diawri issue.

Add ldfevi bit to diawri.F90.
It got past the diawri issue so the problem must be with the equation of state (EOS) pot/con and abs/pra types.

Try adding part of the EOS bit only temperature.
Previous slurm error with dia_wri.

Removed problematic EOS code from diawri.F90 again to continue. Something to come back to.
The abort output with ssh going to infinate persists.

ocean.output error:
  ===>>> : E R R O R

          ===========

   stp_ctl: |ssh| > 20 m  or  |U| > 10 m/s  or  S <= 0  or  S >= 100  or  NaN encounter in the tests
 
 kt 2 |ssh| max    Infinity at i j   3564 1967    found in   24 MPI tasks, spread out among ranks  156 to 1110
 kt 2 |U|   max   9.561     at i j k 2732 2671  1 MPI rank  586
 kt 2 |V|   max   10.33     at i j k 2732 2670  1 MPI rank  586
 kt 2 Sal   min   1.020     at i j k  674 3527  1 MPI rank 1389
 kt 2 Sal   max   39.21     at i j k 3848 2492 30 MPI rank  454
 
        ===> output of last computed fields in output.abort* files

Slurm error:
slurmstepd: error: Detected 1 oom-kill event(s) in StepId=11088795.0+0. Some of your processes may have been killed by the cgroup out-of-memory handler.
srun: error: nid003400: task 1275: Out Of Memory

Tried turning boundaries off in the namelist_cfg_template. Slurm error still related to Out Of Memory but at earlier stage (line 1567 instead of 3094). No error in the ocean.output and no abort files produced. core file produced.

Setting nn_fsbc = 0 made a core file.
Setting ln_blk, ln_traqsr and ln_ssr to false, nn_ice = 0. Error about SI3 needing blk.
Setting ln_blk = true but rest off. Error still out of memory.

Set namelist back to normal. The abort files suggest the blow up is happening in the Arctic.
Trying turning the sea ice (key_si3) off in the cpp compile. Made a core file without any clues.
Turn sea ice back on. Isabella from MO said they had similar issues in the beginning with this error which in our case was found to be related with vectorisation.
They suggest adding the following to the FCFLAGS in arch:
-hvector0 -Ovector0
but this seems to be a cray thing not a mpich and compilation fails because it doesn't recognise the -Ovector0 bit.

Still getting abort files with -hvector0.

James made a working nemo v5 NAARC configuration with some of the MY_SRC.
Apparetly adding -O2 back in compile arch fine.
Take a few lateral boundary conditions bits out of the namelist.
DIA MY_SRC files are for diagnostics and could be added in probably without too much trouble.
Apparently some MY_SRC files don't have any differences.



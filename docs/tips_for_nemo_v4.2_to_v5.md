Copied the NEMO version 4.2.2 files to 5.0 for the folders /arch, /CPP, EXPREF.

The compile command needs "-n NAARC" for correct directory structure so "./makenemo -m archer2-gnu-mpich -r NAARC -n NAARC -j 16". Edited scripts/setup/NAARC_setup

NAARC compiled with v5 but does not run.

INPUTS was empty because a folder for 5.0 did not exist in /work/n01/shared/NAARC/ to link to. The folder /work/n01/shared/NAARC/4.2.2 is now linked to /work/n01/shared/NAARC/5.0 and /work/n01/shared/NAARC/5.0.1.

Changes have been made to translate the namelist_cfg_template and namelist_ice_cfg_template to v5 namelists by removing variables not present in the v5 namelist_ref. This was done by comparing v4.2.2 namelist_ref and v5.0 namelist_ref with tkdiff. Some variables changed name in namelist_cfg_template namsbc_blk, changed csw, csa, cra, crw, cfa, cfw to rn_Cs_io, rn_Cs_ia, rn_Cr_ia, rn_Cr_io, rn_Cf_ia, rn_Cf_io. 

Values of variables have been left unchanged with the exeption of rn_alb_dpnd and nn_fct_imp. rn_alb_dpnd which was set to 0.30 but a new comment suggest that is outside the obs range 0.12 -- 0.25 so the new default 0.18 has be used for rn_alb_dpnd. I've set nn_fct_imp = 2 (default is 1) because nn_fct_h and nn_ft_v = 4 which is non-default. **Not sure about this**

rn_frm_ht0 was left out of the v5 namelist_cfg_template because ice form drag is included in nemo v5. This means the ice\*.F90 files probably won't need to be carried over in MY_SRC.

In the v5 namelist_ref namrun didn't have ln_rstdate=XXX_RSD_XXX and ln_reset_ts=.false. so these have been added to v5 namelist_ref.

MY_SRC has been addapted to v5 from v4.2.2. This has been done by copying v5 src to MY_SRC for the respective file and editing it. The edits have been identified by comparing v4.2.2 MY_SRC files with v4.2 src files with tkdiff.


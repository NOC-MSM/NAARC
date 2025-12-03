Copied the NEMO version 4.2.2 files to 5.0 for the folders /arch, /CPP, EXPREF.

The compile command needs "-n NAARC" for correct directory structure so "./makenemo -m archer2-gnu-mpich -r NAARC -n NAARC -j 16". Edited scripts/setup/NAARC_setup

NAARC compiled with v5 but does not run.

INPUTS was empty because a folder for 5.0 did not exist in /work/n01/shared/NAARC/ to link to. The folder /work/n01/shared/NAARC/4.2.2 is now linked to /work/n01/shared/NAARC/5.0 and /work/n01/shared/NAARC/5.0.1.

Changes have been made to translate the namelist_cfg_template and namelist_ice_cfg_template to v5 namelists by removing variables not present in the v5 namelist_ref. This was done by comparing v4.2.2 namelist_ref and v5.0 namelist_ref with tkdiff. Some variables changed name in namelist_cfg_template namsbc_blk, changed csw, csa, cra, crw, cfa, cfw to rn_Cs_io, rn_Cs_ia, rn_Cr_ia, rn_Cr_io, rn_Cf_ia, rn_Cf_io. 

Values of variables have been left unchanged with the exeption of rn_alb_dpnd and nn_fct_imp. rn_alb_dpnd which was set to 0.30 but a new comment suggest that is outside the obs range 0.12 -- 0.25 so the new default 0.18 has be used for rn_alb_dpnd. I've left nn_fct_imp out of the namelist but perhaps change it in the future. 

rn_frm_ht0 was left out of the v5 namelist_cfg_template because ice form drag is included in nemo v5. This means the ice\*.F90 files probably won't need to be carried over in MY_SRC.

In the v5 namelist_ref namrun didn't have ln_rstdate=XXX_RSD_XXX and ln_reset_ts=.false. so these have been added to v5 namelist_ref.

MY_SRC has been adapted to v5 from v4.2.2. This has been done by copying v5 src to MY_SRC for the respective file and editing it. The edits have been identified by comparing v4.2.2 MY_SRC files with v4.2 src files with tkdiff.

In NAARC, the MY_SRC files were grouped into changes in files that go together. These were 7 groups in branches:
1. nemo_v5_base - domain.F90, dtatsd.F90, fldread.F90, in_out_manager.F90, iom.F90, istate.F90, restart.F90, traqsr.F90, zdftke.F90
2. nemo_v5_diagnostics - diadct.F90, diahth.F90, diapea.F90, diawri.F90, nemogcm.F90
3. nemo_v5_lat_slip - dommsk.F90
4. nemo_v5_triad - traldf_triad.F90
5. nemo_v5_eos_river - diawri.F90, sbc_oce.F90, sbcrnf.F90, sbcssm.F90
6. nemo_v5_iso_neut - ldftra.F90, trdmxl.F90, trdmxl_rst.F90, trdtra.F90
7. nemo_v5_trend - dynspg_ts.F90
NOTE: diawri.F90 had changes from two of these groups.

All MY_SRC files now work and the separate branches for the MY_SRC groups have been combined to nemo_v5. The sea ice MY_SRC (ice\*.F90) has been substantially updated in NEMOv5.0 compared to NEMOv4.2.2 probably including what was done in the MY_SRC but it is not clear if the ice MY_SRC is needed or how it could be integrated. The Momentum budget MY_SRC has been ported from NEMOv4.2.2 to NEMOv5.0 but untested. Instead the Momentum budget code in NEMOv4.2.3 had been ported to NEMO-main and will be relaesed in NEMOv5.1. dynspg_ts.F90 had the main difference between NEMOv5.1 and NEMOv5.0 so I can probably use the draft of dynspg_ts.F90 here to have it working in NEMOv5.0.

Problems Encountered:
A segmentation fault was fixed by removing 	  
<field field_ref="fmmflx"       name="fsitherm"   />
from file_def_nemo-oce.xml

An Out of Memory (OOM) error was fixed by diabling ice output in file_def_nemo-ice.xml.

One bug I found in the sbcrnf.F90 MY_SRC was:
rnf_tsc(:,:,jp_tem) = MAX( sst_m_con(:,:), 0.0_wp ) * rnf(:,:) * r1_rho0 ! SLWA
needed to be
rnf_tsc(:,:,jp_tem) = MAX( sst_m_con(A2D(0)), 0.0_wp ) * rnf(A2D(0)) * r1_rho0



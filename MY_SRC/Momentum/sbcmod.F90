MODULE sbcmod
   !!======================================================================
   !!                       ***  MODULE  sbcmod  ***
   !! Surface module :  provide to the ocean its surface boundary condition
   !!======================================================================
   !! History :  3.0  ! 2006-07  (G. Madec)  Original code
   !!            3.1  ! 2008-08  (S. Masson, A. Caubel, E. Maisonnave, G. Madec) coupled interface
   !!            3.3  ! 2010-04  (M. Leclair, G. Madec)  Forcing averaged over 2 time steps
   !!            3.3  ! 2010-10  (S. Masson)  add diurnal cycle
   !!            3.3  ! 2010-09  (D. Storkey) add ice boundary conditions (BDY)
   !!             -   ! 2010-11  (G. Madec) ice-ocean stress always computed at each ocean time-step
   !!             -   ! 2010-10  (J. Chanut, C. Bricaud, G. Madec)  add the surface pressure forcing
   !!            3.4  ! 2011-11  (C. Harris) CICE added as an option
   !!            3.5  ! 2012-11  (A. Coward, G. Madec) Rethink of heat, mass and salt surface fluxes
   !!            3.6  ! 2014-11  (P. Mathiot, C. Harris) add ice shelves melting
   !!            4.0  ! 2016-06  (L. Brodeau) new general bulk formulation
   !!            4.0  ! 2019-03  (F. Lemarié & G. Samson)  add ABL compatibility (ln_abl=TRUE)
   !!            4.2  ! 2020-12  (G. Madec, E. Clementi) modified wave forcing and coupling
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_init      : read namsbc namelist
   !!   sbc           : surface ocean momentum, heat and freshwater boundary conditions
   !!   sbc_final     : Finalize CICE ice model (if used)
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE closea         ! closed seas
   USE phycst         ! physical constants
   USE sbc_phy, ONLY : pp_cldf
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE trc_oce        ! shared ocean-passive tracers variables
   USE sbc_ice        ! Surface boundary condition: ice fields
   USE sbcdcy         ! surface boundary condition: diurnal cycle
   USE sbcssm         ! surface boundary condition: sea-surface mean variables
   USE sbcflx         ! surface boundary condition: flux formulation
   USE sbcblk         ! surface boundary condition: bulk formulation
   USE sbcabl         ! atmospheric boundary layer
   USE sbcice_if      ! surface boundary condition: ice-if sea-ice model
#if defined key_si3
   USE par_ice        ! SI3 parameters
   USE icestp         ! surface boundary condition: SI3 sea-ice model
#endif
   USE sbcice_cice    ! surface boundary condition: CICE sea-ice model
   USE sbccpl         ! surface boundary condition: coupled formulation
   USE cpl_oasis3     ! OASIS routines for coupling
   USE sbcclo         ! surface boundary condition: closed sea correction
   USE sbcssr         ! surface boundary condition: sea surface restoring
   USE sbcrnf         ! surface boundary condition: runoffs
   USE sbcapr         ! surface boundary condition: atmo pressure
   USE sbcfwb         ! surface boundary condition: freshwater budget
   USE icbstp         ! Icebergs
   USE icb_oce  , ONLY : ln_passive_mode      ! iceberg interaction mode
   USE traqsr         ! active tracers: light penetration
   USE sbcwave        ! Wave module
   USE bdy_oce   , ONLY: ln_bdy
   USE usrdef_sbc     ! user defined: surface boundary condition
   USE closea         ! closed sea
   USE trd_oce, ONLY: l_trddyn ! flag for dynamics trends diagnostics
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   !
   USE prtctl         ! Print control                    (prt_ctl routine)
   USE iom            ! IOM library
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE timing         ! Timing
   USE wet_dry
   USE diu_bulk, ONLY:   ln_diurnal_only   ! diurnal SST diagnostic

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc        ! routine called by step.F90
   PUBLIC   sbc_init   ! routine called by opa.F90

   INTEGER ::   nsbc   ! type of surface boundary condition (deduced from namsbc informations)
   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_init( Kbb, Kmm, Kaa )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_init ***
      !!
      !! ** Purpose :   Initialisation of the ocean surface boundary computation
      !!
      !! ** Method  :   Read the namsbc namelist and set derived parameters
      !!                Call init routines for all other SBC modules that have one
      !!
      !! ** Action  : - read namsbc parameters
      !!              - nsbc: type of sbc
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kbb, Kmm, Kaa         ! ocean time level indices
      INTEGER ::   ios, icpt                         ! local integer
      LOGICAL ::   ll_purecpl, ll_opa, ll_not_nemo   ! local logical
      !!
      NAMELIST/namsbc/ nn_fsbc  ,                                                    &
         &             ln_usr   , ln_flx   , ln_blk   , ln_abl,                      &
         &             ln_cpl   , ln_mixcpl, nn_components,                          &
         &             nn_ice   , ln_ice_embd,                                       &
         &             ln_traqsr, ln_dm2dc ,                                         &
         &             ln_rnf   , nn_fwb     , ln_ssr   , ln_apr_dyn,                &
         &             ln_wave  , nn_lsm
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sbc_init : surface boundary condition setting'
         WRITE(numout,*) '~~~~~~~~ '
      ENDIF
      !
      !                       !**  read Surface Module namelist
      READ_NML_REF(numnam,namsbc)
      READ_NML_CFG(numnam,namsbc)
      IF(lwm) WRITE( numond, namsbc )
      !
#if ! defined key_mpi_off
      ncom_fsbc = nn_fsbc    ! make nn_fsbc available for lib_mpp
#endif
#if ! defined key_si3
      IF( nn_ice == 2 )    nn_ice = 0  ! without key key_si3 you cannot use si3...
#endif
#if defined key_agrif
      ! In case of an agrif zoom, the freshwater water budget is determined by parent:
      IF (.NOT.Agrif_Root()) nn_fwb = Agrif_Parent(nn_fwb) 
#endif
      !
      !
      IF(lwp) THEN                  !* Control print
         WRITE(numout,*) '   Namelist namsbc (partly overwritten with CPP key setting)'
         WRITE(numout,*) '      frequency update of sbc (and ice)             nn_fsbc       = ', nn_fsbc
         WRITE(numout,*) '      Type of air-sea fluxes : '
         WRITE(numout,*) '         user defined formulation                   ln_usr        = ', ln_usr
         WRITE(numout,*) '         flux         formulation                   ln_flx        = ', ln_flx
         WRITE(numout,*) '         bulk         formulation                   ln_blk        = ', ln_blk
         WRITE(numout,*) '         ABL          formulation                   ln_abl        = ', ln_abl
         WRITE(numout,*) '         Surface wave (forced or coupled)           ln_wave       = ', ln_wave
         WRITE(numout,*) '      Type of coupling (Ocean/Ice/Atmosphere) : '
         WRITE(numout,*) '         ocean-atmosphere coupled formulation       ln_cpl        = ', ln_cpl
         WRITE(numout,*) '         mixed forced-coupled     formulation       ln_mixcpl     = ', ln_mixcpl
!!gm  lk_oasis is controlled by key_oasis3  ===>>>  It shoud be removed from the namelist
         WRITE(numout,*) '         OASIS coupling (with atm or sas)           lk_oasis      = ', lk_oasis
         WRITE(numout,*) '         components of your executable              nn_components = ', nn_components
         WRITE(numout,*) '      Sea-ice : '
         WRITE(numout,*) '         ice management in the sbc (=0/1/2/3)       nn_ice        = ', nn_ice
         WRITE(numout,*) '         ice embedded into ocean                    ln_ice_embd   = ', ln_ice_embd
         WRITE(numout,*) '      Misc. options of sbc : '
         WRITE(numout,*) '         Light penetration in temperature Eq.       ln_traqsr     = ', ln_traqsr
         WRITE(numout,*) '            daily mean to diurnal cycle qsr            ln_dm2dc   = ', ln_dm2dc
         WRITE(numout,*) '         Sea Surface Restoring on SST and/or SSS    ln_ssr        = ', ln_ssr
         WRITE(numout,*) '         FreshWater Budget control  (=0/1/2)        nn_fwb        = ', nn_fwb
         WRITE(numout,*) '         Patm gradient added in ocean & ice Eqs.    ln_apr_dyn    = ', ln_apr_dyn
         WRITE(numout,*) '         runoff / runoff mouths                     ln_rnf        = ', ln_rnf
         WRITE(numout,*) '         nb of iterations if land-sea-mask applied  nn_lsm        = ', nn_lsm
      ENDIF
      !
      IF( .NOT.ln_usr ) THEN     ! the model calendar needs some specificities (except in user defined case)
         IF( MOD( rday , rn_Dt ) /= 0. )   CALL ctl_stop( 'the time step must devide the number of second of in a day' )
         IF( MOD( rday , 2.    ) /= 0. )   CALL ctl_stop( 'the number of second of in a day must be an even number'    )
         IF( MOD( rn_Dt, 2.    ) /= 0. )   CALL ctl_stop( 'the time step (in second) must be an even number'           )
      ENDIF
      !                       !**  check option consistency
      !
      IF(lwp) WRITE(numout,*)       !* Single / Multi - executable (NEMO / OCE+SAS)
      SELECT CASE( nn_components )
      CASE( jp_iam_nemo )
         IF(lwp) WRITE(numout,*) '   ==>>>   NEMO configured as a single executable (i.e. including both OCE and Surface module)'
      CASE( jp_iam_oce  )
         IF(lwp) WRITE(numout,*) '   ==>>>   Multi executable configuration. Here, OCE component'
         IF( .NOT.lk_oasis )   CALL ctl_stop( 'sbc_init : OCE-SAS coupled via OASIS, but key_oasis3 disabled' )
         IF( ln_cpl        )   CALL ctl_stop( 'sbc_init : OCE-SAS coupled via OASIS, but ln_cpl = T in OCE'   )
         IF( ln_mixcpl     )   CALL ctl_stop( 'sbc_init : OCE-SAS coupled via OASIS, but ln_mixcpl = T in OCE' )
      CASE( jp_iam_sas  )
         IF(lwp) WRITE(numout,*) '   ==>>>   Multi executable configuration. Here, SAS component'
         IF( .NOT.lk_oasis )   CALL ctl_stop( 'sbc_init : OCE-SAS coupled via OASIS, but key_oasis3 disabled' )
         IF( ln_mixcpl     )   CALL ctl_stop( 'sbc_init : OCE-SAS coupled via OASIS, but ln_mixcpl = T in OCE' )
      CASE DEFAULT
         CALL ctl_stop( 'sbc_init : unsupported value for nn_components' )
      END SELECT
      !                             !* coupled options
      IF( ln_cpl ) THEN
         IF( .NOT. lk_oasis )   CALL ctl_stop( 'sbc_init : coupled mode with an atmosphere model (ln_cpl=T)',   &
            &                                  '           required to defined key_oasis3' )
      ENDIF
      IF( ln_mixcpl ) THEN
         IF( .NOT. lk_oasis )   CALL ctl_stop( 'sbc_init : mixed forced-coupled mode (ln_mixcpl=T) ',   &
            &                                  '           required to defined key_oasis3' )
         IF( .NOT.ln_cpl    )   CALL ctl_stop( 'sbc_init : mixed forced-coupled mode (ln_mixcpl=T) requires ln_cpl = T' )
         IF( nn_components /= jp_iam_nemo )    &
            &                   CALL ctl_stop( 'sbc_init : the mixed forced-coupled mode (ln_mixcpl=T) ',   &
            &                                   '          not yet working with sas-opa coupling via oasis' )
      ENDIF
      !                             !* sea-ice
      SELECT CASE( nn_ice )
      CASE( 0 )                        !- no ice in the domain
      CASE( 1 )                        !- Ice-cover climatology ("Ice-if" model)
      CASE( 2 )                        !- SI3  ice model
         IF( .NOT.( ln_blk .OR. ln_cpl .OR. ln_abl .OR. ln_usr ) )   &
            &                   CALL ctl_stop( 'sbc_init : SI3 sea-ice model requires ln_blk or ln_cpl or ln_abl or ln_usr = T' )
      CASE( 3 )                        !- CICE ice model
         IF( .NOT.( ln_blk .OR. ln_cpl .OR. ln_abl .OR. ln_usr ) )   &
            &                   CALL ctl_stop( 'sbc_init : CICE sea-ice model requires ln_blk or ln_cpl or ln_abl or ln_usr = T' )
         IF( lk_agrif                                )   &
            &                   CALL ctl_stop( 'sbc_init : CICE sea-ice model not currently available with AGRIF' )
      CASE DEFAULT                     !- not supported
      END SELECT
      IF( ln_diurnal .AND. .NOT. (ln_blk.OR.ln_abl) )   CALL ctl_stop( "sbc_init: diurnal flux processing only implemented for bulk forcing" )
      !
      !                       !**  allocate and set required variables
      !
      !                             !* allocate sbc arrays
      IF( sbc_oce_alloc() /= 0 )   CALL ctl_stop( 'sbc_init : unable to allocate sbc_oce arrays' )
#if ! defined key_si3 && ! defined key_cice
      IF( sbc_ice_alloc() /= 0 )   CALL ctl_stop( 'sbc_init : unable to allocate sbc_ice arrays' )
#endif
      !
      !
      IF( sbc_ssr_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'sbc_init : unable to allocate sbc_ssr arrays' )
      IF( .NOT.ln_ssr ) THEN               !* Initialize qrp and erp if no restoring
         qrp(:,:) = 0._wp
         erp(:,:) = 0._wp
      ENDIF
      !
      IF( nn_ice == 0 ) THEN        !* No sea-ice in the domain : ice fraction is always zero
         IF( nn_components /= jp_iam_oce )   fr_i(:,:) = 0._wp    ! except for OCE in SAS-OCE coupled case
      ENDIF
      !
      sfx   (:,:) = 0._wp           !* salt flux due to freezing/melting
      fwfice(:,:) = 0._wp           !* ice-ocean freshwater flux
      cloud_fra(:,:) = pp_cldf      !* cloud fraction over sea ice (used in si3)

      taum(:,:) = 0._wp             !* wind stress module (needed in GLS in case of reduced restart)

      !                          ! Choice of the Surface Boudary Condition (set nsbc)
      nday_qsr = -1   ! allow initialization at the 1st call !LB: now warm-layer of COARE* calls "sbc_dcy_param" of sbcdcy.F90!
      IF( ln_dm2dc ) THEN           !* daily mean to diurnal cycle
         !LB:nday_qsr = -1   ! allow initialization at the 1st call
         IF( .NOT.( ln_flx .OR. ln_blk .OR. ln_abl ) .AND. nn_components /= jp_iam_oce )   &
            &   CALL ctl_stop( 'qsr diurnal cycle from daily values requires flux, bulk or abl formulation' )
      ENDIF
      !                             !* Choice of the Surface Boudary Condition
      !                             (set nsbc)
      !
      ll_purecpl  = ln_cpl .AND. .NOT.ln_mixcpl
      ll_opa      = nn_components == jp_iam_oce
      ll_not_nemo = nn_components /= jp_iam_nemo
      icpt = 0
      !
      IF( ln_usr          ) THEN   ;   nsbc = jp_usr     ; icpt = icpt + 1   ;   ENDIF       ! user defined         formulation
      IF( ln_flx          ) THEN   ;   nsbc = jp_flx     ; icpt = icpt + 1   ;   ENDIF       ! flux                 formulation
      IF( ln_blk          ) THEN   ;   nsbc = jp_blk     ; icpt = icpt + 1   ;   ENDIF       ! bulk                 formulation
      IF( ln_abl          ) THEN   ;   nsbc = jp_abl     ; icpt = icpt + 1   ;   ENDIF       ! ABL                  formulation
      IF( ll_purecpl      ) THEN   ;   nsbc = jp_purecpl ; icpt = icpt + 1   ;   ENDIF       ! Pure Coupled         formulation
      IF( ll_opa          ) THEN   ;   nsbc = jp_none    ; icpt = icpt + 1   ;   ENDIF       ! opa coupling via SAS module
      !
      IF( icpt /= 1 )    CALL ctl_stop( 'sbc_init : choose ONE and only ONE sbc option' )
      !
      IF(lwp) THEN                     !- print the choice of surface flux formulation
         WRITE(numout,*)
         SELECT CASE( nsbc )
         CASE( jp_usr     )   ;   WRITE(numout,*) '   ==>>>   user defined forcing formulation'
         CASE( jp_flx     )   ;   WRITE(numout,*) '   ==>>>   flux formulation'
         CASE( jp_blk     )   ;   WRITE(numout,*) '   ==>>>   bulk formulation'
         CASE( jp_abl     )   ;   WRITE(numout,*) '   ==>>>   ABL  formulation'
         CASE( jp_purecpl )   ;   WRITE(numout,*) '   ==>>>   pure coupled formulation'
!!gm abusive use of jp_none ??   ===>>> need to be check and changed by adding a jp_sas parameter
         CASE( jp_none    )   ;   WRITE(numout,*) '   ==>>>   OCE coupled to SAS via oasis'
            IF( ln_mixcpl )       WRITE(numout,*) '               + forced-coupled mixed formulation'
         END SELECT
         IF( ll_not_nemo  )       WRITE(numout,*) '               + OASIS coupled SAS'
      ENDIF
      !
      !                             !* OASIS initialization
      !
      IF( lk_oasis )   CALL sbc_cpl_init( nn_ice )   ! Must be done before: (1) first time step
      !                                              !                      (2) the use of nn_fsbc
      !     nn_fsbc initialization if OCE-SAS coupling via OASIS
      !     SAS time-step has to be declared in OASIS (mandatory) -> nn_fsbc has to be modified accordingly
      IF( nn_components /= jp_iam_nemo ) THEN
         IF( nn_components == jp_iam_oce )   nn_fsbc = cpl_freq('O_SFLX') / NINT(rn_Dt)
         IF( nn_components == jp_iam_sas )   nn_fsbc = cpl_freq('I_SFLX') / NINT(rn_Dt)
         !
         IF(lwp)THEN
            WRITE(numout,*)
            WRITE(numout,*)"   OCE-SAS coupled via OASIS : nn_fsbc re-defined from OASIS namcouple ", nn_fsbc
            WRITE(numout,*)
         ENDIF
      ENDIF
      !
      !                             !* check consistency between model timeline and nn_fsbc
      IF( ln_rst_list .OR. nn_stock /= -1 ) THEN   ! we will do restart files
         IF( MOD( nitend - nit000 + 1, nn_fsbc) /= 0 ) THEN
            WRITE(ctmp1,*) 'sbc_init : experiment length (', nitend - nit000 + 1, ') is NOT a multiple of nn_fsbc (', nn_fsbc, ')'
            CALL ctl_stop( ctmp1, 'Impossible to properly do model restart' )
         ENDIF
         IF( .NOT. ln_rst_list .AND. MOD( nn_stock, nn_fsbc) /= 0 ) THEN   ! we don't use nn_stock if ln_rst_list
            WRITE(ctmp1,*) 'sbc_init : nn_stock (', nn_stock, ') is NOT a multiple of nn_fsbc (', nn_fsbc, ')'
            CALL ctl_stop( ctmp1, 'Impossible to properly do model restart' )
         ENDIF
      ENDIF
      !
      IF( MOD( rday, REAL(nn_fsbc, wp) * rn_Dt ) /= 0 )   &
         &  CALL ctl_warn( 'sbc_init : nn_fsbc is NOT a multiple of the number of time steps in a day' )
      !
      IF( ln_dm2dc .AND. NINT(rday) / ( nn_fsbc * NINT(rn_Dt) ) < 8  )   &
         &   CALL ctl_warn( 'sbc_init : diurnal cycle for qsr: the sampling of the diurnal cycle is too small...' )
      !

      !                       !**  associated modules : initialization
      !
                          CALL sbc_ssm_init ( Kbb, Kmm ) ! Sea-surface mean fields initialization
      !
      IF( l_sbc_clo   )   CALL sbc_clo_init              ! closed sea surface initialisation
      !
      IF( ln_blk      )   CALL sbc_blk_init              ! bulk formulae initialization

      IF( ln_abl      )   CALL sbc_abl_init              ! Atmospheric Boundary Layer (ABL)

      IF( ln_ssr      )   CALL sbc_ssr_init              ! Sea-Surface Restoring initialization
      !
      !
                          CALL sbc_rnf_init( Kmm )       ! Runof initialization
      !
      IF( ln_apr_dyn )    CALL sbc_apr_init              ! Atmo Pressure Forcing initialization
      !
#if defined key_si3
      IF( nn_ice == 0 ) THEN
         ! allocate ice arrays in case agrif + ice-model + no-ice in child grid
         jpl = 1 ; nlay_i = 1 ; nlay_s = 1
         IF( sbc_ice_alloc() /= 0 )   CALL ctl_stop('STOP', 'sbc_ice_alloc : unable to allocate arrays' )
#if defined key_agrif
         CALL Agrif_Declare_Var_ice  !  "      "   "   "      "  Sea ice
#endif

      ELSEIF( nn_ice == 2 ) THEN
                          CALL ice_init( Kbb, Kmm, Kaa )         ! ICE initialization
      ENDIF
#endif
      IF( nn_ice == 3 )   CALL cice_sbc_init( nsbc, Kbb, Kmm )   ! CICE initialization
      !
      IF( ln_wave     ) THEN
                          CALL sbc_wave_init                     ! surface wave initialisation
      ELSE
                          IF(lwp) WRITE(numout,*)
                          IF(lwp) WRITE(numout,*) '   No surface waves : all wave related logical set to false'
                          ln_sdw       = .false.
                          ln_stcor     = .false.
                          ln_cdgw      = .false.
                          ln_tauoc     = .false.
                          ln_wave_test = .false.
                          ln_charn     = .false.
                          ln_taw       = .false.
                          ln_phioc     = .false.
                          ln_bern_srfc = .false.
                          ln_breivikFV_2016 = .false.
                          ln_vortex_force = .false.
                          ln_stshear  = .false.
      ENDIF
      !
   END SUBROUTINE sbc_init


   SUBROUTINE sbc( kt, Kbb, Kmm )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc  ***
      !!
      !! ** Purpose :   provide at each time-step the ocean surface boundary
      !!                condition (momentum, heat and freshwater fluxes)
      !!
      !! ** Method  :   blah blah  to be written ?????????
      !!                CAUTION : never mask the surface stress field (tke sbc)
      !!
      !! ** Action  : - set the ocean surface boundary condition at before and now
      !!                time step, i.e.
      !!                utau_b, vtau_b, qns_b, qsr_b, emp_b, sfx_b
      !!                utau  , vtau  , qns  , qsr  , emp  , sfx  , qrp  , erp
      !!              - updte the ice fraction : fr_i
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Kmm   ! ocean time level indices
      INTEGER  ::   jj, ji          ! dummy loop argument
      !
      LOGICAL  ::   ll_sas, ll_opa  ! local logical
      !
      REAL(wp) ::  zthscl        ! wd  tanh scale
      REAL(wp) ::  zwdht, zwght  ! wd dep over wd limit, wgt
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('sbc')
      !
#if ! defined key_RK3
      !                                            ! ---------------------------------------- !
      IF( kt /= nit000 ) THEN                      !          Swap of forcing fields          !
         !                                         ! ---------------------------------------- !
         utau_b(:,:) = utauU(:,:)                        ! Swap the ocean forcing fields
         vtau_b(:,:) = vtauV(:,:)                        ! (except at nit000 where before fields
         qns_b (:,:) = qns  (:,:)                        !  are set at the end of the routine)
         emp_b (:,:) = emp  (:,:)
         sfx_b (:,:) = sfx  (:,:)
         IF( ln_rnf ) THEN
            rnf_b    (:,:  ) = rnf    (:,:  )
            rnf_tsc_b(:,:,:) = rnf_tsc(:,:,:)
         ENDIF
         IF( l_trddyn .and. nn_ice == 2 ) THEN
            uiceoc_b(:,:) = uiceoc(:,:)
            viceoc_b(:,:) = viceoc(:,:)
         ENDIF
         !
      ENDIF
#endif
      !                                            ! ---------------------------------------- !
      !                                            !        forcing field computation         !
      !                                            ! ---------------------------------------- !
      ! most of the following routines update fields only in the interior
      ! with the exception of sbcssm, sbcrnf and sbcwave modules
      !
      ll_sas = nn_components == jp_iam_sas               ! component flags
      ll_opa = nn_components == jp_iam_oce
      !
      IF( .NOT.ll_sas )   CALL sbc_ssm ( kt, Kbb, Kmm )  ! mean ocean sea surface variables (sst_m, sss_m, ssu_m, ssv_m)
      !
      !                                            !==  sbc formulation  ==!
      !
      IF( ln_blk .OR. ln_abl ) THEN
         IF( ll_sas  )         CALL sbc_cpl_rcv ( kt, nn_fsbc, nn_ice, Kbb, Kmm )   ! OCE-SAS coupling: SAS receiving fields from OCE
         IF( ln_wave ) THEN
            IF ( lk_oasis )    CALL sbc_cpl_rcv ( kt, nn_fsbc, nn_ice, Kbb, Kmm )   ! OCE-wave coupling
                               CALL sbc_wave ( kt, Kmm )
         ENDIF
      ENDIF
      !
      SELECT CASE( nsbc )                                ! Compute ocean surface boundary condition
      !                                                  ! (i.e. utau,vtau, qns, qsr, emp, sfx)
      CASE( jp_usr     )   ;   CALL usrdef_sbc_oce( kt, Kbb )                        ! user defined formulation
      CASE( jp_flx     )   ;   CALL sbc_flx       ( kt )                             ! flux formulation
      CASE( jp_blk     )   ;   CALL sbc_blk       ( kt )                             ! bulk formulation for the ocean
      CASE( jp_abl     )   ;   CALL sbc_abl       ( kt )                             ! ABL  formulation for the ocean
      CASE( jp_purecpl )   ;   CALL sbc_cpl_rcv   ( kt, nn_fsbc, nn_ice, Kbb, Kmm )  ! pure coupled formulation
      CASE( jp_none    )
         IF( ll_opa    )       CALL sbc_cpl_rcv   ( kt, nn_fsbc, nn_ice, Kbb, Kmm )   ! OCE-SAS coupling: OCE receiving fields from SAS
      END SELECT
      !
      IF( ln_mixcpl )          CALL sbc_cpl_rcv   ( kt, nn_fsbc, nn_ice, Kbb, Kmm )   ! forced-coupled mixed formulation after forcing
      !
      IF( ln_wave .AND. ln_tauoc ) THEN             ! Wave stress reduction
         !
         DO_2D( 0, 0, 0, 0 )
            utau(ji,jj) = utau(ji,jj) * tauoc_wave(ji,jj)
            vtau(ji,jj) = vtau(ji,jj) * tauoc_wave(ji,jj)
            taum(ji,jj) = taum(ji,jj) * tauoc_wave(ji,jj)
         END_2D
         !
         IF( kt == nit000 )   CALL ctl_warn( 'sbc: You are subtracting the wave stress to the ocean.',   &
            &                                'If not requested select ln_tauoc=.false.' )
         !
      ELSEIF( ln_wave .AND. ln_taw ) THEN           ! Wave stress reduction
         DO_2D( 0, 0, 0, 0 )
            utau(ji,jj) = utau(ji,jj) - ( tawx(ji,jj) - twox(ji,jj) )
            vtau(ji,jj) = vtau(ji,jj) - ( tawy(ji,jj) - twoy(ji,jj) )
            taum(ji,jj) = SQRT( utau(ji,jj)*utau(ji,jj) + vtau(ji,jj)*vtau(ji,jj) )
         END_2D
         !
         IF( kt == nit000 )   CALL ctl_warn( 'sbc: You are subtracting the wave stress to the ocean.',   &
            &                                'If not requested select ln_taw=.false.' )
         !
      ENDIF
      !
      !clem: these calls are needed for sbccpl only => only for SAS I think?
      IF( ll_sas .OR. ll_opa )   CALL lbc_lnk( 'sbcmod', sst_m, 'T', 1.0_wp, sss_m, 'T', 1.0_wp, ssh_m, 'T', 1.0_wp, &
         &                                               frq_m, 'T', 1.0_wp, e3t_m, 'T', 1.0_wp, fr_i , 'T', 1.0_wp )
      !clem : these calls are needed for sbccpl => it needs an IF statement but it's complicated
      IF( ln_rnf .AND. l_rnfcpl )     CALL lbc_lnk( 'sbcmod', rnf, 'T', 1.0_wp )
      !
      IF( ln_icebergs ) THEN  ! save pure stresses (with no ice-ocean stress) for use by icebergs
         !     Note the use of 0.5*(2-umask) in order to unmask the stress along coastlines
         !      and the use of MAX(tmask(i,j),tmask(i+1,j) is to mask tau over ice shelves
         ! (PM) cannot be move to icb because we need pure stresses. Why not extract directly wind from sbcblk i
         ! (icb only need wind)!
         CALL lbc_lnk( 'sbcmod', utau, 'T', -1.0_wp, vtau, 'T', -1.0_wp )
         DO_2D( 0, 0, 0, 0 )
            utau_icb(ji,jj) = 0.5_wp * ( utau(ji,jj) + utau(ji+1,jj) ) * &
               &                       ( 2. - umask(ji,jj,1) ) * MAX( tmask(ji,jj,1), tmask(ji+1,jj,1) ) * umask(ji,jj,1)
            vtau_icb(ji,jj) = 0.5_wp * ( vtau(ji,jj) + vtau(ji,jj+1) ) * &
               &                       ( 2. - vmask(ji,jj,1) ) * MAX( tmask(ji,jj,1), tmask(ji,jj+1,1) ) * vmask(ji,jj,1)
         END_2D
         CALL lbc_lnk( 'sbcmod', utau_icb, 'U', -1.0_wp, vtau_icb, 'V', -1.0_wp )
      ENDIF
      !
      !                                            !==  Misc. Options  ==!
      !
      SELECT CASE( nn_ice )                                            ! Update heat and freshwater fluxes over sea-ice areas
      CASE(  1 )   ;         CALL sbc_ice_if   ( kt, Kbb, Kmm )        ! Ice-cover climatology ("Ice-if" model)
#if defined key_si3
      CASE(  2 )   ;         CALL ice_stp  ( kt, Kbb, Kmm, nsbc )      ! SI3 ice model
#endif
      CASE(  3 )   ;         CALL sbc_ice_cice ( kt, nsbc )            ! CICE ice model
      END SELECT
      !==> clem: from here on, the following fields are ok on the halos: snwice_mass, snwice_mass_b, snwice_fmass
      !                        but not utau, vtau, emp (must be done later on)
      
      IF( ln_icebergs    )   CALL icb_stp( kt, Kmm )                   ! compute icebergs

      ! Icebergs do not melt over the haloes.
      ! So emp values over the haloes are no more consistent with the inner domain values.
      ! A lbc_lnk is therefore needed to ensure reproducibility and restartability.
      ! see ticket #2113 for discussion about this lbc_lnk.
      ! (PM) same consideration on qns (no heat flux from iceberg on haloes.
!!$      IF( ln_icebergs .AND. .NOT. ln_passive_mode )   CALL lbc_lnk( 'sbcmod', emp, 'T', 1.0_wp )
      !clem: not needed anymore since lbc is done afterwards
      
      IF( ln_rnf         )   CALL sbc_rnf( kt )                        ! add runoffs to fresh water fluxes

      IF( ln_ssr         )   CALL sbc_ssr( kt )                        ! add SST/SSS damping term

      IF( nn_fwb    /= 0 )   CALL sbc_fwb( kt, nn_fwb, nn_fsbc, Kmm )  ! control the freshwater budget

      ! Special treatment of freshwater fluxes over closed seas in the model domain
      ! Should not be run if ln_diurnal_only
      IF( l_sbc_clo      )   CALL sbc_clo( kt )

      IF( ll_wd ) THEN     ! If near WAD point limit the flux for now
         zthscl = atanh(rn_wd_sbcfra)                     ! taper frac default is .999
         DO_2D( 0, 0, 0, 0 )
            zwdht = ssh(ji,jj,Kmm) + ht_0(ji,jj) - rn_wdmin1   ! do this calc of water depth above wd limit once
            zwght = TANH(zthscl*zwdht)
            IF( zwdht <= 0.0 ) THEN
               taum(ji,jj) = 0.0
               utau(ji,jj) = 0.0
               vtau(ji,jj) = 0.0
               qns (ji,jj) = 0.0
               qsr (ji,jj) = 0.0
               emp (ji,jj) = MIN(emp(ji,jj),0.0) !can allow puddles to grow but not shrink
               sfx (ji,jj) = 0.0
            ELSEIF( zwdht > 0.0  .AND. zwdht < rn_wd_sbcdep ) THEN !  5 m hard limit here is arbitrary
               qsr  (ji,jj) =  qsr(ji,jj)  * zwght
               qns  (ji,jj) =  qns(ji,jj)  * zwght
               taum (ji,jj) =  taum(ji,jj) * zwght
               utau (ji,jj) =  utau(ji,jj) * zwght
               vtau (ji,jj) =  vtau(ji,jj) * zwght
               sfx  (ji,jj) =  sfx(ji,jj)  * zwght
               emp  (ji,jj) =  emp(ji,jj)  * zwght
            ENDIF
         END_2D
      ENDIF

      ! clem: these should be the only fields that are needed over the entire domain
      !       (in addition to snwice_mass)
      ! (PM): ldfull is required for iceberg (need to update all the halo and the inner band of the halo on the north fold)
      !       lbclnk on qns is needed because of the iceberg (need the inner band of the halo on the north fold to be correct)
      IF( ln_rnf ) THEN
         CALL lbc_lnk( 'sbcmod', utau, 'T', -1.0_wp, vtau  , 'T', -1.0_wp, emp, 'T', 1.0_wp, &
            &                    rnf , 'T',  1.0_wp, qns, 'T',  1.0_wp, ldfull=.TRUE. )
      ELSE
         CALL lbc_lnk( 'sbcmod', utau, 'T', -1.0_wp, vtau  , 'T', -1.0_wp, emp, 'T', 1.0_wp, qns, 'T',  1.0_wp, ldfull=.TRUE. )
      ENDIF
      ! --- calculate utau and vtau on U,V-points --- !
      !     Note the use of 0.5*(2-umask) in order to unmask the stress along coastlines
      !      and the use of MAX(tmask(i,j),tmask(i+1,j) is to mask tau over ice shelves
      DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
         utauU  (ji,jj) = 0.5_wp * ( utau(ji,jj) + utau(ji+1,jj) ) * &
            &                      ( 2. - umask(ji,jj,1) ) * MAX( tmask(ji,jj,1), tmask(ji+1,jj,1) )
         vtauV  (ji,jj) = 0.5_wp * ( vtau(ji,jj) + vtau(ji,jj+1) ) * &
            &                      ( 2. - vmask(ji,jj,1) ) * MAX( tmask(ji,jj,1), tmask(ji,jj+1,1) )
      END_2D
      !
      IF( kt == nit000 ) THEN                          !   set the forcing field at nit000 - 1    !
         !                                             ! ---------------------------------------- !
#if defined key_RK3
         IF( ln_rstart .AND. lk_SWE ) THEN                      !* RK3 + SWE: Restart: read in restart file
#else
         IF( ln_rstart .AND. .NOT.l_1st_euler ) THEN            !* MLF: Restart: read in restart file
#endif
            IF(lwp) WRITE(numout,*) '          nit000-1 surface forcing fields read in the restart file'
            CALL iom_get( numror, jpdom_auto, 'utau_b', utau_b, cd_type = 'U', psgn = -1._wp )   ! i-stress
            CALL iom_get( numror, jpdom_auto, 'vtau_b', vtau_b, cd_type = 'V', psgn = -1._wp )   ! j-stress
            CALL iom_get( numror, jpdom_auto,  'qns_b',  qns_b, cd_type = 'T', psgn =  1._wp )   ! non solar heat flux
            CALL iom_get( numror, jpdom_auto,  'emp_b',  emp_b, cd_type = 'T', psgn =  1._wp )   ! freshwater flux
            ! NB: The 3D heat content due to qsr forcing (qsr_hc_b) is treated in traqsr
            ! To ensure restart capability with 3.3x/3.4 restart files    !! to be removed in v3.6
            IF( iom_varid( numror, 'sfx_b', ldstop = .FALSE. ) > 0 ) THEN
               CALL iom_get( numror, jpdom_auto, 'sfx_b', sfx_b, cd_type = 'T', psgn = 1._wp )   ! before salt flux (T-point)
            ELSE
               sfx_b (:,:) = sfx(:,:)
            ENDIF
         ELSE                                                   !* no restart: set from nit000 values
            IF(lwp) WRITE(numout,*) '          nit000-1 surface forcing fields set to nit000'
            utau_b(:,:) = utauU(:,:)
            vtau_b(:,:) = vtauV(:,:)
            qns_b (:,:) = qns (:,:)
            emp_b (:,:) = emp (:,:)
            sfx_b (:,:) = sfx (:,:)
         ENDIF
         !
         IF( l_trddyn .and. nn_ice == 2 ) THEN
            uiceoc_b(:,:) = uiceoc(:,:)
            viceoc_b(:,:) = viceoc(:,:)
         ENDIF
      ENDIF
      !
      !
#if defined key_RK3
      !                                                ! ---------------------------------------- !
      IF( lrst_oce .AND. lk_SWE ) THEN                 !   RK3: Write in the ocean restart file   !
         !                                             ! ---------------------------------------- !
#else
      !                                                ! ---------------------------------------- !
      IF( lrst_oce ) THEN                              !   MLF: Write in the ocean restart file   !
         !                                             ! ---------------------------------------- !
#endif
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'sbc : ocean surface forcing fields written in ocean restart file ',   &
            &                    'at it= ', kt,' date= ', ndastp
         IF(lwp) WRITE(numout,*) '~~~~'
         CALL iom_rstput( kt, nitrst, numrow, 'utau_b' , utauU )
         CALL iom_rstput( kt, nitrst, numrow, 'vtau_b' , vtauV )
         CALL iom_rstput( kt, nitrst, numrow, 'qns_b'  , qns   )
         ! The 3D heat content due to qsr forcing is treated in traqsr
         ! CALL iom_rstput( kt, nitrst, numrow, 'qsr_b'  , qsr  )
         CALL iom_rstput( kt, nitrst, numrow, 'emp_b'  , emp   )
         CALL iom_rstput( kt, nitrst, numrow, 'sfx_b'  , sfx   )
      ENDIF
      !                                                ! ---------------------------------------- !
      !                                                !        Outputs and control print         !
      !                                                ! ---------------------------------------- !
      IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN
         IF( ln_rnf ) THEN
            CALL iom_put( "empmr"  , emp(A2D(0))-rnf(A2D(0))     )                ! upward water flux
            CALL iom_put( "empbmr" , emp_b(A2D(0))-rnf(A2D(0))   )                ! before upward water flux (for ssh in offline )
         ELSE
            CALL iom_put( "empmr"  , emp(A2D(0))    )          ! upward water flux
            CALL iom_put( "empbmr" , emp_b(A2D(0))  )          ! before upward water flux (for ssh in offline )
         ENDIF
         CALL iom_put( "saltflx", sfx         )                ! downward salt flux
         CALL iom_put( "fwfice" , fwfice      )                ! ice-ocean freshwater flux
         CALL iom_put( "qt"     , qns+qsr     )                ! total heat flux
         CALL iom_put( "qns"    , qns         )                ! solar heat flux
         CALL iom_put( "qsr"    , qsr         )                ! solar heat flux
         IF( nn_ice > 0 .OR. ll_opa )   CALL iom_put( "ice_cover", fr_i(:,:) )   ! ice fraction
         CALL iom_put( "taum"   , taum        )                ! wind stress module
         CALL iom_put( "wspd"   , wndm        )                ! wind speed  module over free ocean or leads in presence of sea-ice
         CALL iom_put( "qrp"    , qrp         )                ! heat flux damping
         CALL iom_put( "erp"    , erp         )                ! freshwater flux damping
      ENDIF
      !
      IF(sn_cfctl%l_prtctl) THEN     ! print mean trends (used for debugging)
         CALL prt_ctl(tab2d_1=fr_i                , clinfo1=' fr_i     - : ', mask1=tmask )
         IF( ln_rnf ) THEN
            CALL prt_ctl(tab2d_1=emp-rnf          , clinfo1=' emp-rnf  - : ', mask1=tmask )
            CALL prt_ctl(tab2d_1=sfx-rnf(A2D(0))  , clinfo1=' sfx-rnf  - : ', mask1=tmask )
         ELSE
            CALL prt_ctl(tab2d_1=emp              , clinfo1=' emp      - : ', mask1=tmask )
            CALL prt_ctl(tab2d_1=sfx              , clinfo1=' sfx      - : ', mask1=tmask )
         ENDIF
         CALL prt_ctl(tab2d_1=qns                 , clinfo1=' qns      - : ', mask1=tmask )
         CALL prt_ctl(tab2d_1=qsr                 , clinfo1=' qsr      - : ', mask1=tmask )
         CALL prt_ctl(tab3d_1=tmask               , clinfo1=' tmask    - : ', mask1=tmask, kdim=jpk )
         CALL prt_ctl(tab2d_1=sst_m               , clinfo1=' sst      - : ', mask1=tmask )
         CALL prt_ctl(tab2d_1=sss_m               , clinfo1=' sss      - : ', mask1=tmask )
         CALL prt_ctl(tab2d_1=utau                , clinfo1=' utau     - : ', mask1=tmask,                      &
            &         tab2d_2=vtau                , clinfo2=' vtau     - : ', mask2=tmask )
      ENDIF

      IF( kt == nitend )   CALL sbc_final         ! Close down surface module if necessary
      !
      IF( ln_timing )   CALL timing_stop('sbc')
      !
   END SUBROUTINE sbc


   SUBROUTINE sbc_final
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_final  ***
      !!
      !! ** Purpose :   Finalize CICE (if used)
      !!---------------------------------------------------------------------
      !
      IF( nn_ice == 3 )   CALL cice_sbc_final
      !
   END SUBROUTINE sbc_final

   !!======================================================================
END MODULE sbcmod

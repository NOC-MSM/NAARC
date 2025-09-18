MODULE dynspg_ts

   !! Includes ROMS wd scheme with diagnostic outputs ; puu(:,:,:,Kmm) and puu(:,:,:,Krhs) updates are commented out ! 

   !!======================================================================
   !!                   ***  MODULE  dynspg_ts  ***
   !! Ocean dynamics:  surface pressure gradient trend, split-explicit scheme
   !!======================================================================
   !! History :   1.0  ! 2004-12  (L. Bessieres, G. Madec)  Original code
   !!              -   ! 2005-11  (V. Garnier, G. Madec)  optimization
   !!              -   ! 2006-08  (S. Masson)  distributed restart using iom
   !!             2.0  ! 2007-07  (D. Storkey) calls to BDY routines
   !!              -   ! 2008-01  (R. Benshila)  change averaging method
   !!             3.2  ! 2009-07  (R. Benshila, G. Madec) Complete revisit associated to vvl reactivation
   !!             3.3  ! 2010-09  (D. Storkey, E. O'Dea) update for BDY for Shelf configurations
   !!             3.3  ! 2011-03  (R. Benshila, R. Hordoir, P. Oddo) update calculation of ub_b
   !!             3.5  ! 2013-07  (J. Chanut) Switch to Forward-backward time stepping
   !!             3.6  ! 2013-11  (A. Coward) Update for z-tilde compatibility
   !!             3.7  ! 2015-11  (J. Chanut) free surface simplification
   !!              -   ! 2016-12  (G. Madec, E. Clementi) update for Stoke-Drift divergence
   !!             4.0  ! 2017-05  (G. Madec)  drag coef. defined at t-point (zdfdrg.F90)
   !!---------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_spg_ts     : compute surface pressure gradient trend using a time-splitting scheme 
   !!   dyn_spg_ts_init: initialisation of the time-splitting scheme
   !!   ts_wgt         : set time-splitting weights for temporal averaging (or not)
   !!   ts_rst         : read/write time-splitting fields in restart file
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! surface boundary condition: ocean
   USE isf_oce         ! ice shelf variable (fwfisf)
   USE zdf_oce         ! vertical physics: variables
   USE zdfdrg          ! vertical physics: top/bottom drag coef.
   USE sbcapr          ! surface boundary condition: atmospheric pressure
   USE dynadv    , ONLY: ln_dynadv_vec
   USE dynvor          ! vortivity scheme indicators
   USE phycst          ! physical constants
   USE wet_dry         ! wetting/drying flux limter
   USE bdy_oce         ! open boundary
   USE bdyvol          ! open boundary volume conservation
   USE bdytides        ! open boundary condition data
   USE bdydyn2d        ! open boundary conditions on barotropic variables
   USE tide_mod        !
   USE sbcwave         ! surface wave
   USE daymod, ONLY : ndt05   ! half-length of time step
#if defined key_agrif
   USE agrif_oce_interp ! agrif
   USE agrif_oce
#endif
#if defined key_asminc   
   USE asminc          ! Assimilation increment
#endif
   !
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distributed memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control
   USE iom             ! IOM library
   USE restart         ! only for lrst_oce

   USE timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC dyn_spg_ts        ! called by dyn_spg 
   PUBLIC dyn_spg_ts_init   !    -    - dyn_spg_init
   PUBLIC dyn_drg_init      ! called by stp2d
   !
   INTEGER         ::   icycle         ! Number of barotropic sub-steps for each internal step nn_e <= 2.5 nn_e
   REAL(wp)        ::   rDt_e          ! Barotropic time step
   LOGICAL, PUBLIC ::   ll_cold_start  ! =T : missing restart fields  ==> cold start  =F : ==> hot start 
   REAL(wp) ::   r1_wgt1s       ! Temporal weights sums
   REAL(wp) ::   r1_wgt2s
   !
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:)   ::   wgtbtp1, wgtbtp2   ! 1st & 2nd weights used in time filtering of barotropic fields
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ffu_nw, ffu_ne       ! 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ffu_sw, ffu_se       ! 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ffv_nw, ffv_ne       ! 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ffv_sw, ffv_se       ! 


   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sshe_rhs           ! RHS of ssh Eq.
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   Ue_rhs, Ve_rhs    ! RHS of barotropic velocity Eq.
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   CdU_u, CdU_v       ! top/bottom stress at u- & v-points

   REAL(wp) ::   r1_12 = 1._wp / 12._wp   ! local ratios
   REAL(wp) ::   r1_8  = 0.125_wp         !
   REAL(wp) ::   r1_4  = 0.25_wp          !
   REAL(wp) ::   r1_2  = 0.5_wp           !

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION dyn_spg_ts_alloc()
      !!----------------------------------------------------------------------
      !!                  ***  routine dyn_spg_ts_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(2)
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( wgtbtp1(3*nn_e), wgtbtp2(3*nn_e), STAT=ierr(1) )
         !
      ALLOCATE( ffu_nw(A2D(0)) , ffu_ne(A2D(0)) , ffu_sw(A2D(0)) , ffu_se(A2D(0)), &
         &      ffv_nw(A2D(0)) , ffv_ne(A2D(0)) , ffv_sw(A2D(0)) , ffv_se(A2D(0)), STAT=ierr(2) )
         !
      !
      dyn_spg_ts_alloc = MAXVAL( ierr(:) )
      !
      CALL mpp_sum( 'dynspg_ts', dyn_spg_ts_alloc )
      IF( dyn_spg_ts_alloc /= 0 )   CALL ctl_stop( 'STOP', 'dyn_spg_ts_alloc: failed to allocate arrays' )
      !
   END FUNCTION dyn_spg_ts_alloc

   SUBROUTINE dyn_spg_ts( kt, Kbb, Kmm, pssh, puu_b, pvv_b, Kaa )
      !!----------------------------------------------------------------------
      !!
      !! ** Purpose : - Compute the now trend due to the explicit time stepping
      !!              of the quasi-linear barotropic system, and add it to the
      !!              general momentum trend. 
      !!
      !! ** Method  : - split-explicit schem (time splitting) :
      !!      Barotropic variables are advanced from internal time steps Kbb to Kaa
      !!      thanks to a generalized forward-backward time stepping (see ref. below).
      !!
      !! ** Action :
      !!      -Update the filtered free surface at step "n+1"      : pssh(:,:,Kaa)
      !!      -Update filtered barotropic velocities at step "n+1" : puu_b(:,:,:,Kaa), vv_b(:,:,:,Kaa)
      !!      -Compute barotropic advective fluxes at step "n"     : un_adv, vn_adv
      !!      These are used to advect tracers and are compliant with discrete
      !!      continuity equation taken at the baroclinic time steps. This 
      !!      ensures tracers conservation.
      !!
      !! References : Shchepetkin and McWilliams, Ocean Modelling, 2005. 
      !!---------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::  kt                  ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::  Kbb, Kmm, Kaa       ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpt), INTENT(inout) ::  pssh, puu_b, pvv_b  ! SSH and barotropic velocities at main time levels
      !
      INTEGER  ::   ji, jj, jn                ! dummy loop indices
      LOGICAL  ::   ll_init                   ! =T : special startup of 2d equations
      LOGICAL, SAVE :: ll_bt_av               ! =T : boxcard time averaging   =F : foreward backward dissipation
      REAL(wp) ::   z1_hu , z1_hv             ! local scalars
      REAL(wp) ::   zzsshu, zzsshv            !   -      -
      REAL(wp) ::   za0, za1, za2, za3        !   -      -
      REAL(wp) ::   zldg, zkmax               !   -      -
      REAL(wp) ::   zhu_bck, zhv_bck, zhdiv   !   -      -
      REAL(wp), DIMENSION(jpi,jpj) :: zu_trd, zu_spg
      REAL(wp), DIMENSION(jpi,jpj) :: zv_trd, zv_spg
      REAL(wp), DIMENSION(A2D(0) ) :: zu_frc, zv_frc
      REAL(wp), DIMENSION(jpi,jpj) :: zsshu_a, zhup2_e, zhtp2_e
      REAL(wp), DIMENSION(jpi,jpj) :: zsshv_a, zhvp2_e, zsshp2_e
      REAL(wp), DIMENSION(jpi,jpj) :: zCdU_u, zCdU_v   ! top/bottom stress at u- & v-points
      REAL(wp), DIMENSION(jpi,jpj) :: zhU, zhV         ! fluxes
      REAL(wp), DIMENSION(jpi,jpj) :: zks
      !
      REAL(wp) ::   zwdramp                     ! local scalar - only used if ln_wd_dl = .True. 

      INTEGER  :: iwdg, jwdg, kwdg   ! short-hand values for the indices of the output point

      REAL(wp) ::   zepsilon, zgamma            !   -      -
#if ! defined key_PSYCLONE_2p5p0
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ztwdmask, zuwdmask, zvwdmask ! ROMS wetting and drying masks at t,u,v points
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zuwdav2, zvwdav2    ! averages over the sub-steps of zuwdmask and zvwdmask
#else
      REAL(wp), DIMENSION(jpi,jpj)          :: ztwdmask, zuwdmask, zvwdmask ! ROMS wetting and drying masks at t,u,v points
      REAL(wp), DIMENSION(jpi,jpj)          :: zuwdav2, zvwdav2    ! averages over the sub-steps of zuwdmask and zvwdmask
#endif
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: z2d          ! 2D workspace
      REAL(wp) ::   zt0substep !   Time of day at the beginning of the time substep
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: ztidepotgradu, ztidepotgradv ! temp variables for tide potential gradients
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dyn_spg_ts')
      !
      !                                         !* Allocate temporary arrays
#if ! defined key_PSYCLONE_2p5p0
      IF( ln_wd_dl ) ALLOCATE( ztwdmask(jpi,jpj), zuwdmask(jpi,jpj), zvwdmask(jpi,jpj), zuwdav2(jpi,jpj), zvwdav2(jpi,jpj))
#endif
      !
      zwdramp = r_rn_wdmin1               ! simplest ramp 
!     zwdramp = 1._wp / (rn_wdmin2 - rn_wdmin1) ! more general ramp
      !                                         ! inverse of baroclinic time step 
      !
      ll_bt_av    = .TRUE.                      ! boxcard or no time-filtering if nn_bt_flt != 3
      IF( nn_bt_flt==3 )   ll_bt_av = .FALSE.
      !
      ll_init     = ll_bt_av                    ! if no time averaging, then no specific restart 
      !                                         ! time offset in steps for bdy data update
      !
      IF( kt == nit000 ) THEN                   !* initialisation
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_spg_ts : surface pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~   free surface with time splitting'
         IF(lwp) WRITE(numout,*)
         !                    ! read bb and b field or start from 0
         IF( nn_bt_flt==3 ) THEN    
            IF( ln_rstart ) THEN           ! init bb fields with restart
                IF( ll_cold_start ) THEN
                   ll_init=.TRUE.
                ELSE
                   ll_init=.FALSE.
                ENDIF
            ELSE                           ! init bb fields with 0
                ll_init=.TRUE.
            ENDIF
         ENDIF
         !
         !                    ! Set averaging weights and cycle length:
         CALL ts_wgt( ll_bt_av, icycle, wgtbtp1, wgtbtp2 )
         ! Save weights sums:
         r1_wgt1s = SUM(wgtbtp1)
         r1_wgt2s = SUM(wgtbtp2)
         !
         !
      ENDIF
      !
      !
      ! -----------------------------------------------------------------------------
      !  Phase 1 : Coupling between general trend and barotropic estimates (1st step)
      ! -----------------------------------------------------------------------------
      !
      !                    !====================================!
      !                    !==  Phase 1 for time integration  ==!
      !                    !====================================!
      !
      !                          ! set values computed in dynspg_ts
       ssh_frc(:,:) = sshe_rhs(:,:)
        zu_frc(:,:) =   Ue_rhs(:,:)
        zv_frc(:,:) =   Ve_rhs(:,:)
      zCdU_u  (:,:) = CdU_u   (:,:)
      zCdU_v  (:,:) = CdU_v   (:,:)

!!gm ==>>> !!ts    ISSUe her on en discute    changer les cas ENS ENE  et triad ?
      !
      !                                   !=  remove 2D Coriolis trend  =!
      !                                   !  --------------------------  !
      !
      IF( kt == nit000 .OR. .NOT. lk_linssh )   CALL dyn_cor_2D_init( Kmm )   ! Set zwz, the barotropic Coriolis force coefficient
      !                      ! recompute zwz = f/depth  at every time step for (.NOT.lk_linssh) as the water colomn height changes
      !
      !
      CALL dyn_cor_2D( puu_b(:,:,Kmm), pvv_b(:,:,Kmm), zu_trd, zv_trd )
      !
      DO_2D( 0, 0, 0, 0 )                          ! Remove coriolis term (and possibly spg) from barotropic trend
         zu_frc(ji,jj) = zu_frc(ji,jj) - zu_trd(ji,jj) * ssumask(ji,jj)
         zv_frc(ji,jj) = zv_frc(ji,jj) - zv_trd(ji,jj) * ssvmask(ji,jj)
      END_2D


      !                                   != Fill boundary data arrays for AGRIF
      !                                   ! ------------------------------------
#if defined key_agrif
         IF( .NOT.Agrif_Root() ) CALL agrif_dta_ts( kt )
#endif
      !
      ! -----------------------------------------------------------------------
      !  Phase 2 : Integration of the barotropic equations 
      ! -----------------------------------------------------------------------
      !
      !                                             ! ==================== !
      !                                             !    Initialisations   !
      !                                             ! ==================== !  
      ! Initialize barotropic variables:      
      IF( ll_init )THEN
         sshbb_e(:,:) = 0._wp
         ubb_e  (:,:) = 0._wp
         vbb_e  (:,:) = 0._wp
         sshb_e (:,:) = 0._wp
         ub_e   (:,:) = 0._wp
         vb_e   (:,:) = 0._wp
      ENDIF
      !
      IF( lk_linssh ) THEN    ! mid-step ocean depth is fixed (hup2_e=hu_n=hu_0)
         zhup2_e(:,:) = hu_0(:,:)
         zhvp2_e(:,:) = hv_0(:,:)
         zhtp2_e(:,:) = ht_0(:,:)
      ENDIF
      !
      sshn_e(:,:) =    pssh (:,:,Kmm)            
      un_e  (:,:) =    puu_b(:,:,Kmm)            
      vn_e  (:,:) =    pvv_b(:,:,Kmm)
      !
      hu_e  (:,:) =    hu(:,:,Kmm)       
      hv_e  (:,:) =    hv(:,:,Kmm) 
      hur_e (:,:) = r1_hu(:,:,Kmm)    
      hvr_e (:,:) = r1_hv(:,:,Kmm)
      !
      ! Initialize sums:
      puu_b (:,:,Kaa) = 0._wp       ! After barotropic velocities (or transport if flux form)          
      pvv_b (:,:,Kaa) = 0._wp
      pssh  (:,:,Kaa) = 0._wp       ! Sum for after averaged sea level
      un_adv(:,:)     = 0._wp       ! Sum for now transport issued from ts loop
      vn_adv(:,:)     = 0._wp
      !
      IF( ln_wd_dl ) THEN
         zuwdmask(:,:) = 0._wp  ! set to zero for definiteness (not sure this is necessary) 
         zvwdmask(:,:) = 0._wp  ! 
         zuwdav2 (:,:) = 0._wp 
         zvwdav2 (:,:) = 0._wp   
      END IF 

      !                                             ! ==================== !
      DO jn = 1, icycle                             !  sub-time-step loop  !
         !                                          ! ==================== !
         !
         !                    !==  Update the forcing ==! (BDY and tides)
         !
         IF( ln_bdy      .AND. ln_tide )   CALL bdy_dta_tides( kt, kit=jn, pt_offset= 1._wp )
         ! Update tide potential at the beginning of current time substep
         IF( ln_tide_pot .AND. ln_tide ) THEN
            zt0substep = REAL(nsec_day, wp) - 0.5_wp*rn_Dt + (jn + noffset - 1) * rn_Dt / REAL(nn_e, wp)
            CALL upd_tide(zt0substep, Kmm)
         END IF
         !
         !                    !==  extrapolation at mid-step  ==!   (jn+1/2)
         !
         !                       !* Set extrapolation coefficients for predictor step:
         IF ((jn<3).AND.ll_init) THEN      ! Forward           
           za1 = 1._wp                                          
           za2 = 0._wp                        
           za3 = 0._wp                        
         ELSE                              ! AB3-AM4 Coefficients: bet=0.281105 
           za1 =  1.781105_wp              ! za1 =   3/2 +   bet
           za2 = -1.06221_wp               ! za2 = -(1/2 + 2*bet)
           za3 =  0.281105_wp              ! za3 = bet
         ENDIF
         !
         !                       !* Extrapolate barotropic velocities at mid-step (jn+1/2)
         !--        m+1/2               m                m-1           m-2       --!
         !--       u      = (3/2+beta) u   -(1/2+2beta) u      + beta u          --!
         !-------------------------------------------------------------------------!
         DO_2D( 2, 1, 1, 1 )
            ua_e(ji,jj) = za1 * un_e(ji,jj) + za2 * ub_e(ji,jj) + za3 * ubb_e(ji,jj)
         END_2D
         DO_2D( 1, 1, 2, 1 )
            va_e(ji,jj) = za1 * vn_e(ji,jj) + za2 * vb_e(ji,jj) + za3 * vbb_e(ji,jj)
         END_2D

         IF( .NOT.lk_linssh ) THEN                        !* Update ocean depth (variable volume case only)
            !                                             !  ------------------
            ! Extrapolate Sea Level at step jit+0.5:
            !--         m+1/2                 m                  m-1             m-2       --!
            !--      ssh      = (3/2+beta) ssh   -(1/2+2beta) ssh      + beta ssh          --!
            !--------------------------------------------------------------------------------!
            zsshp2_e(:,:) = za1 * sshn_e(:,:)  + za2 * sshb_e(:,:) + za3 * sshbb_e(:,:)
            
            ! set wetting & drying mask at tracer points for this barotropic mid-step
            IF( ln_wd_dl )   CALL wad_tmsk( zsshp2_e, ztwdmask )
            !
            !                          ! ocean t-depth at mid-step
            DO_2D( 0, 1, 0, 1 )
               zhtp2_e(ji,jj) = ht_0(ji,jj) + zsshp2_e(ji,jj)
            END_2D
            !
            !                          ! ocean u- and v-depth at mid-step   (separate DO-loops remove the need of a lbc_lnk)
#if defined key_qcoTest_FluxForm
            !                                ! 'key_qcoTest_FluxForm' : simple ssh average
            DO_2D( 2, 1, 1, 1 )   ! not jpi-column
               zhup2_e(ji,jj) = hu_0(ji,jj) + r1_2 * (  zsshp2_e(ji,jj) + zsshp2_e(ji+1,jj  )  ) * ssumask(ji,jj)
            END_2D
            DO_2D( 1, 1, 2, 1 )
               zhvp2_e(ji,jj) = hv_0(ji,jj) + r1_2 * (  zsshp2_e(ji,jj) + zsshp2_e(ji  ,jj+1)  ) * ssvmask(ji,jj)
            END_2D
#else
            !                                ! no 'key_qcoTest_FluxForm' : surface weighted ssh average
            DO_2D( 2, 1, 1, 1 )   ! not jpi-column
               zhup2_e(ji,jj) = hu_0(ji,jj) + r1_2 * r1_e1e2u(ji,jj)                        &
                    &                              * (  e1e2t(ji  ,jj) * zsshp2_e(ji  ,jj)  &
                    &                                 + e1e2t(ji+1,jj) * zsshp2_e(ji+1,jj)  ) * ssumask(ji,jj)
            END_2D
            DO_2D( 1, 1, 2, 1 )   ! not jpj-row
               zhvp2_e(ji,jj) = hv_0(ji,jj) + r1_2 * r1_e1e2v(ji,jj)                        &
                    &                              * (  e1e2t(ji,jj  ) * zsshp2_e(ji,jj  )  &
                    &                                 + e1e2t(ji,jj+1) * zsshp2_e(ji,jj+1)  ) * ssvmask(ji,jj)
            END_2D
#endif               
            !
         ENDIF
         !
         !                    !==  after SSH  ==!   (jn+1)
         !
         !                             ! update (ua_e,va_e) to enforce volume conservation at open boundaries
         !                             ! values of zhup2_e and zhvp2_e on the halo are not needed in bdy_vol2d
         IF( ln_bdy .AND. ln_vol ) CALL bdy_vol2d( kt, jn, ua_e, va_e, zhup2_e, zhvp2_e )
         !      
         !                             ! resulting flux at mid-step (not over the full domain)
         DO_2D( 2, 1, 1, 1 )   ! not jpi-column
            zhU(ji,jj) = e2u(ji,jj) * ua_e(ji,jj) * zhup2_e(ji,jj)
         END_2D
         DO_2D( 1, 1, 2, 1 )   ! not jpj-row
            zhV(ji,jj) = e1v(ji,jj) * va_e(ji,jj) * zhvp2_e(ji,jj)
         END_2D
         !
#if defined key_agrif
         ! Set fluxes during predictor step to ensure volume conservation
         CALL agrif_dyn_ts_flux( jn, zhU, zhV )
#endif

         IF( ln_wd_dl ) THEN           ! un_e and vn_e are set to zero at faces where 
            !                          ! the direction of the flow is from dry cells
            CALL wad_Umsk( ztwdmask, zhU, zhV, un_e, vn_e, zuwdmask, zvwdmask )   ! not jpi colomn for U, not jpj row for V
            !
         ENDIF    
         !
         !
         !     Compute Sea Level at step jit+1
         !--           m+1        m                               m+1/2          --!
         !--        ssh    =  ssh   - delta_t' * [ frc + div( flux      ) ]      --!
         !-------------------------------------------------------------------------!
         DO_2D( 1, 1, 1, 1 )
            zhdiv = (   ( zhU(ji,jj) - zhU(ji-1,jj) ) + ( zhV(ji,jj) - zhV(ji,jj-1) )  ) * r1_e1e2t(ji,jj)
            ssha_e(ji,jj) = (  sshn_e(ji,jj) - rDt_e * ( ssh_frc(ji,jj) + zhdiv )  ) * ssmask(ji,jj)
         END_2D
         !
         ! Duplicate sea level across open boundaries (this is only cosmetic if linssh=T)
         IF( ln_bdy )   CALL bdy_ssh( ssha_e )
#if defined key_agrif
         CALL agrif_ssh_ts( jn )
#endif
         !
         !                             ! Sum dver sub-time-steps to compute advective velocities
         za2 = wgtbtp2(jn)             ! zhU, zhV hold fluxes extrapolated at jn+0.5
         DO_2D( 0, 0, 0, 0 )
            un_adv(ji,jj) = un_adv(ji,jj) + za2 * zhU(ji,jj) * r1_e2u(ji,jj)
            vn_adv(ji,jj) = vn_adv(ji,jj) + za2 * zhV(ji,jj) * r1_e1v(ji,jj)
         END_2D
         ! sum over sub-time-steps to decide which baroclinic velocities to set to zero (zuwdav2 is only used when ln_wd_dl_bc=True) 
         IF ( ln_wd_dl_bc ) THEN
            DO_2D( 0, 0, 0, 0 )
               zuwdav2(ji,jj) = zuwdav2(ji,jj) + za2 * zuwdmask(ji,jj)
               zvwdav2(ji,jj) = zvwdav2(ji,jj) + za2 * zvwdmask(ji,jj)
            END_2D
         END IF
         !
         !  
         ! Sea Surface Height at u-,v-points (vvl case only)
         IF( .NOT.lk_linssh ) THEN
#if defined key_qcoTest_FluxForm
            !                                ! 'key_qcoTest_FluxForm' : simple ssh average
            DO_2D( 0, 0, 0, 0 )
               zsshu_a(ji,jj) = r1_2 * (  ssha_e(ji,jj) + ssha_e(ji+1,jj  )  ) * ssumask(ji,jj)
               zsshv_a(ji,jj) = r1_2 * (  ssha_e(ji,jj) + ssha_e(ji  ,jj+1)  ) * ssvmask(ji,jj)
            END_2D
#else
            DO_2D( 0, 0, 0, 0 )
               zsshu_a(ji,jj) = r1_2 * r1_e1e2u(ji,jj) * ( e1e2t(ji  ,jj  ) * ssha_e(ji  ,jj  )   &
                  &                                      + e1e2t(ji+1,jj  ) * ssha_e(ji+1,jj  ) ) * ssumask(ji,jj)
               zsshv_a(ji,jj) = r1_2 * r1_e1e2v(ji,jj) * ( e1e2t(ji  ,jj  ) * ssha_e(ji  ,jj  )   &
                  &                                      + e1e2t(ji  ,jj+1) * ssha_e(ji  ,jj+1) ) * ssvmask(ji,jj)
            END_2D
#endif
         ENDIF
         !         
         ! Half-step back interpolation of SSH for surface pressure computation at step jit+1/2
         !--            m+1/2           m+1              m               m-1              m-2     --!
         !--        ssh'    =  za0 * ssh     +  za1 * ssh   +  za2 * ssh      +  za3 * ssh        --!
         !------------------------------------------------------------------------------------------!
         CALL ts_bck_interp( jn, ll_init, za0, za1, za2, za3 )   ! coeficients of the interpolation
         DO_2D( 1, 1, 1, 1 )
            zsshp2_e(ji,jj) = za0 *  ssha_e(ji,jj) + za1 *  sshn_e (ji,jj)   &
               &            + za2 *  sshb_e(ji,jj) + za3 *  sshbb_e(ji,jj)
         END_2D
         !
         !                             ! Surface pressure gradient
         zldg = ( 1._wp - rn_scal_load ) * grav    ! local factor
         DO_2D( 0, 0, 0, 0 )
            zu_spg(ji,jj) = - zldg * ( zsshp2_e(ji+1,jj) - zsshp2_e(ji,jj) ) * r1_e1u(ji,jj)
            zv_spg(ji,jj) = - zldg * ( zsshp2_e(ji,jj+1) - zsshp2_e(ji,jj) ) * r1_e2v(ji,jj)
         END_2D
         !
         ! Add Coriolis trend:
         ! zwz array below or triads normally depend on sea level with lk_linssh=F and should be updated
         ! at each time step. We however keep them constant here for optimization.
         ! Recall that zhU and zhV hold fluxes at jn+0.5 (extrapolated not backward interpolated)
         CALL dyn_cor_2D( ua_e, va_e, zu_trd, zv_trd )
         !
         ! Add tidal astronomical forcing if defined
         IF ( ln_tide .AND. ln_tide_pot ) THEN
            IF(l_trddyn) THEN
               ALLOCATE( ztidepotgradu(jpi,jpj), ztidepotgradv(jpi,jpj) )            
               ztidepotgradu(:,:) = zu_trd(:,:)
               ztidepotgradv(:,:) = zv_trd(:,:) 
            ENDIF
            !
            DO_2D( 0, 0, 0, 0 )
               zu_trd(ji,jj) = zu_trd(ji,jj) + grav * ( pot_astro(ji+1,jj) - pot_astro(ji,jj) ) * r1_e1u(ji,jj)
               zv_trd(ji,jj) = zv_trd(ji,jj) + grav * ( pot_astro(ji,jj+1) - pot_astro(ji,jj) ) * r1_e2v(ji,jj)
            END_2D
            !
            IF( l_trddyn ) THEN
               za2 = wgtbtp2(jn)
               ztidepotgradu(:,:) = zu_trd(:,:) - ztidepotgradu(:,:)
               ztidepotgradv(:,:) = zv_trd(:,:) - ztidepotgradv(:,:)
               zspgtrdu(:,:) = zspgtrdu + za2 * ztidepotgradu(:,:) * ssumask(:,:)
               zspgtrdv(:,:) = zspgtrdv + za2 * ztidepotgradv(:,:) * ssvmask(:,:)
               DEALLOCATE( ztidepotgradu, ztidepotgradv )
            ENDIF 
         ENDIF
         !
         ! Add bottom stresses:
!jth do implicitly instead
         IF ( .NOT. ll_wd ) THEN ! Revert to explicit for bit comparison tests in non wad runs
            DO_2D( 0, 0, 0, 0 )
               zu_trd(ji,jj) = zu_trd(ji,jj) + zCdU_u(ji,jj) * un_e(ji,jj) * hur_e(ji,jj)
               zv_trd(ji,jj) = zv_trd(ji,jj) + zCdU_v(ji,jj) * vn_e(ji,jj) * hvr_e(ji,jj)
            END_2D
         ENDIF
         !
         ! Set next velocities:
         !     Compute barotropic speeds at step jit+1    (h : total height of the water colomn)
         !--                              VECTOR FORM
         !--   m+1                 m               /                                                       m+1/2           \    --!
         !--  u     =             u   + delta_t' * \         (1-r)*g * grad_x( ssh') -         f * k vect u      +     frc /    --!
         !--                                                                                                                    --!
         !--                             FLUX FORM                                                                              --!
         !--  m+1   __1__  /  m    m               /  m+1/2                             m+1/2              m+1/2    n      \ \  --!
         !-- u    =   m+1 |  h  * u   + delta_t' * \ h     * (1-r)*g * grad_x( ssh') - h     * f * k vect u      + h * frc /  | --!
         !--         h     \                                                                                                 /  --!
         !------------------------------------------------------------------------------------------------------------------------!
         IF( ln_dynadv_vec .OR. lk_linssh ) THEN      !* Vector form
            DO_2D( 0, 0, 0, 0 )
               ua_e(ji,jj) = (                                 un_e(ji,jj)   & 
                         &     + rDt_e * (                   zu_spg(ji,jj)   &
                         &                                 + zu_trd(ji,jj)   &
                         &                                 + zu_frc(ji,jj) ) & 
                         &   ) * ssumask(ji,jj)

               va_e(ji,jj) = (                                 vn_e(ji,jj)   &
                         &     + rDt_e * (                   zv_spg(ji,jj)   &
                         &                                 + zv_trd(ji,jj)   &
                         &                                 + zv_frc(ji,jj) ) &
                         &   ) * ssvmask(ji,jj)
            END_2D
            !
         ELSE                           !* Flux form
            DO_2D( 0, 0, 0, 0 )
               !                    ! hu_e, hv_e hold depth at jn,  zhup2_e, zhvp2_e hold extrapolated depth at jn+1/2
               !                    ! backward interpolated depth used in spg terms at jn+1/2
#if defined key_qcoTest_FluxForm
            !                                ! 'key_qcoTest_FluxForm' : simple ssh average
               zhu_bck = hu_0(ji,jj) + r1_2 * (  zsshp2_e(ji,jj) + zsshp2_e(ji+1,jj  )  ) * ssumask(ji,jj)
               zhv_bck = hv_0(ji,jj) + r1_2 * (  zsshp2_e(ji,jj) + zsshp2_e(ji  ,jj+1)  ) * ssvmask(ji,jj)
#else
               zhu_bck = hu_0(ji,jj) + r1_2*r1_e1e2u(ji,jj) * (  e1e2t(ji  ,jj) * zsshp2_e(ji  ,jj)    &
                    &                                          + e1e2t(ji+1,jj) * zsshp2_e(ji+1,jj)  ) * ssumask(ji,jj)
               zhv_bck = hv_0(ji,jj) + r1_2*r1_e1e2v(ji,jj) * (  e1e2t(ji,jj  ) * zsshp2_e(ji,jj  )    &
                    &                                          + e1e2t(ji,jj+1) * zsshp2_e(ji,jj+1)  ) * ssvmask(ji,jj)
#endif
               !                    ! inverse depth at jn+1
               z1_hu = ssumask(ji,jj) / ( hu_0(ji,jj) + zsshu_a(ji,jj) + 1._wp - ssumask(ji,jj) )
               z1_hv = ssvmask(ji,jj) / ( hv_0(ji,jj) + zsshv_a(ji,jj) + 1._wp - ssvmask(ji,jj) )
               !
               ua_e(ji,jj) = (               hu_e  (ji,jj) *   un_e (ji,jj)      & 
                    &            + rDt_e * (  zhu_bck        * zu_spg (ji,jj)  &   !
                    &                       + zhup2_e(ji,jj) * zu_trd (ji,jj)  &   !
                    &                       +  hu(ji,jj,Kmm) * zu_frc (ji,jj)  )   ) * z1_hu
               !
               va_e(ji,jj) = (               hv_e  (ji,jj) *   vn_e (ji,jj)      &
                    &            + rDt_e * (  zhv_bck        * zv_spg (ji,jj)  &   !
                    &                       + zhvp2_e(ji,jj) * zv_trd (ji,jj)  &   !
                    &                       +  hv(ji,jj,Kmm) * zv_frc (ji,jj)  )   ) * z1_hv
            END_2D
         ENDIF
!jth implicit bottom friction:
         IF ( ll_wd ) THEN ! revert to explicit for bit comparison tests in non wad runs
            DO_2D( 0, 0, 0, 0 )
               ua_e(ji,jj) =  ua_e(ji,jj) / ( 1._wp - rDt_e * zCdU_u(ji,jj) * hur_e(ji,jj) )
               va_e(ji,jj) =  va_e(ji,jj) / ( 1._wp - rDt_e * zCdU_v(ji,jj) * hvr_e(ji,jj) )
            END_2D
         ENDIF
       
         IF( .NOT.lk_linssh ) THEN !* Update ocean depth (variable volume case only)
            DO_2D( 0, 0, 0, 0 )
               hu_e (ji,jj) =    hu_0(ji,jj) + zsshu_a(ji,jj)
               hur_e(ji,jj) = ssumask(ji,jj) / (  hu_e(ji,jj) + 1._wp - ssumask(ji,jj)  )
               hv_e (ji,jj) =    hv_0(ji,jj) + zsshv_a(ji,jj)
               hvr_e(ji,jj) = ssvmask(ji,jj) / (  hv_e(ji,jj) + 1._wp - ssvmask(ji,jj)  )
            END_2D
         ENDIF
         !
         IF( .NOT.lk_linssh ) THEN   !* Update ocean depth (variable volume case only)
            IF( ln_wd_dl ) THEN
               CALL lbc_lnk( 'dynspg_ts', ua_e , 'U', -1._wp, va_e , 'V', -1._wp  &
                    &                   , hu_e , 'U',  1._wp, hv_e , 'V',  1._wp  &
                    &                   , hur_e, 'U',  1._wp, hvr_e, 'V',  1._wp  , ssha_e, 'T', 1._wp  & 
                    &                   , zuwdmask, 'U', -1._wp, zvwdmask, 'V', -1._wp, ldfull=.TRUE. )
            ELSE
               CALL lbc_lnk( 'dynspg_ts', ua_e , 'U', -1._wp, va_e , 'V', -1._wp  &
                    &                   , hu_e , 'U',  1._wp, hv_e , 'V',  1._wp  &
                    &                   , hur_e, 'U',  1._wp, hvr_e, 'V',  1._wp  , ssha_e, 'T', 1._wp, ldfull=.TRUE. )
            ENDIF
         ELSE
            IF( ln_wd_dl ) THEN
               CALL lbc_lnk( 'dynspg_ts', ua_e , 'U', -1._wp, va_e , 'V', -1._wp  , ssha_e, 'T', 1._wp  &
                    &                   , zuwdmask, 'U', -1._wp, zvwdmask, 'V', -1._wp, ldfull=.TRUE. )
            ELSE
               CALL lbc_lnk( 'dynspg_ts', ua_e , 'U', -1._wp, va_e , 'V', -1._wp  , ssha_e, 'T', 1._wp, ldfull=.TRUE. )
            ENDIF
         ENDIF
         !                                                 ! open boundaries
         IF( ln_bdy )   CALL bdy_dyn2d( jn, ua_e, va_e, un_e, vn_e, hur_e, hvr_e, ssha_e )
#if defined key_agrif                                                           
         CALL agrif_dyn_ts( jn )  ! Agrif
#endif
         !                                             !* Swap
         !                                             !  ----
         ubb_e  (:,:) = ub_e  (:,:)
         ub_e   (:,:) = un_e  (:,:)
         un_e   (:,:) = ua_e  (:,:)
         !
         vbb_e  (:,:) = vb_e  (:,:)
         vb_e   (:,:) = vn_e  (:,:)
         vn_e   (:,:) = va_e  (:,:)
         !
         sshbb_e(:,:) = sshb_e(:,:)
         sshb_e (:,:) = sshn_e(:,:)
         sshn_e (:,:) = ssha_e(:,:)
         !
         !                                             !* Sum over whole bt loop (except in weight average)
         !                                             !  ----------------------
         IF( ll_bt_av ) THEN
            za1 = wgtbtp1(jn)                                    
            IF( ln_dynadv_vec .OR. lk_linssh ) THEN    ! Sum velocities
               puu_b  (:,:,Kaa) = puu_b  (:,:,Kaa) + za1 * ua_e  (:,:) 
               pvv_b  (:,:,Kaa) = pvv_b  (:,:,Kaa) + za1 * va_e  (:,:) 
            ELSE                                       ! Sum transports
               IF ( .NOT.ln_wd_dl ) THEN  
                  puu_b  (:,:,Kaa) = puu_b  (:,:,Kaa) + za1 * ua_e  (:,:) * hu_e (:,:)
                  pvv_b  (:,:,Kaa) = pvv_b  (:,:,Kaa) + za1 * va_e  (:,:) * hv_e (:,:)
               ELSE 
                  puu_b  (:,:,Kaa) = puu_b  (:,:,Kaa) + za1 * ua_e  (:,:) * hu_e (:,:) * zuwdmask(:,:)
                  pvv_b  (:,:,Kaa) = pvv_b  (:,:,Kaa) + za1 * va_e  (:,:) * hv_e (:,:) * zvwdmask(:,:)
               ENDIF
            ENDIF
            !                                          ! Sum sea level
            pssh(:,:,Kaa) = pssh(:,:,Kaa) + za1 * ssha_e(:,:)
         ENDIF
         !
         !                                                 ! ==================== !
      END DO                                               !        end loop      !
      !                                                    ! ==================== !
      ! Finalize sums:
      ! (Divide by weights sums here to prevent from truncation errors)
      un_adv(:,:) = un_adv(:,:) / r1_wgt2s 
      vn_adv(:,:) = vn_adv(:,:) / r1_wgt2s 
      puu_b  (:,:,Kaa) = puu_b  (:,:,Kaa) / r1_wgt1s 
      pvv_b  (:,:,Kaa) = pvv_b  (:,:,Kaa) / r1_wgt1s 
      pssh   (:,:,Kaa) = pssh   (:,:,Kaa) / r1_wgt1s

      IF( ln_wd_dl .AND. ln_wd_dl_bc ) THEN
         zuwdav2(:,:) = zuwdav2(:,:) / r1_wgt2s
         zvwdav2(:,:) = zvwdav2(:,:) / r1_wgt2s
         CALL lbc_lnk( 'dynspg_ts', un_adv, 'U', -1._wp, vn_adv, 'V', -1._wp, zuwdav2, 'U', -1._wp, zvwdav2, 'V', -1._wp ) ! Boundary conditions
      ELSE
         CALL lbc_lnk( 'dynspg_ts', un_adv, 'U', -1._wp, vn_adv, 'V', -1._wp ) ! Boundary conditions
      ENDIF
      
      ! -----------------------------------------------------------------------------
      ! Phase 3. update the general trend with the barotropic trend
      ! -----------------------------------------------------------------------------
      !
      IF(.NOT.ll_bt_av ) THEN                          !* Update Kaa barotropic external mode 
         puu_b(:,:,Kaa) = ua_e  (:,:)
         pvv_b(:,:,Kaa) = va_e  (:,:)
         pssh (:,:,Kaa) = ssha_e(:,:)
      ENDIF
      !
      IF( (.NOT.(ln_dynadv_vec .OR. lk_linssh)) .AND. ll_bt_av ) THEN                  ! at this stage, pssh(:,:,:,Krhs) has been corrected: compute new depths at velocity points
         !
# if defined key_qcoTest_FluxForm
         !                                       ! 'key_qcoTest_FluxForm' : simple ssh average
         DO_2D( 0, 0, 0, 0 )
            zzsshu = r1_2 * ( pssh(ji,jj,Kaa) + pssh(ji+1,jj  ,Kaa) ) * ssumask(ji,jj)
            zzsshv = r1_2 * ( pssh(ji,jj,Kaa) + pssh(ji  ,jj+1,Kaa) ) * ssvmask(ji,jj)
            !
            !                              ! Save barotropic velocities (not transport)
            puu_b(ji,jj,Kaa) = puu_b(ji,jj,Kaa) / ( hu_0(ji,jj) + zzsshu + 1._wp - ssumask(ji,jj) )
            pvv_b(ji,jj,Kaa) = pvv_b(ji,jj,Kaa) / ( hv_0(ji,jj) + zzsshv + 1._wp - ssvmask(ji,jj) )
         END_2D
# else
         DO_2D( 0, 0, 0, 0 )
            zzsshu = r1_2 * r1_e1e2u(ji,jj) * ( e1e2t(ji  ,jj) * pssh(ji  ,jj,Kaa)   &
               &                              + e1e2t(ji+1,jj) * pssh(ji+1,jj,Kaa) ) * ssumask(ji,jj)
            zzsshv = r1_2 * r1_e1e2v(ji,jj) * ( e1e2t(ji,jj  ) * pssh(ji,jj  ,Kaa)   &
               &                              + e1e2t(ji,jj+1) * pssh(ji,jj+1,Kaa) ) * ssvmask(ji,jj)
            !
            !                             ! Save barotropic velocities (not transport)
            puu_b(ji,jj,Kaa) = puu_b(ji,jj,Kaa) / ( hu_0(ji,jj) + zzsshu + 1._wp - ssumask(ji,jj) )
            pvv_b(ji,jj,Kaa) = pvv_b(ji,jj,Kaa) / ( hv_0(ji,jj) + zzsshv + 1._wp - ssvmask(ji,jj) )
         END_2D
# endif   
         !
         CALL lbc_lnk( 'dynspg_ts', puu_b, 'U', -1._wp, pvv_b, 'V', -1._wp ) ! Boundary conditions
         !
      ENDIF
      !
      IF( iom_use("ubar") ) THEN    ! single first: hu[N+1/2] = 1/2 ( hu[N] + hu[N+1] ) 
         ALLOCATE( z2d(jpi,jpj) )
         z2d(:,:) = 2._wp / ( hu_e(:,:) + hu(:,:,Kbb) + 1._wp - ssumask(:,:) ) 
         CALL iom_put(  "ubar", un_adv(:,:)*z2d(:,:) )    ! barotropic i-current
         z2d(:,:) = 2._wp / ( hv_e(:,:) + hv(:,:,Kbb) + 1._wp - ssvmask(:,:) )
         CALL iom_put(  "vbar", vn_adv(:,:)*z2d(:,:) )    ! barotropic i-current
         DEALLOCATE( z2d )
      ENDIF
      !
      !                    !==  END Phase 3 (forward mode) ==!
      !
      !
#if defined key_agrif
      !
      ! Save time integrated fluxes during child grid integration
      ! (used to update coarse grid transports at next time step)
      !
      IF( .NOT.Agrif_Root() .AND. ln_agrif_2way ) THEN
         IF( Agrif_NbStepint() == 0 ) THEN
            ub2_i_b(:,:) = 0._wp
            vb2_i_b(:,:) = 0._wp
         END IF
         !
         za1 = 1._wp / REAL(Agrif_rhot(), wp)
         ub2_i_b(:,:) = ub2_i_b(:,:) + za1 * un_adv(:,:)
         vb2_i_b(:,:) = vb2_i_b(:,:) + za1 * vn_adv(:,:)
      ENDIF
#endif
      !                                   !: write time-spliting arrays in the restart
      IF( lrst_oce )   CALL ts_rst( kt, 'WRITE' )
      !
#if ! defined key_PSYCLONE_2p5p0
      IF( ln_wd_dl )   DEALLOCATE( ztwdmask, zuwdmask, zvwdmask, zuwdav2, zvwdav2 )
#endif
      !
      CALL iom_put( "baro_u" , puu_b(:,:,Kmm) )  ! Barotropic  U Velocity
      CALL iom_put( "baro_v" , pvv_b(:,:,Kmm) )  ! Barotropic  V Velocity
      !
      IF( ln_timing )   CALL timing_stop('dyn_spg_ts')
      !
   END SUBROUTINE dyn_spg_ts

   
   SUBROUTINE ts_wgt( ll_av, Kpit, zwgt1, zwgt2)
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE ts_wgt  ***
      !!
      !! ** Purpose : Set time-splitting weights for temporal averaging (or not)
      !!----------------------------------------------------------------------
      LOGICAL, INTENT(in   ) ::   ll_av     ! temporal averaging=.true.
      INTEGER, INTENT(inout) ::   Kpit      ! cycle length
      !!
      INTEGER ::  jic, jn, ji   ! local integers
      REAL(wp) :: za1, za2      ! loca scalars
      REAL(wp), DIMENSION(3*nn_e), INTENT(inout) ::   zwgt1, zwgt2   ! Primary & Secondary weights
      !!----------------------------------------------------------------------
      !
      zwgt1(:) = 0._wp
      zwgt2(:) = 0._wp
      !
      !                          !==  Set time index when averaged value is requested  ==!
      jic = nn_e
      !
      !                          !==  Set primary weights  ==!
      !
      IF (ll_av) THEN               != Define simple boxcar window for primary weights 
         !                                       ! (width = nn_e, centered around jic)     
         SELECT CASE( nn_bt_flt )
         !
         CASE( 0 )                  ! No averaging
            zwgt1(jic) = 1._wp
            Kpit = jic
            !
         CASE( 1 )                  ! Boxcar, width = nn_e
            DO jn = 1, 3*nn_e
               za1 = ABS( REAL( jn-jic, wp) ) / REAL( nn_e, wp ) 
               IF( za1 < 0.5_wp ) THEN
                  zwgt1(jn) = 1._wp
                  Kpit = jn
               ENDIF
            END DO
            !
         CASE( 2 )                  ! Boxcar, width = 2 * nn_e
            DO jn = 1, 3*nn_e
               za1 = ABS(REAL( jn-jic, wp) ) / REAL( nn_e, wp ) 
               IF( za1 < 1._wp ) THEN
                  zwgt1(jn) = 1._wp
                  Kpit = jn
               ENDIF
            END DO
            !
         CASE DEFAULT   ;   CALL ctl_stop( 'unrecognised value for nn_bt_flt' )
         !
         END SELECT

      ELSE                          !=  No time averaging
         zwgt1(jic) = 1._wp
         Kpit = jic
      ENDIF
      !
      !                          !==  Set secondary weights  ==!
      DO jn = 1, Kpit
         DO ji = jn, Kpit
            zwgt2(jn) = zwgt2(jn) + zwgt1(ji)
         END DO
      END DO
      !
   END SUBROUTINE ts_wgt


   SUBROUTINE ts_rst( kt, cdrw )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE ts_rst  ***
      !!
      !! ** Purpose : Read or write time-splitting arrays in restart file
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt     ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      !!----------------------------------------------------------------------
      !
      IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialise 
         !                                   ! ---------------
         IF( ln_rstart ) THEN                           !* Read the restart file
            IF( nn_bt_flt == 3 ) THEN
               IF( iom_varid( numror, 'sshbb_e', ldstop = .FALSE. ) > 0 ) THEN
                  CALL iom_get( numror, jpdom_auto, 'sshbb_e'  , sshbb_e(:,:), cd_type = 'T', psgn =  1._wp )
                  CALL iom_get( numror, jpdom_auto, 'ubb_e'    ,   ubb_e(:,:), cd_type = 'U', psgn = -1._wp )
                  CALL iom_get( numror, jpdom_auto, 'vbb_e'    ,   vbb_e(:,:), cd_type = 'V', psgn = -1._wp )
                  CALL iom_get( numror, jpdom_auto, 'sshb_e'   ,  sshb_e(:,:), cd_type = 'T', psgn =  1._wp )
                  CALL iom_get( numror, jpdom_auto, 'ub_e'     ,    ub_e(:,:), cd_type = 'U', psgn = -1._wp )
                  CALL iom_get( numror, jpdom_auto, 'vb_e'     ,    vb_e(:,:), cd_type = 'V', psgn = -1._wp )
               ELSE    ! if restart does not contain these variables, cold barotropic start
                  ll_cold_start = .TRUE. 
               ENDIF
            ENDIF
         ENDIF
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
         !                                   ! -------------------
         IF(lwp) WRITE(numout,*) '---- ts_rst ----'
         !
         IF( nn_bt_flt == 3 ) THEN
            CALL iom_rstput( kt, nitrst, numrow, 'sshbb_e'  , sshbb_e(:,:) ) 
            CALL iom_rstput( kt, nitrst, numrow, 'ubb_e'    ,   ubb_e(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'vbb_e'    ,   vbb_e(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'sshb_e'   ,  sshb_e(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'ub_e'     ,    ub_e(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'vb_e'     ,    vb_e(:,:) )
         ENDIF
      ENDIF
      !
   END SUBROUTINE ts_rst


   SUBROUTINE dyn_spg_ts_init
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE dyn_spg_ts_init  ***
      !!
      !! ** Purpose : Set time splitting options
      !!----------------------------------------------------------------------
      INTEGER  ::   ji ,jj              ! dummy loop indices
      REAL(wp) ::   zxr2, zyr2, zcmax   ! local scalar
      REAL(wp) ::   zc0max, zzc0, zzc1, zzd   ! local scalar
      REAL(wp), DIMENSION(jpi,jpj) ::   zcu
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'dyn_spg_ts_init : split-explicit free surface'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~'
      !
      IF( ( nn_bt_flt==1 .OR. nn_bt_flt==2 ) .AND. ( rn_bt_alpha>0._wp ) ) THEN
         rn_bt_alpha = 0._wp
         CALL ctl_warn( 'dynspg_ts: remove temporal averaging when nn_bt_flt=1 or 2 ==> rn_bt_alpha = 0' )
      ENDIF
      IF( nn_bt_flt==3 .AND. rn_bt_alpha<=0._wp ) &
         &   CALL ctl_stop( 'dynspg_ts ERROR: if nn_bt_flt=3, then rn_bt_alpha must be /=0 (typical value=0.07)' )
      !
      !
      ! ---------------------------------------
      ! Max courant number for ext. grav. waves
      ! ---------------------------------------
      DO_2D( 0, 0, 0, 0 )
         zxr2 = r1_e1t(ji,jj) * r1_e1t(ji,jj)
         zyr2 = r1_e2t(ji,jj) * r1_e2t(ji,jj)
         zcu(ji,jj) = SQRT( grav * MAX( ht_0(ji,jj), 0._wp ) * ( zxr2 + zyr2 ) )
      END_2D
      !
#if defined key_agrif
      ! Discard points that are not stepped by 2d mode:
      zcu(Nis0:Nie0,Njs0:Nje0) = zcu(Nis0:Nie0,Njs0:Nje0) * ( 1._wp - tmask_upd(Nis0:Nie0,Njs0:Nje0) )
#endif
      !
      zcmax = MAXVAL( zcu(Nis0:Nie0,Njs0:Nje0) )
      CALL mpp_max( 'dynspg_ts', zcmax )
      !
      ! -----------------------------
      ! Estimate number of iterations   ( to satisfy a max courant number = rn_bt_cmax )
      ! -----------------------------
      IF( ln_bt_auto )   nn_e = CEILING( rn_Dt / rn_bt_cmax * zcmax)

      ! Estimate number of iterations and FB dissipative parameter (ongoing work) => not working for now
      zc0max = rn_Dt * zcmax
      IF( nn_bt_flt==3 .AND. rn_bt_alpha==0._wp ) THEN
         zzc0 = zc0max*zc0max
         zzc1 = rpi * rpi * 3._wp / 4._wp
         zzd  = SQRT(49284._wp * zzc0 - 375._wp * zzc1*zzc1) - 35._wp * zzc1
         nn_e = 240._wp * zzc0 / zzd
         rn_bt_alpha = zzc1 / zzc0 * nn_e / 6._wp
      ENDIF
      !
      ! --------------
      ! Courant Number
      ! --------------
      rDt_e = rn_Dt / REAL( nn_e , wp )
      zcmax = zcmax * rDt_e
      !
      IF( zcmax>0.9_wp )   CALL ctl_stop( 'dynspg_ts ERROR: Maximum Courant number is greater than 0.9 => increase nn_e' )          
      !
      ! -------------
      ! Print results
      ! -------------
      IF(lwp) WRITE(numout,*) '     Barotropic time filter => nn_bt_flt = ', nn_bt_flt
      SELECT CASE ( nn_bt_flt )
      CASE( 0 )      ;   IF(lwp) WRITE(numout,*) '        Dirac'
      CASE( 1 )      ;   IF(lwp) WRITE(numout,*) '        Boxcar: width = nn_e'
      CASE( 2 )      ;   IF(lwp) WRITE(numout,*) '        Boxcar: width = 2*nn_e' 
      CASE( 3 )      ;   IF(lwp) WRITE(numout,*) '        Demange time filter'
      CASE DEFAULT   ;   CALL ctl_stop( 'unrecognised value for nn_bt_flt: should 0,1,2 or 3' )
      END SELECT
      IF( lwp .AND. nn_bt_flt==3 )   WRITE(numout,*) '           rn_bt_alpha: ', rn_bt_alpha

      IF( lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) '     Barotropic time steps => in seconds         = ', rDt_e
         WRITE(numout,*) '                              in iterations nn_e = ', nn_e
         IF( ln_bt_auto ) THEN ; WRITE(numout,*) '        set auto (ln_bt_auto=T) with max courant number = ', zcmax
         ELSE                  ; WRITE(numout,*) '        set      (ln_bt_auto=F) with the namelist parameter nn_e '
         ENDIF
      ENDIF
      !
      !                             ! Allocate time-splitting arrays
      IF( dyn_spg_ts_alloc() /= 0    )   CALL ctl_stop('STOP', 'dyn_spg_init: failed to allocate dynspg_ts  arrays' )
      !
      ! init some arrays for debug sette
      ssha_e(:,:) = 0._wp
      !
      !                      !: restart/initialise
      ll_cold_start = .FALSE.        ! flag for restart issue
      !
      CALL ts_rst( nit000, 'READ' )
      !
   END SUBROUTINE dyn_spg_ts_init

   
   SUBROUTINE dyn_cor_2D_init( Kmm )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE dyn_cor_2D_init  ***
      !!
      !! ** Purpose : Set time splitting options
      !! Set arrays to remove/compute coriolis trend.
      !! Do it once during initialization if volume is fixed, else at each long time step.
      !! Note that these arrays are also used during barotropic loop. These are however frozen
      !! although they should be updated in the variable volume case. Not a big approximation.
      !! To remove this approximation, copy lines below inside barotropic loop
      !! and update depths at T- points (ht) at each barotropic time step
      !!
      !!----------------------------------------------------------------------
      INTEGER,  INTENT(in)         ::  Kmm  ! Time index
      INTEGER  ::   ji ,jj, jk              ! dummy loop indices
      REAL(wp) ::   z1_ht, zr1_e3f
      REAL(wp) ::   zpvo_nw, zpvo_ne, zpvo_sw, zpvo_se
      !!----------------------------------------------------------------------
      !
      ffu_nw(:,:) = 0._wp   ;   ffu_ne(:,:) = 0._wp   ;   ffu_sw(:,:) = 0._wp   ;   ffu_se(:,:) = 0._wp
      ffv_nw(:,:) = 0._wp   ;   ffv_ne(:,:) = 0._wp   ;   ffv_sw(:,:) = 0._wp   ;   ffv_se(:,:) = 0._wp
      !
      SELECT CASE( nvor_scheme )
      CASE( np_EEN )
         !
         DO_2D( 0, 0, 0, 0)
            DO jk = 1, mbku(ji,jj)
               zpvo_nw = ff_f(ji-1,jj  ) / e3f_vor(ji-1,jj  ,jk) + & 
                       & ff_f(ji  ,jj  ) / e3f_vor(ji  ,jj  ,jk) + & 
                       & ff_f(ji  ,jj-1) / e3f_vor(ji  ,jj-1,jk)
               zpvo_ne = ff_f(ji  ,jj-1) / e3f_vor(ji  ,jj-1,jk) + & 
                       & ff_f(ji  ,jj  ) / e3f_vor(ji  ,jj  ,jk) + & 
                       & ff_f(ji+1,jj  ) / e3f_vor(ji+1,jj  ,jk)
               zpvo_sw = ff_f(ji  ,jj  ) / e3f_vor(ji  ,jj  ,jk) + & 
                       & ff_f(ji  ,jj-1) / e3f_vor(ji  ,jj-1,jk) + & 
                       & ff_f(ji-1,jj-1) / e3f_vor(ji-1,jj-1,jk)
               zpvo_se = ff_f(ji+1,jj-1) / e3f_vor(ji+1,jj-1,jk) + & 
                       & ff_f(ji  ,jj-1) / e3f_vor(ji  ,jj-1,jk) + & 
                       & ff_f(ji  ,jj  ) / e3f_vor(ji  ,jj  ,jk)
               !
               ffu_nw(ji,jj) = ffu_nw(ji,jj) + e3u(ji  ,jj, jk, Kmm) * e3v(ji  ,jj  , jk, Kmm) * vmask(ji  ,jj  ,jk) * zpvo_nw
               ffu_ne(ji,jj) = ffu_ne(ji,jj) + e3u(ji  ,jj, jk, Kmm) * e3v(ji+1,jj  , jk, Kmm) * vmask(ji+1,jj  ,jk) * zpvo_ne
               ffu_sw(ji,jj) = ffu_sw(ji,jj) + e3u(ji  ,jj, jk, Kmm) * e3v(ji  ,jj-1, jk, Kmm) * vmask(ji  ,jj-1,jk) * zpvo_sw
               ffu_se(ji,jj) = ffu_se(ji,jj) + e3u(ji  ,jj, jk, Kmm) * e3v(ji+1,jj-1, jk, Kmm) * vmask(ji+1,jj-1,jk) * zpvo_se
            END DO
            ffu_nw(ji,jj) = r1_12 * r1_e1u(ji,jj) * r1_hu(ji,jj,Kmm) * e1v(ji  ,jj  ) * ffu_nw(ji,jj)
            ffu_ne(ji,jj) = r1_12 * r1_e1u(ji,jj) * r1_hu(ji,jj,Kmm) * e1v(ji+1,jj  ) * ffu_ne(ji,jj)
            ffu_sw(ji,jj) = r1_12 * r1_e1u(ji,jj) * r1_hu(ji,jj,Kmm) * e1v(ji  ,jj-1) * ffu_sw(ji,jj)
            ffu_se(ji,jj) = r1_12 * r1_e1u(ji,jj) * r1_hu(ji,jj,Kmm) * e1v(ji+1,jj-1) * ffu_se(ji,jj)
            !
         END_2D
         !
         DO_2D( 0, 0, 0, 0)
            DO jk = 1, mbkv(ji,jj)
               zpvo_se = ff_f(ji-1,jj  ) / e3f_vor(ji-1,jj  ,jk) + & 
                       & ff_f(ji  ,jj  ) / e3f_vor(ji  ,jj  ,jk) + & 
                       & ff_f(ji  ,jj-1) / e3f_vor(ji  ,jj-1,jk)
               zpvo_sw = ff_f(ji-1,jj-1) / e3f_vor(ji-1,jj-1,jk) + & 
                       & ff_f(ji-1,jj  ) / e3f_vor(ji-1,jj  ,jk) + & 
                       & ff_f(ji  ,jj  ) / e3f_vor(ji  ,jj  ,jk)
               zpvo_ne = ff_f(ji  ,jj+1) / e3f_vor(ji  ,jj+1,jk) + & 
                       & ff_f(ji  ,jj  ) / e3f_vor(ji  ,jj  ,jk) + & 
                       & ff_f(ji-1,jj  ) / e3f_vor(ji-1,jj  ,jk)
               zpvo_nw = ff_f(ji  ,jj  ) / e3f_vor(ji  ,jj  ,jk) + & 
                       & ff_f(ji-1,jj  ) / e3f_vor(ji-1,jj  ,jk) + & 
                       & ff_f(ji-1,jj+1) / e3f_vor(ji-1,jj+1,jk)
               !
               ffv_nw(ji,jj) = ffv_nw(ji,jj) + e3v(ji  ,jj, jk, Kmm) * e3u(ji-1,jj+1, jk, Kmm) * umask(ji-1,jj+1, jk) * zpvo_nw
               ffv_ne(ji,jj) = ffv_ne(ji,jj) + e3v(ji  ,jj, jk, Kmm) * e3u(ji  ,jj+1, jk, Kmm) * umask(ji  ,jj+1, jk) * zpvo_ne
               ffv_sw(ji,jj) = ffv_sw(ji,jj) + e3v(ji  ,jj, jk, Kmm) * e3u(ji-1,jj  , jk, Kmm) * umask(ji-1,jj  , jk) * zpvo_sw
               ffv_se(ji,jj) = ffv_se(ji,jj) + e3v(ji  ,jj, jk, Kmm) * e3u(ji  ,jj  , jk, Kmm) * umask(ji  ,jj  , jk) * zpvo_se
            END DO
            ffv_nw(ji,jj) = r1_12 * r1_e2v(ji,jj) * r1_hv(ji,jj,Kmm) * e2u(ji-1,jj+1) * ffv_nw(ji,jj)
            ffv_ne(ji,jj) = r1_12 * r1_e2v(ji,jj) * r1_hv(ji,jj,Kmm) * e2u(ji  ,jj+1) * ffv_ne(ji,jj)
            ffv_sw(ji,jj) = r1_12 * r1_e2v(ji,jj) * r1_hv(ji,jj,Kmm) * e2u(ji-1,jj  ) * ffv_sw(ji,jj)
            ffv_se(ji,jj) = r1_12 * r1_e2v(ji,jj) * r1_hv(ji,jj,Kmm) * e2u(ji  ,jj  ) * ffv_se(ji,jj)
            !
         END_2D
         !
      CASE( np_ENE )
         !
         DO_2D( 0, 0, 0, 0)
            DO jk = 1, mbku(ji,jj)
               ffu_nw(ji,jj) = ffu_nw(ji,jj) + e3u(ji  ,jj, jk, Kmm) * e3v(ji  ,jj  , jk, Kmm) * vmask(ji  ,jj  ,jk) / e3f_vor(ji,jj  ,jk)
               ffu_ne(ji,jj) = ffu_ne(ji,jj) + e3u(ji  ,jj, jk, Kmm) * e3v(ji+1,jj  , jk, Kmm) * vmask(ji+1,jj  ,jk) / e3f_vor(ji,jj  ,jk)
               ffu_sw(ji,jj) = ffu_sw(ji,jj) + e3u(ji  ,jj, jk, Kmm) * e3v(ji  ,jj-1, jk, Kmm) * vmask(ji  ,jj-1,jk) / e3f_vor(ji,jj-1,jk)
               ffu_se(ji,jj) = ffu_se(ji,jj) + e3u(ji  ,jj, jk, Kmm) * e3v(ji+1,jj-1, jk, Kmm) * vmask(ji+1,jj-1,jk) / e3f_vor(ji,jj-1,jk)
            END DO
            ffu_nw(ji,jj) = r1_4 * r1_e1u(ji,jj) * r1_hu(ji,jj,Kmm) * e1v(ji  ,jj  ) * ff_f(ji  ,jj  ) * ffu_nw(ji,jj)
            ffu_ne(ji,jj) = r1_4 * r1_e1u(ji,jj) * r1_hu(ji,jj,Kmm) * e1v(ji+1,jj  ) * ff_f(ji  ,jj  ) * ffu_ne(ji,jj)
            ffu_sw(ji,jj) = r1_4 * r1_e1u(ji,jj) * r1_hu(ji,jj,Kmm) * e1v(ji  ,jj-1) * ff_f(ji  ,jj-1) * ffu_sw(ji,jj)
            ffu_se(ji,jj) = r1_4 * r1_e1u(ji,jj) * r1_hu(ji,jj,Kmm) * e1v(ji+1,jj-1) * ff_f(ji  ,jj-1) * ffu_se(ji,jj)
            !
         END_2D
         !
         DO_2D( 0, 0, 0, 0)
            DO jk = 1, mbkv(ji,jj)
               ffv_nw(ji,jj) = ffv_nw(ji,jj) + e3v(ji  ,jj, jk, Kmm) * e3u(ji-1,jj+1, jk, Kmm) * umask(ji-1,jj+1, jk) / e3f_vor(ji-1,jj,jk)
               ffv_ne(ji,jj) = ffv_ne(ji,jj) + e3v(ji  ,jj, jk, Kmm) * e3u(ji  ,jj+1, jk, Kmm) * umask(ji  ,jj+1, jk) / e3f_vor(ji  ,jj,jk)
               ffv_sw(ji,jj) = ffv_sw(ji,jj) + e3v(ji  ,jj, jk, Kmm) * e3u(ji-1,jj  , jk, Kmm) * umask(ji-1,jj  , jk) / e3f_vor(ji-1,jj,jk)
               ffv_se(ji,jj) = ffv_se(ji,jj) + e3v(ji  ,jj, jk, Kmm) * e3u(ji  ,jj  , jk, Kmm) * umask(ji  ,jj  , jk) / e3f_vor(ji  ,jj,jk)
            END DO
            ffv_nw(ji,jj) = r1_4 * r1_e2v(ji,jj) * r1_hv(ji,jj,Kmm) * e2u(ji-1,jj+1) * ff_f(ji-1,jj  ) * ffv_nw(ji,jj)
            ffv_ne(ji,jj) = r1_4 * r1_e2v(ji,jj) * r1_hv(ji,jj,Kmm) * e2u(ji  ,jj+1) * ff_f(ji  ,jj  ) * ffv_ne(ji,jj)
            ffv_sw(ji,jj) = r1_4 * r1_e2v(ji,jj) * r1_hv(ji,jj,Kmm) * e2u(ji-1,jj  ) * ff_f(ji-1,jj  ) * ffv_sw(ji,jj)
            ffv_se(ji,jj) = r1_4 * r1_e2v(ji,jj) * r1_hv(ji,jj,Kmm) * e2u(ji  ,jj  ) * ff_f(ji  ,jj  ) * ffv_se(ji,jj)
            !
         END_2D
         !
      CASE( np_ENS )
         !
         DO_2D( 0, 0, 0, 0)
            DO jk = 1, mbku(ji,jj)
               zr1_e3f = ff_f(ji,jj) / e3f_vor(ji,jj  ,jk) + ff_f(ji,jj-1) / e3f_vor(ji,jj-1,jk)
               ffu_nw(ji,jj) = ffu_nw(ji,jj) + e3u(ji  ,jj, jk, Kmm) * e3v(ji  ,jj  , jk, Kmm) * vmask(ji  ,jj  ,jk) * zr1_e3f
               ffu_ne(ji,jj) = ffu_ne(ji,jj) + e3u(ji  ,jj, jk, Kmm) * e3v(ji+1,jj  , jk, Kmm) * vmask(ji+1,jj  ,jk) * zr1_e3f
               ffu_sw(ji,jj) = ffu_sw(ji,jj) + e3u(ji  ,jj, jk, Kmm) * e3v(ji  ,jj-1, jk, Kmm) * vmask(ji  ,jj-1,jk) * zr1_e3f
               ffu_se(ji,jj) = ffu_se(ji,jj) + e3u(ji  ,jj, jk, Kmm) * e3v(ji+1,jj-1, jk, Kmm) * vmask(ji+1,jj-1,jk) * zr1_e3f
            END DO
            ffu_nw(ji,jj) = r1_8 * r1_e1u(ji,jj) * r1_hu(ji,jj,Kmm) * e1v(ji  ,jj  ) * ffu_nw(ji,jj)
            ffu_ne(ji,jj) = r1_8 * r1_e1u(ji,jj) * r1_hu(ji,jj,Kmm) * e1v(ji+1,jj  ) * ffu_ne(ji,jj)
            ffu_sw(ji,jj) = r1_8 * r1_e1u(ji,jj) * r1_hu(ji,jj,Kmm) * e1v(ji  ,jj-1) * ffu_sw(ji,jj)
            ffu_se(ji,jj) = r1_8 * r1_e1u(ji,jj) * r1_hu(ji,jj,Kmm) * e1v(ji+1,jj-1) * ffu_se(ji,jj)
            !
         END_2D
         !
         DO_2D( 0, 0, 0, 0)
            DO jk = 1, mbkv(ji,jj)
               zr1_e3f = ff_f(ji,jj) / e3f_vor(ji,jj  ,jk) + ff_f(ji-1,jj) / e3f_vor(ji-1,jj,jk)
               ffv_nw(ji,jj) = ffv_nw(ji,jj) + e3v(ji  ,jj, jk, Kmm) * e3u(ji-1,jj+1, jk, Kmm) * umask(ji-1,jj+1, jk) * zr1_e3f
               ffv_ne(ji,jj) = ffv_ne(ji,jj) + e3v(ji  ,jj, jk, Kmm) * e3u(ji  ,jj+1, jk, Kmm) * umask(ji  ,jj+1, jk) * zr1_e3f
               ffv_sw(ji,jj) = ffv_sw(ji,jj) + e3v(ji  ,jj, jk, Kmm) * e3u(ji-1,jj  , jk, Kmm) * umask(ji-1,jj  , jk) * zr1_e3f
               ffv_se(ji,jj) = ffv_se(ji,jj) + e3v(ji  ,jj, jk, Kmm) * e3u(ji  ,jj  , jk, Kmm) * umask(ji  ,jj  , jk) * zr1_e3f
            END DO
            ffv_nw(ji,jj) = r1_8 * r1_e2v(ji,jj) * r1_hv(ji,jj,Kmm) * e2u(ji-1,jj+1) * ffv_nw(ji,jj)
            ffv_ne(ji,jj) = r1_8 * r1_e2v(ji,jj) * r1_hv(ji,jj,Kmm) * e2u(ji  ,jj+1) * ffv_ne(ji,jj)
            ffv_sw(ji,jj) = r1_8 * r1_e2v(ji,jj) * r1_hv(ji,jj,Kmm) * e2u(ji-1,jj  ) * ffv_sw(ji,jj)
            ffv_se(ji,jj) = r1_8 * r1_e2v(ji,jj) * r1_hv(ji,jj,Kmm) * e2u(ji  ,jj  ) * ffv_se(ji,jj)
            !
         END_2D
         !

      CASE( np_ENT )
         !
         DO_2D( 0, 0, 0, 0)
            !
            DO jk = 1, mbku(ji,jj)
               ffu_nw(ji,jj) = ffu_nw(ji,jj) + e3t(ji  ,jj, jk, Kmm) * vmask(ji  , jj  , jk)
               ffu_ne(ji,jj) = ffu_ne(ji,jj) + e3t(ji+1,jj, jk, Kmm) * vmask(ji+1, jj  , jk)
               ffu_sw(ji,jj) = ffu_sw(ji,jj) + e3t(ji  ,jj, jk, Kmm) * vmask(ji  , jj-1, jk)
               ffu_se(ji,jj) = ffu_se(ji,jj) + e3t(ji+1,jj, jk, Kmm) * vmask(ji+1, jj-1, jk)
            END DO
            ffu_nw(ji,jj) = r1_4 * r1_e1e2u(ji,jj) * r1_hu(ji,jj,Kmm) * e1e2t(ji  ,jj) * ff_t(ji  ,jj) * ffu_nw(ji,jj)
            ffu_ne(ji,jj) = r1_4 * r1_e1e2u(ji,jj) * r1_hu(ji,jj,Kmm) * e1e2t(ji+1,jj) * ff_t(ji+1,jj) * ffu_ne(ji,jj)
            ffu_sw(ji,jj) = r1_4 * r1_e1e2u(ji,jj) * r1_hu(ji,jj,Kmm) * e1e2t(ji  ,jj) * ff_t(ji  ,jj) * ffu_sw(ji,jj)
            ffu_se(ji,jj) = r1_4 * r1_e1e2u(ji,jj) * r1_hu(ji,jj,Kmm) * e1e2t(ji+1,jj) * ff_t(ji+1,jj) * ffu_se(ji,jj)
         END_2D
         !
         DO_2D( 0, 0, 0, 0)
            !
            DO jk = 1, mbkv(ji,jj)
               ffv_nw(ji,jj) = ffv_nw(ji,jj) + e3t(ji,jj+1, jk, Kmm) * umask(ji-1, jj+1, jk)
               ffv_sw(ji,jj) = ffv_sw(ji,jj) + e3t(ji,jj  , jk, Kmm) * umask(ji-1, jj  , jk)
               ffv_ne(ji,jj) = ffv_ne(ji,jj) + e3t(ji,jj+1, jk, Kmm) * umask(ji  , jj+1, jk)
               ffv_se(ji,jj) = ffv_se(ji,jj) + e3t(ji,jj  , jk, Kmm) * umask(ji  , jj  , jk)
            END DO
            ffv_nw(ji,jj) = r1_4 * r1_e1e2v(ji,jj) * r1_hv(ji,jj,Kmm) * e1e2t(ji,jj+1) * ff_t(ji,jj+1) * ffv_nw(ji,jj)
            ffv_sw(ji,jj) = r1_4 * r1_e1e2v(ji,jj) * r1_hv(ji,jj,Kmm) * e1e2t(ji,jj  ) * ff_t(ji,jj  ) * ffv_sw(ji,jj)
            ffv_ne(ji,jj) = r1_4 * r1_e1e2v(ji,jj) * r1_hv(ji,jj,Kmm) * e1e2t(ji,jj+1) * ff_t(ji,jj+1) * ffv_ne(ji,jj)
            ffv_se(ji,jj) = r1_4 * r1_e1e2v(ji,jj) * r1_hv(ji,jj,Kmm) * e1e2t(ji,jj  ) * ff_t(ji,jj  ) * ffv_se(ji,jj)
            !
            !
         END_2D
         !
      END SELECT
      !
   END SUBROUTINE dyn_cor_2D_init


   SUBROUTINE dyn_cor_2D( punb, pvnb, zu_trd, zv_trd   )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE dyn_cor_2D  ***
      !!
      !! ** Purpose : Compute u and v coriolis trends
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj                             ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: punb, pvnb
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) :: zu_trd, zv_trd
      !!----------------------------------------------------------------------
      !
      DO_2D( 0, 0, 0, 0 )
         zu_trd(ji,jj) = + (  (  ffu_nw(ji,jj) * pvnb(ji  ,jj  )   & 
            &                  + ffu_se(ji,jj) * pvnb(ji+1,jj-1) ) &
            &               + (  ffu_sw(ji,jj) * pvnb(ji  ,jj-1)   & 
            &                  + ffu_ne(ji,jj) * pvnb(ji+1,jj  ) ) ) 
            !
         zv_trd(ji,jj) = - (  (  ffv_sw(ji,jj) * punb(ji-1,jj  )   & 
            &                  + ffv_ne(ji,jj) * punb(ji  ,jj+1) ) &
            &               + (  ffv_nw(ji,jj) * punb(ji-1,jj+1)   &
            &                  + ffv_se(ji,jj) * punb(ji  ,jj  ) ) )
      END_2D
      !
   END SUBROUTINE dyn_cor_2D


   SUBROUTINE wad_tmsk( pssh, ptmsk )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE wad_tmsk  ***
      !!                    
      !! ** Purpose :   set wetting & drying mask at tracer points 
      !!              for the current barotropic sub-step 
      !!
      !! ** Method  :   ??? 
      !!
      !! ** Action  :  ptmsk : wetting & drying t-mask
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pssh    !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   ptmsk   !
      !
      INTEGER  ::   ji, jj   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF( ln_wd_dl_rmp ) THEN     
         DO_2D( 0, 1, 0, 1 )
            IF    ( pssh(ji,jj) + ht_0(ji,jj) >  2._wp * rn_wdmin1 ) THEN 
               !           IF    ( pssh(ji,jj) + ht_0(ji,jj) >          rn_wdmin2 ) THEN 
               ptmsk(ji,jj) = 1._wp
            ELSEIF( pssh(ji,jj) + ht_0(ji,jj) >          rn_wdmin1 ) THEN
               ptmsk(ji,jj) = TANH( 50._wp*( ( pssh(ji,jj) + ht_0(ji,jj) -  rn_wdmin1 )*r_rn_wdmin1) )
            ELSE 
               ptmsk(ji,jj) = 0._wp
            ENDIF
         END_2D
      ELSE  
         DO_2D( 0, 1, 0, 1 )
            IF ( pssh(ji,jj) + ht_0(ji,jj) >  rn_wdmin1 ) THEN   ;   ptmsk(ji,jj) = 1._wp
            ELSE                                                 ;   ptmsk(ji,jj) = 0._wp
            ENDIF
         END_2D
      ENDIF
      !
   END SUBROUTINE wad_tmsk


   SUBROUTINE wad_Umsk( pTmsk, phU, phV, pu, pv, pUmsk, pVmsk )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE wad_Umsk  ***
      !!                    
      !! ** Purpose :   set wetting & drying mask at tracer points 
      !!              for the current barotropic sub-step 
      !!
      !! ** Method  :   ??? 
      !!
      !! ** Action  :  ptmsk : wetting & drying t-mask
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pTmsk              ! W & D t-mask
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   phU, phV, pu, pv   ! ocean velocities and transports
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   pUmsk, pVmsk       ! W & D u- and v-mask
      !
      INTEGER  ::   ji, jj   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      DO_2D( 0, 0, 0, 0 )   ! not jpi-column
         IF ( phU(ji,jj) > 0._wp ) THEN   ;   pUmsk(ji,jj) = pTmsk(ji  ,jj) 
         ELSE                             ;   pUmsk(ji,jj) = pTmsk(ji+1,jj)  
         ENDIF
         phU(ji,jj) = pUmsk(ji,jj)*phU(ji,jj)
         pu (ji,jj) = pUmsk(ji,jj)*pu (ji,jj)
      END_2D
      !
      DO_2D( 0, 0, 0, 0 )   ! not jpj-row
         IF ( phV(ji,jj) > 0._wp ) THEN   ;   pVmsk(ji,jj) = pTmsk(ji,jj  )
         ELSE                             ;   pVmsk(ji,jj) = pTmsk(ji,jj+1)  
         ENDIF
         phV(ji,jj) = pVmsk(ji,jj)*phV(ji,jj) 
         pv (ji,jj) = pVmsk(ji,jj)*pv (ji,jj)
      END_2D
      !
   END SUBROUTINE wad_Umsk

   SUBROUTINE dyn_drg_init( Kbb, puu, pvv, puu_b ,pvv_b, pu_RHSi, pv_RHSi, pCdU_u, pCdU_v )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_drg_init  ***
      !!                    
      !! ** Purpose : - add the baroclinic top/bottom drag contribution to 
      !!              the baroclinic part of the barotropic RHS
      !!              - compute the barotropic drag coefficients
      !!
      !! ** Method  :   computation done over the INNER domain only 
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT(in   ) ::  Kbb                ! ocean time level indice
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(in   ) ::  puu, pvv           ! ocean velocities and RHS of momentum equation
      REAL(wp), DIMENSION(jpi,jpj,jpt)    , INTENT(in   ) ::  puu_b, pvv_b       ! barotropic velocities at main time levels
      REAL(wp), DIMENSION(A2D(0))         , INTENT(inout) ::  pu_RHSi, pv_RHSi   ! baroclinic part of the barotropic RHS
      REAL(wp), DIMENSION(jpi,jpj)        , INTENT(  out) ::  pCdU_u , pCdU_v    ! barotropic drag coefficients
      !
      INTEGER  ::   ji, jj   ! dummy loop indices
      INTEGER  ::   ikbu, ikbv, iktu, iktv
      REAL(wp) ::   zztmp
      REAL(wp), DIMENSION(jpi,jpj) ::   zu_i, zv_i
      !!----------------------------------------------------------------------
      !
      !                    !==  Set the barotropic drag coef.  ==!
      !
      IF( ln_isfcav.OR.ln_drgice_imp ) THEN          ! top+bottom friction (ocean cavities)
         
         DO_2D( 0, 0, 0, 0 )
            pCdU_u(ji,jj) = r1_2*( ( rCdU_bot(ji+1,jj) + rCdU_bot(ji,jj) ) + ( rCdU_top(ji+1,jj) + rCdU_top(ji,jj) ) )
            pCdU_v(ji,jj) = r1_2*( ( rCdU_bot(ji,jj+1) + rCdU_bot(ji,jj) ) + ( rCdU_top(ji,jj+1) + rCdU_top(ji,jj) ) )
         END_2D
      ELSE                          ! bottom friction only
         DO_2D( 0, 0, 0, 0 )
            pCdU_u(ji,jj) = r1_2*( rCdU_bot(ji+1,jj) + rCdU_bot(ji,jj) )
            pCdU_v(ji,jj) = r1_2*( rCdU_bot(ji,jj+1) + rCdU_bot(ji,jj) )
         END_2D
      ENDIF
      !
      !                    !==  BOTTOM stress contribution from baroclinic velocities  ==!
      !
      DO_2D( 0, 0, 0, 0 )
         ikbu = mbku(ji,jj)       
         ikbv = mbkv(ji,jj)    
         zu_i(ji,jj) = puu(ji,jj,ikbu,Kbb) - puu_b(ji,jj,Kbb)
         zv_i(ji,jj) = pvv(ji,jj,ikbv,Kbb) - pvv_b(ji,jj,Kbb)
      END_2D
      !
      ! use "unclipped" drag (even if explicit friction is used in 3D calculation)
         
      DO_2D( 0, 0, 0, 0 )
         pu_RHSi(ji,jj) = pu_RHSi(ji,jj) + r1_hu(ji,jj,Kbb) * r1_2*( rCdU_bot(ji+1,jj)+rCdU_bot(ji,jj) ) * zu_i(ji,jj)
         pv_RHSi(ji,jj) = pv_RHSi(ji,jj) + r1_hv(ji,jj,Kbb) * r1_2*( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj) ) * zv_i(ji,jj)
      END_2D
      !
      !                    !==  TOP stress contribution from baroclinic velocities  ==!   (no W/D case)
      !
      IF( ln_isfcav.OR.ln_drgice_imp ) THEN
         !
         DO_2D( 0, 0, 0, 0 )
            iktu = miku(ji,jj)
            iktv = mikv(ji,jj)
            zu_i(ji,jj) = puu(ji,jj,iktu,Kbb) - puu_b(ji,jj,Kbb)
            zv_i(ji,jj) = pvv(ji,jj,iktv,Kbb) - pvv_b(ji,jj,Kbb)
         END_2D
         !                    ! use "unclipped" top drag (even if explicit friction is used in 3D calculation)       
         DO_2D( 0, 0, 0, 0 )
            pu_RHSi(ji,jj) = pu_RHSi(ji,jj) + r1_hu(ji,jj,Kbb) * r1_2*( rCdU_top(ji+1,jj)+rCdU_top(ji,jj) ) * zu_i(ji,jj)
            pv_RHSi(ji,jj) = pv_RHSi(ji,jj) + r1_hv(ji,jj,Kbb) * r1_2*( rCdU_top(ji,jj+1)+rCdU_top(ji,jj) ) * zv_i(ji,jj)
         END_2D
         !
      ENDIF
      !
   END SUBROUTINE dyn_drg_init


   SUBROUTINE ts_bck_interp( jn, ll_init,       &   ! <== in
      &                      za0, za1, za2, za3 )   ! ==> out
      !!----------------------------------------------------------------------
      INTEGER ,INTENT(in   ) ::   jn                   ! index of sub time step
      LOGICAL ,INTENT(in   ) ::   ll_init              !
      REAL(wp),INTENT(  out) ::   za0, za1, za2, za3   ! Half-step back interpolation coefficient
      !
      REAL(wp) ::   zepsilon, zgamma                   !   -      -
      !!----------------------------------------------------------------------
      !                             ! set Half-step back interpolation coefficient
      IF    ( jn==1 .AND. ll_init ) THEN   !* Forward-backward
         za0 = 1._wp                        
         za1 = 0._wp                           
         za2 = 0._wp
         za3 = 0._wp
      ELSEIF( jn==2 .AND. ll_init ) THEN   !* AB2-AM3 Coefficients; bet=0 ; gam=-1/6 ; eps=1/12
         za0 = 1.0833333333333_wp                 ! za0 = 1-gam-eps
         za1 =-0.1666666666666_wp                 ! za1 = gam
         za2 = 0.0833333333333_wp                 ! za2 = eps
         za3 = 0._wp              
      ELSE                                 !* AB3-AM4 Coefficients; bet=0.281105 ; eps=0.013 ; gam=0.0880 
         IF( rn_bt_alpha == 0._wp ) THEN      ! Time diffusion  
            za0 = 0.614_wp                        ! za0 = 1/2 +   gam + 2*eps
            za1 = 0.285_wp                        ! za1 = 1/2 - 2*gam - 3*eps
            za2 = 0.088_wp                        ! za2 = gam
            za3 = 0.013_wp                        ! za3 = eps
         ELSE                                 ! no time diffusion
            zepsilon = 0.00976186_wp - 0.13451357_wp * rn_bt_alpha
            zgamma   = 0.08344500_wp - 0.51358400_wp * rn_bt_alpha
            za0 = 0.5_wp + zgamma + 2._wp * rn_bt_alpha + 2._wp * zepsilon
            za1 = 1._wp - za0 - zgamma - zepsilon
            za2 = zgamma
            za3 = zepsilon
         ENDIF 
      ENDIF
   END SUBROUTINE ts_bck_interp

   !!======================================================================
END MODULE dynspg_ts

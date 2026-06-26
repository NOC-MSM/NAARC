MODULE dynspg
   !!======================================================================
   !!                       ***  MODULE  dynspg  ***
   !! Ocean dynamics:  surface pressure gradient control
   !!======================================================================
   !! History :  1.0  ! 2005-12  (C. Talandier, G. Madec, V. Garnier)  Original code
   !!            3.2  ! 2009-07  (R. Benshila)  Suppression of rigid-lid option
   !!            4.2  ! 2020-12  (G. Madec, E. Clementi) add Bernoulli Head for
   !!                            wave coupling
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_spg     : update the dynamics trend with surface pressure gradient
   !!   dyn_spg_init: initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain variables
   USE c1d            ! 1D vertical configuration
   USE phycst         ! physical constants
   USE sbc_oce        ! surface boundary condition: ocean
   USE sbc_ice , ONLY : snwice_mass, snwice_mass_b
   USE sbcapr         ! surface boundary condition: atmospheric pressure
   USE sbcwave,  ONLY : bhd_wave
   USE dynspg_exp     ! surface pressure gradient     (dyn_spg_exp routine)
   USE dynspg_ts      ! surface pressure gradient     (dyn_spg_ts  routine)
   USE tide_mod       !
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   !
   USE prtctl         ! Print control                     (prt_ctl routine)
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_spg        ! routine called by step module
   PUBLIC   dyn_spg_init   ! routine called by opa module

   INTEGER ::   nspg = 0   ! type of surface pressure gradient scheme defined from lk_dynspg_...

   !                       ! Parameter to control the surface pressure gradient scheme
   INTEGER, PARAMETER ::   np_TS  = 1   ! split-explicit time stepping (Time-Splitting)
   INTEGER, PARAMETER ::   np_EXP = 0   !       explicit time stepping
   INTEGER, PARAMETER ::   np_NO  =-1   ! no surface pressure gradient, no scheme
   !
   REAL(wp) ::   zt0step !   Time of day at the beginning of the time step

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_spg( kt, Kbb, Kmm, Krhs, puu, pvv, pssh, puu_b, pvv_b, Kaa )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_spg  ***
      !!
      !! ** Purpose :   compute surface pressure gradient including the
      !!              atmospheric pressure forcing (ln_apr_dyn=T).
      !!
      !! ** Method  :   Two schemes:
      !!              - explicit       : the spg is evaluated at now
      !!              - split-explicit : a time splitting technique is used
      !!
      !!              ln_apr_dyn=T : the atmospheric pressure forcing is applied
      !!             as the gradient of the inverse barometer ssh:
      !!                apgu = - 1/rho0 di[apr] = 0.5*grav di[ssh_ib+ssh_ibb]
      !!                apgv = - 1/rho0 dj[apr] = 0.5*grav dj[ssh_ib+ssh_ibb]
      !!             Note that as all external forcing a time averaging over a two rn_Dt
      !!             period is used to prevent the divergence of odd and even time step.
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT( in )  ::  kt                  ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  Kbb, Kmm, Krhs, Kaa ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv            ! ocean velocities and RHS of momentum equation
      REAL(wp), DIMENSION(jpi,jpj,jpt)    , INTENT(inout) ::  pssh, puu_b, pvv_b  ! SSH and barotropic velocities at main time levels
      !
      INTEGER  ::   ji, jj, jk                   ! dummy loop indices
      REAL(wp) ::   zg_2, zintp, zgrho0r, zld    ! local scalars
      REAL(wp)             , DIMENSION(jpi,jpj) ::   zpgu, zpgv   ! 2D workspace
#if ! defined key_PSYCLONE_2p5p0
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)     ::   zpice
#else
      REAL(wp), DIMENSION(jpi,jpj)              ::   zpice
#endif
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   ::   ztrdu, ztrdv
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dyn_spg')
      !
      IF( l_trddyn .AND. nspg == np_EXP )   THEN                      ! temporary save of ta and sa trends
         ALLOCATE( ztrdu(T2D(0),jpk), ztrdv(T2D(0),jpk) )
         ztrdu(:,:,:) = puu(T2D(0),:,Krhs)
         ztrdv(:,:,:) = pvv(T2D(0),:,Krhs)
      ENDIF
      !
      IF(      ln_apr_dyn                                                &   ! atmos. pressure
         .OR.  ( .NOT.ln_dynspg_ts .AND. (ln_tide_pot .AND. ln_tide) )   &   ! tide potential (no time slitting)
         .OR.  ln_ice_embd                                               &   ! embedded sea-ice
         .OR.  ( ln_wave .and. ln_bern_srfc ) ) THEN                         ! depth-independent Bernoulli head
         !
         DO_2D( 0, 0, 0, 0 )
            zpgu(ji,jj) = 0._wp
            zpgv(ji,jj) = 0._wp
         END_2D
         !
         IF( ln_apr_dyn .AND. .NOT.ln_dynspg_ts ) THEN   !==  Atmospheric pressure gradient (added later in time-split case) ==!
            zg_2 = grav * 0.5
            DO_2D( 0, 0, 0, 0 )                       ! gradient of Patm using inverse barometer ssh
               zpgu(ji,jj) = zpgu(ji,jj) + zg_2 * (  ( ssh_ib (ji+1,jj) - ssh_ib (ji,jj) )    &   ! add () for NP repro
                  &                                + ( ssh_ibb(ji+1,jj) - ssh_ibb(ji,jj) )  ) * r1_e1u(ji,jj)
               zpgv(ji,jj) = zpgv(ji,jj) + zg_2 * (  ( ssh_ib (ji,jj+1) - ssh_ib (ji,jj) )    &   ! add () for NP repro
                  &                                + ( ssh_ibb(ji,jj+1) - ssh_ibb(ji,jj) )  ) * r1_e2v(ji,jj)
            END_2D
         ENDIF
         !
         !                                    !==  tide potential forcing term  ==!
         IF( .NOT.ln_dynspg_ts .AND. ( ln_tide_pot .AND. ln_tide )  ) THEN   ! N.B. added directly at sub-time-step in ts-case
            !
            ! Update tide potential at the beginning of current time step
            zt0step = REAL(nsec_day, wp)-0.5_wp*rn_Dt
            CALL upd_tide(zt0step, Kmm)
            !
            DO_2D( 0, 0, 0, 0 )                      ! add tide potential forcing
               zpgu(ji,jj) = zpgu(ji,jj) + grav * ( pot_astro(ji+1,jj) - pot_astro(ji,jj) ) * r1_e1u(ji,jj)
               zpgv(ji,jj) = zpgv(ji,jj) + grav * ( pot_astro(ji,jj+1) - pot_astro(ji,jj) ) * r1_e2v(ji,jj)
            END_2D
            !
            IF (ln_scal_load) THEN
               zld = rn_scal_load * grav
               DO_2D( 0, 0, 0, 0 )                   ! add scalar approximation for load potential
                  zpgu(ji,jj) = zpgu(ji,jj) + zld * ( pssh(ji+1,jj,Kmm) - pssh(ji,jj,Kmm) ) * r1_e1u(ji,jj)
                  zpgv(ji,jj) = zpgv(ji,jj) + zld * ( pssh(ji,jj+1,Kmm) - pssh(ji,jj,Kmm) ) * r1_e2v(ji,jj)
               END_2D
            ENDIF
         ENDIF
         !
         IF( ln_ice_embd ) THEN              !== embedded sea ice: Pressure gradient due to snow-ice mass ==!
#if ! defined key_PSYCLONE_2p5p0
            ALLOCATE( zpice(jpi,jpj) )
#endif
            zintp = REAL( MOD( kt-1, nn_fsbc ) ) / REAL( nn_fsbc )
            zgrho0r     = - grav * r1_rho0
            zpice(:,:) = (  zintp * snwice_mass(:,:) + ( 1.- zintp ) * snwice_mass_b(:,:)  ) * zgrho0r
            DO_2D( 0, 0, 0, 0 )
               zpgu(ji,jj) = zpgu(ji,jj) + ( zpice(ji+1,jj) - zpice(ji,jj) ) * r1_e1u(ji,jj)
               zpgv(ji,jj) = zpgv(ji,jj) + ( zpice(ji,jj+1) - zpice(ji,jj) ) * r1_e2v(ji,jj)
            END_2D
#if ! defined key_PSYCLONE_2p5p0
            DEALLOCATE( zpice )
#endif
         ENDIF
         !
         IF( ln_wave .and. ln_bern_srfc ) THEN          !== Add J terms: depth-independent Bernoulli head
            DO_2D( 0, 0, 0, 0 )
               zpgu(ji,jj) = zpgu(ji,jj) + ( bhd_wave(ji+1,jj) - bhd_wave(ji,jj) ) * r1_e1u(ji,jj)   !++ bhd_wave from wave model in m2/s2 [BHD parameters in WW3]
               zpgv(ji,jj) = zpgv(ji,jj) + ( bhd_wave(ji,jj+1) - bhd_wave(ji,jj) ) * r1_e2v(ji,jj)
            END_2D
         ENDIF
         !
         DO_3D( 0, 0, 0, 0, 1, jpkm1 )       !== Add all terms to the general trend
            puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) + zpgu(ji,jj)
            pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) + zpgv(ji,jj)
         END_3D
         !
!!gm add here a call to dyn_trd for ice pressure gradient, the surf pressure trends ????
         !
      ENDIF
      !
      SELECT CASE ( nspg )                   !== surface pressure gradient computed and add to the general trend ==!
      CASE ( np_EXP )   ;   CALL dyn_spg_exp( kt,      Kmm,       puu, pvv, Krhs )                    ! explicit
      CASE ( np_TS  )   ;   CALL dyn_spg_ts ( kt, Kbb, Kmm, Krhs, puu, pvv, pssh, puu_b, pvv_b, Kaa ) ! time-splitting
      END SELECT
      !
      IF( l_trddyn .AND. nspg == np_EXP )   THEN                  ! save the surface pressure gradient trends for further diagnostics
         ztrdu(:,:,:) = puu(T2D(0),:,Krhs) - ztrdu(:,:,:)
         ztrdv(:,:,:) = pvv(T2D(0),:,Krhs) - ztrdv(:,:,:)
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_spg, kt, Kmm )
         ! In this case we also need to finalise the write-out of the PVO 3D diagnostic with zero correction
         ztrdu(:,:,1) = 0._wp
         ztrdv(:,:,1) = 0._wp
         CALL trd_dyn( ztrdu(:,:,1), ztrdv(:,:,1), jpdyn_pvo, kt, Kmm )
         DEALLOCATE( ztrdu , ztrdv )
      ENDIF
      !                                      ! print mean trends (used for debugging)
      IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab3d_1=puu(:,:,:,Krhs), clinfo1=' spg  - Ua: ', mask1=umask, &
         &                                  tab3d_2=pvv(:,:,:,Krhs), clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( ln_timing )   CALL timing_stop('dyn_spg')
      !
   END SUBROUTINE dyn_spg


   SUBROUTINE dyn_spg_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_spg_init  ***
      !!
      !! ** Purpose :   Control the consistency between namelist options for
      !!              surface pressure gradient schemes
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio, ios   ! local integers
      !
      NAMELIST/namdyn_spg/ ln_dynspg_exp    , ln_dynspg_ts,   &
      &                    ln_bt_fw         , ln_bt_auto  ,   &
      &                    nn_e , rn_bt_cmax, nn_bt_flt, rn_bt_alpha
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_spg_init : choice of the surface pressure gradient scheme'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      READ_NML_REF(numnam,namdyn_spg)
      READ_NML_CFG(numnam,namdyn_spg)
      IF(lwm) WRITE ( numond, namdyn_spg )
      !
      IF(lwp) THEN             ! Namelist print
         WRITE(numout,*) '   Namelist : namdyn_spg                    '
         WRITE(numout,*) '      Explicit free surface                  ln_dynspg_exp = ', ln_dynspg_exp
         WRITE(numout,*) '      Free surface with time splitting       ln_dynspg_ts  = ', ln_dynspg_ts
      ENDIF
      !                          ! Control of surface pressure gradient scheme options
                                     nspg =  np_NO    ;   ioptio = 0
      IF( ln_dynspg_exp ) THEN   ;   nspg =  np_EXP   ;   ioptio = ioptio + 1   ;   ENDIF
      IF( ln_dynspg_ts  ) THEN   ;   nspg =  np_TS    ;   ioptio = ioptio + 1   ;   ENDIF
      !
      IF( ioptio  > 1 )   CALL ctl_stop( 'Choose only one surface pressure gradient scheme' )
      IF( ioptio == 0 )   CALL ctl_warn( 'NO surface pressure gradient trend in momentum Eqs.' )
      IF( ln_dynspg_exp .AND. ln_isfcav )   &
           &   CALL ctl_stop( ' dynspg_exp not tested with ice shelf cavity ' )
      !
      IF( ln_dynspg_exp .AND. lk_RK3 )   &
           &   CALL ctl_stop( ' dynspg_exp not implemented with RK3 time-stepping ' )
      !
      IF(lwp) THEN
         WRITE(numout,*)
         IF( nspg == np_EXP )   WRITE(numout,*) '   ==>>>   explicit free surface'
         IF( nspg == np_TS  )   WRITE(numout,*) '   ==>>>   free surface with time splitting scheme'
         IF( nspg == np_NO  )   WRITE(numout,*) '   ==>>>   No surface surface pressure gradient trend in momentum Eqs.'
      ENDIF
      !
      SELECT CASE( nspg )
      CASE ( np_TS )             ! split-explicit scheme initialisation
         CALL dyn_spg_ts_init          ! do it first: set nn_e used to allocate some arrays later on
      CASE ( np_NO )             ! set barotropic fields to zeros
         IF( lk_RK3 ) THEN 
            uu_b(:,:,:) = 0._wp
            vv_b(:,:,:) = 0._wp
            un_adv(:,:) = 0._wp
            vn_adv(:,:) = 0._wp
         ENDIF
      END SELECT
!!      IF( nspg == np_TS ) THEN   ! split-explicit scheme initialisation
!!         CALL dyn_spg_ts_init          ! do it first: set nn_e used to allocate some arrays later on
!!      ENDIF
      !
   END SUBROUTINE dyn_spg_init

  !!======================================================================
END MODULE dynspg

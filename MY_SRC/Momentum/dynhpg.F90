MODULE dynhpg
   !!======================================================================
   !!                       ***  MODULE  dynhpg  ***
   !! Ocean dynamics:  hydrostatic pressure gradient trend
   !!======================================================================
   !! History :  OPA  !  1987-09  (P. Andrich, M.-A. Foujols)  hpg_zco: Original code
   !!            5.0  !  1991-11  (G. Madec)
   !!            7.0  !  1996-01  (G. Madec)  hpg_sco: Original code for s-coordinates
   !!            8.0  !  1997-05  (G. Madec)  split dynber into dynkeg and dynhpg
   !!            8.5  !  2002-07  (G. Madec)  F90: Free form and module
   !!            8.5  !  2002-08  (A. Bozec)  hpg_zps: Original code
   !!   NEMO     1.0  !  2005-10  (A. Beckmann, B.W. An)  various s-coordinate options
   !!                 !         Original code for hpg_ctl, hpg_hel hpg_wdj, hpg_djc, hpg_rot
   !!             -   !  2005-11  (G. Madec) style & small optimisation
   !!            3.3  !  2010-10  (C. Ethe, G. Madec)  reorganisation of initialisation phase
   !!            3.4  !  2011-11  (H. Liu)  hpg_prj: Original code for s-coordinates
   !!                 !           (A. Coward)  suppression of hel, wdj and rot options
   !!            3.6  !  2014-11  (P. Mathiot)  hpg_isf: original code for ice shelf cavity
   !!            4.2  !  2020-12  (M. Bell, A. Young)  hpg_djc: revised djc scheme
   !!            5.0  !  2022-12  (S. Techene, G. Madec)  remove hpg_zps
   !!             -   !  2024-02  (S. Techene, G. Madec)  RK3: dyn_hpg always called first
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_hpg       : update the momentum trend with the now horizontal
   !!                   gradient of the hydrostatic pressure
   !!   dyn_hpg_init  : initialisation and control of options
   !!       hpg_zco   : z-coordinate scheme
   !!       hpg_sco   : s-coordinate (standard jacobian formulation)
   !!       hpg_isf   : s-coordinate (sco formulation) adapted to ice shelf
   !!       hpg_djc   : s-coordinate (Density Jacobian with Cubic polynomial)
   !!       hpg_prj   : s-coordinate (Pressure Jacobian with Cubic polynomial)
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE sbc_oce        ! surface variable (only for the flag with ice shelf)
   USE dom_oce        ! ocean space and time domain
   USE eosbn2         ! compute density
   USE isf_oce , ONLY : risfload  ! ice shelf  (risfload variable)
   USE isfload , ONLY : isf_load  ! ice shelf  (isf_load routine )
   USE phycst         ! physical constants
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   USE wet_dry        ! wetting and drying
   !
   USE in_out_manager ! I/O manager
   USE prtctl         ! Print control
   USE lib_mpp        ! MPP library
   USE timing         ! Timing
   USE iom            ! I/O manager
   USE lib_fortran    ! to use sign with key_nosignedzero

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_hpg        ! routine called by step module
   PUBLIC   dyn_hpg_init   ! routine called by opa module

   !                                !!* Namelist namdyn_hpg : hydrostatic pressure gradient
   LOGICAL, PUBLIC ::   ln_hpg_zco   !: z-coordinate - full steps
   LOGICAL, PUBLIC ::   ln_hpg_sco   !: s-coordinate (standard jacobian formulation)
   LOGICAL, PUBLIC ::   ln_hpg_djc   !: s-coordinate (Density Jacobian with Cubic polynomial)
   LOGICAL, PUBLIC ::   ln_hpg_prj   !: s-coordinate (Pressure Jacobian scheme)
   LOGICAL, PUBLIC ::   ln_hpg_isf   !: s-coordinate similar to sco modify for isf

   !                                !! Flag to control the type of hydrostatic pressure gradient
   INTEGER, PARAMETER ::   np_ERROR  =-10   ! error in specification of lateral diffusion
   INTEGER, PARAMETER ::   np_zco    =  0   ! z-coordinate - full steps
   INTEGER, PARAMETER ::   np_sco    =  1   ! s-coordinate (standard jacobian formulation)
   INTEGER, PARAMETER ::   np_djc    =  2   ! s-coordinate (Density Jacobian with Cubic polynomial)
   INTEGER, PARAMETER ::   np_prj    =  3   ! s-coordinate (Pressure Jacobian scheme)
   INTEGER, PARAMETER ::   np_isf    =  4   ! s-coordinate similar to sco modify for isf
   !
   INTEGER, PUBLIC  ::   nhpg         !: type of pressure gradient scheme used ! (deduced from ln_hpg_... flags) (PUBLIC for TAM)
   !
   LOGICAL          ::   ln_hpg_djc_vnh, ln_hpg_djc_vnv                 ! flag to specify hpg_djc boundary condition type
   REAL(wp), PUBLIC ::   aco_bc_hor, bco_bc_hor, aco_bc_vrt, bco_bc_vrt !: coefficients for hpg_djc hor and vert boundary conditions

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
#  include "domzgr_substitute.h90"

   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_hpg( kt, Kmm, puu, pvv, Krhs )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_hpg  ***
      !!
      !! ** Method  :   Call the hydrostatic pressure gradient routine
      !!              using the scheme defined in the namelist
      !!
      !! ** Action : - set (puu,pvv)_Krhs with the Kmm HPG trend
      !!             as in RK3, dyn_hpg always the first computed trend.
      !!             - send trends to trd_dyn for diagnostics (l_trddyn=T)
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT(in   ) ::  kt, Kmm, Krhs   ! ocean time-step and -level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv        ! ocean velocities and RHS of momentum Eq.
#if ! defined key_RK3
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   ztrdu, ztrdv
#endif
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dyn_hpg')
      !
#if ! defined key_RK3
      IF( l_trddyn ) THEN                    ! Temporary saving of puu(:,:,:,Krhs) and pvv(:,:,:,Krhs) trends (l_trddyn)
         ALLOCATE( ztrdu(T2D(0),jpk), ztrdv(T2D(0),jpk) )
         ztrdu(:,:,:) = puu(T2D(0),:,Krhs)
         ztrdv(:,:,:) = pvv(T2D(0),:,Krhs)
      ENDIF
#endif
      !
      SELECT CASE ( nhpg )      ! Hydrostatic pressure gradient computation
      CASE ( np_zco )   ;   CALL hpg_zco    ( kt, Kmm, puu, pvv, Krhs )  ! z-coordinate
      CASE ( np_sco )   ;   CALL hpg_sco    ( kt, Kmm, puu, pvv, Krhs )  ! s-coordinate (standard jacobian formulation)
      CASE ( np_djc )   ;   CALL hpg_djc    ( kt, Kmm, puu, pvv, Krhs )  ! s-coordinate (Density Jacobian with Cubic polynomial)
      CASE ( np_prj )   ;   CALL hpg_prj    ( kt, Kmm, puu, pvv, Krhs )  ! s-coordinate (Pressure Jacobian scheme)
      CASE ( np_isf )   ;   CALL hpg_isf    ( kt, Kmm, puu, pvv, Krhs )  ! s-coordinate similar to sco modify for ice shelf
      END SELECT
      !
#if defined key_RK3
      IF( l_trddyn ) THEN      ! RK3 save the HPG trends for diagnostics directly with Krhs field
         CALL trd_dyn( puu(:,:,:,Krhs), pvv(:,:,:,Krhs), jpdyn_hpg, kt, Kmm )
      ENDIF
#else
      IF( l_trddyn ) THEN      ! save the hydrostatic pressure gradient trends for momentum trend diagnostics
         ztrdu(:,:,:) = puu(T2D(0),:,Krhs) - ztrdu(:,:,:)
         ztrdv(:,:,:) = pvv(T2D(0),:,Krhs) - ztrdv(:,:,:)
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_hpg_save, kt, Kmm )
         DEALLOCATE( ztrdu , ztrdv )
      ENDIF
#endif
      !
      IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab3d_1=puu(:,:,:,Krhs), clinfo1=' hpg  - Ua: ', mask1=umask,   &
         &                                  tab3d_2=pvv(:,:,:,Krhs), clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( ln_timing )   CALL timing_stop('dyn_hpg')
      !
   END SUBROUTINE dyn_hpg


   SUBROUTINE dyn_hpg_init( Kmm )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dyn_hpg_init  ***
      !!
      !! ** Purpose :   initializations for the hydrostatic pressure gradient
      !!              computation and consistency control
      !!
      !! ** Action  :   Read the namelist namdyn_hpg and check the consistency
      !!      with the type of vertical coordinate used (zco, zps, sco)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   Kmm   ! ocean time level index
      !!
      INTEGER ::   ioptio, ios   ! local integer
      !!
      NAMELIST/namdyn_hpg/ ln_hpg_zco    , ln_hpg_sco    ,               &
         &                 ln_hpg_djc    , ln_hpg_prj    , ln_hpg_isf,   &
         &                 ln_hpg_djc_vnh, ln_hpg_djc_vnv
      !!----------------------------------------------------------------------
      !
      READ_NML_REF(numnam,namdyn_hpg)
      READ_NML_CFG(numnam,namdyn_hpg)
      IF(lwm) WRITE ( numond, namdyn_hpg )
      !
      IF(lwp) THEN                   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_hpg_init : hydrostatic pressure gradient initialisation'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namdyn_hpg : choice of hpg scheme'
         WRITE(numout,*) '      z-coord. - full steps                             ln_hpg_zco    = ', ln_hpg_zco
         WRITE(numout,*) '      s-coord. (standard jacobian formulation)          ln_hpg_sco    = ', ln_hpg_sco
         WRITE(numout,*) '      s-coord. (standard jacobian formulation) for isf  ln_hpg_isf    = ', ln_hpg_isf
         WRITE(numout,*) '      s-coord. (Density Jacobian: Cubic polynomial)     ln_hpg_djc    = ', ln_hpg_djc
         WRITE(numout,*) '      s-coord. (Pressure Jacobian: Cubic polynomial)    ln_hpg_prj    = ', ln_hpg_prj
      ENDIF
      !
      IF( .NOT.lk_linssh .AND. ln_hpg_zco )   &
         &   CALL ctl_stop( 'dyn_hpg_init : non-linear free surface incompatible with hpg_zco' )
      !
      IF( (.NOT.ln_hpg_isf .AND. ln_isfcav) .OR. (ln_hpg_isf .AND. .NOT.ln_isfcav) )                  &
         &   CALL ctl_stop( 'dyn_hpg_init : ln_hpg_isf=T requires ln_isfcav=T and vice versa' )  
      !
      !
      !                               ! Set nhpg from ln_hpg_... flags & consistency check
      nhpg   = np_ERROR
      ioptio = 0
      IF( ln_hpg_zco ) THEN   ;   nhpg = np_zco   ;   ioptio = ioptio +1   ;   ENDIF
      IF( ln_hpg_sco ) THEN   ;   nhpg = np_sco   ;   ioptio = ioptio +1   ;   ENDIF
      IF( ln_hpg_djc ) THEN   ;   nhpg = np_djc   ;   ioptio = ioptio +1   ;   ENDIF
      IF( ln_hpg_prj ) THEN   ;   nhpg = np_prj   ;   ioptio = ioptio +1   ;   ENDIF
      IF( ln_hpg_isf ) THEN   ;   nhpg = np_isf   ;   ioptio = ioptio +1   ;   ENDIF
      !
      IF( ioptio /= 1 )   CALL ctl_stop( 'NO or several hydrostatic pressure gradient options used' )
      ! 
      IF(lwp) THEN
         WRITE(numout,*)
         SELECT CASE( nhpg )
         CASE( np_zco )   ;   WRITE(numout,*) '   ==>>>   z-coord. - full steps '
         CASE( np_sco )   ;   WRITE(numout,*) '   ==>>>   s-coord. (standard jacobian formulation)'
         CASE( np_djc )   ;   WRITE(numout,*) '   ==>>>   s-coord. (Density Jacobian: Cubic polynomial)'
         CASE( np_prj )   ;   WRITE(numout,*) '   ==>>>   s-coord. (Pressure Jacobian: Cubic polynomial)'
         CASE( np_isf )   ;   WRITE(numout,*) '   ==>>>   s-coord. (standard jacobian formulation) for isf'
         END SELECT
         WRITE(numout,*)
      ENDIF
      !                          
      IF ( ln_hpg_djc ) THEN
         IF (ln_hpg_djc_vnh) THEN ! Von Neumann boundary condition
           IF(lwp) WRITE(numout,*) '           horizontal bc: von Neumann '
           aco_bc_hor = 6.0_wp/5.0_wp
           bco_bc_hor = 7.0_wp/15.0_wp
         ELSE ! Linear extrapolation
           IF(lwp) WRITE(numout,*) '           horizontal bc: linear extrapolation'
           aco_bc_hor = 3.0_wp/2.0_wp
           bco_bc_hor = 1.0_wp/2.0_wp
         END IF
         IF (ln_hpg_djc_vnv) THEN ! Von Neumann boundary condition
           IF(lwp) WRITE(numout,*) '           vertical bc: von Neumann '
           aco_bc_vrt = 6.0_wp/5.0_wp
           bco_bc_vrt = 7.0_wp/15.0_wp
         ELSE ! Linear extrapolation
           IF(lwp) WRITE(numout,*) '           vertical bc: linear extrapolation'
           aco_bc_vrt = 3.0_wp/2.0_wp
           bco_bc_vrt = 1.0_wp/2.0_wp
         END IF
      END IF
      !
   END SUBROUTINE dyn_hpg_init


   SUBROUTINE hpg_zco( kt, Kmm, puu, pvv, Krhs )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_zco  ***
      !!
      !! ** Method  :   z-coordinate case, levels are horizontal surfaces.
      !!      The now hydrostatic pressure gradient at a given level, jk,
      !!      is computed by taking the vertical integral of the in-situ
      !!      density gradient along the model level from the suface to that
      !!      level:    zhpi = grav .....
      !!                zhpj = grav .....
      !!      set it to the general momentum trend (puu,pvv)_Krhs.
      !!            puu_Krhs = puu_Krhs - 1/e1u * zhpi
      !!            pvv_Krhs = pvv_Krhs - 1/e2v * zhpj
      !!
      !! ** Action : - set (puu,pvv)_Krh with the Kmm HPG trend
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT(in   ) ::  kt, Kmm, Krhs   ! ocean time-step and level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv        ! ocean velocities and RHS of momentum Eq.
      !
      INTEGER  ::   ji, jj, jk       ! dummy loop indices
      REAL(wp) ::   zcoef0, zcoef1   ! local scalars
      REAL(wp), DIMENSION(T2D(nn_hls)) ::   zhpi, zhpj   ! 2D workspace
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF( .NOT. l_istiled .OR. ntile == 1 )  THEN   ! Do only on the first tile
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn:hpg_zco : hydrostatic pressure gradient trend'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   z-coordinate case '
         ENDIF
      ENDIF
      !
      zcoef0 = - grav * 0.5_wp         ! Local constant initialization
      !
      DO_2D( 0, 0, 0, 0 )              ! Surface value
         zcoef1 = zcoef0 * e3w(ji,jj,1,Kmm)
         !                                ! hydrostatic pressure gradient
         zhpi(ji,jj) = zcoef1 * ( rhd(ji+1,jj,1) - rhd(ji,jj,1) ) * r1_e1u(ji,jj)
         zhpj(ji,jj) = zcoef1 * ( rhd(ji,jj+1,1) - rhd(ji,jj,1) ) * r1_e2v(ji,jj)
         !                                ! add to the general momentum trend
#if defined key_RK3
         puu(ji,jj,1,Krhs) = zhpi(ji,jj)     ! RK3 case: dyn_hpg always called first
         pvv(ji,jj,1,Krhs) = zhpj(ji,jj)
#else
         puu(ji,jj,1,Krhs) = puu(ji,jj,1,Krhs) + zhpi(ji,jj)
         pvv(ji,jj,1,Krhs) = pvv(ji,jj,1,Krhs) + zhpj(ji,jj)
#endif
      END_2D
      !
      DO_3D( 0, 0, 0, 0, 2, jpkm1 )     ! interior value (2=<jk=<jpkm1)
         zcoef1 = zcoef0 * e3w(ji,jj,jk,Kmm)
         !                                ! hydrostatic pressure gradient
         zhpi(ji,jj) = zhpi(ji,jj) + zcoef1 * (  ( rhd(ji+1,jj,jk)+rhd(ji+1,jj,jk-1) )  &
            &                                  - ( rhd(ji  ,jj,jk)+rhd(ji  ,jj,jk-1) )  ) * r1_e1u(ji,jj)

         zhpj(ji,jj) = zhpj(ji,jj) + zcoef1 * (  ( rhd(ji,jj+1,jk)+rhd(ji,jj+1,jk-1) )  &
            &                                  - ( rhd(ji,jj,  jk)+rhd(ji,jj  ,jk-1) )  ) * r1_e2v(ji,jj)
         !                                ! add to the general momentum trend
#if defined key_RK3
         puu(ji,jj,jk,Krhs) = zhpi(ji,jj)    ! RK3 case: dyn_hpg always called first
         pvv(ji,jj,jk,Krhs) = zhpj(ji,jj)
#else
         puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) + zhpi(ji,jj)
         pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) + zhpj(ji,jj)
#endif
      END_3D
      !
   END SUBROUTINE hpg_zco


   SUBROUTINE hpg_sco( kt, Kmm, puu, pvv, Krhs )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_sco  ***
      !!
      !! ** Method  :   s-coordinate case. Jacobian scheme.
      !!      The now hydrostatic pressure gradient at a given level, jk,
      !!      is computed by taking the vertical integral of the in-situ
      !!      density gradient along the model level from the suface to that
      !!      level. s-coordinates (l_sco): a corrective term is added
      !!      to the horizontal pressure gradient :
      !!         zhpi = grav .....  + 1/e1u mi(rhd) di[ grav dep3w ]
      !!         zhpj = grav .....  + 1/e2v mj(rhd) dj[ grav dep3w ]
      !!      add it to the general momentum trend (puu,pvv)_Krhs.
      !!         puu_Krhs = puu_Krhs - 1/e1u * zhpi
      !!         pvv_Krhs = pvv_Krhs - 1/e2v * zhpj
      !!
      !! ** Action : - set (puu,pvv)_Krh with the Kmm HPG trend
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT(in   ) ::  kt          ! ocean time-step index
      INTEGER                             , INTENT(in   ) ::  Kmm, Krhs   ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv    ! ocean velocities and RHS of momentum equation
      !!
      INTEGER  ::   ji, jj, jk, jii, jjj           ! dummy loop indices
      REAL(wp) ::   zcoef0, zuap, zvap, ztmp       ! local scalars
      REAL(wp), DIMENSION(T2D(0))           ::   zhpi, zhpj
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF( .NOT. l_istiled .OR. ntile == 1 )  THEN   ! Do only on the first tile
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn:hpg_sco : hydrostatic pressure gradient trend'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   s-coordinate case, OCE original scheme used'
         ENDIF
      ENDIF
      !
      zcoef0 = - grav * 0.5_wp
      !
      !
      DO_2D( 0, 0, 0, 0 )              ! Surface value
         !                                   ! hydrostatic pressure gradient along s-surfaces
         zhpi(ji,jj) = zcoef0 * r1_e1u(ji,jj)                      &
            &        * (  e3w(ji+1,jj  ,1,Kmm) * rhd(ji+1,jj  ,1)  &
            &           - e3w(ji  ,jj  ,1,Kmm) * rhd(ji  ,jj  ,1)  )
         zhpj(ji,jj) = zcoef0 * r1_e2v(ji,jj)                      &
            &        * (  e3w(ji  ,jj+1,1,Kmm) * rhd(ji  ,jj+1,1)  &
            &           - e3w(ji  ,jj  ,1,Kmm) * rhd(ji  ,jj  ,1)  )
         !                                   ! s-coordinate pressure gradient correction
         zuap = -zcoef0 * ( rhd     (ji+1,jj,1)     + rhd     (ji,jj,1) )   &
            &           * ( gdept_z0(ji+1,jj,1,Kmm) - gdept_z0(ji,jj,1,Kmm) ) * r1_e1u(ji,jj)
         zvap = -zcoef0 * ( rhd     (ji,jj+1,1)     + rhd     (ji,jj,1) )   &
            &           * ( gdept_z0(ji,jj+1,1,Kmm) - gdept_z0(ji,jj,1,Kmm) ) * r1_e2v(ji,jj)
         !
         !                                   ! add to the general momentum trend
#if defined key_RK3
         puu(ji,jj,1,Krhs) = zhpi(ji,jj) + zuap   ! RK3 case: dyn_hpg always called first
         pvv(ji,jj,1,Krhs) = zhpj(ji,jj) + zvap
#else
         puu(ji,jj,1,Krhs) = puu(ji,jj,1,Krhs) + zhpi(ji,jj) + zuap
         pvv(ji,jj,1,Krhs) = pvv(ji,jj,1,Krhs) + zhpj(ji,jj) + zvap
#endif
      END_2D
      !
      DO jk= 2, jpkm1
         DO_2D( 0, 0, 0, 0 )    ! interior value (2=<jk=<jpkm1)
            !                          ! hydrostatic pressure gradient along s-surfaces
            zhpi(ji,jj) = zhpi(ji,jj) + zcoef0 * r1_e1u(ji,jj)                         &
               &        * (  e3w(ji+1,jj,jk,Kmm) * ( rhd(ji+1,jj,jk) + rhd(ji+1,jj,jk-1) )  &
               &           - e3w(ji  ,jj,jk,Kmm) * ( rhd(ji  ,jj,jk) + rhd(ji  ,jj,jk-1) )  )
            zhpj(ji,jj) = zhpj(ji,jj) + zcoef0 * r1_e2v(ji,jj)                         &
               &        * (  e3w(ji,jj+1,jk,Kmm) * ( rhd(ji,jj+1,jk) + rhd(ji,jj+1,jk-1) )  &
               &           - e3w(ji,jj  ,jk,Kmm) * ( rhd(ji,jj,  jk) + rhd(ji,jj  ,jk-1) )  )
            !                         ! s-coordinate pressure gradient correction
            zuap = -zcoef0 * ( rhd     (ji+1,jj  ,jk)     + rhd     (ji,jj,jk)     ) &
               &           * ( gdept_z0(ji+1,jj  ,jk,Kmm) - gdept_z0(ji,jj,jk,Kmm) ) * r1_e1u(ji,jj)
            zvap = -zcoef0 * ( rhd     (ji  ,jj+1,jk)     + rhd     (ji,jj,jk)     ) &
               &           * ( gdept_z0(ji  ,jj+1,jk,Kmm) - gdept_z0(ji,jj,jk,Kmm) ) * r1_e2v(ji,jj)
            !                          ! add to the general momentum trend
#if defined key_RK3
            puu(ji,jj,jk,Krhs) = zhpi(ji,jj) + zuap   ! RK3 case: dyn_hpg always called first
            pvv(ji,jj,jk,Krhs) = zhpj(ji,jj) + zvap
#else
            puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) + zhpi(ji,jj) + zuap
            pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) + zhpj(ji,jj) + zvap
#endif
         END_2D
      END DO
      !
      !
   END SUBROUTINE hpg_sco


   SUBROUTINE hpg_isf( kt, Kmm, puu, pvv, Krhs )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_isf  ***
      !!
      !! ** Method  :   s-coordinate case. Jacobian scheme.
      !!      The now hydrostatic pressure gradient at a given level, jk,
      !!      is computed by taking the vertical integral of the in-situ
      !!      density gradient along the model level from the suface to that
      !!      level. s-coordinates (l_sco): a corrective term is added
      !!      to the horizontal pressure gradient :
      !!         zhpi = grav .....  + 1/e1u mi(rhd) di[ grav dep3w ]
      !!         zhpj = grav .....  + 1/e2v mj(rhd) dj[ grav dep3w ]
      !!      add it to the general momentum trend (puu,pvv)_Krhs.
      !!         puu_Krhs = puu_Krhs - 1/e1u * zhpi
      !!         pvv_Krhs = pvv_Krhs - 1/e2v * zhpj
      !!      iceload is added
      !!      
      !! ** Action : - set (puu,pvv)_Krhs with the Kmm HPG trend
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT(in   ) ::  kt          ! ocean time-step index
      INTEGER                             , INTENT(in   ) ::  Kmm, Krhs   ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv    ! ocean velocities and RHS of momentum equation
      !!
      INTEGER  ::   ji, jj, jk             ! dummy loop indices
      INTEGER  ::   ikt ,  ikti1,  iktj1   ! local integer
      REAL(wp) ::   ze3w, ze3wi1, ze3wj1   ! local scalars
      REAL(wp) ::   zcoef0, zuap, zvap     !   -      -
      REAL(wp), DIMENSION(T2D(1),jpts) ::  zts_top
      REAL(wp), DIMENSION(T2D(1))      ::  zrhd_top, zdep_top
      REAL(wp), DIMENSION(T2D(0))      ::  zhpi, zhpj
      !!----------------------------------------------------------------------
      !
      zcoef0 = - grav * 0.5_wp   ! Local constant initialization
      !
      !                          ! iniitialised to 0. zhpi zhpi   !!st is it necessary ???
      zhpi(:,:) = 0._wp   ;   zhpj(:,:) = 0._wp

      ! compute rhd at the ice/oce interface (ocean side)
      ! usefull to reduce residual current in the test case ISOMIP with no melting
      DO_2D( 1, 1, 1, 1 )
         ikt = mikt(ji,jj)
         zts_top(ji,jj,1) = ts(ji,jj,ikt,1,Kmm)
         zts_top(ji,jj,2) = ts(ji,jj,ikt,2,Kmm)
         zdep_top(ji,jj)  = MAX( risfdep(ji,jj) , gdept(ji,jj,1,Kmm) )
      END_2D
      ! eos 2D interface
      CALL eos( zts_top, zdep_top, zrhd_top, kbnd=1 )
      !
      !                     !===========================!
      !                     !=====  surface value  =====!
      !                     !===========================!
      DO_2D( 0, 0, 0, 0 )
         ikt   = mikt(ji  ,jj  )   ;   ze3w   = e3w(ji  ,jj  ,ikt  ,Kmm)
         ikti1 = mikt(ji+1,jj  )   ;   ze3wi1 = e3w(ji+1,jj  ,ikti1,Kmm)
         iktj1 = mikt(ji  ,jj+1)   ;   ze3wj1 = e3w(ji  ,jj+1,iktj1,Kmm)
         !                          ! hydrostatic pressure gradient along s-surfaces and ice shelf pressure
         !                          ! we assume ISF is in isostatic equilibrium
         zhpi(ji,jj) = zcoef0 * r1_e1u(ji,jj) * (   risfload(ji+1,jj) - risfload(ji,jj)  &
            &                                     + 0.5_wp * ( ze3wi1 * ( rhd(ji+1,jj,ikti1) + zrhd_top(ji+1,jj) )     &
            &                                                - ze3w   * ( rhd(ji  ,jj,ikt  ) + zrhd_top(ji  ,jj) ) )   )
         zhpj(ji,jj) = zcoef0 * r1_e2v(ji,jj) * (   risfload(ji,jj+1) - risfload(ji,jj)  &
            &                                     + 0.5_wp * ( ze3wj1 * ( rhd(ji,jj+1,iktj1) + zrhd_top(ji,jj+1) )      &
            &                                                - ze3w   * ( rhd(ji,jj  ,ikt  ) + zrhd_top(ji,jj  ) ) )   )
         !                          ! s-coordinate pressure gradient correction (=0 if z coordinate)
         zuap = -zcoef0 * ( rhd     (ji+1,jj,1)     + rhd     (ji,jj,1)     )   &
            &           * ( gdept_z0(ji+1,jj,1,Kmm) - gdept_z0(ji,jj,1,Kmm) ) * r1_e1u(ji,jj)
         zvap = -zcoef0 * ( rhd     (ji,jj+1,1)     + rhd     (ji,jj,1)     )   &
            &           * ( gdept_z0(ji,jj+1,1,Kmm) - gdept_z0(ji,jj,1,Kmm) ) * r1_e2v(ji,jj)
         !                          ! add to the general momentum trend
#if defined key_RK3
         puu(ji,jj,1,Krhs) = (zhpi(ji,jj) + zuap) * umask(ji,jj,1)   ! RK3 case: dyn_hpg always called first
         pvv(ji,jj,1,Krhs) = (zhpj(ji,jj) + zvap) * vmask(ji,jj,1)
#else
         puu(ji,jj,1,Krhs) = puu(ji,jj,1,Krhs) + (zhpi(ji,jj) + zuap) * umask(ji,jj,1)
         pvv(ji,jj,1,Krhs) = pvv(ji,jj,1,Krhs) + (zhpj(ji,jj) + zvap) * vmask(ji,jj,1)
#endif
      END_2D
      !   
      !                     !=============================!
      !                     !=====  interior values  =====!
      !                     !=============================!
      DO jk=2, jpkm1
         DO_2D( 0, 0, 0, 0 )
            ze3w   = e3w(ji  ,jj  ,jk,Kmm)
            ze3wi1 = e3w(ji+1,jj  ,jk,Kmm)
            ze3wj1 = e3w(ji  ,jj+1,jk,Kmm)
            !                          ! hydrostatic pressure gradient along s-surfaces
            zhpi(ji,jj) = zhpi(ji,jj) + zcoef0 / e1u(ji,jj)   &
               &        * (  ze3wi1 * ( rhd(ji+1,jj,jk) + rhd(ji+1,jj,jk-1) ) * wmask(ji+1,jj,jk)   &
               &           - ze3w   * ( rhd(ji  ,jj,jk) + rhd(ji  ,jj,jk-1) ) * wmask(ji  ,jj,jk)   )
            zhpj(ji,jj) = zhpj(ji,jj) + zcoef0 / e2v(ji,jj)   &
               &        * (  ze3wj1 * ( rhd(ji,jj+1,jk) + rhd(ji,jj+1,jk-1) ) * wmask(ji,jj+1,jk)   &
               &           - ze3w   * ( rhd(ji,jj,  jk) + rhd(ji,jj  ,jk-1) ) * wmask(ji,jj  ,jk)   )
            !                          ! s-coordinate pressure gradient correction
            zuap = -zcoef0 * ( rhd     (ji+1,jj  ,jk)     + rhd     (ji,jj,jk) )   &
               &           * ( gdept_z0(ji+1,jj  ,jk,Kmm) - gdept_z0(ji,jj,jk,Kmm) ) / e1u(ji,jj)
            zvap = -zcoef0 * ( rhd     (ji  ,jj+1,jk)     + rhd     (ji,jj,jk) )   &
               &           * ( gdept_z0(ji  ,jj+1,jk,Kmm) - gdept_z0(ji,jj,jk,Kmm) ) / e2v(ji,jj)
            !                          ! add to the general momentum trend
#if defined key_RK3
            puu(ji,jj,jk,Krhs) = (zhpi(ji,jj) + zuap) * umask(ji,jj,jk)   ! RK3 case: dyn_hpg always called first
            pvv(ji,jj,jk,Krhs) = (zhpj(ji,jj) + zvap) * vmask(ji,jj,jk)
#else
            puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) + (zhpi(ji,jj) + zuap) * umask(ji,jj,jk)
            pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) + (zhpj(ji,jj) + zvap) * vmask(ji,jj,jk)
#endif
         END_2D
      END DO
      !
   END SUBROUTINE hpg_isf


   SUBROUTINE hpg_djc( kt, Kmm, puu, pvv, Krhs )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_djc  ***
      !!
      !! ** Method  :   Density Jacobian with Cubic polynomial scheme
      !!
      !! ** Action : - set (puu,pvv)_Krhs with the Kmm HPG trend
      !!
      !! Reference: Shchepetkin and McWilliams, J. Geophys. Res. 2003
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT( in )  ::  kt          ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  Kmm, Krhs   ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv    ! ocean velocities and RHS of momentum equation
      !!
      INTEGER  ::   ji, jj, jk          ! dummy loop indices
      INTEGER  ::   iktb, iktt          ! jk indices at tracer points for top and bottom points 
      REAL(wp) ::   zcoef0, zep, cffw   ! temporary scalars
      REAL(wp) ::   z_grav_10, z1_12, z1_cff
      REAL(wp) ::   cffu, cffx          !    "         "
      REAL(wp) ::   cffv, cffy          !    "         "
      REAL(wp), DIMENSION(T2D(nn_hls),jpk) ::   zhpi, zhpj

      REAL(wp), DIMENSION(T2D(nn_hls),jpk) ::   zdzx, zdzy, zdzz                          ! Primitive grid differences ('delta_xyz')
      REAL(wp), DIMENSION(T2D(nn_hls),jpk) ::   zdz_i, zdz_j, zdz_k                       ! Harmonic average of primitive grid differences ('d_xyz')
      REAL(wp), DIMENSION(T2D(nn_hls),jpk) ::   zdrhox, zdrhoy, zdrhoz                    ! Primitive rho differences ('delta_rho')
      REAL(wp), DIMENSION(T2D(nn_hls),jpk) ::   zdrho_i, zdrho_j, zdrho_k                 ! Harmonic average of primitive rho differences ('d_rho')
      REAL(wp), DIMENSION(T2D(nn_hls),jpk) ::   z_rho_i, z_rho_j, z_rho_k                 ! Face intergrals
      REAL(wp), DIMENSION(T2D(nn_hls))     ::   zz_dz_i, zz_dz_j, zz_drho_i, zz_drho_j    ! temporary arrays
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF( .NOT. l_istiled .OR. ntile == 1 )  THEN   ! Do only on the first tile
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn:hpg_djc : hydrostatic pressure gradient trend'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   s-coordinate case, density Jacobian with cubic polynomial scheme'
         ENDIF
      ENDIF
      !
      ! Local constant initialization
      zcoef0     = - grav * 0.5_wp
      z_grav_10  =   grav / 10._wp
      z1_12      = 1.0_wp / 12._wp

      !----------------------------------------------------------------------------------------
      !  1. compute and store elementary vertical differences in provisional arrays 
      !----------------------------------------------------------------------------------------

!!bug gm   Not a true bug, but... zdzz=e3w  for zdzx, zdzy verify what it is really

      DO_3D( 1, 1, 1, 1, 2, jpkm1 )  
         zdrhoz(ji,jj,jk) =   rhd     (ji,jj,jk    ) - rhd     (ji,jj,jk-1    )
         zdzz  (ji,jj,jk) = - gdept_z0(ji,jj,jk,Kmm) + gdept_z0(ji,jj,jk-1,Kmm)
      END_3D

      !-------------------------------------------------------------------------
      ! 2. compute harmonic averages for vertical differences using eq. 5.18
      !-------------------------------------------------------------------------
      zep = 1.e-15

!! mb zdrho_k, zdz_k, zdrho_i, zdz_i, zdrho_j, zdz_j re-centred about the point (ji,jj,jk) 
      zdrho_k(:,:,:) = 0._wp
      zdz_k  (:,:,:) = 0._wp

      DO_3D( 1, 1, 1, 1, 2, jpk-2 )
         cffw = MAX( 2._wp * zdrhoz(ji,jj,jk) * zdrhoz(ji,jj,jk+1), 0._wp )
         z1_cff = zdrhoz(ji,jj,jk) + zdrhoz(ji,jj,jk+1)
         zdrho_k(ji,jj,jk) = cffw / SIGN( MAX( ABS(z1_cff), zep ), z1_cff )
         zdz_k(ji,jj,jk) = 2._wp *   zdzz(ji,jj,jk) * zdzz(ji,jj,jk+1)   &
            &                    / ( zdzz(ji,jj,jk) + zdzz(ji,jj,jk+1) )
      END_3D

      !----------------------------------------------------------------------------------
      ! 3. apply boundary conditions at top and bottom using 5.36-5.37
      !----------------------------------------------------------------------------------

! mb for sea-ice shelves we will need to re-write this upper boundary condition in the same form as the lower boundary condition
      DO_2D( 1, 1, 1, 1 )
         zdrho_k(ji,jj,1) = aco_bc_vrt * ( rhd     (ji,jj,2    ) - rhd     (ji,jj,1    ) ) - bco_bc_vrt * zdrho_k(ji,jj,2)
         zdz_k  (ji,jj,1) = aco_bc_vrt * (-gdept_z0(ji,jj,2,Kmm) + gdept_z0(ji,jj,1,Kmm) ) - bco_bc_vrt * zdz_k  (ji,jj,2)
      END_2D

      DO_2D( 1, 1, 1, 1 )
         IF ( mbkt(ji,jj)>1 ) THEN
            iktb = mbkt(ji,jj)
            zdrho_k(ji,jj,iktb) = aco_bc_vrt * (      rhd(ji,jj,iktb    ) -      rhd(ji,jj,iktb-1    ) ) - bco_bc_vrt * zdrho_k(ji,jj,iktb-1)
            zdz_k  (ji,jj,iktb) = aco_bc_vrt * (-gdept_z0(ji,jj,iktb,Kmm) + gdept_z0(ji,jj,iktb-1,Kmm) ) - bco_bc_vrt * zdz_k  (ji,jj,iktb-1)
         END IF
      END_2D

      !--------------------------------------------------------------
      ! 4. Compute side face integrals
      !-------------------------------------------------------------

!! ssh replaces e3w_n ; gde3w is a depth; the formulae involve heights  
!! rho_k stores grav * FX / rho_0  

      !--------------------------------------------------------------
      ! 4. a) Upper half of top-most grid box, compute and store
      !-------------------------------------------------------------
! *** AY note: ssh(ji,jj,Kmm) + gde3w(ji,jj,1) = e3w(ji,jj,1,Kmm)
      DO_2D( 0, 1, 0, 1)
         z_rho_k(ji,jj,1) =  grav * gdept(ji,jj,1,Kmm)                             & 
            &                     * (             rhd(ji,jj,1)                     &
            &                         -0.5_wp * ( rhd(ji,jj,2) - rhd(ji,jj,1) )    &
            &                              * gdept(ji,jj,1,Kmm) / e3w(ji,jj,2,Kmm) &
            &                       )
      END_2D

      !--------------------------------------------------------------
      ! 4. b) Interior faces, compute and store
      !-------------------------------------------------------------

      DO_3D( 0, 1, 0, 1, 2, jpkm1 )
         z_rho_k(ji,jj,jk) = zcoef0 * (   rhd     (ji,jj,jk    ) + rhd     (ji,jj,jk-1    ) )                                 &
            &                       * ( - gdept_z0(ji,jj,jk,Kmm) + gdept_z0(ji,jj,jk-1,Kmm) )                                 &
            &                       + z_grav_10 * (                                                                           &
            &     (     zdrho_k (ji,jj,jk    ) - zdrho_k (ji,jj,jk-1    ) )                                                   &
            &   * ( - ( gdept_z0(ji,jj,jk,Kmm) - gdept_z0(ji,jj,jk-1,Kmm) ) - z1_12 * ( zdz_k  (ji,jj,jk) + zdz_k  (ji,jj,jk-1) ) )   &
            &   - (     zdz_k   (ji,jj,jk)     - zdz_k   (ji,jj,jk-1)     )                                                           &
            &   * (   ( rhd     (ji,jj,jk)     - rhd     (ji,jj,jk-1)     ) - z1_12 * ( zdrho_k(ji,jj,jk) + zdrho_k(ji,jj,jk-1) ) )   &
            &                                      )
      END_3D

      !----------------------------------------------------------------------------------------
      !  5. compute and store elementary horizontal differences in provisional arrays 
      !----------------------------------------------------------------------------------------
      zdrhox(:,:,:) = 0._wp
      zdzx  (:,:,:) = 0._wp
      zdrhoy(:,:,:) = 0._wp
      zdzy  (:,:,:) = 0._wp

      DO_3D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 1, jpkm1 )
         zdrhox(ji,jj,jk) = rhd     (ji+1,jj  ,jk    ) - rhd     (ji  ,jj  ,jk)
         zdzx  (ji,jj,jk) = gdept_z0(ji  ,jj  ,jk,Kmm) - gdept_z0(ji+1,jj  ,jk,Kmm)
         zdrhoy(ji,jj,jk) = rhd     (ji  ,jj+1,jk    ) - rhd     (ji  ,jj  ,jk)
         zdzy  (ji,jj,jk) = gdept_z0(ji  ,jj  ,jk,Kmm) - gdept_z0(ji  ,jj+1,jk,Kmm)
      END_3D

      !-------------------------------------------------------------------------
      ! 6. compute harmonic averages using eq. 5.18
      !-------------------------------------------------------------------------

      DO_3D( 0, 1, 0, 1, 1, jpkm1 )
         cffu = MAX( 2._wp * zdrhox(ji-1,jj,jk) * zdrhox(ji,jj,jk), 0._wp )
         z1_cff = zdrhox(ji-1,jj,jk) + zdrhox(ji,jj,jk)
         zdrho_i(ji,jj,jk) = cffu / SIGN( MAX( ABS(z1_cff), zep ), z1_cff )

         cffx = MAX( 2._wp * zdzx(ji-1,jj,jk)   * zdzx(ji,jj,jk), 0._wp )
         z1_cff = zdzx(ji-1,jj,jk)   + zdzx(ji,jj,jk)
         zdz_i(ji,jj,jk)   = cffx / SIGN( MAX( ABS(z1_cff), zep ), z1_cff )

         cffv = MAX( 2._wp * zdrhoy(ji,jj-1,jk) * zdrhoy(ji,jj,jk), 0._wp )
         z1_cff = zdrhoy(ji,jj-1,jk) + zdrhoy(ji,jj,jk)
         zdrho_j(ji,jj,jk) = cffv / SIGN( MAX( ABS(z1_cff), zep ), z1_cff )

         cffy = MAX( 2._wp * zdzy(ji,jj-1,jk)   * zdzy(ji,jj,jk), 0._wp )
         z1_cff = zdzy(ji,jj-1,jk)   + zdzy(ji,jj,jk)
         zdz_j(ji,jj,jk)   = cffy / SIGN( MAX( ABS(z1_cff), zep ), z1_cff )
      END_3D
      
!!! Note that zdzx, zdzy, zdzz, zdrhox, zdrhoy and zdrhoz should NOT be used beyond this point      

      !----------------------------------------------------------------------------------
      ! 6B. apply boundary conditions at side boundaries using 5.36-5.37
      !----------------------------------------------------------------------------------

      DO jk = 1, jpkm1
         zz_drho_i(:,:) = zdrho_i(:,:,jk)
         zz_dz_i  (:,:) = zdz_i  (:,:,jk)
         zz_drho_j(:,:) = zdrho_j(:,:,jk)
         zz_dz_j  (:,:) = zdz_j  (:,:,jk)
         ! Walls coming from left: should check from 2 to jpi-1 (and jpj=2-jpj)
         DO_2D( 0, 0, 0, 1 )
            IF ( umask(ji,jj,jk) > 0.5_wp .AND. umask(ji-1,jj,jk) < 0.5_wp .AND. umask(ji+1,jj,jk) > 0.5_wp)  THEN
               zz_drho_i(ji,jj) = aco_bc_hor * ( rhd    (ji+1,jj,jk) - rhd    (ji,jj,jk) ) - bco_bc_hor * zdrho_i(ji+1,jj,jk)
               zz_dz_i  (ji,jj) = aco_bc_hor * (-gdept_z0(ji+1,jj,jk,Kmm) + gdept_z0(ji,jj,jk,Kmm) ) - bco_bc_hor * zdz_i  (ji+1,jj,jk)
            END IF
         END_2D
         ! Walls coming from right: should check from 3 to jpi (and jpj=2-jpj)
         DO_2D( -1, 1, 0, 1 )
            IF ( umask(ji,jj,jk) < 0.5_wp .AND. umask(ji-1,jj,jk) > 0.5_wp .AND. umask(ji-2,jj,jk) > 0.5_wp) THEN
               zz_drho_i(ji,jj) = aco_bc_hor * ( rhd     (ji,jj,jk    ) - rhd     (ji-1,jj,jk    ) ) - bco_bc_hor * zdrho_i(ji-1,jj,jk)
               zz_dz_i  (ji,jj) = aco_bc_hor * (-gdept_z0(ji,jj,jk,Kmm) + gdept_z0(ji-1,jj,jk,Kmm) ) - bco_bc_hor * zdz_i  (ji-1,jj,jk)
            END IF
         END_2D
         ! Walls coming from left: should check from 2 to jpj-1 (and jpi=2-jpi)
         DO_2D( 0, 1, 0, 0 )
            IF ( vmask(ji,jj,jk) > 0.5_wp .AND. vmask(ji,jj-1,jk) < 0.5_wp .AND. vmask(ji,jj+1,jk) > 0.5_wp)  THEN
               zz_drho_j(ji,jj) = aco_bc_hor * ( rhd     (ji,jj+1,jk    ) - rhd     (ji,jj,jk    ) ) - bco_bc_hor * zdrho_j(ji,jj+1,jk)
               zz_dz_j  (ji,jj) = aco_bc_hor * (-gdept_z0(ji,jj+1,jk,Kmm) + gdept_z0(ji,jj,jk,Kmm) ) - bco_bc_hor * zdz_j  (ji,jj+1,jk)
            END IF
         END_2D
         ! Walls coming from right: should check from 3 to jpj (and jpi=2-jpi)
         DO_2D( 0, 1, -1, 1 )
            IF ( vmask(ji,jj,jk) < 0.5_wp .AND. vmask(ji,jj-1,jk) > 0.5_wp .AND. vmask(ji,jj-2,jk) > 0.5_wp) THEN
               zz_drho_j(ji,jj) = aco_bc_hor * ( rhd     (ji,jj,jk    ) - rhd     (ji,jj-1,jk    ) ) - bco_bc_hor * zdrho_j(ji,jj-1,jk)
               zz_dz_j  (ji,jj) = aco_bc_hor * (-gdept_z0(ji,jj,jk,Kmm) + gdept_z0(ji,jj-1,jk,Kmm) ) - bco_bc_hor * zdz_j  (ji,jj-1,jk)
            END IF
         END_2D
         zdrho_i(:,:,jk) = zz_drho_i(:,:)
         zdz_i  (:,:,jk) = zz_dz_i  (:,:)
         zdrho_j(:,:,jk) = zz_drho_j(:,:)
         zdz_j  (:,:,jk) = zz_dz_j  (:,:)
      END DO

      !--------------------------------------------------------------
      ! 7. Calculate integrals on side faces  
      !-------------------------------------------------------------

      DO_3D( 0, 0, 0, 0, 1, jpkm1 )
! two -ve signs cancel in next two lines (within zcoef0 and because gde3w is a depth not a height)
         z_rho_i(ji,jj,jk) = zcoef0 * ( rhd   (ji+1,jj,jk    ) + rhd     (ji,jj,jk    ) )                                       &
             &                    * ( gdept_z0(ji+1,jj,jk,Kmm) - gdept_z0(ji,jj,jk,Kmm) )
         IF ( umask(ji-1, jj, jk) > 0.5 .OR. umask(ji+1, jj, jk) > 0.5 ) THEN
            z_rho_i(ji,jj,jk) = z_rho_i(ji,jj,jk) - z_grav_10 * (                                                               &
             &     (     zdrho_i  (ji+1,jj,jk)    - zdrho_i  (ji,jj,jk) )                                                            &
             &   * ( - ( gdept_z0(ji+1,jj,jk,Kmm) - gdept_z0(ji,jj,jk,Kmm) ) - z1_12 * ( zdz_i  (ji+1,jj,jk) + zdz_i  (ji,jj,jk) ) ) &
             &   - (     zdz_i    (ji+1,jj,jk)    - zdz_i    (ji,jj,jk) )                                                            &
             &   * (   ( rhd    (ji+1,jj,jk)      - rhd    (ji,jj,jk)      ) - z1_12 * ( zdrho_i(ji+1,jj,jk) + zdrho_i(ji,jj,jk) ) ) &
             &                                                  )
         END IF
  
         z_rho_j(ji,jj,jk) = zcoef0 * ( rhd   (ji,jj+1,jk    ) + rhd     (ji,jj,jk    ) )                                       &
             &                    * ( gdept_z0(ji,jj+1,jk,Kmm) - gdept_z0(ji,jj,jk,Kmm) )
         IF ( vmask(ji, jj-1, jk) > 0.5 .OR. vmask(ji, jj+1, jk) > 0.5 ) THEN
            z_rho_j(ji,jj,jk) = z_rho_j(ji,jj,jk) - z_grav_10 * (                                                               &
             &     (     zdrho_j (ji,jj+1,jk)     - zdrho_j (ji,jj,jk)       )                                                        &
             &   * ( - ( gdept_z0(ji,jj+1,jk,Kmm) - gdept_z0(ji,jj,jk,Kmm) ) - z1_12 * ( zdz_j  (ji,jj+1,jk) + zdz_j  (ji,jj,jk) ) )  &
             &   - (        zdz_j(ji,jj+1,jk)     -    zdz_j(ji,jj,jk)       )                                                        &
             &   * (   (      rhd(ji,jj+1,jk)     -      rhd(ji,jj,jk)     ) - z1_12 * ( zdrho_j(ji,jj+1,jk) + zdrho_j(ji,jj,jk) ) )  &
             &                                                  )
         END IF
      END_3D

      !--------------------------------------------------------------
      ! 8. Integrate in the vertical   
      !-------------------------------------------------------------
      !
      ! ---------------
      !  Surface value
      ! ---------------
      DO_2D( 0, 0, 0, 0 )
         zhpi(ji,jj,1) = ( ( z_rho_k(ji,jj,1) - z_rho_k(ji+1,jj  ,1) ) - z_rho_i(ji,jj,1) ) * r1_e1u(ji,jj)   ! add () for NP repro
         zhpj(ji,jj,1) = ( ( z_rho_k(ji,jj,1) - z_rho_k(ji  ,jj+1,1) ) - z_rho_j(ji,jj,1) ) * r1_e2v(ji,jj)
         ! add to the general momentum trend
#if defined key_RK3
         puu(ji,jj,1,Krhs) = zhpi(ji,jj,1)   ! RK3 case: dyn_hpg always called first
         pvv(ji,jj,1,Krhs) = zhpj(ji,jj,1)
#else
         puu(ji,jj,1,Krhs) = puu(ji,jj,1,Krhs) + zhpi(ji,jj,1)
         pvv(ji,jj,1,Krhs) = pvv(ji,jj,1,Krhs) + zhpj(ji,jj,1)
#endif
      END_2D

      ! ----------------
      !  interior value   (2=<jk=<jpkm1)
      ! ----------------
      DO_3D( 0, 0, 0, 0, 2, jpkm1 )
         ! hydrostatic pressure gradient along s-surfaces
         zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1)                                                     &
            &           + (  ( z_rho_k(ji,jj,jk) - z_rho_k(ji+1,jj,jk  ) )                     &
            &              - ( z_rho_i(ji,jj,jk) - z_rho_i(ji  ,jj,jk-1) )  ) * r1_e1u(ji,jj)
         zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1)                                                     &
            &           + (  ( z_rho_k(ji,jj,jk) - z_rho_k(ji,jj+1,jk  ) )                     &
            &               -( z_rho_j(ji,jj,jk) - z_rho_j(ji,jj  ,jk-1) )  ) * r1_e2v(ji,jj)
         ! add to the general momentum trend
#if defined key_RK3
         puu(ji,jj,jk,Krhs) = zhpi(ji,jj,jk)   ! RK3 case: dyn_hpg always called first
         pvv(ji,jj,jk,Krhs) = zhpj(ji,jj,jk)
#else
         puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) + zhpi(ji,jj,jk)
         pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) + zhpj(ji,jj,jk)
#endif
      END_3D
      !
      !
   END SUBROUTINE hpg_djc


   SUBROUTINE hpg_prj( kt, Kmm, puu, pvv, Krhs )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_prj  ***
      !!
      !! ** Method  :   s-coordinate case.
      !!      A Pressure-Jacobian horizontal pressure gradient method
      !!      based on the constrained cubic-spline interpolation for
      !!      all vertical coordinate systems
      !!
      !! ** Action : - set (puu,pvv)_Krhs with the Kmm HPG trend
      !!----------------------------------------------------------------------
      INTEGER, PARAMETER  :: polynomial_type = 1    ! 1: cubic spline, 2: linear
      INTEGER                             , INTENT( in )  ::  kt          ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  Kmm, Krhs   ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv    ! ocean velocities and RHS of momentum equation
      !!
      INTEGER  ::   ji, jj, jk, jkk                 ! dummy loop indices
      REAL(wp) ::   zcoef0, znad                    ! local scalars
      !
      !! The local variables for the correction term
      INTEGER  :: jk1, jis, jid, jjs, jjd
      REAL(wp) :: zuijk, zvijk, zpwes, zpwed, zpnss, zpnsd, zdeps
      REAL(wp) :: zrhdt1
      REAL(wp) :: zdpdx1, zdpdx2, zdpdy1, zdpdy2
      REAL(wp), DIMENSION(T2D(nn_hls))     ::   zpgu, zpgv   ! 2D workspace
      REAL(wp), DIMENSION(T2D(nn_hls))     ::   zsshu_n, zsshv_n
      REAL(wp), DIMENSION(T2D(nn_hls),jpk) ::   zdept, zrhh
      REAL(wp), DIMENSION(T2D(nn_hls),jpk) ::   zhpi, zu, zv, fsp, xsp, asp, bsp, csp, dsp
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF( .NOT. l_istiled .OR. ntile == 1 )  THEN   ! Do only on the first tile
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn:hpg_prj : hydrostatic pressure gradient trend'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   s-coordinate case, cubic spline pressure Jacobian'
         ENDIF
      ENDIF

      ! Local constant initialization
      zcoef0 = - grav
      znad = 1._wp
      IF( lk_linssh )   znad = 1._wp
      !
      ! ---------------
      !  Surface pressure gradient to be removed
      ! ---------------
      DO_2D( 0, 0, 0, 0 )
         zpgu(ji,jj) = - grav * ( ssh(ji+1,jj,Kmm) - ssh(ji,jj,Kmm) ) * r1_e1u(ji,jj)
         zpgv(ji,jj) = - grav * ( ssh(ji,jj+1,Kmm) - ssh(ji,jj,Kmm) ) * r1_e2v(ji,jj)
      END_2D
      !

      ! Clean 3-D work arrays
      zhpi(:,:,:) = 0._wp
      zrhh(:,:,:) = rhd(T2D(nn_hls),:)

      ! Preparing vertical density profile "zrhh(:,:,:)" for hybrid-sco coordinate
      DO_2D( 1, 1, 1, 1 )
         jk = mbkt(ji,jj)
         IF(     jk <=  1   ) THEN   ;   zrhh(ji,jj,    :   ) = 0._wp
         ELSEIF( jk ==  2   ) THEN   ;   zrhh(ji,jj,jk+1:jpk) = rhd(ji,jj,jk)
         ELSEIF( jk < jpkm1 ) THEN
            DO jkk = jk+1, jpk
               zrhh(ji,jj,jkk) = interp1(gdept_z0(ji,jj,jkk  ,Kmm), gdept_z0(ji,jj,jkk-1,Kmm),   &
                  &                      gdept_z0(ji,jj,jkk-2,Kmm), zrhh (ji,jj,jkk-1), zrhh(ji,jj,jkk-2))
            END DO
         ENDIF
      END_2D

      ! Transfer the depth of "T(:,:,:)" to vertical coordinate "zdept(:,:,:)"
      DO_2D( 1, 1, 1, 1 )
         zdept(ji,jj,1) = 0.5_wp * e3w(ji,jj,1,Kmm) - ssh(ji,jj,Kmm)
      END_2D

      DO_3D( 1, 1, 1, 1, 2, jpk )
         zdept(ji,jj,jk) = zdept(ji,jj,jk-1) + e3w(ji,jj,jk,Kmm)
      END_3D

      fsp(:,:,:) = zrhh (:,:,:)
      xsp(:,:,:) = zdept(:,:,:)

      ! Construct the vertical density profile with the
      ! constrained cubic spline interpolation
      ! rho(z) = asp + bsp*z + csp*z^2 + dsp*z^3
      CALL cspline( fsp, xsp, asp, bsp, csp, dsp, polynomial_type )

      ! Integrate the hydrostatic pressure "zhpi(:,:,:)" at "T(ji,jj,1)"
      DO_2D( 0, 1, 0, 1 )
         zrhdt1 = zrhh(ji,jj,1) - interp3( zdept(ji,jj,1), asp(ji,jj,1), bsp(ji,jj,1),  &
            &                                              csp(ji,jj,1), dsp(ji,jj,1) ) * 0.25_wp * e3w(ji,jj,1,Kmm)

         ! assuming linear profile across the top half surface layer
         zhpi(ji,jj,1) =  0.5_wp * e3w(ji,jj,1,Kmm) * zrhdt1
      END_2D

      ! Calculate the pressure "zhpi(:,:,:)" at "T(ji,jj,2:jpkm1)"
      DO_3D( 0, 1, 0, 1, 2, jpkm1 )
         zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1) +                                  &
            &             integ_spline( zdept(ji,jj,jk-1), zdept(ji,jj,jk),   &
            &                           asp  (ji,jj,jk-1), bsp  (ji,jj,jk-1), &
            &                           csp  (ji,jj,jk-1), dsp  (ji,jj,jk-1)  )
      END_3D

      ! Z coordinate of U(ji,jj,1:jpkm1) and V(ji,jj,1:jpkm1)

      ! Prepare zsshu_n and zsshv_n
      DO_2D( 0, 0, 0, 0 )
!!gm BUG ?    if it is ssh at u- & v-point then it should be:
!          zsshu_n(ji,jj) = (e1e2t(ji,jj) * ssh(ji,jj,Kmm) + e1e2t(ji+1,jj) * ssh(ji+1,jj,Kmm)) * &
!                         & r1_e1e2u(ji,jj) * umask(ji,jj,1) * 0.5_wp 
!          zsshv_n(ji,jj) = (e1e2t(ji,jj) * ssh(ji,jj,Kmm) + e1e2t(ji,jj+1) * ssh(ji,jj+1,Kmm)) * &
!                         & r1_e1e2v(ji,jj) * vmask(ji,jj,1) * 0.5_wp 
!!gm not this:
         zsshu_n(ji,jj) = (e1e2u(ji,jj) * ssh(ji,jj,Kmm) + e1e2u(ji+1, jj) * ssh(ji+1,jj,Kmm)) * &
                        & r1_e1e2u(ji,jj) * umask(ji,jj,1) * 0.5_wp
         zsshv_n(ji,jj) = (e1e2v(ji,jj) * ssh(ji,jj,Kmm) + e1e2v(ji+1, jj) * ssh(ji,jj+1,Kmm)) * &
                        & r1_e1e2v(ji,jj) * vmask(ji,jj,1) * 0.5_wp
      END_2D

      DO_2D( 0, 0, 0, 0 )
         zu(ji,jj,1) = - ( e3u(ji,jj,1,Kmm) - zsshu_n(ji,jj) )
         zv(ji,jj,1) = - ( e3v(ji,jj,1,Kmm) - zsshv_n(ji,jj) )
      END_2D

      DO_3D( 0, 0, 0, 0, 2, jpkm1 )
         zu(ji,jj,jk) = zu(ji,jj,jk-1) - e3u(ji,jj,jk,Kmm)
         zv(ji,jj,jk) = zv(ji,jj,jk-1) - e3v(ji,jj,jk,Kmm)
      END_3D

      DO_3D( 0, 0, 0, 0, 1, jpkm1 )
         zu(ji,jj,jk) = zu(ji,jj,jk) + 0.5_wp * e3u(ji,jj,jk,Kmm)
         zv(ji,jj,jk) = zv(ji,jj,jk) + 0.5_wp * e3v(ji,jj,jk,Kmm)
      END_3D

      DO_3D( 0, 0, 0, 0, 1, jpkm1 )
         zu(ji,jj,jk) = MIN(  zu(ji,jj,jk) , MAX( -zdept(ji,jj,jk) , -zdept(ji+1,jj,jk) )  )
         zu(ji,jj,jk) = MAX(  zu(ji,jj,jk) , MIN( -zdept(ji,jj,jk) , -zdept(ji+1,jj,jk) )  )
         zv(ji,jj,jk) = MIN(  zv(ji,jj,jk) , MAX( -zdept(ji,jj,jk) , -zdept(ji,jj+1,jk) )  )
         zv(ji,jj,jk) = MAX(  zv(ji,jj,jk) , MIN( -zdept(ji,jj,jk) , -zdept(ji,jj+1,jk) )  )
      END_3D


      DO_3D( 0, 0, 0, 0, 1, jpkm1 )
         zpwes = 0._wp; zpwed = 0._wp
         zpnss = 0._wp; zpnsd = 0._wp
         zuijk = zu(ji,jj,jk)
         zvijk = zv(ji,jj,jk)

         !!!!!     for u equation
         IF( jk <= mbku(ji,jj) ) THEN
            IF( -zdept(ji+1,jj,jk) >= -zdept(ji,jj,jk) ) THEN
              jis = ji + 1; jid = ji
            ELSE
              jis = ji;     jid = ji +1
            ENDIF

            ! integrate the pressure on the shallow side
            jk1 = jk
            DO WHILE ( -zdept(jis,jj,jk1) > zuijk )
               IF( jk1 == mbku(ji,jj) ) THEN
                  zuijk = -zdept(jis,jj,jk1)
                  EXIT
               ENDIF
               zdeps = MIN(zdept(jis,jj,jk1+1), -zuijk)
               zpwes = zpwes +                                      &
                  integ_spline(zdept(jis,jj,jk1), zdeps,            &
                                 asp(jis,jj,jk1), bsp(jis,jj,jk1),  &
                                 csp(jis,jj,jk1), dsp(jis,jj,jk1))
               jk1 = jk1 + 1
            END DO

            ! integrate the pressure on the deep side
            jk1 = jk
            DO WHILE ( -zdept(jid,jj,jk1) < zuijk )
               IF( jk1 == 1 ) THEN
                  zdeps = zdept(jid,jj,1) + MIN(zuijk, ssh(jid,jj,Kmm)*znad)
                  zrhdt1 = zrhh(jid,jj,1) - interp3(zdept(jid,jj,1), asp(jid,jj,1), &
                                                    bsp(jid,jj,1)  , csp(jid,jj,1), &
                                                    dsp(jid,jj,1)) * zdeps
                  zpwed  = zpwed + 0.5_wp * (zrhh(jid,jj,1) + zrhdt1) * zdeps
                  EXIT
               ENDIF
               zdeps = MAX(zdept(jid,jj,jk1-1), -zuijk)
               zpwed = zpwed +                                        &
                  integ_spline(zdeps,             zdept(jid,jj,jk1),  &
                               asp(jid,jj,jk1-1), bsp(jid,jj,jk1-1),  &
                               csp(jid,jj,jk1-1), dsp(jid,jj,jk1-1) )
               jk1 = jk1 - 1
            END DO

            ! update the momentum trends in u direction
            zdpdx1 = zcoef0 * r1_e1u(ji,jj) * ( zhpi(ji+1,jj,jk) - zhpi(ji,jj,jk) )
            IF( .NOT.lk_linssh ) THEN
               zdpdx2 = zcoef0 * r1_e1u(ji,jj) * &
                  &    ( REAL(jis-jid, wp) * (zpwes + zpwed) + (ssh(ji+1,jj,Kmm)-ssh(ji,jj,Kmm)) )
            ELSE
               zdpdx2 = zcoef0 * r1_e1u(ji,jj) * REAL(jis-jid, wp) * (zpwes + zpwed)
            ENDIF
#if defined key_RK3
            puu(ji,jj,jk,Krhs) = ( zdpdx1 + zdpdx2 - zpgu(ji,jj) ) * umask(ji,jj,jk)   ! RK3 case: dyn_hpg always called first
#else
            puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) + (zdpdx1 + zdpdx2 - zpgu(ji,jj)) * umask(ji,jj,jk)
#endif
         ENDIF

         !!!!!     for v equation
         IF( jk <= mbkv(ji,jj) ) THEN
            IF( -zdept(ji,jj+1,jk) >= -zdept(ji,jj,jk) ) THEN
               jjs = jj + 1; jjd = jj
            ELSE
               jjs = jj    ; jjd = jj + 1
            ENDIF

            ! integrate the pressure on the shallow side
            jk1 = jk
            DO WHILE ( -zdept(ji,jjs,jk1) > zvijk )
               IF( jk1 == mbkv(ji,jj) ) THEN
                  zvijk = -zdept(ji,jjs,jk1)
                  EXIT
               ENDIF
               zdeps = MIN(zdept(ji,jjs,jk1+1), -zvijk)
               zpnss = zpnss +                                       &
                  integ_spline(zdept(ji,jjs,jk1), zdeps,             &
                               asp(ji,jjs,jk1),   bsp(ji,jjs,jk1),   &
                               csp(ji,jjs,jk1),   dsp(ji,jjs,jk1) )
              jk1 = jk1 + 1
            END DO

            ! integrate the pressure on the deep side
            jk1 = jk
            DO WHILE ( -zdept(ji,jjd,jk1) < zvijk )
               IF( jk1 == 1 ) THEN
                  zdeps = zdept(ji,jjd,1) + MIN(zvijk, ssh(ji,jjd,Kmm)*znad)
                  zrhdt1 = zrhh(ji,jjd,1) - interp3(zdept(ji,jjd,1), asp(ji,jjd,1), &
                                                    bsp(ji,jjd,1)  , csp(ji,jjd,1), &
                                                    dsp(ji,jjd,1) ) * zdeps
                  zpnsd  = zpnsd + 0.5_wp * (zrhh(ji,jjd,1) + zrhdt1) * zdeps
                  EXIT
               ENDIF
               zdeps = MAX(zdept(ji,jjd,jk1-1), -zvijk)
               zpnsd = zpnsd +                                        &
                  integ_spline(zdeps,             zdept(ji,jjd,jk1),  &
                               asp(ji,jjd,jk1-1), bsp(ji,jjd,jk1-1),  &
                               csp(ji,jjd,jk1-1), dsp(ji,jjd,jk1-1) )
               jk1 = jk1 - 1
            END DO

            ! update the momentum trends in v direction
            zdpdy1 = zcoef0 * r1_e2v(ji,jj) * ( zhpi(ji,jj+1,jk) - zhpi(ji,jj,jk) )
            IF( .NOT.lk_linssh ) THEN
               zdpdy2 = zcoef0 * r1_e2v(ji,jj) * &
                       ( REAL(jjs-jjd, wp) * (zpnss + zpnsd) + (ssh(ji,jj+1,Kmm)-ssh(ji,jj,Kmm)) )
            ELSE
               zdpdy2 = zcoef0 * r1_e2v(ji,jj) * REAL(jjs-jjd, wp) * (zpnss + zpnsd )
            ENDIF

#if defined key_RK3
            pvv(ji,jj,jk,Krhs) = ( zdpdy1 + zdpdy2 - zpgv(ji,jj) ) * vmask(ji,jj,jk)   ! RK3 case: dyn_hpg always called first
#else
            pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) + (zdpdy1 + zdpdy2 - zpgv(ji,jj)) * vmask(ji,jj,jk)
#endif
         ENDIF
         !
      END_3D
      !
   END SUBROUTINE hpg_prj


   SUBROUTINE cspline( fsp, xsp, asp, bsp, csp, dsp, polynomial_type )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE cspline  ***
      !!
      !! ** Purpose :   constrained cubic spline interpolation
      !!
      !! ** Method  :   f(x) = asp + bsp*x + csp*x^2 + dsp*x^3
      !!
      !! Reference: CJC Kruger, Constrained Cubic Spline Interpoltation
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(T2D(nn_hls),jpk), INTENT(in   ) ::   fsp, xsp           ! value and coordinate
      REAL(wp), DIMENSION(T2D(nn_hls),jpk), INTENT(  out) ::   asp, bsp, csp, dsp ! coefficients of the interpoated function
      INTEGER                             , INTENT(in   ) ::   polynomial_type    ! 1: cubic spline   ;   2: Linear
      !
      INTEGER  ::   ji, jj, jk                 ! dummy loop indices
      REAL(wp) ::   zdf1, zdf2, zddf1, zddf2, ztmp1, ztmp2, zdxtmp
      REAL(wp) ::   zdxtmp1, zdxtmp2, zalpha
      REAL(wp) ::   zdf(jpk)
      !!----------------------------------------------------------------------
      !
      IF (polynomial_type == 1) THEN     ! Constrained Cubic Spline
         DO_2D( 1, 1, 1, 1 )
            !!Fritsch&Butland's method, 1984 (preferred, but more computation)
            !    DO jk = 2, jpkm1-1
            !       zdxtmp1 = xsp(ji,jj,jk)   - xsp(ji,jj,jk-1)
            !       zdxtmp2 = xsp(ji,jj,jk+1) - xsp(ji,jj,jk)
            !       zdf1    = ( fsp(ji,jj,jk)   - fsp(ji,jj,jk-1) ) / zdxtmp1
            !       zdf2    = ( fsp(ji,jj,jk+1) - fsp(ji,jj,jk)   ) / zdxtmp2
            !
            !       zalpha = ( zdxtmp1 + 2._wp * zdxtmp2 ) / ( zdxtmp1 + zdxtmp2 ) / 3._wp
            !
            !       IF(zdf1 * zdf2 <= 0._wp) THEN
            !           zdf(jk) = 0._wp
            !       ELSE
            !         zdf(jk) = zdf1 * zdf2 / ( ( 1._wp - zalpha ) * zdf1 + zalpha * zdf2 )
            !       ENDIF
            !    END DO

            !!Simply geometric average
            DO jk = 2, jpk-2
               zdf1 = (fsp(ji,jj,jk  ) - fsp(ji,jj,jk-1)) / (xsp(ji,jj,jk  ) - xsp(ji,jj,jk-1))
               zdf2 = (fsp(ji,jj,jk+1) - fsp(ji,jj,jk  )) / (xsp(ji,jj,jk+1) - xsp(ji,jj,jk  ))

               IF(zdf1 * zdf2 <= 0._wp) THEN
                  zdf(jk) = 0._wp
               ELSE
                  zdf(jk) = 2._wp * zdf1 * zdf2 / (zdf1 + zdf2)
               ENDIF
            END DO

            zdf(1)     = 1.5_wp * ( fsp(ji,jj,2) - fsp(ji,jj,1) ) / &
                       &          ( xsp(ji,jj,2) - xsp(ji,jj,1) )           -  0.5_wp * zdf(2)
            zdf(jpkm1) = 1.5_wp * ( fsp(ji,jj,jpkm1) - fsp(ji,jj,jpkm1-1) ) / &
                       &          ( xsp(ji,jj,jpkm1) - xsp(ji,jj,jpkm1-1) ) - 0.5_wp * zdf(jpk - 2)

            DO jk = 1, jpk-2
               zdxtmp = xsp(ji,jj,jk+1) - xsp(ji,jj,jk)
               ztmp1  = (zdf(jk+1) + 2._wp * zdf(jk)) / zdxtmp
               ztmp2  =  6._wp * (fsp(ji,jj,jk+1) - fsp(ji,jj,jk)) / zdxtmp / zdxtmp
               zddf1  = -2._wp * ztmp1 + ztmp2
               ztmp1  = (2._wp * zdf(jk+1) + zdf(jk)) / zdxtmp
               zddf2  =  2._wp * ztmp1 - ztmp2

               dsp(ji,jj,jk) = (zddf2 - zddf1) / 6._wp / zdxtmp
               csp(ji,jj,jk) = ( xsp(ji,jj,jk+1) * zddf1 - xsp(ji,jj,jk)*zddf2 ) / 2._wp / zdxtmp
               bsp(ji,jj,jk) = ( fsp(ji,jj,jk+1) - fsp(ji,jj,jk) ) / zdxtmp - &
                             & csp(ji,jj,jk) * ( xsp(ji,jj,jk+1) + xsp(ji,jj,jk) ) - &
                             & dsp(ji,jj,jk) * ((xsp(ji,jj,jk+1) + xsp(ji,jj,jk))**2 - &
                             &                   xsp(ji,jj,jk+1) * xsp(ji,jj,jk))
               asp(ji,jj,jk) = fsp(ji,jj,jk) - xsp(ji,jj,jk) * (bsp(ji,jj,jk) + &
                             &                (xsp(ji,jj,jk) * (csp(ji,jj,jk) + &
                             &                 dsp(ji,jj,jk) * xsp(ji,jj,jk))))
            END DO
         END_2D

      ELSEIF ( polynomial_type == 2 ) THEN     ! Linear
         DO_3D( 1, 1, 1, 1, 1, jpk-2 )
            zdxtmp =xsp(ji,jj,jk+1) - xsp(ji,jj,jk)
            ztmp1 = fsp(ji,jj,jk+1) - fsp(ji,jj,jk)

            dsp(ji,jj,jk) = 0._wp
            csp(ji,jj,jk) = 0._wp
            bsp(ji,jj,jk) = ztmp1 / zdxtmp
            asp(ji,jj,jk) = fsp(ji,jj,jk) - bsp(ji,jj,jk) * xsp(ji,jj,jk)
         END_3D
         !
      ELSE
         CALL ctl_stop( 'invalid polynomial type in cspline' )
      ENDIF
      !
   END SUBROUTINE cspline


   FUNCTION interp1(x, xl, xr, fl, fr)  RESULT(f)
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE interp1  ***
      !!
      !! ** Purpose :   1-d linear interpolation
      !!
      !! ** Method  :   interpolation is straight forward
      !!                extrapolation is also permitted (no value limit)
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in) ::  x, xl, xr, fl, fr
      REAL(wp)             ::  f ! result of the interpolation (extrapolation)
      REAL(wp)             ::  zdeltx
      !!----------------------------------------------------------------------
      !
      zdeltx = xr - xl
      IF( abs(zdeltx) <= 10._wp * EPSILON(x) ) THEN
         f = 0.5_wp * (fl + fr)
      ELSE
         f = ( (x - xl ) * fr - ( x - xr ) * fl ) / zdeltx
      ENDIF
      !
   END FUNCTION interp1


   FUNCTION interp2( x, a, b, c, d )  RESULT(f)
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE interp1  ***
      !!
      !! ** Purpose :   1-d constrained cubic spline interpolation
      !!
      !! ** Method  :  cubic spline interpolation
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in) ::  x, a, b, c, d
      REAL(wp)             ::  f ! value from the interpolation
      !!----------------------------------------------------------------------
      !
      f = a + x* ( b + x * ( c + d * x ) )
      !
   END FUNCTION interp2


   FUNCTION interp3( x, a, b, c, d )  RESULT(f)
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE interp1  ***
      !!
      !! ** Purpose :   Calculate the first order of derivative of
      !!                a cubic spline function y=a+b*x+c*x^2+d*x^3
      !!
      !! ** Method  :   f=dy/dx=b+2*c*x+3*d*x^2
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in) ::  x, a, b, c, d
      REAL(wp)             ::  f ! value from the interpolation
      !!----------------------------------------------------------------------
      !
      f = b + x * ( 2._wp * c + 3._wp * d * x)
      !
   END FUNCTION interp3


   FUNCTION integ_spline( xl, xr, a, b, c, d )  RESULT(f)
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE interp1  ***
      !!
      !! ** Purpose :   1-d constrained cubic spline integration
      !!
      !! ** Method  :  integrate polynomial a+bx+cx^2+dx^3 from xl to xr
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in) ::  xl, xr, a, b, c, d
      REAL(wp)             ::  za1, za2, za3
      REAL(wp)             ::  f                   ! integration result
      !!----------------------------------------------------------------------
      !
      za1 = 0.5_wp * b
      za2 = c / 3.0_wp
      za3 = 0.25_wp * d
      !
      f  = xr * ( a + xr * ( za1 + xr * ( za2 + za3 * xr ) ) ) - &
         & xl * ( a + xl * ( za1 + xl * ( za2 + za3 * xl ) ) )
      !
   END FUNCTION integ_spline

   !!======================================================================
END MODULE dynhpg

MODULE trddyn
   !!======================================================================
   !!                       ***  MODULE  trddyn  ***
   !! Ocean diagnostics:  ocean dynamic trends
   !!=====================================================================
   !! History :  3.5  !  2012-02  (G. Madec) creation from trdmod: split DYN and TRA trends
   !!                                        and manage  3D trends output for U, V, and KE
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   trd_dyn       : manage the type of momentum trend diagnostics (3D I/O, domain averaged, KE)
   !!   trd_dyn_iom   : output 3D momentum and/or tracer trends using IOM
   !!   trd_dyn_init  : initialization step
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain variables
   USE phycst         ! physical constants
   USE sbc_oce        ! surface boundary condition: ocean
   USE zdf_oce        ! ocean vertical physics: variables
!!gm   USE zdfdrg         ! ocean vertical physics: bottom friction
   USE trd_oce        ! trends: ocean variables
   USE trdken         ! trends: Kinetic ENergy 
   USE trdglo         ! trends: global domain averaged
   USE trdvor         ! trends: vertical averaged vorticity 
   USE trdmxl         ! trends: mixed layer averaged 
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! lateral boundary condition 
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC trd_dyn        ! called by all dynXXX modules

   INTERFACE trd_dyn
      module procedure trd_dyn_3d, trd_dyn_2d
   END INTERFACE

   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: sutrd_hpg, svtrd_hpg
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: sutrd_pvo, svtrd_pvo
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: sutrd_tfre, svtrd_tfre
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: sutrd_bfre, svtrd_bfre
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: sutrd_tfr, svtrd_tfr
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: sutrd_bfr, svtrd_bfr
   REAL(wp), ALLOCATABLE, DIMENSION(:,:)  , SAVE :: sutrd_iceoc, svtrd_iceoc
   REAL(wp), ALLOCATABLE, DIMENSION(:,:)  , SAVE :: sutrd_tau2d, svtrd_tau2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:)  , SAVE :: sutrd_iceoc2d, svtrd_iceoc2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:)  , SAVE :: sutrd_tfr2d, svtrd_tfr2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:)  , SAVE :: sutrd_bfr2d, svtrd_bfr2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:)  , SAVE :: sutrd_atm2d, svtrd_atm2d
   REAL(wp), ALLOCATABLE, DIMENSION(:,:)  , SAVE :: sutrd_bdy2d, svtrd_bdy2d

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trd_dyn_3d( putrd, pvtrd, ktrd, kt, Kmm, Kaa )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_dyn_3d  ***
      !! 
      !! ** Purpose :   Dispatch momentum trend computation, e.g. 3D output, 
      !!              integral constraints, barotropic vorticity, kinetic enrgy, 
      !!              and/or mixed layer budget.
      !!----------------------------------------------------------------------
      REAL(dp), DIMENSION(:,:,:), INTENT(inout)           ::   putrd, pvtrd   ! U and V trends 
      INTEGER                   , INTENT(in   )           ::   ktrd           ! trend index
      INTEGER                   , INTENT(in   )           ::   kt             ! time step
      INTEGER                   , INTENT(in   ), OPTIONAL ::   Kmm            ! time level index
      INTEGER                   , INTENT(in   ), OPTIONAL ::   Kaa            ! time level index
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   zue, zve                     ! temporary 2D arrays
      INTEGER                                 ::   ji, jj, jk     ! loop indices
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   zutrd, zvtrd   ! working arrays
      !!----------------------------------------------------------------------
      !
      DO_3D( 0, 0, 0, 0, 1, jpkm1 )
         putrd(ji,jj,jk) = putrd(ji,jj,jk) * umask(ji,jj,jk)                       ! mask the trends
         pvtrd(ji,jj,jk) = pvtrd(ji,jj,jk) * vmask(ji,jj,jk)
      END_3D

      ! KE and vorticity trends require data on 1st halo points
      IF( ln_KE_trd .OR. ln_vor_trd ) THEN
         ALLOCATE( zutrd(T2D(1),jpk), zvtrd(T2D(1),jpk) )
         DO_3D( 0, 0, 0, 0, 1, jpk )
            zutrd(ji,jj,jk) = putrd(ji,jj,jk)
            zvtrd(ji,jj,jk) = pvtrd(ji,jj,jk)
         END_3D
         CALL lbc_lnk( 'trddyn', zutrd, 'U', -1.0_wp , zvtrd, 'V', -1.0_wp )
      ENDIF
      !

!!gm NB : here a lbc_lnk should probably be added

      SELECT CASE( ktrd )
      CASE( jpdyn_hpg_save ) 
         !
         ! save 3D HPG trends to possibly have barotropic part corrected later before writing out
         ALLOCATE( sutrd_hpg(jpi,jpj,jpk), svtrd_hpg(jpi,jpj,jpk) )
         sutrd_hpg(:,:,:) = putrd(:,:,:)
         svtrd_hpg(:,:,:) = pvtrd(:,:,:)

      CASE( jpdyn_pvo_save ) 
         !
         ! save 3D coriolis trends to possibly have barotropic part corrected later before writing out
         ALLOCATE( sutrd_pvo(jpi,jpj,jpk), svtrd_pvo(jpi,jpj,jpk) )
         sutrd_pvo(:,:,:) = putrd(:,:,:)
         svtrd_pvo(:,:,:) = pvtrd(:,:,:)

      CASE( jpdyn_spg ) 
         ! For explicit scheme SPG trends come here as 3D fields
         ! Add SPG trend to 3D HPG trend and also output as 2D diagnostic in own right.
         CALL trd_dyn_iom_2d( putrd(:,:,1), pvtrd(:,:,1), jpdyn_spg, kt ) 
         sutrd_hpg(:,:,:) = sutrd_hpg(:,:,:) + putrd(:,:,:)  
         svtrd_hpg(:,:,:) = svtrd_hpg(:,:,:) + pvtrd(:,:,:)  
         CALL trd_dyn_iom_3d( svtrd_hpg, svtrd_hpg, jpdyn_hpg, kt ) 
         DEALLOCATE( sutrd_hpg, svtrd_hpg )

      CASE( jpdyn_tfre )
         !
         ! Explicit top drag trend calculated in zdf_drg. Save to add to 
         ! ZDF trend later and add to 3D TFR trend. 
         IF( .NOT. ALLOCATED(sutrd_tfre) ) THEN 
            ALLOCATE( sutrd_tfre(jpi,jpj,jpk), svtrd_tfre(jpi,jpj,jpk) )
            sutrd_tfre(:,:,:) = putrd(:,:,:)
            svtrd_tfre(:,:,:) = pvtrd(:,:,:)
         ENDIF
         IF( .NOT. ALLOCATED(sutrd_tfr) ) THEN 
            ALLOCATE( sutrd_tfr(jpi,jpj,jpk), svtrd_tfr(jpi,jpj,jpk) )
            sutrd_tfr(:,:,:) = 0.0_wp
            svtrd_tfr(:,:,:) = 0.0_wp
         ENDIF
         sutrd_tfr(:,:,:) = sutrd_tfr(:,:,:) + putrd(:,:,:) 
         svtrd_tfr(:,:,:) = svtrd_tfr(:,:,:) + pvtrd(:,:,:)

      CASE( jpdyn_tfre_bt, jpdyn_tfri )
         !
         ! Add various top friction terms for baroclinic trend to saved quantity.
         ! Any depth-mean component removed later when TFR trend written out. 
         IF( .NOT. ALLOCATED(sutrd_tfr) ) THEN 
            ALLOCATE( sutrd_tfr(jpi,jpj,jpk), svtrd_tfr(jpi,jpj,jpk) )
            sutrd_tfr(:,:,:) = 0.0_wp
            svtrd_tfr(:,:,:) = 0.0_wp
         ENDIF
         sutrd_tfr(:,:,:) = sutrd_tfr(:,:,:) + putrd(:,:,:) 
         svtrd_tfr(:,:,:) = svtrd_tfr(:,:,:) + pvtrd(:,:,:)

      CASE( jpdyn_bfre )
         !
         ! Explicit bottom drag trend calculated in zdf_drg. Save to add to 
         ! ZDF trend later and add to 3D BFR trend. 
         IF( .NOT. ALLOCATED(sutrd_bfre) ) THEN 
            ALLOCATE( sutrd_bfre(jpi,jpj,jpk), svtrd_bfre(jpi,jpj,jpk) )
            sutrd_bfre(:,:,:) = putrd(:,:,:)
            svtrd_bfre(:,:,:) = pvtrd(:,:,:)
         ENDIF
         IF( .NOT. ALLOCATED(sutrd_bfr) ) THEN 
            ALLOCATE( sutrd_bfr(jpi,jpj,jpk), svtrd_bfr(jpi,jpj,jpk) )
            sutrd_bfr(:,:,:) = 0.0_wp
            svtrd_bfr(:,:,:) = 0.0_wp
         ENDIF
         sutrd_bfr(:,:,:) = sutrd_bfr(:,:,:) + putrd(:,:,:) 
         svtrd_bfr(:,:,:) = svtrd_bfr(:,:,:) + pvtrd(:,:,:)

      CASE( jpdyn_bfre_bt, jpdyn_bfri )
         !
         ! Add various bottom friction terms for baroclinic trend to saved quantity.
         ! Any depth-mean component removed later when BFR trend written out. 
         IF( .NOT. ALLOCATED(sutrd_bfr) ) THEN 
            ALLOCATE( sutrd_bfr(jpi,jpj,jpk), svtrd_bfr(jpi,jpj,jpk) )
            sutrd_bfr(:,:,:) = 0.0_wp
            svtrd_bfr(:,:,:) = 0.0_wp
         ENDIF
         sutrd_bfr(:,:,:) = sutrd_bfr(:,:,:) + putrd(:,:,:) 
         svtrd_bfr(:,:,:) = svtrd_bfr(:,:,:) + pvtrd(:,:,:)

      CASE( jpdyn_zdf ) 
         ! ZDF trend: Add explicit top/bottom friction if necessary. If ln_dynspg_ts, remove barotropic component 
         !            and add wind stress, and top and bottom friction trends from dynspg_ts.
         !
         ! If TFRE or BFRE arrays allocated at this stage then they will contain trends due
         ! to explicit top or bottom drag components which need to be added to the ZDF trend. 
         IF( ALLOCATED( sutrd_tfre ) ) THEN
            DO jk = 1, jpkm1
               putrd(:,:,jk) = ( putrd(:,:,jk) + sutrd_tfre(:,:,jk) ) * umask(:,:,jk)
               pvtrd(:,:,jk) = ( pvtrd(:,:,jk) + svtrd_tfre(:,:,jk) ) * vmask(:,:,jk)
            END DO
            DEALLOCATE( sutrd_tfre, svtrd_tfre )
         ENDIF
         IF( ALLOCATED( sutrd_bfre ) ) THEN
            DO jk = 1, jpkm1
               putrd(:,:,jk) = ( putrd(:,:,jk) + sutrd_bfre(:,:,jk) ) * umask(:,:,jk)
               pvtrd(:,:,jk) = ( pvtrd(:,:,jk) + svtrd_bfre(:,:,jk) ) * vmask(:,:,jk)
            END DO
            DEALLOCATE( sutrd_bfre, svtrd_bfre )
         ENDIF
         IF( ln_dynspg_ts ) THEN
            ALLOCATE( zue(A2D(0)), zve(A2D(0)) )
            DO_2D( 0, 0, 0, 0 )
               zue(ji,jj) = e3u(ji,jj,1,Kaa) * putrd(ji,jj,1) * umask(ji,jj,1)
               zve(ji,jj) = e3v(ji,jj,1,Kaa) * pvtrd(ji,jj,1) * vmask(ji,jj,1)
            END_2D
            DO_3D( 0, 0, 0, 0, 2, jpkm1 )
               zue(ji,jj) = zue(ji,jj) + e3u(ji,jj,jk,Kaa) * putrd(ji,jj,jk) * umask(ji,jj,jk)
               zve(ji,jj) = zve(ji,jj) + e3v(ji,jj,jk,Kaa) * pvtrd(ji,jj,jk) * vmask(ji,jj,jk)
            END_3D
            DO_3D( 0, 0, 0, 0, 1, jpkm1 )
               putrd(ji,jj,jk) = ( sutrd_tau2d(ji,jj) + sutrd_bfr2d(ji,jj) + putrd(ji,jj,jk) - &
                  &                zue(ji,jj) * r1_hu(ji,jj,Kaa) ) * umask(ji,jj,jk)
               pvtrd(ji,jj,jk) = ( svtrd_tau2d(ji,jj) + svtrd_bfr2d(ji,jj) + pvtrd(ji,jj,jk) - &
                  &                zve(ji,jj) * r1_hv(ji,jj,Kaa) ) * vmask(ji,jj,jk)
            END_3D
            DEALLOCATE( zue, zve, sutrd_tau2d, svtrd_tau2d, sutrd_bfr2d, svtrd_bfr2d)
            IF( ALLOCATED( sutrd_tfr2d ) ) THEN
               DO jk = 1, jpkm1
                  putrd(:,:,jk) = ( putrd(:,:,jk) + sutrd_tfr2d(:,:) ) * umask(:,:,jk)
                  pvtrd(:,:,jk) = ( pvtrd(:,:,jk) + svtrd_tfr2d(:,:) ) * vmask(:,:,jk)
               END DO
               DEALLOCATE( sutrd_tfr2d, svtrd_tfr2d )
            ENDIF
            !
         ENDIF
         CASE( jpdyn_bdy )
            !
            IF( ALLOCATED( sutrd_bdy2d ) ) THEN
               DO jk = 1, jpkm1
                  putrd(:,:,jk) = ( putrd(:,:,jk) + sutrd_bdy2d(:,:) ) * umask(:,:,jk)
                  pvtrd(:,:,jk) = ( pvtrd(:,:,jk) + svtrd_bdy2d(:,:) ) * vmask(:,:,jk)
               END DO
               DEALLOCATE( sutrd_bdy2d, svtrd_bdy2d )
            ENDIF
         !
      CASE( jpdyn_tot ) 
         ! Don't need to do anything special for TOT trends, but we are at the end of the timestep, so
         ! write out total top and bottom friction "trends" for the surface / bottom layers after 
         ! removing any depth-mean component.
         IF( ALLOCATED( sutrd_tfr ) .OR. ALLOCATED( sutrd_iceoc ) ) THEN
            ! With explicit top and bottom friction, the top friction diagnostic
            ! is initialised here.
            IF( .NOT. ALLOCATED( sutrd_tfr ) ) THEN
               ALLOCATE( sutrd_tfr(jpi,jpj,jpk), svtrd_tfr(jpi,jpj,jpk) )
               sutrd_tfr(:,:,:) = 0.0_wp
               svtrd_tfr(:,:,:) = 0.0_wp
            ENDIF   
            IF( ALLOCATED( sutrd_iceoc ) ) THEN
               ! Add trend due to ice-ocean stress at the surface
               sutrd_tfr(:,:,1) = sutrd_tfr(:,:,1) + sutrd_iceoc(:,:)
               svtrd_tfr(:,:,1) = svtrd_tfr(:,:,1) + svtrd_iceoc(:,:)
               DEALLOCATE( sutrd_iceoc, svtrd_iceoc )
            ENDIF
            ALLOCATE( zue(A2D(0)), zve(A2D(0)) )
            DO_2D( 0, 0, 0, 0 )
               zue(ji,jj) = e3u(ji,jj,1,Kaa) * sutrd_tfr(ji,jj,1) * umask(ji,jj,1)
               zve(ji,jj) = e3v(ji,jj,1,Kaa) * svtrd_tfr(ji,jj,1) * vmask(ji,jj,1)
            END_2D
            DO_3D( 0, 0, 0, 0, 2, jpkm1 )
               zue(ji,jj) = zue(ji,jj) + e3u(ji,jj,jk,Kaa) * sutrd_tfr(ji,jj,jk) * umask(ji,jj,jk)
               zve(ji,jj) = zve(ji,jj) + e3v(ji,jj,jk,Kaa) * svtrd_tfr(ji,jj,jk) * vmask(ji,jj,jk)
            END_3D
            DO_3D( 0, 0, 0, 0, 1, jpkm1 )
               sutrd_tfr(ji,jj,jk) = ( sutrd_tfr(ji,jj,jk) - zue(ji,jj) * r1_hu(ji,jj,Kaa) ) * umask(ji,jj,jk)
               svtrd_tfr(ji,jj,jk) = ( svtrd_tfr(ji,jj,jk) - zve(ji,jj) * r1_hv(ji,jj,Kaa) ) * vmask(ji,jj,jk)
            END_3D
            CALL trd_dyn_iom_3d( sutrd_tfr, svtrd_tfr, jpdyn_tfr, kt )
            DEALLOCATE( zue, zve, sutrd_tfr, svtrd_tfr )
         ENDIF
         IF( ALLOCATED( sutrd_bfr ) ) THEN
            ALLOCATE( zue(A2D(0)), zve(A2D(0)) )
            DO_2D( 0, 0, 0, 0 )
               zue(ji,jj) = e3u(ji,jj,1,Kaa) * sutrd_bfr(ji,jj,1) * umask(ji,jj,1)
               zve(ji,jj) = e3v(ji,jj,1,Kaa) * svtrd_bfr(ji,jj,1) * vmask(ji,jj,1)
            END_2D
            DO_3D( 0, 0, 0, 0, 2, jpkm1 )
               zue(ji,jj) = zue(ji,jj) + e3u(ji,jj,jk,Kaa) * sutrd_bfr(ji,jj,jk) * umask(ji,jj,jk)
               zve(ji,jj) = zve(ji,jj) + e3v(ji,jj,jk,Kaa) * svtrd_bfr(ji,jj,jk) * vmask(ji,jj,jk)
            END_3D
            DO_3D( 0, 0, 0, 0, 1, jpkm1 )
               sutrd_bfr(ji,jj,jk) = ( sutrd_bfr(ji,jj,jk) - zue(ji,jj) * r1_hu(ji,jj,Kaa) ) * umask(ji,jj,jk)
               svtrd_bfr(ji,jj,jk) = ( svtrd_bfr(ji,jj,jk) - zve(ji,jj) * r1_hv(ji,jj,Kaa) ) * vmask(ji,jj,jk)
            END_3D
            CALL trd_dyn_iom_3d( sutrd_bfr, svtrd_bfr, jpdyn_bfr, kt )
            DEALLOCATE( zue, zve, sutrd_bfr, svtrd_bfr )
         ENDIF

      END SELECT

      IF ( ktrd <= jptot_dyn ) THEN  ! output of 3D trends and use for other diagnostics
         !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         !   3D output of momentum and/or tracers trends using IOM interface
         !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         IF( ln_dyn_trd )   CALL trd_dyn_iom( putrd, pvtrd, ktrd, kt, Kmm )
         
         !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         !  Integral Constraints Properties for momentum and/or tracers trends
         !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         IF( ln_glo_trd )   CALL trd_glo( putrd, pvtrd, ktrd, 'DYN', kt, Kmm )

         !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         !  Kinetic Energy trends
         !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         IF( ln_KE_trd  )   CALL trd_ken( zutrd, zvtrd, ktrd, kt, Kmm )

         !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         !  Vorticity trends
         !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         IF( ln_vor_trd )   CALL trd_vor( zutrd, zvtrd, ktrd, kt, Kmm )

         !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         !  Mixed layer trends for active tracers
         !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         !!gm      IF( ln_dyn_mxl )   CALL trd_mxl_dyn   
         !
      ENDIF
      IF( ln_KE_trd .OR. ln_vor_trd ) DEALLOCATE( zutrd, zvtrd )
      !
   END SUBROUTINE trd_dyn_3d


   SUBROUTINE trd_dyn_2d( putrd, pvtrd, ktrd, kt, Kmm )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_dyn_2d  ***
      !! 
      !! ** Purpose :   Dispatch momentum trend computation, e.g. 2D output, 
      !!              integral constraints, barotropic vorticity, kinetic enrgy, 
      !!              and/or mixed layer budget.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(inout)           ::   putrd, pvtrd   ! U and V trends 
      INTEGER                 , INTENT(in   )           ::   ktrd           ! trend index
      INTEGER                 , INTENT(in   )           ::   kt             ! time step
      INTEGER                 , INTENT(in   ), OPTIONAL ::   Kmm            ! time level index
      INTEGER                                 ::   jk
      !!----------------------------------------------------------------------
      !
      putrd(:,:) = putrd(:,:) * umask(:,:,1)                       ! mask the trends
      pvtrd(:,:) = pvtrd(:,:) * vmask(:,:,1)
      !

!!gm NB : here a lbc_lnk should probably be added

      SELECT CASE(ktrd)

      CASE( jpdyn_pvo_corr )
         !
         ! Remove "first-guess" barotropic coriolis trend from 3D PVO trend. 
         DO jk = 1, jpkm1
            sutrd_pvo(:,:,jk) = sutrd_pvo(:,:,jk) - putrd(:,:)
            svtrd_pvo(:,:,jk) = svtrd_pvo(:,:,jk) - pvtrd(:,:)
         ENDDO

      CASE( jpdyn_spg )
          !
          ! For split-explicit scheme SPG trends come here as 2D fields
          ! Add SPG trend to 3D HPG trend and also output as 2D diagnostic in own right.
          IF( ALLOCATED(sutrd_atm2d) ) THEN
             DO jk = 1, jpkm1
                sutrd_hpg(:,:,jk) = sutrd_hpg(:,:,jk) + sutrd_atm2d(:,:) + putrd(:,:)
                svtrd_hpg(:,:,jk) = svtrd_hpg(:,:,jk) + svtrd_atm2d(:,:) + pvtrd(:,:)
             ENDDO
             DEALLOCATE( sutrd_atm2d, svtrd_atm2d )
          ELSE
             DO jk = 1, jpkm1
                sutrd_hpg(:,:,jk) = sutrd_hpg(:,:,jk) + putrd(:,:)
                svtrd_hpg(:,:,jk) = svtrd_hpg(:,:,jk) + pvtrd(:,:)
             ENDDO
          ENDIF
          CALL trd_dyn_3d( sutrd_hpg, svtrd_hpg, jpdyn_hpg, kt, Kmm )
          DEALLOCATE( sutrd_hpg, svtrd_hpg )

      CASE( jpdyn_pvo )
          !
          ! Add 2D PVO trend to 3D PVO trend and also output as diagnostic in own right.
          DO jk = 1, jpkm1
             sutrd_pvo(:,:,jk) = sutrd_pvo(:,:,jk) + putrd(:,:)
             svtrd_pvo(:,:,jk) = svtrd_pvo(:,:,jk) + pvtrd(:,:)
          ENDDO
          CALL trd_dyn_3d( sutrd_pvo, svtrd_pvo, jpdyn_pvo, kt, Kmm )
          DEALLOCATE( sutrd_pvo, svtrd_pvo )

      CASE( jpdyn_iceoc )
          !
          ! Save surface ice-ocean stress trend locally to be subtracted from
          ! surface wind stress trend and added to 3D top friction trend. 
          IF( .NOT. ALLOCATED(sutrd_iceoc) ) ALLOCATE( sutrd_iceoc(jpi,jpj), svtrd_iceoc(jpi,jpj) )
          sutrd_iceoc(:,:) = putrd(:,:)
          svtrd_iceoc(:,:) = pvtrd(:,:)

      CASE( jpdyn_tau )
          !
          ! Subtract ice-ocean stress from surface wind forcing
          IF( ALLOCATED(sutrd_iceoc) ) THEN
             putrd(:,:) = putrd(:,:) - sutrd_iceoc(:,:) 
             pvtrd(:,:) = pvtrd(:,:) - svtrd_iceoc(:,:) 
          ENDIF

      CASE( jpdyn_iceoc2d )
          !
          ! Save 2D ice-ocean stress trend locally as the first installment of top friction.
          ! Subtracted from 2D wind stress trend later. 
          IF( .NOT. ALLOCATED(sutrd_tfr2d) ) ALLOCATE( sutrd_tfr2d(jpi,jpj), svtrd_tfr2d(jpi,jpj) )
          sutrd_tfr2d(:,:) = putrd(:,:)
          svtrd_tfr2d(:,:) = pvtrd(:,:)

      CASE( jpdyn_atm )
          !
          ! Save 2D atmospheric pressure effect on pressure
          ALLOCATE( sutrd_atm2d(jpi,jpj), svtrd_atm2d(jpi,jpj) )
          sutrd_atm2d(:,:) = putrd(:,:)
          svtrd_atm2d(:,:) = pvtrd(:,:)

      CASE( jpdyn_bdy2d )
          !
          ! Save 2D lateral boundary forcing
          ALLOCATE( sutrd_bdy2d(jpi,jpj), svtrd_bdy2d(jpi,jpj) )
          sutrd_bdy2d(:,:) = putrd(:,:)
          svtrd_bdy2d(:,:) = pvtrd(:,:)

      CASE( jpdyn_tau2d )
          !
          ! Subtract ice-ocean stress from depth-mean trend due to wind forcing
          ! and save to be added to ZDF trend later. Output as a trend in its own right (below).
          ! Note at this stage, sutrd_tfr2d should only contain the contribution to top friction
          ! from (partial) ice-ocean stress.
          ALLOCATE( sutrd_tau2d(jpi,jpj), svtrd_tau2d(jpi,jpj) )
          IF( ALLOCATED(sutrd_tfr2d) ) THEN
             putrd(:,:) = putrd(:,:) - sutrd_tfr2d(:,:) 
             pvtrd(:,:) = pvtrd(:,:) - svtrd_tfr2d(:,:) 
          ENDIF
          sutrd_tau2d(:,:) = putrd(:,:)
          svtrd_tau2d(:,:) = pvtrd(:,:)

      CASE( jpdyn_tfr )
          !
          ! Add ice-ocean stress from depth-mean trend due to top friction
          ! and save to be added to ZDF trend later. Output as a trend in its own right (below).
          IF( .NOT. ALLOCATED(sutrd_tfr2d) ) THEN
             ALLOCATE( sutrd_tfr2d(jpi,jpj), svtrd_tfr2d(jpi,jpj) )
             sutrd_tfr2d(:,:) = 0._wp ; svtrd_tfr2d(:,:) = 0._wp 
          ENDIF
          sutrd_tfr2d(:,:) = sutrd_tfr2d(:,:) + putrd(:,:)
          svtrd_tfr2d(:,:) = svtrd_tfr2d(:,:) + pvtrd(:,:)
          ! update (putrd,pvtrd) so that total tfr2d trend is output by call to trd_dyn_iom_2d
          putrd(:,:) = sutrd_tfr2d(:,:)
          pvtrd(:,:) = svtrd_tfr2d(:,:)

      CASE( jpdyn_bfr )
          !
          !  Save 2D field to add to ZDF trend  and also output 2D field as diagnostic in own right (below).
          ALLOCATE( sutrd_bfr2d(jpi,jpj), svtrd_bfr2d(jpi,jpj) )
          sutrd_bfr2d(:,:) = putrd(:,:)
          svtrd_bfr2d(:,:) = pvtrd(:,:)

      END SELECT

      IF( ktrd <= jptot_dyn ) THEN ! output of 2D trends and use for other diagnostics

         !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         !   2D output of momentum and/or tracers trends using IOM interface
         !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         IF( ln_dyn_trd )   CALL trd_dyn_iom_2d( putrd, pvtrd, ktrd, kt )

!!$   CALLS TO THESE ROUTINES FOR 2D DIAGOSTICS NOT CODED YET
!!$         !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!$         !  Integral Constraints Properties for momentum and/or tracers trends
!!$         !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!$         IF( ln_glo_trd )   CALL trd_glo( putrd, pvtrd, ktrd, 'DYN', kt )
!!$
!!$         !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!$         !  Kinetic Energy trends
!!$         !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!!$         IF( ln_KE_trd  )   CALL trd_ken( putrd, pvtrd, ktrd, kt )
!!$
!!$         !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!$         !  Vorticity trends
!!$         !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!!$         IF( ln_vor_trd )   CALL trd_vor( putrd, pvtrd, ktrd, kt )
!!$
!!$         !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!!$         !  Mixed layer trends for active tracers
!!$         !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!$         IF( ln_dyn_mxl )   CALL trd_mxl_dyn   

      ENDIF
      !

   END SUBROUTINE trd_dyn_2d


   SUBROUTINE trd_dyn_iom_3d( putrd, pvtrd, ktrd, kt, Kmm )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_dyn_iom_3d  ***
      !! 
      !! ** Purpose :   output 3D trends using IOM
      !!----------------------------------------------------------------------
      REAL(dp), DIMENSION(:,:,:), INTENT(inout)           ::   putrd, pvtrd   ! U and V trends
      INTEGER                   , INTENT(in   )           ::   ktrd           ! trend index
      INTEGER                   , INTENT(in   )           ::   kt             ! time step
      INTEGER                   , INTENT(in   ), OPTIONAL ::   Kmm            ! time level index
      !
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   z3dx, z3dy   ! 3D workspace 
      !!----------------------------------------------------------------------
      !
      SELECT CASE( ktrd )
      CASE( jpdyn_hpg )   ;   CALL iom_put( "utrd_hpg", putrd )    ! hydrostatic pressure gradient
                              CALL iom_put( "vtrd_hpg", pvtrd )
      CASE( jpdyn_pvo )   ;   CALL iom_put( "utrd_pvo", putrd )    ! planetary vorticity
                              CALL iom_put( "vtrd_pvo", pvtrd )
      CASE( jpdyn_rvo )   ;   CALL iom_put( "utrd_rvo", putrd )    ! relative  vorticity     (or metric term)
                              CALL iom_put( "vtrd_rvo", pvtrd )
      CASE( jpdyn_keg )   ;   CALL iom_put( "utrd_keg", putrd )    ! Kinetic Energy gradient (or had)
                              CALL iom_put( "vtrd_keg", pvtrd )
                              ALLOCATE( z3dx(T2D(0),jpk) , z3dy(T2D(0),jpk) )   ! U.dxU & V.dyV (approximation)
                              z3dx(:,:,jpk) = 0._wp
                              z3dy(:,:,jpk) = 0._wp
                              DO_3D( 0, 0, 0, 0, 1, jpkm1 )                          ! no mask as un,vn are masked
                                 z3dx(ji,jj,jk) = uu(ji,jj,jk,Kmm) * ( uu(ji+1,jj,jk,Kmm) - uu(ji-1,jj,jk,Kmm) ) / ( 2._wp * e1u(ji,jj) )
                                 z3dy(ji,jj,jk) = vv(ji,jj,jk,Kmm) * ( vv(ji,jj+1,jk,Kmm) - vv(ji,jj-1,jk,Kmm) ) / ( 2._wp * e2v(ji,jj) )
                              END_3D
                              CALL iom_put( "utrd_udx", z3dx  )
                              CALL iom_put( "vtrd_vdy", z3dy  )
                              DEALLOCATE( z3dx , z3dy )
      CASE( jpdyn_zad )   ;   CALL iom_put( "utrd_zad", putrd )    ! vertical advection
                              CALL iom_put( "vtrd_zad", pvtrd )
      CASE( jpdyn_bdy )   ;   CALL iom_put( "utrd_bdy", putrd )    ! lateral boundary forcing trend
                              CALL iom_put( "vtrd_bdy", pvtrd )
      CASE( jpdyn_ldf )   ;   CALL iom_put( "utrd_ldf", putrd )    ! lateral  diffusion
                              CALL iom_put( "vtrd_ldf", pvtrd )
      CASE( jpdyn_zdf )   ;   CALL iom_put( "utrd_zdf", putrd )    ! vertical diffusion 
                              CALL iom_put( "vtrd_zdf", pvtrd )
      CASE( jpdyn_bfr )   ;   CALL iom_put( "utrd_bfr", putrd )    ! bottom friction for bottom layer
                              CALL iom_put( "vtrd_bfr", pvtrd )
      CASE( jpdyn_tfr )   ;   CALL iom_put( "utrd_tfr", putrd )    ! total top friction for top layer
                              CALL iom_put( "vtrd_tfr", pvtrd )
      CASE( jpdyn_tot )   ;   CALL iom_put( "utrd_tot", putrd )    ! total trends excluding asselin filter
                              CALL iom_put( "vtrd_tot", pvtrd )
      CASE( jpdyn_atf )   ;   CALL iom_put( "utrd_atf", putrd )        ! asselin filter trends 
                              CALL iom_put( "vtrd_atf", pvtrd )
      END SELECT
      !
   END SUBROUTINE trd_dyn_iom_3d


   SUBROUTINE trd_dyn_iom_2d( putrd, pvtrd, ktrd, kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_dyn_iom_2d  ***
      !! 
      !! ** Purpose :   output 2D trends using IOM
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   putrd, pvtrd   ! U and V trends
      INTEGER                 , INTENT(in   ) ::   ktrd           ! trend index
      INTEGER                 , INTENT(in   ) ::   kt             ! time step
      !
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      SELECT CASE( ktrd )
      CASE( jpdyn_spg )      ;   CALL iom_put( "utrd_spg2d", putrd )      ! surface pressure gradient
                                 CALL iom_put( "vtrd_spg2d", pvtrd )
      CASE( jpdyn_pvo )      ;   CALL iom_put( "utrd_pvo2d", putrd )      ! planetary vorticity (barotropic part)
                                 CALL iom_put( "vtrd_pvo2d", pvtrd )
      CASE( jpdyn_atm )      ;   CALL iom_put( "utrd_atm2d", putrd )      ! atmospheric pressure gradient trend
                                 CALL iom_put( "vtrd_atm2d", pvtrd )
      CASE( jpdyn_bdy2d )    ;   CALL iom_put( "utrd_bdy2d", putrd )      ! lateral boundary forcing trend
                                 CALL iom_put( "vtrd_bdy2d", pvtrd )
      CASE( jpdyn_frc2d )    ;   CALL iom_put( "utrd_frc2d", putrd )      ! constant forcing term from barotropic calcn.
                                 CALL iom_put( "vtrd_frc2d", pvtrd ) 
      CASE( jpdyn_tau )      ;   CALL iom_put( "utrd_tau", putrd )        ! surface wind stress trend
                                 CALL iom_put( "vtrd_tau", pvtrd )
      CASE( jpdyn_tau2d )    ;   CALL iom_put( "utrd_tau2d", putrd )      ! wind stress depth-mean trend
                                 CALL iom_put( "vtrd_tau2d", pvtrd )
      CASE( jpdyn_bfr )      ;   CALL iom_put( "utrd_bfr2d", putrd )      ! bottom friction depth-mean trend
                                 CALL iom_put( "vtrd_bfr2d", pvtrd )
      CASE( jpdyn_tfr )      ;   CALL iom_put( "utrd_tfr2d", putrd )      ! top friction depth-mean trend
                                 CALL iom_put( "vtrd_tfr2d", pvtrd )
      CASE( jpdyn_tot )      ;   CALL iom_put( "utrd_tot2d", putrd )      ! total 2D trend, excluding time filter
                                 CALL iom_put( "vtrd_tot2d", pvtrd )
      END SELECT
      !
   END SUBROUTINE trd_dyn_iom_2d

   !!======================================================================
END MODULE trddyn

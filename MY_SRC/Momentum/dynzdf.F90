MODULE dynzdf
   !!==============================================================================
   !!                 ***  MODULE  dynzdf  ***
   !! Ocean dynamics :  vertical component of the momentum mixing trend
   !!                   (optionally, Courant-number-dependent Implicit Vertical Advection
   !!                    and vertical component of rotated lateral mixing)
   !!==============================================================================
   !! History :  1.0  !  2005-11  (G. Madec)  Original code
   !!            3.3  !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!            4.0  !  2017-06  (G. Madec) remove the explicit time-stepping option + avm at t-point
   !!            4.5  !  2022-06  (S. Techene, G, Madec) refactorization to reduce local memory usage
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_zdf       : compute the after velocity through implicit calculation of vertical mixing
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE phycst         ! physical constants
   USE dom_oce        ! ocean space and time domain variables 
   USE sbc_oce        ! surface boundary condition: ocean
   USE zdf_oce        ! ocean vertical physics variables
   USE zdfdrg         ! vertical physics: top/bottom drag coef.
   USE dynadv    ,ONLY: ln_dynadv_vec    ! dynamics: advection form
   USE dynldf_iso,ONLY: akzu, akzv       ! dynamics: vertical component of rotated lateral mixing 
   USE ldfdyn         ! lateral diffusion: eddy viscosity coef. and type of operator
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   !
   USE in_out_manager ! I/O manager
   USE iom            ! IOM library
   USE lib_mpp        ! MPP library
   USE prtctl         ! Print control
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_zdf   !  routine called by step.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
   
   SUBROUTINE dyn_zdf( kt, Kbb, Kmm, Krhs, puu, pvv, Kaa )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_zdf  ***
      !!
      !! ** Purpose :   compute the trend due to the vert. momentum diffusion
      !!              (optionally, Courant-number-dependent Implicit Vertical 
      !!              Advection and vertical component of rotated lateral mixing)
      !!              together with the time stepping using an implicit scheme.
      !!
      !! ** Method  :  - Leap-Frog time stepping on all trends but the vertical mixing
      !!         u(after) =         u(before) + 2*dt *       u(rhs)                vector form or linear free surf.
      !!         u(after) = ( e3u_b*u(before) + 2*dt * e3u_n*u(rhs) ) / e3u_after   otherwise
      !!               - update the after velocity with the implicit vertical mixing.
      !!      This requires to solver the following system: 
      !!         u(after) = u(after) + 1/e3u_after  dk+1[ mi(avm) / e3uw_after dk[ua] ]
      !!      with the following surface/top/bottom boundary condition:
      !!      surface: wind stress input (averaged over kt-1/2 & kt+1/2)
      !!      top & bottom : top stress (iceshelf-ocean) & bottom stress (cf zdfdrg.F90)
      !!
      !! ** Action :   (puu(:,:,:,Kaa),pvv(:,:,:,Kaa))   after velocity 
      !!---------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::  kt, Kbb, Kmm, Krhs, Kaa ! ocean time-step and -level indices
      REAL(wp), DIMENSION(:,:,:,:), INTENT(inout) ::  puu, pvv                ! ocean velocities and their RHS
      !
      INTEGER  ::   ji, jj, jk           ! dummy loop indices
      INTEGER  ::   iku, ikv             ! local integers
      REAL(wp) ::   zzwi, ze3ua          ! local scalars
      REAL(wp) ::   zDt_2, zDt_4         !   -      -
      REAL(wp) ::   zDt_2bua, zDt_2bva   !   -      -
      REAL(wp) ::   zzWp_e3, zzWn_e3     !   -      -
      REAL(wp) ::   zzws, ze3va          !   -      -
      REAL(wp) ::   z1_e3ua, z1_e3va     !   -      -
      REAL(wp) ::   zWu , zWv            !   -      -
      REAL(wp) ::   zWui, zWvi           !   -      -
      REAL(wp) ::   zWus, zWvs           !   -      -
      REAL(wp), DIMENSION(T1Di(0),jpkm1)      ::   zwi, zwd, zws  ! 2D workspace
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   ztrdu, ztrdv   !  -      -
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   ztrdu_fr, ztrdv_fr   !  -      -
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dyn_zdf')
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN       !* initialization
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn_zdf_imp : vertical momentum diffusion implicit operator'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
         ENDIF
      ENDIF
      !
      zDt_2 = rDt * 0.5_wp
      zDt_4 = rDt * 0.25_wp
      !
      !                             !* explicit top/bottom drag case
      IF( .NOT.ln_drgimp )   CALL zdf_drg_exp( kt, Kmm, puu(:,:,:,Kbb), pvv(:,:,:,Kbb), puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )  ! add top/bottom friction trend to (puu(Kaa),pvv(Kaa))
      !
      !
      !                                               ! ================= !
      DO_1Dj( 0, 0 )                                  !  i-k slices loop  !
         !                                            ! ================= !
         !
         !
         !              !--------------------------------------------------------------!
         !              !==  RHS : time-stepping of all trends but the implicit one  ==!   ( put in (uu,vv)(Kaa) )
         !              !--------------------------------------------------------------!
         !
         IF( ln_dynadv_vec .OR. lk_linssh ) THEN   ! applied on velocity
            DO_2Dik( 0, 0,    1, jpkm1, 1 )
               puu(ji,jj,jk,Kaa) = ( puu(ji,jj,jk,Kbb) + rDt * puu(ji,jj,jk,Krhs) ) * umask(ji,jj,jk)
               pvv(ji,jj,jk,Kaa) = ( pvv(ji,jj,jk,Kbb) + rDt * pvv(ji,jj,jk,Krhs) ) * vmask(ji,jj,jk)
            END_2D
         ELSE                                      ! applied on thickness weighted velocity
            DO_2Dik( 0, 0,    1, jpkm1, 1 )
#if defined key_qco
               puu(ji,jj,jk,Kaa) = (        ( 1._wp + r3u(ji,jj,Kbb) ) * puu(ji,jj,jk,Kbb )  &
                  &                 + rDt * ( 1._wp + r3u(ji,jj,Kmm) ) * puu(ji,jj,jk,Krhs)  )   &
                  &              /          ( 1._wp + r3u(ji,jj,Kaa) ) * umask(ji,jj,jk)
               pvv(ji,jj,jk,Kaa) = (        ( 1._wp + r3v(ji,jj,Kbb) ) * pvv(ji,jj,jk,Kbb )  &
                  &                 + rDt * ( 1._wp + r3v(ji,jj,Kmm) ) * pvv(ji,jj,jk,Krhs)  )   &
                  &              /          ( 1._wp + r3v(ji,jj,Kaa) ) * vmask(ji,jj,jk)
#else
               puu(ji,jj,jk,Kaa) = (        e3u(ji,jj,jk,Kbb) * puu(ji,jj,jk,Kbb )  &
                  &                 + rDt * e3u(ji,jj,jk,Kmm) * puu(ji,jj,jk,Krhs)  )   &
                  &              /          e3u(ji,jj,jk,Kaa) * umask(ji,jj,jk)
               pvv(ji,jj,jk,Kaa) = (        e3v(ji,jj,jk,Kbb) * pvv(ji,jj,jk,Kbb )  &
                  &                 + rDt * e3v(ji,jj,jk,Kmm) * pvv(ji,jj,jk,Krhs)  )   &
                  &              /          e3v(ji,jj,jk,Kaa) * vmask(ji,jj,jk)
#endif
            END_2D
         ENDIF
         !
         IF( l_trddyn )   THEN         !* temporary save of after velocities
            ALLOCATE( ztrdu(jpi,jpj,jpk), ztrdv(jpi,jpj,jpk), ztrdu_fr(jpi,jpj,jpk), ztrdv_fr(jpi,jpj,jpk) ) 
            ztrdu(:,:,:) = puu(:,:,:,Kaa)
            ztrdv(:,:,:) = pvv(:,:,:,Kaa)
         ENDIF
         !
         !                    ! add top/bottom friction
         !     With split-explicit free surface, barotropic stress is treated explicitly Update velocities at the bottom.
         !     J. Chanut: The bottom stress is computed considering after barotropic velocities, which does
         !                not lead to the effective stress seen over the whole barotropic loop.
         !     G. Madec : in linear free surface, e3u(:,:,:,Kaa) = e3u_0, so systematic use of e3u(:,:,:,Kaa)
         IF( ln_drgimp .AND. ln_dynspg_ts ) THEN
            DO_2Dik( 0, 0,     1, jpkm1, 1 )      ! remove barotropic velocities
               puu(ji,jj,jk,Kaa) = ( puu(ji,jj,jk,Kaa) - uu_b(ji,jj,Kaa) ) * umask(ji,jj,jk)
               pvv(ji,jj,jk,Kaa) = ( pvv(ji,jj,jk,Kaa) - vv_b(ji,jj,Kaa) ) * vmask(ji,jj,jk)
            END_2D
            IF( l_trddyn )   THEN         !* temporary save of ta and sa trends
               ztrdu_fr(:,:,:) = puu(:,:,:,Kaa)
               ztrdv_fr(:,:,:) = pvv(:,:,:,Kaa)
            ENDIF
            DO_1Di( 0, 0 )      ! Add bottom/top stress due to barotropic component only
               iku = mbku(ji,jj)         ! ocean bottom level at u- and v-points
               ikv = mbkv(ji,jj)         ! (deepest ocean u- and v-points)
               puu(ji,jj,iku,Kaa) = puu(ji,jj,iku,Kaa) + zDt_2 * ( rCdU_bot(ji+1,jj)+rCdU_bot(ji,jj) ) * uu_b(ji,jj,Kaa)   &
                  &                                            / e3u(ji,jj,iku,Kaa)
               pvv(ji,jj,ikv,Kaa) = pvv(ji,jj,ikv,Kaa) + zDt_2 * ( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj) ) * vv_b(ji,jj,Kaa)   &
                  &                                            / e3v(ji,jj,ikv,Kaa)
            END_1D
            IF( l_trddyn )   THEN    ! save baroclinic bottom friction trends due to barotropic currents  
               ztrdu_fr(:,:,:) = ( puu(:,:,:,Kaa) - ztrdu_fr(:,:,:) ) * r1_Dt 
               ztrdv_fr(:,:,:) = ( pvv(:,:,:,Kaa) - ztrdv_fr(:,:,:) ) * r1_Dt 
               CALL trd_dyn( ztrdu_fr, ztrdv_fr, jpdyn_bfre_bt, kt, Kmm )
            ENDIF
            IF( ln_isfcav .OR. ln_drgice_imp ) THEN    ! Ocean cavities (ISF)
               IF( l_trddyn )   THEN         !* temporary save of ta and sa trends
                  ztrdu_fr(:,:,:) = puu(:,:,:,Kaa)
                  ztrdv_fr(:,:,:) = pvv(:,:,:,Kaa)
               ENDIF
               DO_1Di( 0, 0 )
                  iku = miku(ji,jj)         ! top ocean level at u- and v-points
                  ikv = mikv(ji,jj)         ! (first wet ocean u- and v-points)
                  puu(ji,jj,iku,Kaa) = puu(ji,jj,iku,Kaa) + zDt_2 * ( rCdU_top(ji+1,jj)+rCdU_top(ji,jj) ) * uu_b(ji,jj,Kaa)   &
                     &                                            / e3u(ji,jj,iku,Kaa)
                  pvv(ji,jj,ikv,Kaa) = pvv(ji,jj,ikv,Kaa) + zDt_2 * ( rCdU_top(ji,jj+1)+rCdU_top(ji,jj) ) * vv_b(ji,jj,Kaa)   &
                     &                                            / e3v(ji,jj,ikv,Kaa)
               END_1D
               IF( l_trddyn )   THEN   ! save baroclinic top friction trends due to barotropic currents
                  ztrdu_fr(:,:,:) = ( puu(:,:,:,Kaa) - ztrdu_fr(:,:,:) ) * r1_Dt
                  ztrdv_fr(:,:,:) = ( pvv(:,:,:,Kaa) - ztrdv_fr(:,:,:) ) * r1_Dt 
                  CALL trd_dyn( ztrdu_fr, ztrdv_fr, jpdyn_tfre_bt, kt, Kmm )
               ENDIF
            ENDIF
         ENDIF
         !
         !              !-------------------------------!
         !              !==  Vertical diffusion on u  ==!
         !              !-------------------------------!
         !
         !                    !**  tridiagonal matrix construction  **!    diagonal (zwd), lower (zwi), upper (zws)
         !
         !                                   !*  ZDF contribution  *! 
         !
         DO_2Dik( 0, 0,   2, jpkm1, 1 )            ! inner values
            zzwi = - zDt_2 * ( avm(ji+1,jj,jk  )     +  avm(ji,jj,jk  )     ) &
               &           / ( e3u(ji  ,jj,jk  ,Kaa) * e3uw(ji,jj,jk  ,Kmm) ) * wumask(ji,jj,jk  )
            zzws = - zDt_2 * ( avm(ji+1,jj,jk+1)     +  avm(ji,jj,jk+1)     ) &
               &           / ( e3u(ji  ,jj,jk  ,Kaa) * e3uw(ji,jj,jk+1,Kmm) ) * wumask(ji,jj,jk+1)
            zwi(ji,jk) = zzwi
            zws(ji,jk) = zzws
            zwd(ji,jk) = 1._wp - zzwi - zzws
         END_2D
         DO_1Di( 0, 0 )                            ! surface boundary conditions
            zzws = - zDt_2 * ( avm(ji+1,jj,2)     +  avm(ji,jj,2)     ) &
               &           / ( e3u(ji  ,jj,1,Kaa) * e3uw(ji,jj,2,Kmm) ) * wumask(ji,jj,2)
            zwi(ji,1) = 0._wp
            zws(ji,1) = zzws
            zwd(ji,1) = 1._wp - zzws
         END_1D

! !! gm  NOUVEAU CAS implicite :
!        
!        
!          IF( ln_zad_Aimp ) THEN               !*  add Courant-number-dependent Implicit Vertical Advection  *!
!             !
!             !                                !  Flux form  !   dk[ WI_uw * upstream u ] 
!             DO_2Dik( 0, 0,   2, jpkm1, 1 )
!                zDt_2bua =  zDt_2 * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kaa)
!                zWui = ( e1e2t(ji,jj)*wi(ji,jj,jk  ) + e1e2t(ji+1,jj)*wi(ji+1,jj,jk  ) ) * zDt_2bua 
!                zWus = ( e1e2t(ji,jj)*wi(ji,jj,jk+1) + e1e2t(ji+1,jj)*wi(ji+1,jj,jk+1) ) * zDt_2bua
!                zwi(ji,jk) = zwi(ji,jk) + MIN( zWui, 0._wp )
!                zws(ji,jk) = zws(ji,jk)                      - MAX( zWus, 0._wp )
!                zwd(ji,jk) = zwd(ji,jk) + MAX( zWui, 0._wp ) - MIN( zWus, 0._wp )
!             END_2D
!             !
!             DO_1Di( 0, 0 )          ! surface boundary conditions   ( wi(jk=1) always 0 )
!                zDt_2bua = zDt_2 * r1_e1e2u(ji,jj) / e3u(ji,jj,1,Kaa)
!                zWus = ( e1e2t(ji,jj)*wi(ji,jj,2) + e1e2t(ji+1,jj)*wi(ji+1,jj,2) ) * zDt_2bua
!                zwi(ji,1) = 0._wp
!                zws(ji,1) = zws(ji,1)                        - MAX( zWus, 0._wp )
!                zwd(ji,1) = zwd(ji,1)                        - MIN( zWus, 0._wp )
!             END_1D
!             IF( ln_dynadv_vec ) THEN            !  Vector Invariant Form  !   remove 
!                DO_2Dik( 0, 0,   1, jpkm1, 1 )
!                   zDt_2bua =  zDt_2 * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kaa)
!                   zwd(ji,jk) = zwd(ji,jk) - zDt_2bua * (   e1e2t(ji  ,jj)*( wi(ji  ,jj,jk) - wi(ji  ,jj,jk+1) )   &
!                      &                                   + e1e2t(ji+1,jj)*( wi(ji+1,jj,jk) - wi(ji+1,jj,jk+1) )   )
!                END_2D
!             ENDIF
!             !
!          ENDIF
! 
! !!gm end

         !
         IF( ln_zad_Aimp ) THEN              !*  add Courant-number-dependent Implicit Vertical Advection  *!
            !
            IF( ln_dynadv_vec ) THEN            !  Vector Invariant Form  !   Wu * upstream dk[U]
               DO_2Dik( 0, 0,   2, jpkm1, 1 )
                  zWu  = zDt_4 * (  ( e1e2t(ji,jj)*wi(ji,jj,jk  ) + e1e2t(ji+1,jj)*wi(ji+1,jj,jk  ) )  &
                     &            + ( e1e2t(ji,jj)*wi(ji,jj,jk+1) + e1e2t(ji+1,jj)*wi(ji+1,jj,jk+1) )  ) * r1_e1e2u(ji,jj)
                  zzWn_e3 = MIN( zWu, 0._wp ) / e3uw(ji,jj,jk  ,Kaa)
                  zzWp_e3 = MAX( zWu, 0._wp ) / e3uw(ji,jj,jk+1,Kaa)
                  zwi(ji,jk) = zwi(ji,jk) + zzWn_e3
                  zws(ji,jk) = zws(ji,jk)           - zzWp_e3
                  zwd(ji,jk) = zwd(ji,jk) - zzWn_e3 + zzWp_e3
               END_2D
               !
               DO_1Di( 0, 0 )          ! surface boundary conditions  ( wi(jk=1) always 0 )
                  zWu  = zDt_4 * ( e1e2t(ji,jj)*wi(ji,jj,2) + e1e2t(ji+1,jj)*wi(ji+1,jj,2) ) * r1_e1e2u(ji,jj)
                  zzWp_e3 = MAX( zWu, 0._wp ) / e3uw(ji,jj,2,Kaa)
                  zws(ji,1) = zws(ji,1)             - zzWp_e3
                  zwd(ji,1) = zwd(ji,1)             + zzWp_e3
               END_1D
            ELSE                                !  Flux form  !   dk[ WI_uw * upstream u ] 
               DO_2Dik( 0, 0,   2, jpkm1, 1 )
                  zDt_2bua =  zDt_2 * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kaa)
                  zWui = ( e1e2t(ji,jj)*wi(ji,jj,jk  ) + e1e2t(ji+1,jj)*wi(ji+1,jj,jk  ) ) * zDt_2bua 
                  zWus = ( e1e2t(ji,jj)*wi(ji,jj,jk+1) + e1e2t(ji+1,jj)*wi(ji+1,jj,jk+1) ) * zDt_2bua
                  zwi(ji,jk) = zwi(ji,jk) + MIN( zWui, 0._wp )
                  zws(ji,jk) = zws(ji,jk)                      - MAX( zWus, 0._wp )
                  zwd(ji,jk) = zwd(ji,jk) + MAX( zWui, 0._wp ) - MIN( zWus, 0._wp )
               END_2D
               !
               DO_1Di( 0, 0 )          ! surface boundary conditions  ( wi(jk=1) always 0 )
                  zDt_2bua = zDt_2 * r1_e1e2u(ji,jj) / e3u(ji,jj,1,Kaa)
                  zWus = ( e1e2t(ji,jj)*wi(ji,jj,2) + e1e2t(ji+1,jj)*wi(ji+1,jj,2) ) * zDt_2bua
                  zws(ji,1) = zws(ji,1)                        - MAX( zWus, 0._wp )
                  zwd(ji,1) = zwd(ji,1)                        - MIN( zWus, 0._wp )
               END_1D
            ENDIF
            !
         ENDIF
         !
         IF( nldf_dyn == np_lap_i ) THEN     !*  add LDF contribution  *!   (rotated lateral diffusion case)
            DO_2Dik( 0, 0,   2, jpkm1, 1 )            ! inner values
               zzwi = - rDt * akzu(ji,jj,jk  ) / ( e3u(ji,jj,jk,Kaa) * e3uw(ji,jj,jk  ,Kmm) ) * wumask(ji,jj,jk  )
               zzws = - rDt * akzu(ji,jj,jk+1) / ( e3u(ji,jj,jk,Kaa) * e3uw(ji,jj,jk+1,Kmm) ) * wumask(ji,jj,jk+1)
               zwi(ji,jk) = zwi(ji,jk) + zzwi
               zws(ji,jk) = zws(ji,jk)        + zzws
               zwd(ji,jk) = zwd(ji,jk) - zzwi - zzws
            END_2D
            DO_1Di( 0, 0 )                            ! surface boundary conditions (zwi already set to 0)
               zzws = - rDt * akzu(ji,jj,2) / ( e3u(ji,jj,1,Kaa) * e3uw(ji,jj,2,Kmm) ) * wumask(ji,jj,2)
               zws(ji,1) = zws(ji,1)          + zzws
               zwd(ji,1) = zwd(ji,1)          - zzws
            END_1D
         ENDIF
         !
         !                                   !*  add semi-implicit bottom friction  *!
         !
         !     Only needed for semi-implicit bottom friction setup. The explicit
         !     bottom friction has been included in "u(v)a" which act as the R.H.S
         !     column vector of the tri-diagonal matrix equation
         !
         IF( ln_drgimp ) THEN      ! implicit bottom friction
            DO_1Di( 0, 0 )
               iku = mbku(ji,jj)       ! ocean bottom level at u- and v-points
               zwd(ji,iku) = zwd(ji,iku) - zDt_2 *( rCdU_bot(ji+1,jj)+rCdU_bot(ji,jj) ) / e3u(ji,jj,iku,Kaa)
            END_1D
            IF( ln_isfcav .OR. ln_drgice_imp ) THEN   ! top friction (always implicit)
               DO_1Di( 0, 0 )
               !!gm   top Cd is masked (=0 outside cavities) no need of test on mik>=2  ==>> it has been suppressed
                  iku = miku(ji,jj)       ! ocean top level at u- and v-points
                  zwd(ji,iku) = zwd(ji,iku) - zDt_2 *( rCdU_top(ji+1,jj)+rCdU_top(ji,jj) ) / e3u(ji,jj,iku,Kaa)
               END_1D
            ENDIF
         ENDIF
         !
         ! Matrix inversion starting from the first level
         !-----------------------------------------------------------------------
         !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
         !
         !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
         !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
         !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
         !        (        ...               )( ...  ) ( ...  )
         !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
         !
         !   m is decomposed in the product of an upper and a lower triangular matrix
         !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
         !   The solution (the after velocity) is in puu(:,:,:,Kaa)
         !-----------------------------------------------------------------------
         !
         DO_2Dik( 0, 0,    2, jpkm1, 1 )   !==  First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1   (increasing k)  ==
            zwd(ji,jk) = zwd(ji,jk) - zwi(ji,jk) * zws(ji,jk-1) / zwd(ji,jk-1)
         END_2D
         !
         DO_1Di( 0, 0 )                    !==  second recurrence:    SOLk = RHSk - Lk / Dk-1  Lk-1  ==!
#if defined key_RK3
            !                                  ! RK3: use only utau (not utau_b)
            puu(ji,jj,1,Kaa) = puu(ji,jj,1,Kaa) + rDt * utauU(ji,jj)   &
                 &                                    / ( e3u(ji,jj,1,Kaa) * rho0 ) * umask(ji,jj,1)
#else
            !                                  ! MLF: average of utau and utau_b
            puu(ji,jj,1,Kaa) = puu(ji,jj,1,Kaa) + zDt_2 * ( utau_b(ji,jj) + utauU(ji,jj) )   &
                 &                                      / ( e3u(ji,jj,1,Kaa) * rho0 ) * umask(ji,jj,1)
#endif
         END_1D
         DO_2Dik( 0, 0,     2, jpkm1, 1 )
            puu(ji,jj,jk,Kaa) = puu(ji,jj,jk,Kaa) - zwi(ji,jk) / zwd(ji,jk-1) * puu(ji,jj,jk-1,Kaa)
         END_2D
         !
         DO_1Di( 0, 0 )                    !==  third recurrence : SOLk = ( Lk - Uk * Ek+1 ) / Dk  ==!
            puu(ji,jj,jpkm1,Kaa) = puu(ji,jj,jpkm1,Kaa) / zwd(ji,jpkm1)
         END_1D
         DO_2Dik( 0, 0,    jpk-2, 1, -1 )
            puu(ji,jj,jk,Kaa) = ( puu(ji,jj,jk,Kaa) - zws(ji,jk) * puu(ji,jj,jk+1,Kaa) ) / zwd(ji,jk)
         END_2D
         !
         !              !-------------------------------!
         !              !==  Vertical diffusion on v  ==!
         !              !-------------------------------!
         !
         !                    !**  tridiagonal matrix construction  **!    diagonal (zwd), lower (zwi), upper (zws)
         !
         !                                   !*  ZDF contribution  *! 
         !
         DO_2Dik( 0, 0,    2, jpkm1, 1 )
            zzwi = - zDt_2 * ( avm(ji,jj+1,jk  )     +  avm(ji,jj,jk  )     ) &
               &           / ( e3v(ji,jj  ,jk  ,Kaa) * e3vw(ji,jj,jk  ,Kmm) ) * wvmask(ji,jj,jk  )
            zzws = - zDt_2 * ( avm(ji,jj+1,jk+1)     +  avm(ji,jj,jk+1)     ) &
               &           / ( e3v(ji,jj  ,jk  ,Kaa) * e3vw(ji,jj,jk+1,Kmm) ) * wvmask(ji,jj,jk+1)
            zwi(ji,jk) = zzwi
            zws(ji,jk) = zzws
            zwd(ji,jk) = 1._wp - zwi(ji,jk) - zws(ji,jk)
         END_2D
         !
         DO_1Di( 0, 0 )          ! surface boundary conditions
            zzws = - zDt_2 * ( avm(ji,jj+1,2)     +  avm(ji,jj,2)     ) &
               &           / ( e3v(ji,jj  ,1,Kaa) * e3vw(ji,jj,2,Kmm) ) * wvmask(ji,jj,2)
            zwi(ji,1) = 0._wp
            zws(ji,1) = zzws
            zwd(ji,1) = 1._wp - zzws
         END_1D

! !! gm  NOUVEAU CAS implicite :                 <<<<<<<<<<=======     A FAIRE !!!!
!        
!        
!          IF( ln_zad_Aimp ) THEN               !*  add Courant-number-dependent Implicit Vertical Advection  *!
!             !
!             !                                !  Flux form  !   dk[ WI_vw * upstream v ] 
!             DO_2Dik( 0, 0,    2, jpkm1, 1 )
!                zDt_2bva = zDt_2 * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kaa)
!                zWvi = ( e1e2t(ji,jj)*wi(ji,jj,jk  ) + e1e2t(ji,jj+1)*wi(ji,jj+1,jk  ) ) * zDt_2bva
!                zWvs = ( e1e2t(ji,jj)*wi(ji,jj,jk+1) + e1e2t(ji,jj+1)*wi(ji,jj+1,jk+1) ) * zDt_2bva
!                zwi(ji,jk) = zwi(ji,jk) + MIN( zWvi, 0._wp )
!                zws(ji,jk) = zws(ji,jk)                      - MAX( zWvs, 0._wp )
!                zwd(ji,jk) = zwd(ji,jk) + MAX( zWvi, 0._wp ) - MIN( zWvs, 0._wp )
!             END_2D
!             DO_1Di( 0, 0 )                         ! surface boundary conditions   (zwi already set to 0)
!                zDt_2bva = zDt_2 * r1_e1e2v(ji,jj) / e3v(ji,jj,1,Kaa)
!                zWvs = ( e1e2t(ji,jj  )*wi(ji,jj  ,2) +  e1e2t(ji,jj+1)*wi(ji,jj+1,2) ) * zDt_2bva
!                zws(ji,1) = zws(ji,1)                        - MAX( zWvs, 0._wp )
!                zwd(ji,1) = zwd(ji,1)                        - MIN( zWvs, 0._wp )
!             END_1D
!             IF( ln_dynadv_vec ) THEN            !  Vector Invariant Form  !   remove 
!                DO_2Dik( 0, 0,   1, jpkm1, 1 )
!                   zDt_2bva =  zDt_2 * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kaa)
!                   zwd(ji,jk) = zwd(ji,jk) - zDt_2bva * (   e1e2t(ji,jj  )*( wi(ji,jj  ,jk) - wi(ji,jj  ,jk+1) )   &
!                      &                                   + e1e2t(ji,jj+1)*( wi(ji,jj+1,jk) - wi(ji,jj+1,jk+1) )   )
!                END_2D
!             ENDIF
!             !
!          ENDIF
! 
! !!gm end
         !
         IF( ln_zad_Aimp ) THEN              !*  add Courant-number-dependent Implicit Vertical Advection  *!
            !
            IF( ln_dynadv_vec ) THEN            !  Vector Invariant Form  !   Wv * upstream dk[V]
               DO_2Dik( 0, 0,    2, jpkm1, 1 )
                  zWv  = zDt_4 * (  ( e1e2t(ji,jj)*wi(ji,jj,jk  ) + e1e2t(ji,jj+1)*wi(ji,jj+1,jk  ) )  &
                     &            + ( e1e2t(ji,jj)*wi(ji,jj,jk+1) + e1e2t(ji,jj+1)*wi(ji,jj+1,jk+1) )  ) * r1_e1e2v(ji,jj)
                  zzWn_e3 = MIN( zWv, 0._wp ) / e3vw(ji,jj,jk  ,Kaa)
                  zzWp_e3 = MAX( zWv, 0._wp ) / e3vw(ji,jj,jk+1,Kaa)
                  zwi(ji,jk) = zwi(ji,jk) + zzWn_e3
                  zws(ji,jk) = zws(ji,jk)           - zzWp_e3
                  zwd(ji,jk) = zwd(ji,jk) - zzWn_e3 + zzWp_e3
               END_2D
               DO_1Di( 0, 0 )                         ! surface boundary conditions   (zwi already set to 0)
                  zWv  = zDt_4 * ( e1e2t(ji,jj)*wi(ji,jj,2) + e1e2t(ji,jj+1)*wi(ji,jj+1,2) ) * r1_e1e2v(ji,jj)
                  zzWp_e3 = MAX( zWv, 0._wp ) / e3vw(ji,jj,2,Kaa)
                  zws(ji,1) = zws(ji,1)             - zzWp_e3
                  zwd(ji,1) = zwd(ji,1)             + zzWp_e3
               END_1D
               !
            ELSE                       !---- Flux form ----!
               DO_2Dik( 0, 0,    2, jpkm1, 1 )
                  zDt_2bva = zDt_2 * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kaa)
                  zWvi = ( e1e2t(ji,jj)*wi(ji,jj,jk  ) + e1e2t(ji,jj+1)*wi(ji,jj+1,jk  ) ) * zDt_2bva
                  zWvs = ( e1e2t(ji,jj)*wi(ji,jj,jk+1) + e1e2t(ji,jj+1)*wi(ji,jj+1,jk+1) ) * zDt_2bva
                  zwi(ji,jk) = zwi(ji,jk) + MIN( zWvi, 0._wp )
                  zws(ji,jk) = zws(ji,jk)                      - MAX( zWvs, 0._wp )
                  zwd(ji,jk) = zwd(ji,jk) + MAX( zWvi, 0._wp ) - MIN( zWvs, 0._wp )
               END_2D
               DO_1Di( 0, 0 )                         ! surface boundary conditions   (zwi already set to 0)
                  zDt_2bva = zDt_2 * r1_e1e2v(ji,jj) / e3v(ji,jj,1,Kaa)
                  zWvs = ( e1e2t(ji,jj  )*wi(ji,jj  ,2) +  e1e2t(ji,jj+1)*wi(ji,jj+1,2) ) * zDt_2bva
                  zws(ji,1) = zws(ji,1)                        - MAX( zWvs, 0._wp )
                  zwd(ji,1) = zwd(ji,1)                        - MIN( zWvs, 0._wp )
               END_1D
            ENDIF
            !
         ENDIF
         !
         IF( nldf_dyn == np_lap_i ) THEN     !*  add LDF contribution  *!   (rotated lateral diffusion case)
            DO_2Dik( 0, 0,    2, jpkm1, 1 )           ! inner values
               zzwi = - rDt * akzv(ji,jj,jk  ) / ( e3v(ji,jj,jk,Kaa) * e3vw(ji,jj,jk  ,Kmm) ) * wvmask(ji,jj,jk  )
               zzws = - rDt * akzv(ji,jj,jk+1) / ( e3v(ji,jj,jk,Kaa) * e3vw(ji,jj,jk+1,Kmm) ) * wvmask(ji,jj,jk+1)
               zwi(ji,jk) = zwi(ji,jk) + zzwi
               zws(ji,jk) = zws(ji,jk)        + zzws
               zwd(ji,jk) = zwd(ji,jk) - zzwi - zzws
            END_2D
            DO_1Di( 0, 0 )                            ! surface boundary conditions   (zwi already set to 0)
               zwi(ji,1) = 0._wp
               zzws = - rDt * akzv(ji,jj,2) / ( e3v(ji,jj,1,Kaa) * e3vw(ji,jj,2,Kmm) ) * wvmask(ji,jj,2)
               zws(ji,1) = zws(ji,1)          + zzws
               zwd(ji,1) = zwd(ji,1)          - zzws
            END_1D
         ENDIF
         !
         !                                   !*  add semi-implicit bottom friction  *!
         !
         !     Only needed for semi-implicit bottom friction setup. The explicit
         !     bottom friction has been included in "u(v)a" which act as the R.H.S
         !     column vector of the tri-diagonal matrix equation
         !
         IF ( ln_drgimp ) THEN      ! implicit bottom friction
            DO_1Di( 0, 0 )
               ikv = mbkv(ji,jj)       ! (deepest ocean u- and v-points)
               zwd(ji,ikv) = zwd(ji,ikv) - zDt_2*( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj) )   &
                  &                                   / e3v(ji,jj,ikv,Kaa)
            END_1D
            IF ( ln_isfcav.OR.ln_drgice_imp ) THEN   ! top friction (always implicit)
               DO_1Di( 0, 0 )
                  ikv = mikv(ji,jj)       ! (first wet ocean u- and v-points)
                  zwd(ji,ikv) = zwd(ji,ikv) - zDt_2*( rCdU_top(ji,jj+1)+rCdU_top(ji,jj) )   &
                     &                                   / e3v(ji,jj,ikv,Kaa)
               END_1D
            ENDIF
         ENDIF
         !
         ! Matrix inversion starting from the first level
         !-----------------------------------------------------------------------
         !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
         !
         !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
         !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
         !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
         !        (        ...               )( ...  ) ( ...  )
         !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
         !
         !   m is decomposed in the product of an upper and lower triangular matrix
         !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
         !   The solution (after velocity) is in 2d array va
         !-----------------------------------------------------------------------
         !
         DO_2Dik( 0, 0,    2, jpkm1, 1 )   !==  First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1   (increasing k)  ==
            zwd(ji,jk) = zwd(ji,jk) - zwi(ji,jk) * zws(ji,jk-1) / zwd(ji,jk-1)
         END_2D
         !
         DO_1Di( 0, 0 )                    !==  second recurrence:    SOLk = RHSk - Lk / Dk-1  Lk-1  ==!
#if defined key_RK3
            !                                  ! RK3: use only vtau (not vtau_b)
            pvv(ji,jj,1,Kaa) = pvv(ji,jj,1,Kaa) + rDt * vtauV(ji,jj)   &
               &                                   / ( e3v(ji,jj,1,Kaa) * rho0 ) * vmask(ji,jj,1)
#else
            !                                  ! MLF: average of vtau and vtau_b
            pvv(ji,jj,1,Kaa) = pvv(ji,jj,1,Kaa) + zDt_2*( vtau_b(ji,jj) + vtauV(ji,jj) )   &
                 &                                 / ( e3v(ji,jj,1,Kaa) * rho0 ) * vmask(ji,jj,1)
#endif
         END_1D
         DO_2Dik( 0, 0,    2, jpkm1, 1 )
            pvv(ji,jj,jk,Kaa) = pvv(ji,jj,jk,Kaa) - zwi(ji,jk) / zwd(ji,jk-1) * pvv(ji,jj,jk-1,Kaa)
         END_2D
         !
         DO_1Di( 0, 0 )                    !==  third recurrence : SOLk = ( Lk - Uk * SOLk+1 ) / Dk  ==!
            pvv(ji,jj,jpkm1,Kaa) = pvv(ji,jj,jpkm1,Kaa) / zwd(ji,jpkm1)
         END_1D
         DO_2Dik( 0, 0,   jpk-2, 1, -1 )
            pvv(ji,jj,jk,Kaa) = ( pvv(ji,jj,jk,Kaa) - zws(ji,jk) * pvv(ji,jj,jk+1,Kaa) ) / zwd(ji,jk)
         END_2D
         !
         !                                            ! ================= !
      END_1D                                          !  i-k slices loop  !
      !                                               ! ================= !
      !
      IF( l_trddyn )   THEN                      
         ! save the vertical diffusive trends for further diagnostics
         ! then calculate implicit drag trends
         ! (The order of operation is important here). 
         !
         ztrdu(:,:,:) = ( puu(:,:,:,Kaa) - ztrdu(:,:,:) ) * r1_Dt 
         ztrdv(:,:,:) = ( pvv(:,:,:,Kaa) - ztrdv(:,:,:) ) * r1_Dt 
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_zdf, kt, Kmm, Kaa )
         !
         ztrdu_fr(:,:,:) = 0._wp    ;   ztrdv_fr(:,:,:) = 0._wp
         IF( nn_ice == 2 ) THEN
            !                             ! trends due to ice-ocean drag in surface layer
            ztrdu_fr(:,:,1) = ( uiceoc_b(:,:) + uiceoc(:,:) ) / ( e3u(:,:,1,Kmm) * rho0 )
            ztrdv_fr(:,:,1) = ( viceoc_b(:,:) + viceoc(:,:) ) / ( e3v(:,:,1,Kmm) * rho0 )
            CALL trd_dyn( ztrdu_fr(:,:,1), ztrdv_fr(:,:,1), jpdyn_iceoc, kt, Kmm )
         ENDIF
         !
         !                             ! wind stress trends in surface layer
         ztrdu_fr(:,:,1) = ( utau_b(:,:) + utau(:,:) ) / ( e3u(:,:,1,Kmm) * rho0 )
         ztrdv_fr(:,:,1) = ( vtau_b(:,:) + vtau(:,:) ) / ( e3v(:,:,1,Kmm) * rho0 )
         CALL trd_dyn( ztrdu_fr(:,:,1), ztrdv_fr(:,:,1), jpdyn_tau, kt, Kmm)
         !
         IF( ln_drgimp ) THEN
            !                          ! trends due to implicit bottom friction
            ztrdu_fr(:,:,1) = 0._wp    ;   ztrdv_fr(:,:,1) = 0._wp
            DO_2D( 0, 0, 0, 0 )
               iku = mbku(ji,jj)          ! deepest ocean u- & v-levels
               ikv = mbkv(ji,jj)
               ztrdu_fr(ji,jj,iku) = 0.5 * ( rCdU_bot(ji+1,jj) + rCdU_bot(ji,jj) )  &
                  &                * puu(ji,jj,iku,Kmm) / e3u(ji,jj,iku,Kmm)
               ztrdv_fr(ji,jj,ikv) = 0.5 * ( rCdU_bot(ji,jj+1) + rCdU_bot(ji,jj) )  &
                  &                * pvv(ji,jj,ikv,Kmm) / e3v(ji,jj,ikv,Kmm)
            END_2D
            CALL trd_dyn( ztrdu_fr, ztrdv_fr, jpdyn_bfri, kt, Kmm )
            !                          ! trends due to implicit top friction
            IF( ln_isfcav.OR.ln_drgice_imp ) THEN    ! Ocean cavities (ISF) or implicit ice-ocean drag
               ztrdu_fr(:,:,:) = 0._wp    ;   ztrdv_fr(:,:,:) = 0._wp
               DO_2D( 0, 0, 0, 0 )
                  iku = miku(ji,jj)       ! ocean top level at u-points
                  ikv = mikv(ji,jj)       ! ocean top level at v-points
                  ztrdu_fr(ji,jj,iku) = 0.5 * ( rCdU_top(ji+1,jj) + rCdU_top(ji,jj) )  &
                     &                * puu(ji,jj,iku,Kmm) / e3u(ji,jj,iku,Kmm)
                  ztrdv_fr(ji,jj,ikv) = 0.5 * ( rCdU_top(ji,jj+1) + rCdU_top(ji,jj) )  &
                     &                * pvv(ji,jj,ikv,Kmm) / e3v(ji,jj,ikv,Kmm)
               END_2D
               CALL trd_dyn( ztrdu_fr, ztrdv_fr, jpdyn_tfri, kt, Kmm )
            ENDIF
         ENDIF
         !
         DEALLOCATE( ztrdu, ztrdv, ztrdu_fr, ztrdv_fr ) 
      ENDIF
      !                                          ! print mean trends (used for debugging)
      IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab3d_1=puu(:,:,:,Kaa), clinfo1=' zdf  - Ua: ', mask1=umask,               &
         &                                  tab3d_2=pvv(:,:,:,Kaa), clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
         !
      IF( ln_timing )   CALL timing_stop('dyn_zdf')
      !
   END SUBROUTINE dyn_zdf

   !!==============================================================================
END MODULE dynzdf


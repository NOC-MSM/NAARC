MODULE traqsr
   !!======================================================================
   !!                       ***  MODULE  traqsr  ***
   !! Ocean physics:   solar radiation penetration in the top ocean levels
   !!======================================================================
   !! History :  OPA  !  1990-10  (B. Blanke)  Original code
   !!            7.0  !  1991-11  (G. Madec)
   !!                 !  1996-01  (G. Madec)  s-coordinates
   !!   NEMO     1.0  !  2002-06  (G. Madec)  F90: Free form and module
   !!             -   !  2005-11  (G. Madec) zco, zps, sco coordinate
   !!            3.2  !  2009-04  (G. Madec & NEMO team)
   !!            3.6  !  2012-05  (C. Rousset) store attenuation coef for use in ice model
   !!            3.6  !  2015-12  (O. Aumont, J. Jouanno, C. Ethe) use vertical profile of chlorophyll
   !!            3.7  !  2015-11  (G. Madec, A. Coward)  remove optimisation for fix volume
   !!            4.0  !  2020-11  (A. Coward)  optimisation
   !!            4.5  !  2021-03  (G. Madec)  further optimisation + adaptation for RK3
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_qsr       : temperature trend due to the penetration of solar radiation
   !!       qsr_RGBc  : IR + RGB      light penetration with Chlorophyll data case
   !!       qsr_RGB   : IR + RGB      light penetration with constant Chlorophyll case
   !!       qsr_2BD   : 2 bands (InfraRed + Visible light) case
   !!       qsr_5BDc  : IR + RGB + UV light penetration with Chlorophyll data case
   !!       qsr_5BD   : IR + RGB + UV light penetration with constant Chlorophyll case
   !!       qsr_ext_lev : level of extinction for each bands
   !!   tra_qsr_init  : initialization of the qsr penetration
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers
   USE phycst         ! physical constants
   USE dom_oce        ! ocean space and time domain
   USE domtile
   USE sbc_oce        ! surface boundary condition: ocean
   USE trc_oce        ! share SMS/Ocean variables
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! trends manager: tracers
   !
   USE in_out_manager ! I/O manager
   USE prtctl         ! Print control
   USE iom            ! I/O library
   USE fldread        ! read input fields
   USE restart        ! ocean restart
   USE lib_mpp        ! MPP library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_qsr       ! routine called by step.F90 (ln_traqsr=T)
   PUBLIC   tra_qsr_init  ! routine called by nemogcm.F90

   !                                 !!* Namelist namtra_qsr: penetrative solar radiation
   LOGICAL , PUBLIC ::   ln_traqsr    !: light absorption (qsr) flag
   LOGICAL , PUBLIC ::   ln_qsr_rgb   !: Red-Green-Blue light absorption flag
   LOGICAL , PUBLIC ::   ln_qsr_2bd   !: 2 band         light absorption flag
   LOGICAL , PUBLIC ::   ln_qsr_5bd   !: 5 bands        light absorption flag
   LOGICAL , PUBLIC ::   ln_qsr_bio   !: bio-model      light absorption flag
   INTEGER , PUBLIC ::   nn_chldta    !: use Chlorophyll data 3D/Surface (=2/1) or not (=0)
   REAL(wp), PUBLIC ::   rn_abs       !: fraction absorbed in the very near surface (RGB, 2 bands & 5 bands)
   REAL(wp), PUBLIC ::   rn_si0       !: very near surface depth of extinction      (RGB, 2 bands & 5 bands)
   REAL(wp), PUBLIC ::   rn_si1       !: deepest depth of extinction (water type I)            (2 bands)
   REAL(wp), PUBLIC ::   rn_par       !: photosynthetically active radiation fraction absorbed (5 bands)
   !
   INTEGER, PARAMETER ::   np_RGB  = 1   ! R-G-B     light penetration with constant Chlorophyll
   INTEGER, PARAMETER ::   np_RGBc = 2   ! R-G-B     light penetration with Chlorophyll data
   INTEGER, PARAMETER ::   np_2BD  = 3   ! 2 bands   light penetration
   INTEGER, PARAMETER ::   np_5BD  = 4   ! R-G-B-UV  light penetration with constant Chlorophyll
   INTEGER, PARAMETER ::   np_5BDc = 5   ! R-G-B-UV  light penetration with Chlorophyll data
   INTEGER, PARAMETER ::   np_BIO  = 6   ! bio-model light penetration
   !
   INTEGER  ::   nqsr     ! user choice of the type of light penetration
   INTEGER  ::   nc_tab   ! RGB with cst Chlorophyll: index associated with the chosen Chl value
   !
   !                       ! extinction level 
   INTEGER  ::   nk0             !: IR (depth larger ~12 m)
   INTEGER  ::   nkV             !: Visible light (depth larger than ~840 m) 
   INTEGER  ::   nkR, nkG, nkB, nkU   !: RGB (depth larger than ~100 m, ~470 m, ~1700 m, resp.) & UV
   !
   INTEGER, PUBLIC  ::   nksr    !: =nkV, i.e. maximum level of light extinction (used in traatf(_qco).F90)
   !
   !                  ! inverse of attenuation length
   REAL(wp) ::   r1_si0                     ! all schemes : infrared  = 1/rn_si0 
   REAL(wp) ::   r1_si1                     ! 2 band      : mean RGB  = 1/rn_si1   
   REAL(wp) ::   r1_LR, r1_LG, r1_LB, r1_LU ! RGB & UV with constant Chl
   REAL(wp) ::   zz0
   !
   REAL(wp) , PUBLIC, DIMENSION(4,61)   ::   rktab    ! tabulated attenuation coefficients for RGB-UV absorption
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_chl   ! structure of input Chl (file informations, fields read)

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_qsr( kt, Kmm, pts, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_qsr  ***
      !!
      !! ** Purpose :   Compute the temperature trend due to the solar radiation
      !!              penetration and add it to the general temperature trend.
      !!
      !! ** Method  : The profile of the solar radiation within the ocean is defined
      !!      through 2 wavebands (rn_si0,rn_si1), 4 wavebands (IR-RGB), 
      !!         5 wavebands (IR-RGB-UV) or computed by the biogeochemical model
      !!         The computation is only done down to the level where
      !!      I(k) < 1.e-15 W/m2 (i.e. over the top nk levels) .
      !!
      !! ** Action  : - update ts(jp_tem) with the penetrative solar radiation trend
      !!              - send  trend for further diagnostics (l_trdtra=T)
      !!----------------------------------------------------------------------
      INTEGER,                                   INTENT(in   ) ::   kt, Kmm, Krhs   ! ocean time-step and time-level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts,jpt), INTENT(inout) ::   pts             ! active tracers and RHS of tracer equation
      !
      INTEGER  ::   ji, jj, jk               ! dummy loop indices
      REAL(wp) ::   z1_2, ze3t                     ! local scalars
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   ztrdt, zetot
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('tra_qsr')
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'tra_qsr : penetration of the surface solar radiation'
            IF(lwp) WRITE(numout,*) '~~~~~~~'
         ENDIF
      ENDIF
      !
      IF( l_trdtra ) THEN     ! trends diagnostic: save the input temperature trend
         ALLOCATE( ztrdt(T2D(0),jpk) )
         ztrdt(:,:,:) = pts(T2D(0),:,jp_tem,Krhs)
      ENDIF
      ! 
#if ! defined key_RK3
      !                        ! MLF only : heat content trend due to Qsr flux (qsr_hc)
      !
      !                         !-----------------------------------!
      !                         !  before qsr induced heat content  !
      !                         !-----------------------------------!
      IF( kt == nit000 ) THEN          !==  1st time step  ==!
         IF( ln_rstart .AND. .NOT.l_1st_euler ) THEN    ! read in restart
            z1_2 = 0.5_wp
            IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                        ! Do only on the first tile
               IF(lwp) WRITE(numout,*) '          nit000-1 qsr tracer content forcing field read in the restart file'
               CALL iom_get( numror, jpdom_auto, 'qsr_hc_b', qsr_hc_b )   ! before heat content trend due to Qsr flux
            ENDIF
         ELSE                                           ! No restart or Euler forward at 1st time step
            z1_2 = 1._wp
            DO_3D( 0, 0, 0, 0, 1, jpk )
               qsr_hc_b(ji,jj,jk) = 0._wp
            END_3D
         ENDIF
      ELSE                             !==  Swap of qsr heat content  ==!
         z1_2 = 0.5_wp
         DO_3D( 0, 0, 0, 0, 1, jpk )
            qsr_hc_b(ji,jj,jk) = qsr_hc(ji,jj,jk)
         END_3D
      ENDIF
#endif

      !                       !----------------------------!
      SELECT CASE( nqsr )     !  qsr induced heat content  !
      !                       !----------------------------!
      !
      CASE( np_RGBc )   ;   CALL qsr_RGBc( kt, Kmm, pts, Krhs )  !==  R-G-B fluxes using chlorophyll data     ==!    with Morel &Berthon (1989) vertical profile
         !
      CASE( np_RGB  )   ;   CALL qsr_RGB ( kt, Kmm, pts, Krhs )  !==  R-G-B fluxes with constant chlorophyll  ==!   
         !
      CASE( np_2BD  )   ;   CALL qsr_2BD (     Kmm, pts, Krhs )  !==  2-bands fluxes                          ==!
         !
      CASE( np_5BDc )   ;   CALL qsr_5BDc( kt, Kmm, pts, Krhs )  !== IR-R-G-B-UV fluxes using chlorophyll data     ==!    with Morel &Berthon (1989) vertical profile
         !
      CASE( np_5BD  )   ;   CALL qsr_5BD ( kt, Kmm, pts, Krhs )  !== IR-R-G-B-UV fluxes with constant chlorophyll  ==!   
         !
      CASE( np_BIO )                                     !==  bio-model fluxes                        ==!
         DO_3D( 0, 0, 0, 0, 1, nkV )
#if defined key_RK3
            !                                                  !- RK3 : temperature trend at jk t-level
            ze3t   = e3t(ji,jj,jk,Kmm)
            pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs) + r1_rho0_rcp * ( etot3(ji,jj,jk) - etot3(ji,jj,jk+1) ) / ze3t
#else
            !                                                  !- MLF : heat content trend due to Qsr flux (qsr_hc)
            qsr_hc(ji,jj,jk) = r1_rho0_rcp * ( etot3(ji,jj,jk) - etot3(ji,jj,jk+1) )
#endif
         END_3D
         !                                                     !- sea-ice : store the 1st level attenuation coefficient
         WHERE( etot3(T2D(0),1) /= 0._wp )   ;   fraqsr_1lev(T2D(0)) = 1._wp - etot3(T2D(0),2) / etot3(T2D(0),1)
         ELSEWHERE                           ;   fraqsr_1lev(T2D(0)) = 1._wp
         END WHERE
         !
      END SELECT
      !
#if ! defined key_RK3
      !                             ! MLF : add the temperature trend
      DO_3D( 0, 0, 0, 0, 1, nksr )
         pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs)   &
            &                      + z1_2 * ( qsr_hc_b(ji,jj,jk) + qsr_hc(ji,jj,jk) )   &
            &                             / e3t(ji,jj,jk,Kmm)
      END_3D
      !
      ! sea-ice: store the 1st ocean level attenuation coefficient
      DO_2D( 0, 0, 0, 0 )
         zz0 = r1_rho0_rcp * qsr(ji,jj)   ! test zz0 and not qsr for rounding errors in single precision
         IF( zz0 /= 0._wp ) THEN   ;   fraqsr_1lev(ji,jj) = qsr_hc(ji,jj,1) / zz0
         ELSE                      ;   fraqsr_1lev(ji,jj) = 1._wp
         ENDIF
      END_2D
#endif
      !
      IF( l_trdtra ) THEN     ! qsr tracers trends saved for diagnostics
         ztrdt(:,:,:) = pts(T2D(0),:,jp_tem,Krhs) - ztrdt(:,:,:)
         CALL trd_tra( kt, Kmm, Krhs, 'TRA', jp_tem, jptra_qsr, ztrdt )
         DEALLOCATE( ztrdt )
      ENDIF
      !
#if ! defined key_RK3
      IF( iom_use('qsr3d') ) THEN      ! output the shortwave Radiation distribution
         ALLOCATE( zetot(T2D(0),jpk) )
         zetot(:,:,nksr+1:jpk) = 0._wp     ! below ~400m set to zero
         DO_3DS(0, 0, 0, 0, nksr, 1, -1)
            zetot(ji,jj,jk) = zetot(ji,jj,jk+1) + qsr_hc(ji,jj,jk) * rho0_rcp
         END_3D
         CALL iom_put( 'qsr3d', zetot )   ! 3D distribution of shortwave Radiation
         DEALLOCATE( zetot )
      ENDIF
#endif
      !
      IF( .NOT. l_istiled .OR. ntile == nijtile )  THEN                ! Do only on the last tile
         IF( lrst_oce ) THEN     ! write in the ocean restart file
#if ! defined key_RK3
            CALL iom_rstput( kt, nitrst, numrow, 'qsr_hc_b'   , qsr_hc      )
#endif
            CALL iom_rstput( kt, nitrst, numrow, 'fraqsr_1lev', fraqsr_1lev )
         ENDIF
      ENDIF
      !
      !                       ! print mean trends (used for debugging)
      IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab3d_1=pts(:,:,:,jp_tem,Krhs), clinfo1=' qsr  - Ta: ', mask1=tmask, clinfo3='tra-ta' )
      !
      IF( ln_timing )   CALL timing_stop('tra_qsr')
      !
   END SUBROUTINE tra_qsr


   SUBROUTINE qsr_RGBc( kt, Kmm, pts, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE qsr_RGBc  ***
      !!
      !! ** Purpose :   Red-Green-Blue solar radiation using chlorophyll data
      !!
      !! ** Method  : The profile of the solar radiation within the ocean is defined
      !!      through 2 wavebands (rn_si0,rn_si1) or 3 wavebands (RGB) and a ratio rn_abs
      !!      Considering the 2 wavebands case:
      !!         I(k) = Qsr*( rn_abs*EXP(z(k)/rn_si0) + (1.-rn_abs)*EXP(z(k)/rn_si1) )
      !!         The temperature trend associated with the solar radiation penetration
      !!         is given by : zta = 1/e3t dk[ I ] / (rho0*Cp)
      !!         At the bottom, boudary condition for the radiation is no flux :
      !!      all heat which has not been absorbed in the above levels is put
      !!      in the last ocean level.
      !!         The computation is only done down to the level where
      !!      I(k) < 1.e-15 W/m2 (i.e. over the top nk levels) .
      !!
      !! ** Action  : - update ta with the penetrative solar radiation trend
      !!              - send  trend for further diagnostics (l_trdtra=T)
      !!
      !! Reference  : Lengaigne et al. 2007, Clim. Dyn., V28, 5, 503-516.
      !!              Morel, A. et Berthon, JF, 1989, Limnol Oceanogr 34(8), 1545-1562
      !!----------------------------------------------------------------------
      INTEGER,                                   INTENT(in   ) ::   kt, Kmm, Krhs   ! ocean time-step and time-level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts,jpt), INTENT(inout) ::   pts         ! active tracers and RHS of tracer equation
      !!
      INTEGER  ::   ji, jj, jk, ik           ! dummy loop indices
      INTEGER  ::   ipk, itab                ! local integer
      REAL(wp) ::   zc1 , zc2 , zc3, zchl    ! local scalars
      REAL(wp) ::   zze0, zzeR, zzeG, zzeB, zzeT              !    -         -
      REAL(wp) ::   zz0 , zz1 , ze3t                          !    -         -
      REAL(wp) ::   zCb, zCmax, zpsi, zpsimax, zrdpsi, zCze   !    -         -
      REAL(wp) ::   zlogc, zlogze, zlogCtot, zlogCze          !    -         -
      !!
      REAL(wp), DIMENSION(T2D(0)) ::   ze0, zeR, zeG, zeB, zeT
#if ! defined key_PSYCLONE_2p5p0
      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::   zc
#else
      REAL(wp), DIMENSION(T2D(0),jpk,0:3) ::   zc
#endif
      !!----------------------------------------------------------------------
      !
      IF    ( nn_chldta == 1 ) THEN   ;   ipk=1
      ELSEIF( nn_chldta == 2 ) THEN   ;   ipk=jpk   ;   ENDIF
#if ! defined key_PSYCLONE_2p5p0
      ALLOCATE( zc(T2D(0),ipk,0:3) )
#endif
      !
      !                       !===========================================!
      !                       !==  R-G-B fluxes using chlorophyll data  ==!    with Morel &Berthon (1989) vertical profile
      !                       !===================================****====!
      !
      !                             !=  Chlorophyll data  =!
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                             ! Do only for the full domain
         IF( ln_tile ) CALL dom_tile_stop( ldhold=.TRUE. )                       ! Use full domain
         CALL fld_read( kt, 1, sf_chl )                                          ! Read Chl data and provides it at the current time step
         IF( ln_tile ) CALL dom_tile_start( ldhold=.TRUE. )                      ! Revert to tile domain
      ENDIF
      !
      DO_3D( 0, 0, 0, 0, 1, ipk )                          ! pre-calculated expensive coefficient
         zlogc = LOG(  MAX( 0.03_wp, MIN( sf_chl(1)%fnow(ji,jj,jk) ,10._wp ) )  ) ! zlogc = log(zchl)   with 0.03 <= Chl >= 10. 
         zc1   = 0.113328685307 + 0.803 * zlogc                                   ! zc1 : log(zCze)  = log (1.12  * zchl**0.803)
         zc2   = 3.703768066608 + 0.459 * zlogc                                   ! zc2 : log(zCtot) = log(40.6  * zchl**0.459)
         zc3   = 6.34247346942  - 0.746 * zc2                                     ! zc3 : log(zze)   = log(568.2 * zCtot**(-0.746))
         IF( zc3 > 4.62497281328 )   zc3 = 5.298317366548 - 0.293 * zc2           ! IF(log(zze)>log(102)) log(zze) = log(200*zCtot**(-0.293))
         !
         zc(ji,jj,jk,0) = zlogc                                                   ! ze(0) = log(zchl)
         zc(ji,jj,jk,1) = EXP( zc1 )                                              ! ze(1) = zCze
         zc(ji,jj,jk,2) = 1._wp / ( 0.710 + zlogc * ( 0.159 + zlogc * 0.021 ) )   ! ze(2) = 1/zdelpsi
         zc(ji,jj,jk,3) = EXP( - zc3 )                                            ! ze(3) = 1/zze
      END_3D
      !
      !                             !=  surface light  =!
      !
      zz0 =           rn_abs              ! Infrared absorption
      zz1 = ( 1._wp - rn_abs ) / 3._wp    ! R-G-B equi-partition
      !
      DO_2D( 0, 0, 0, 0 )                 ! surface light
         ze0(ji,jj) = zz0 * qsr(ji,jj)   ;   zeR(ji,jj) = zz1 * qsr(ji,jj)    ! IR    ; Red
         zeG(ji,jj) = zz1 * qsr(ji,jj)   ;   zeB(ji,jj) = zz1 * qsr(ji,jj)    ! Green ; Blue
         zeT(ji,jj) =       qsr(ji,jj)                                        ! Total
      END_2D
      !              
      !                             !=  interior light  =!
      !
      DO jk = 1, nk0                      !* near surface layers *!   (< ~12 meters : IR + RGB )
         ik = MIN( jk , ipk )
         DO_2D( 0, 0, 0, 0 )
            !                                      !- inverse of RGB attenuation lengths
            zlogc     = zc(ji,jj,ik,0)
            zCb       = 0.768 + zlogc * ( 0.087 - zlogc * ( 0.179 + zlogc * 0.025 ) )
            zCmax     = 0.299 - zlogc * ( 0.289 - zlogc * 0.579 )
            zpsimax   = 0.6   - zlogc * ( 0.640 - zlogc * ( 0.021 + zlogc * 0.115 ) )
            ! zdelpsi = 0.710 + zlogc * ( 0.159 + zlogc * 0.021 )
            zCze   = zc(ji,jj,ik,1)
            zrdpsi = zc(ji,jj,ik,2)                                     ! 1/zdelpsi
!!st05            zpsi   = zc(ji,jj,ik,3) * gdepw(ji,jj,jk,Kmm)               ! gdepw/zze
            zpsi   = zc(ji,jj,ik,3) * gdepw(ji,jj,jk+1,Kmm)               ! gdepw/zze
            !                                                        ! make sure zchl value is such that: 0.03 < zchl < 10. 
            zchl = MAX(  0.03_wp , MIN( zCze * ( zCb + zCmax * EXP( -( (zpsi - zpsimax) * zrdpsi )**2 ) ) , 10._wp )  )
            !                                                        ! Convert chlorophyll value to attenuation coefficient
            itab = NINT( 41 + 20.*LOG10(zchl) + 1.e-15 )             ! look-up table index
            !       Red             !         Green              !         Blue
            r1_LR = rktab(3,itab)   ;   r1_LG = rktab(2,itab)    ;   r1_LB = rktab(1,itab)
            !
            !                                      !- fluxes at jk+1 w-level
            ze3t = e3t(ji,jj,jk,Kmm)
            zze0 = ze0(ji,jj) * EXP( - ze3t*r1_si0 )   ;   zzeR = zeR(ji,jj) * EXP( - ze3t*r1_LR )   ! IR    ; Red  at jk+1 w-level
            zzeG = zeG(ji,jj) * EXP( - ze3t*r1_LG  )   ;   zzeB = zeB(ji,jj) * EXP( - ze3t*r1_LB )   ! Green ; Blue      -      -
            zzeT = ( zze0 + zzeB + zzeG + zzeR ) * wmask(ji,jj,jk+1)                                 ! Total             -      -
!!st01            zzeT = ( zze0 + zzeR + zzeG + zzeB ) * wmask(ji,jj,jk+1)                                 ! Total             -      -
            !
#if defined key_RK3
            !                                      !- RK3 : temperature trend at jk t-level
            pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs) + r1_rho0_rcp * ( zeT(ji,jj) - zzeT ) / ze3t
#else
            !                                      !- MLF : heat content trend due to Qsr flux (qsr_hc)
            qsr_hc(ji,jj,jk) = r1_rho0_rcp * ( zeT(ji,jj) - zzeT )
#endif
            ze0(ji,jj) = zze0   ;   zeR(ji,jj) = zzeR           ! IR    ; Red  store at jk+1 w-level
            zeG(ji,jj) = zzeG   ;   zeB(ji,jj) = zzeB           ! Green ; Blue   -        -      -
            zeT(ji,jj) = zzeT                                   ! total          -        -      -
         END_2D
         !
      END DO
      !
      DO jk = nk0+1, nkR                  !* down to Red extinction *!   (< ~71 meters : RGB , IR removed from calculation)
         ik = MIN( jk , ipk )
         DO_2D( 0, 0, 0, 0 )
            !                                      !- inverse of RGB attenuation lengths
            zlogc     = zc(ji,jj,ik,0)
            zCb       = 0.768 + zlogc * ( 0.087 - zlogc * ( 0.179 + zlogc * 0.025 ) )
            zCmax     = 0.299 - zlogc * ( 0.289 - zlogc * 0.579 )
            zpsimax   = 0.6   - zlogc * ( 0.640 - zlogc * ( 0.021 + zlogc * 0.115 ) )
            ! zdelpsi = 0.710 + zlogc * ( 0.159 + zlogc * 0.021 )
            zCze   = zc(ji,jj,ik,1)
            zrdpsi = zc(ji,jj,ik,2)                               ! 1/zdelpsi
            zpsi   = zc(ji,jj,ik,3) * gdepw(ji,jj,jk+1,Kmm)         ! gdepw/zze
!!st05            zpsi   = zc(ji,jj,ik,3) * gdepw(ji,jj,jk,Kmm)         ! gdepw/zze
            !                                                  ! make sure zchl value is such that: 0.03 < zchl < 10. 
            zchl = MAX(  0.03_wp , MIN( zCze * ( zCb + zCmax * EXP( -( (zpsi - zpsimax) * zrdpsi )**2 ) ) , 10._wp )  )
            !                                                  ! Convert chlorophyll value to attenuation coefficient
            itab = NINT( 41 + 20.*LOG10(zchl) + 1.e-15 )       ! look-up table index
            !       Red             !         Green              !         Blue
            r1_LR = rktab(3,itab)   ;   r1_LG = rktab(2,itab)    ;   r1_LB = rktab(1,itab)
            !
            !                                      !- fluxes at jk+1 w-level
            ze3t = e3t(ji,jj,jk,Kmm)
            zzeR = zeR(ji,jj) * EXP( - ze3t*r1_LR )                                                 ! Red          at jk+1 w-level
            zzeG = zeG(ji,jj) * EXP( - ze3t*r1_LG )   ;   zzeB = zeB(ji,jj) * EXP( - ze3t*r1_LB )   ! Green ; Blue      -      -
            zzeT = ( zzeR + zzeG + zzeB ) * wmask(ji,jj,jk+1)                                       ! Total             -      -
            !
#if defined key_RK3
            !                                      !- RK3 : temperature trend at jk t-level
            pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs) + r1_rho0_rcp * ( zeT(ji,jj) - zzeT ) / ze3t
#else
            !                                      !- MLF : heat content trend due to Qsr flux (qsr_hc)
            qsr_hc(ji,jj,jk) = r1_rho0_rcp * ( zeT(ji,jj) - zzeT )
#endif
            zeR(ji,jj) = zzeR                                  ! Red          store at jk+1 w-level
            zeG(ji,jj) = zzeG   ;   zeB(ji,jj) = zzeB          ! Green ; Blue   -        -      -
            zeT(ji,jj) = zzeT                                  ! total          -        -      -
         END_2D
      END DO
      !
      DO jk = nkR+1, nkG                  !* down to Green extinction *!   (< ~350 m : GB , IR+R removed from calculation)
         ik = MIN( jk , ipk )
         DO_2D( 0, 0, 0, 0 )
            !                                      !- inverse of RGB attenuation lengths
            zlogc     = zc(ji,jj,ik,0)
            zCb       = 0.768 + zlogc * ( 0.087 - zlogc * ( 0.179 + zlogc * 0.025 ) )
            zCmax     = 0.299 - zlogc * ( 0.289 - zlogc * 0.579 )
            zpsimax   = 0.6   - zlogc * ( 0.640 - zlogc * ( 0.021 + zlogc * 0.115 ) )
            ! zdelpsi = 0.710 + zlogc * ( 0.159 + zlogc * 0.021 )
            zCze   = zc(ji,jj,ik,1)
            zrdpsi = zc(ji,jj,ik,2)                               ! 1/zdelpsi
            zpsi   = zc(ji,jj,ik,3) * gdepw(ji,jj,jk+1,Kmm)         ! gdepw/zze
!!st05            zpsi   = zc(ji,jj,ik,3) * gdepw(ji,jj,jk,Kmm)         ! gdepw/zze
            !                                                  ! make sure zchl value is such that: 0.03 < zchl < 10. 
            zchl = MAX(  0.03_wp , MIN( zCze * ( zCb + zCmax * EXP( -( (zpsi - zpsimax) * zrdpsi )**2 ) ) , 10._wp )  )
            !                                                  ! Convert chlorophyll value to attenuation coefficient
            itab = NINT( 41 + 20.*LOG10(zchl) + 1.e-15 )       ! look-up table index
            !     Green              !         Blue
            r1_LG = rktab(2,itab)    ;   r1_LB = rktab(1,itab)
            !
            !                                      !- fluxes at jk+1 w-level
            ze3t = e3t(ji,jj,jk,Kmm)
            zzeG = zeG(ji,jj) * EXP( - ze3t * r1_LG )   ;   zzeB = zeB(ji,jj) * EXP( - ze3t * r1_LB ) ! Green ; Blue
            zzeT = ( zzeG + zzeB ) * wmask(ji,jj,jk+1)                                                ! Total             -      -
#if defined key_RK3
            !                                      !- RK3 : temperature trend at jk t-level
            pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs) + r1_rho0_rcp * ( zeT(ji,jj) - zzeT ) / ze3t
#else
            !                                      !- MLF : heat content trend due to Qsr flux (qsr_hc)
            qsr_hc(ji,jj,jk) = r1_rho0_rcp * ( zeT(ji,jj) - zzeT )
#endif
            zeG(ji,jj) = zzeG   ;   zeB(ji,jj) = zzeB          ! Green ; Blue store at jk+1 w-level
            zeT(ji,jj) = zzeT                                  ! total          -        -      -
         END_2D
      END DO
      !
      DO jk = nkG+1, nkB                  !* down to Blue extinction *!   (< ~1300 m : B , IR+RG removed from calculation)
         ik = MIN( jk , ipk )
         DO_2D( 0, 0, 0, 0 )
            !                                      !- inverse of RGB attenuation lengths
            zlogc     = zc(ji,jj,ik,0)
            zCb       = 0.768 + zlogc * ( 0.087 - zlogc * ( 0.179 + zlogc * 0.025 ) )
            zCmax     = 0.299 - zlogc * ( 0.289 - zlogc * 0.579 )
            zpsimax   = 0.6   - zlogc * ( 0.640 - zlogc * ( 0.021 + zlogc * 0.115 ) )
            ! zdelpsi = 0.710 + zlogc * ( 0.159 + zlogc * 0.021 )
            zCze   = zc(ji,jj,ik,1)
            zrdpsi = zc(ji,jj,ik,2)                               ! 1/zdelpsi
            zpsi   = zc(ji,jj,ik,3) * gdepw(ji,jj,jk+1,Kmm)         ! gdepw/zze
!!st05            zpsi   = zc(ji,jj,ik,3) * gdepw(ji,jj,jk,Kmm)         ! gdepw/zze
            !                                                  ! make sure zchl value is such that: 0.03 < zchl < 10. 
            zchl = MAX(  0.03_wp , MIN( zCze * ( zCb + zCmax * EXP( -( (zpsi - zpsimax) * zrdpsi )**2 ) ) , 10._wp )  )
            !                                                  ! Convert chlorophyll value to attenuation coefficient
            itab = NINT( 41 + 20.*LOG10(zchl) + 1.e-15 )       ! look-up table index
            r1_LB = rktab(1,itab)                              ! Blue
            !
            !                                      !- fluxes at jk+1 w-level
            ze3t = e3t(ji,jj,jk,Kmm)
            zzeB = zeB(ji,jj) * EXP( - ze3t * r1_LB )          ! Blue
            zzeT = ( zzeB ) * wmask(ji,jj,jk+1)                ! Total             -      -
#if defined key_RK3
            !                                      !- RK3 : temperature trend at jk t-level
            pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs) + r1_rho0_rcp * ( zeT(ji,jj) - zzeT ) / ze3t
#else
            !                                      !- MLF : heat content trend due to Qsr flux (qsr_hc)
            qsr_hc(ji,jj,jk) = r1_rho0_rcp * ( zeT(ji,jj) - zzeT )
#endif
            zeB(ji,jj) = zzeB                                  ! Blue store at jk+1 w-level
            zeT(ji,jj) = zzeT                                  ! total  -        -      -
         END_2D
      END DO
      !
#if ! defined key_PSYCLONE_2p5p0
      DEALLOCATE( zc )
#endif
      !
   END SUBROUTINE qsr_RGBc


   SUBROUTINE qsr_RGB( kt, Kmm, pts, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE qsr_RGB  ***
      !!
      !! ** Purpose :   Red-Green-Blue solar radiation with constant chlorophyll
      !!
      !! ** Method  : The profile of the solar radiation within the ocean is defined
      !!      through 2 wavebands (rn_si0,rn_si1) or 1 (rn_si0,rn_abs) + 3 wavebands (RGB) 
      !!         At the bottom, boudary condition for the radiation is no flux :
      !!      all heat which has not been absorbed in the above levels is put
      !!      in the last ocean level.
      !!         For each band, the computation is only done down to the level where
      !!      I(k) < 1.e-15 W/m2 (i.e. over the top nk levels) .
      !!
      !! ** Action  : - update ta with the penetrative solar radiation trend
      !!              - send  trend for further diagnostics (l_trdtra=T)
      !!
      !! Reference  : Lengaigne et al. 2007, Clim. Dyn., V28, 5, 503-516.
      !!              Morel, A. et Berthon, JF, 1989, Limnol Oceanogr 34(8), 1545-1562
      !!----------------------------------------------------------------------
      INTEGER,                                   INTENT(in   ) ::   kt, Kmm, Krhs   ! ocean time-step and time-level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts,jpt), INTENT(inout) ::   pts             ! active tracers and RHS of tracer equation
      !!
      INTEGER  ::   ji, jj, jk               ! dummy loop indices
      REAL(wp) ::   zze0, zzeR, zzeG, zzeB, zzeT              !    -         -
      REAL(wp) ::   zz0 , zz1 , ze3t                          !    -         -
      REAL(wp), DIMENSION(T2D(0))   ::   ze0, zeR, zeG, zeB, zeT
      !!----------------------------------------------------------------------
      !      
      !
      !                       !==============================================!
      !                       !==  R-G-B fluxes with constant chlorophyll  ==!   
      !                       !======================********================!
      !
      !                             !=  surface light  =!
      !
      zz0 =           rn_abs              ! Infrared absorption
      zz1 = ( 1._wp - rn_abs ) / 3._wp    ! surface equi-partition in R-G-B
      !
      DO_2D( 0, 0, 0, 0 )                 ! surface light
         ze0(ji,jj) = zz0 * qsr(ji,jj)   ;   zeR(ji,jj) = zz1 * qsr(ji,jj)    ! IR    ; Red
         zeG(ji,jj) = zz1 * qsr(ji,jj)   ;   zeB(ji,jj) = zz1 * qsr(ji,jj)    ! Green ; Blue
         zeT(ji,jj) =       qsr(ji,jj)                                        ! Total
      END_2D
      !
      !                             !=  interior light  =!
      !
      DO jk = 1, nk0                      !* near surface layers *!   (< ~12 meters : IR + RGB )
          DO_2D( 0, 0, 0, 0 )
            ze3t = e3t(ji,jj,jk,Kmm)
            zze0 = ze0(ji,jj) * EXP( - ze3t * r1_si0 )   ;   zzeR = zeR(ji,jj) * EXP( - ze3t * r1_LR )   ! IR    ; Red  at jk+1 w-level
            zzeG = zeG(ji,jj) * EXP( - ze3t * r1_LG  )   ;   zzeB = zeB(ji,jj) * EXP( - ze3t * r1_LB )   ! Green ; Blue      -      -
            zzeT = ( zze0 + zzeB + zzeG + zzeR ) * wmask(ji,jj,jk+1)                                     ! Total             -      -
!!st7-9            zzeT = ( zze0 + zzeR + zzeG + zzeB ) * wmask(ji,jj,jk+1)                                     ! Total             -      -
#if defined key_RK3
            !                                               ! RK3 : temperature trend at jk t-level
            pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs) + r1_rho0_rcp * ( zeT(ji,jj) - zzeT ) / ze3t
#else
            !                                               ! MLF : heat content trend due to Qsr flux (qsr_hc)
            qsr_hc(ji,jj,jk) = r1_rho0_rcp * ( zeT(ji,jj) - zzeT )
#endif
            ze0(ji,jj) = zze0   ;   zeR(ji,jj) = zzeR           ! IR    ; Red  store at jk+1 w-level
            zeG(ji,jj) = zzeG   ;   zeB(ji,jj) = zzeB           ! Green ; Blue   -        -      -
            zeT(ji,jj) = zzeT                                   ! total          -        -      -
         END_2D
!!stbug         IF( jk == 1 ) THEN               !* sea-ice *!   store the 1st level attenuation coeff.
!!stbug            WHERE( qsr(T2D(0)) /= 0._wp )   ;   fraqsr_1lev(T2D(0)) = 1._wp - zeT(T2D(0)) / qsr(T2D(0))
!!stbug            ELSEWHERE                       ;   fraqsr_1lev(T2D(0)) = 1._wp
!!stbug            END WHERE
!!stbug         ENDIF
      END DO
      !
      DO jk = nk0+1, nkR                  !* down to Red extinction *!   (< ~71 meters : RGB , IR removed from calculation)
          DO_2D( 0, 0, 0, 0 )
            ze3t = e3t(ji,jj,jk,Kmm)
            zzeR = zeR(ji,jj) * EXP( - ze3t * r1_LR )                                                 ! Red          at jk+1 w-level
            zzeG = zeG(ji,jj) * EXP( - ze3t * r1_LG )   ;   zzeB = zeB(ji,jj) * EXP( - ze3t * r1_LB ) ! Green ; Blue      -      -
            zzeT = ( zzeB + zzeG + zzeR ) * wmask(ji,jj,jk+1)                                         ! Total             -      -
!!st7-11            zzeT = ( zzeR + zzeG + zzeB ) * wmask(ji,jj,jk+1)                                         ! Total             -      -
#if defined key_RK3
            !                                               ! RK3 : temperature trend at jk t-level
            pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs) + r1_rho0_rcp * ( zeT(ji,jj) - zzeT ) / ze3t
#else
            !                                               ! MLF : heat content trend due to Qsr flux (qsr_hc)
            qsr_hc(ji,jj,jk) = r1_rho0_rcp * ( zeT(ji,jj) - zzeT )
#endif
            zeR(ji,jj) = zzeR                                   ! Red          store at jk+1 w-level
            zeG(ji,jj) = zzeG   ;   zeB(ji,jj) = zzeB           ! Green ; Blue   -        -      -
            zeT(ji,jj) = zzeT                                   ! total          -        -      -
         END_2D
      END DO
      !
      DO jk = nkR+1, nkG                  !* down to Green extinction *!   (< ~350 m : GB , IR+R removed from calculation)
         DO_2D( 0, 0, 0, 0 )
            ze3t = e3t(ji,jj,jk,Kmm)
            zzeG = zeG(ji,jj) * EXP( - ze3t * r1_LG )   ;   zzeB = zeB(ji,jj) * EXP( - ze3t * r1_LB ) ! Green ; Blue at jk+1 w-level
            zzeT = ( zzeG + zzeB ) * wmask(ji,jj,jk+1)                                                ! Total             -      -
#if defined key_RK3
            !                                               ! RK3 : temperature trend at jk t-level
            pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs) + r1_rho0_rcp * ( zeT(ji,jj) - zzeT ) / ze3t
#else
            !                                               ! MLF : heat content trend due to Qsr flux (qsr_hc)
            qsr_hc(ji,jj,jk) = r1_rho0_rcp * ( zeT(ji,jj) - zzeT )
#endif
            zeG(ji,jj) = zzeG   ;   zeB(ji,jj) = zzeB             ! Green ; Blue store at jk+1 w-level
            zeT(ji,jj) = zzeT                                     ! total          -        -      -
         END_2D
      END DO
      !
      DO jk = nkG+1, nkB                  !* down to Blue extinction *!   (< ~1300 m : B , IR+RG removed from calculation)
         DO_2D( 0, 0, 0, 0 )
            ze3t = e3t(ji,jj,jk,Kmm)
            zzeB = zeB(ji,jj) * EXP( - ze3t * r1_LB )             ! Blue at jk+1 w-level
            zzeT = ( zzeB ) * wmask(ji,jj,jk+1)                   ! Total     -      -
#if defined key_RK3
            !                                               ! RK3 : temperature trend at jk t-level
            pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs) + r1_rho0_rcp * ( zeT(ji,jj) - zzeT ) / ze3t
#else
            !                                               ! MLF : heat content trend due to Qsr flux (qsr_hc)
            qsr_hc(ji,jj,jk) = r1_rho0_rcp * ( zeT(ji,jj) - zzeT )
#endif
            zeB(ji,jj) = zzeB                                     ! Blue store at jk+1 w-level
            zeT(ji,jj) = zzeT                                     ! total  -        -      -
         END_2D
      END DO
      !
   END SUBROUTINE qsr_RGB


   SUBROUTINE qsr_2BD( Kmm, pts, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE qsr_2BD  ***
      !!
      !! ** Purpose :   2 bands (IR+visible) solar radiation with constant chlorophyll
      !!
      !! ** Method  : The profile of the solar radiation within the ocean is defined
      !!      through 2 wavebands (rn_si0,rn_si1) a ratio rn_abs for IR absorbtion.
      !!      Considering the 2 wavebands case:
      !!         I(k) = Qsr*( rn_abs*EXP(z(k)/rn_si0) + (1.-rn_abs)*EXP(z(k)/rn_si1) )
      !!         The temperature trend associated with the solar radiation penetration
      !!         is given by : zta = 1/e3t dk[ I ] / (rho0*Cp)
      !!         At the bottom, boudary condition for the radiation is no flux :
      !!      all heat which has not been absorbed in the above levels is put
      !!      in the last ocean level.
      !!         The computation is only done down to the level where
      !!      I(k) < 1.e-15 W/m2 (i.e. over the top nk levels) .
      !!
      !! ** Action  : - update ta with the penetrative solar radiation trend
      !!              - send  trend for further diagnostics (l_trdtra=T)
      !!
      !! Reference  : Jerlov, N. G., 1968 Optical Oceanography, Elsevier, 194pp.
      !!              Lengaigne et al. 2007, Clim. Dyn., V28, 5, 503-516.
      !!              Morel, A. et Berthon, JF, 1989, Limnol Oceanogr 34(8), 1545-1562
      !!----------------------------------------------------------------------
      INTEGER,                                   INTENT(in   ) ::   Kmm, Krhs   ! ocean time-step and time-level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts,jpt), INTENT(inout) ::   pts         ! active tracers and RHS of tracer equation
      !!
      INTEGER  ::   ji, jj, jk               ! dummy loop indices
      REAL(wp) ::   zzatt                    !    -         -
      REAL(wp) ::   zz0 , zz1 , ze3t         !    -         -
      REAL(wp), DIMENSION(T2D(0)) ::   zatt
      !!----------------------------------------------------------------------
      !      
      !                       !======================!
      !                       !==  2-bands fluxes  ==!
      !                       !======================!
      !
      zz0 =           rn_abs   * r1_rho0_rcp       ! surface equi-partition in 2-bands
      zz1 = ( 1._wp - rn_abs ) * r1_rho0_rcp
      !
      zatt(T2D(0)) = r1_rho0_rcp                   !* surface value *!
      !
      DO_2D( 0, 0, 0, 0 )
         zatt(ji,jj) = (  zz0 * EXP( -gdepw(ji,jj,1,Kmm)*r1_si0 ) + zz1 * EXP( -gdepw(ji,jj,1,Kmm)*r1_si1 )  )
      END_2D
      !
!!st      IF(lwp) WRITE(numout,*) 'level = ', 1, ' qsr max = ' , MAXVAL(zatt)*rho0_rcp, ' W/m2', ' qsr min = ' , MINVAL(zatt)*rho0_rcp, ' W/m2' 
      !
      DO jk = 1, nk0                               !* near surface layers *!   (< ~14 meters : IR + visible light )
         DO_2D( 0, 0, 0, 0 )
            ze3t  = e3t(ji,jj,jk,Kmm)                    ! light attenuation at jk+1 w-level (divided by rho0_rcp)
            zzatt = (   zz0 * EXP( -gdepw(ji,jj,jk+1,Kmm)*r1_si0 )     &
               &      + zz1 * EXP( -gdepw(ji,jj,jk+1,Kmm)*r1_si1 )   ) * wmask(ji,jj,jk+1)
#if defined key_RK3
            !                                            ! RK3 : temperature trend at jk t-level
            pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs) + qsr(ji,jj) * ( zatt(ji,jj) - zzatt ) / ze3t
#else
            !                                            ! MLF : heat content trend due to Qsr flux (qsr_hc)
            qsr_hc(ji,jj,jk) =  qsr(ji,jj) * ( zatt(ji,jj) - zzatt )
#endif
            zatt(ji,jj) = zzatt                          ! save for the next level computation
         END_2D
!!stbug         !                                         !* sea-ice *!   store the 1st level attenuation coeff.
!!stbug         IF( jk == 1 )   fraqsr_1lev(T2D(0)) = 1._wp - zatt(T2D(0)) * rho0_rcp
      END DO
!!st      IF(lwp) WRITE(numout,*) 'nk0+1= ', nk0+1, ' qsr max = ' , MAXVAL(zatt*qsr)*rho0_rcp, ' W/m2' , MAXVAL(zatt*qsr/e3t(:,:,nk0+1,Kmm)), ' K/s' 
      !
      DO jk = nk0+1, nkV                           !* deeper layers *!   (visible light only)
         DO_2D( 0, 0, 0, 0 )
            ze3t  = e3t(ji,jj,jk,Kmm)                    ! light attenuation at jk+1 w-level (divided by rho0_rcp)
            zzatt = (   zz1 * EXP( -gdepw(ji,jj,jk+1,Kmm)*r1_si1 )   ) * wmask(ji,jj,jk+1)
#if defined key_RK3
            !                                            ! RK3 : temperature trend at jk t-level
            pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs) + qsr(ji,jj) * ( zatt(ji,jj) - zzatt ) / ze3t
#else
            !                                            ! MLF : heat content trend due to Qsr flux (qsr_hc)
            qsr_hc(ji,jj,jk) = qsr(ji,jj) * ( zatt(ji,jj) - zzatt )
#endif
            zatt(ji,jj) = zzatt                       ! save for the next level computation
         END_2D
      END DO      
      !
!!st      IF(lwp) WRITE(numout,*) 'nkV+1= ', nkV+1, ' qsr max = ' , MAXVAL(zatt*qsr)*rho0_rcp, ' W/m2' , MAXVAL(zatt*qsr/e3t(:,:,nkV+1,Kmm)), ' K/s' 
   END SUBROUTINE qsr_2bd


   SUBROUTINE qsr_5BDc( kt, Kmm, pts, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE qsr_5BDc  ***
      !!
      !! ** Purpose :   IR-Red-Green-Blue-UV solar radiation using chlorophyll data
      !!
      !! ** Method  : The profile of the solar radiation within the ocean is defined
      !!      through 2 wavebands (rn_si0,rn_si1) or 3 wavebands (RGB) and a ratio rn_abs
      !!      Considering the 2 wavebands case:
      !!         I(k) = Qsr*( rn_abs*EXP(z(k)/rn_si0) + (1.-rn_abs)*EXP(z(k)/rn_si1) )
      !!         The temperature trend associated with the solar radiation penetration
      !!         is given by : zta = 1/e3t dk[ I ] / (rho0*Cp)
      !!         At the bottom, boudary condition for the radiation is no flux :
      !!      all heat which has not been absorbed in the above levels is put
      !!      in the last ocean level.
      !!         The computation is only done down to the level where
      !!      I(k) < 1.e-15 W/m2 (i.e. over the top nk levels) .
      !!
      !! ** Action  : - update ta with the penetrative solar radiation trend
      !!              - send  trend for further diagnostics (l_trdtra=T)
      !!
      !! Reference  : Lengaigne et al. 2007, Clim. Dyn., V28, 5, 503-516.
      !!              Morel, A. et Berthon, JF, 1989, Limnol Oceanogr 34(8), 1545-1562
      !!----------------------------------------------------------------------
      INTEGER,                                   INTENT(in   ) ::   kt, Kmm, Krhs   ! ocean time-step and time-level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts,jpt), INTENT(inout) ::   pts         ! active tracers and RHS of tracer equation
      !!
      INTEGER  ::   ji, jj, jk, ik           ! dummy loop indices
      INTEGER  ::   ipk, itab                ! local integer
      REAL(wp) ::   zc1 , zc2 , zc3, zchl    ! local scalars
      REAL(wp) ::   zze0, zzeR, zzeG, zzeB, zzeU, zzeT        !    -         -
      REAL(wp) ::   zz0 , zz1 , zz2, ze3t                     !    -         -
      REAL(wp) ::   zCb, zCmax, zpsi, zpsimax, zrdpsi, zCze   !    -         -
      REAL(wp) ::   zlogc, zlogze, zlogCtot, zlogCze          !    -         -
      !!
      REAL(wp), DIMENSION(T2D(0)) ::   ze0, zeR, zeG, zeB, zeU, zeT
      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::   zc
      !!----------------------------------------------------------------------
      !
      IF    ( nn_chldta == 1 ) THEN   ;   ipk=1
      ELSEIF( nn_chldta == 2 ) THEN   ;   ipk=jpk   ;   ENDIF
      ALLOCATE( zc(T2D(0),ipk,0:3) )
      !
      !                       !===========================================!
      !                       !==  R-G-B fluxes using chlorophyll data  ==!    with Morel &Berthon (1989) vertical profile
      !                       !===================================****====!
      !
      !                             !=  Chlorophyll data  =!
      !
      IF( ntile == 0 .OR. ntile == 1 )  THEN       ! Do only for the full domain
         IF( ln_tile )   CALL dom_tile( ntsi, ntsj, ntei, ntej, ktile = 0 )   ! Use full domain
         CALL fld_read( kt, 1, sf_chl )                                       ! Read Chl data and provides it at the current time step
         IF( ln_tile )   CALL dom_tile( ntsi, ntsj, ntei, ntej, ktile = 1 )   ! Revert to tile domain
      ENDIF
      !
      DO_3D( 0, 0, 0, 0, 1, ipk )                          ! pre-calculated expensive coefficient
         zlogc = LOG(  MAX( 0.03_wp, MIN( sf_chl(1)%fnow(ji,jj,jk) ,10._wp ) )  ) ! zlogc = log(zchl)   with 0.03 <= Chl >= 10. 
         zc1   = 0.113328685307 + 0.803 * zlogc                                   ! zc1 : log(zCze)  = log (1.12  * zchl**0.803)
         zc2   = 3.703768066608 + 0.459 * zlogc                                   ! zc2 : log(zCtot) = log(40.6  * zchl**0.459)
         zc3   = 6.34247346942  - 0.746 * zc2                                     ! zc3 : log(zze)   = log(568.2 * zCtot**(-0.746))
         IF( zc3 > 4.62497281328 )   zc3 = 5.298317366548 - 0.293 * zc2           ! IF(log(zze)>log(102)) log(zze) = log(200*zCtot**(-0.293))
         !
         zc(ji,jj,jk,0) = zlogc                                                   ! ze(0) = log(zchl)
         zc(ji,jj,jk,1) = EXP( zc1 )                                              ! ze(1) = zCze
         zc(ji,jj,jk,2) = 1._wp / ( 0.710 + zlogc * ( 0.159 + zlogc * 0.021 ) )   ! ze(2) = 1/zdelpsi
         zc(ji,jj,jk,3) = EXP( - zc3 )                                            ! ze(3) = 1/zze
      END_3D
      !
      !                             !=  surface light  =!
      !
      zz0 =           rn_abs              ! Infrared absorption
      zz1 =         ( rn_par ) / 3._wp    ! R-G-B equi-partition
      zz2 =    ( 1._wp - rn_abs - rn_par) ! uv partition
      !
      DO_2D( 0, 0, 0, 0 )                 ! surface light
         ze0(ji,jj) = zz0 * qsr(ji,jj)   ;   zeR(ji,jj) = zz1 * qsr(ji,jj)    ! IR    ; Red
         zeG(ji,jj) = zz1 * qsr(ji,jj)   ;   zeB(ji,jj) = zz1 * qsr(ji,jj)    ! Green ; Blue
         zeU(ji,jj) = zz2 * qsr(ji,jj)   ;   zeT(ji,jj) =       qsr(ji,jj)    ! UV    ; Total
      END_2D
      !              
      !                             !=  interior light  =!
      !
      DO jk = 1, nk0                      !* near surface layers *!   (< ~12 meters : IR + RGB + UV)
         ik = MIN( jk , ipk )
         DO_2D( 0, 0, 0, 0 )
            !                                      !- inverse of RGB attenuation lengths
            zlogc     = zc(ji,jj,ik,0)
            zCb       = 0.768 + zlogc * ( 0.087 - zlogc * ( 0.179 + zlogc * 0.025 ) )
            zCmax     = 0.299 - zlogc * ( 0.289 - zlogc * 0.579 )
            zpsimax   = 0.6   - zlogc * ( 0.640 - zlogc * ( 0.021 + zlogc * 0.115 ) )
            ! zdelpsi = 0.710 + zlogc * ( 0.159 + zlogc * 0.021 )
            zCze   = zc(ji,jj,ik,1)
            zrdpsi = zc(ji,jj,ik,2)                                     ! 1/zdelpsi
!!st05            zpsi   = zc(ji,jj,ik,3) * gdepw(ji,jj,jk,Kmm)               ! gdepw/zze
            zpsi   = zc(ji,jj,ik,3) * gdepw(ji,jj,jk+1,Kmm)               ! gdepw/zze
            !                                                        ! make sure zchl value is such that: 0.03 < zchl < 10. 
            zchl = MAX(  0.03_wp , MIN( zCze * ( zCb + zCmax * EXP( -( (zpsi - zpsimax) * zrdpsi )**2 ) ) , 10._wp )  )
            !                                                        ! Convert chlorophyll value to attenuation coefficient
            itab = NINT( 41 + 20.*LOG10(zchl) + 1.e-15 )             ! look-up table index
            !       Red             !         Green              !         Blue
            r1_LR = rktab(3,itab)   ;   r1_LG = rktab(2,itab)    ;   r1_LB = rktab(1,itab)    ;   r1_LU = rktab(4,itab)
            !
            !                                      !- fluxes at jk+1 w-level
            ze3t = e3t(ji,jj,jk,Kmm)
            zze0 = ze0(ji,jj) * EXP( - ze3t*r1_si0 )   ;   zzeR = zeR(ji,jj) * EXP( - ze3t*r1_LR )   ! IR    ; Red  at jk+1 w-level
            zzeG = zeG(ji,jj) * EXP( - ze3t*r1_LG  )   ;   zzeB = zeB(ji,jj) * EXP( - ze3t*r1_LB )   ! Green ; Blue      -      -
            zzeU = zeU(ji,jj) * EXP( - ze3t*r1_LU  )                                                 !    UV             -      -
            zzeT = ( zze0 + zzeU + zzeB + zzeG + zzeR ) * wmask(ji,jj,jk+1)                          ! Total             -      -
            !
#if defined key_RK3
            !                                      !- RK3 : temperature trend at jk t-level
            pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs) + r1_rho0_rcp * ( zeT(ji,jj) - zzeT ) / ze3t
#else
            !                                      !- MLF : heat content trend due to Qsr flux (qsr_hc)
            qsr_hc(ji,jj,jk) = r1_rho0_rcp * ( zeT(ji,jj) - zzeT )
#endif
            ze0(ji,jj) = zze0   ;   zeR(ji,jj) = zzeR           ! IR    ; Red  store at jk+1 w-level
            zeG(ji,jj) = zzeG   ;   zeB(ji,jj) = zzeB           ! Green ; Blue   -        -      -
            zeU(ji,jj) = zzeU   ;   zeT(ji,jj) = zzeT           ! UV    ; total  -        -      -
         END_2D
         !
      END DO
      !
      DO jk = nk0+1, nkR                  !* down to Red extinction *!   (< ~71 meters : RGB + UV, IR removed from calculation)
         ik = MIN( jk , ipk )
         DO_2D( 0, 0, 0, 0 )
            !                                      !- inverse of RGB attenuation lengths
            zlogc     = zc(ji,jj,ik,0)
            zCb       = 0.768 + zlogc * ( 0.087 - zlogc * ( 0.179 + zlogc * 0.025 ) )
            zCmax     = 0.299 - zlogc * ( 0.289 - zlogc * 0.579 )
            zpsimax   = 0.6   - zlogc * ( 0.640 - zlogc * ( 0.021 + zlogc * 0.115 ) )
            ! zdelpsi = 0.710 + zlogc * ( 0.159 + zlogc * 0.021 )
            zCze   = zc(ji,jj,ik,1)
            zrdpsi = zc(ji,jj,ik,2)                               ! 1/zdelpsi
            zpsi   = zc(ji,jj,ik,3) * gdepw(ji,jj,jk+1,Kmm)         ! gdepw/zze
!!st05            zpsi   = zc(ji,jj,ik,3) * gdepw(ji,jj,jk,Kmm)         ! gdepw/zze
            !                                                  ! make sure zchl value is such that: 0.03 < zchl < 10. 
            zchl = MAX(  0.03_wp , MIN( zCze * ( zCb + zCmax * EXP( -( (zpsi - zpsimax) * zrdpsi )**2 ) ) , 10._wp )  )
            !                                                  ! Convert chlorophyll value to attenuation coefficient
            itab = NINT( 41 + 20.*LOG10(zchl) + 1.e-15 )       ! look-up table index
            !       Red             !         Green              !         Blue
            r1_LR = rktab(3,itab)   ;   r1_LG = rktab(2,itab)    ;   r1_LB = rktab(1,itab)    ;   r1_LU = rktab(4,itab)
            !
            !                                      !- fluxes at jk+1 w-level
            ze3t = e3t(ji,jj,jk,Kmm)
            zzeR = zeR(ji,jj) * EXP( - ze3t*r1_LR )    ;  zzeU = zeU(ji,jj) * EXP( - ze3t*r1_LU )   ! Red   ;   UV at jk+1 w-level
            zzeG = zeG(ji,jj) * EXP( - ze3t*r1_LG )   ;   zzeB = zeB(ji,jj) * EXP( - ze3t*r1_LB )   ! Green ; Blue      -      -
            zzeT = ( zzeR + zzeG + zzeB + zzeU ) * wmask(ji,jj,jk+1)                                ! Total             -      -
            !
#if defined key_RK3
            !                                      !- RK3 : temperature trend at jk t-level
            pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs) + r1_rho0_rcp * ( zeT(ji,jj) - zzeT ) / ze3t
#else
            !                                      !- MLF : heat content trend due to Qsr flux (qsr_hc)
            qsr_hc(ji,jj,jk) = r1_rho0_rcp * ( zeT(ji,jj) - zzeT )
#endif
            zeR(ji,jj) = zzeR   ;   zeU(ji,jj) = zzeU          ! Red   ;   UV store at jk+1 w-level
            zeG(ji,jj) = zzeG   ;   zeB(ji,jj) = zzeB          ! Green ; Blue   -        -      -
            zeT(ji,jj) = zzeT                                  ! total          -        -      -
         END_2D
      END DO
      !
      DO jk = nkR+1, nkG                  !* down to Green extinction *!   (< ~350 m : GB + UV , IR+R removed from calculation)
         ik = MIN( jk , ipk )
         DO_2D( 0, 0, 0, 0 )
            !                                      !- inverse of RGB attenuation lengths
            zlogc     = zc(ji,jj,ik,0)
            zCb       = 0.768 + zlogc * ( 0.087 - zlogc * ( 0.179 + zlogc * 0.025 ) )
            zCmax     = 0.299 - zlogc * ( 0.289 - zlogc * 0.579 )
            zpsimax   = 0.6   - zlogc * ( 0.640 - zlogc * ( 0.021 + zlogc * 0.115 ) )
            ! zdelpsi = 0.710 + zlogc * ( 0.159 + zlogc * 0.021 )
            zCze   = zc(ji,jj,ik,1)
            zrdpsi = zc(ji,jj,ik,2)                               ! 1/zdelpsi
            zpsi   = zc(ji,jj,ik,3) * gdepw(ji,jj,jk+1,Kmm)         ! gdepw/zze
!!st05            zpsi   = zc(ji,jj,ik,3) * gdepw(ji,jj,jk,Kmm)         ! gdepw/zze
            !                                                  ! make sure zchl value is such that: 0.03 < zchl < 10. 
            zchl = MAX(  0.03_wp , MIN( zCze * ( zCb + zCmax * EXP( -( (zpsi - zpsimax) * zrdpsi )**2 ) ) , 10._wp )  )
            !                                                  ! Convert chlorophyll value to attenuation coefficient
            itab = NINT( 41 + 20.*LOG10(zchl) + 1.e-15 )       ! look-up table index
            !     Green              !         Blue
            r1_LG = rktab(2,itab)    ;   r1_LB = rktab(1,itab)    ;   r1_LU = rktab(4,itab)
            !
            !                                      !- fluxes at jk+1 w-level
            ze3t = e3t(ji,jj,jk,Kmm)
            zzeG = zeG(ji,jj) * EXP( - ze3t * r1_LG )   ;   zzeB = zeB(ji,jj) * EXP( - ze3t * r1_LB ) ! Green ; Blue
            zzeU = zeU(ji,jj) * EXP( - ze3t * r1_LU )                                                 !    UV
            zzeT = ( zzeG + zzeB + zzeU ) * wmask(ji,jj,jk+1)                                         ! Total             -      -
#if defined key_RK3
            !                                      !- RK3 : temperature trend at jk t-level
            pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs) + r1_rho0_rcp * ( zeT(ji,jj) - zzeT ) / ze3t
#else
            !                                      !- MLF : heat content trend due to Qsr flux (qsr_hc)
            qsr_hc(ji,jj,jk) = r1_rho0_rcp * ( zeT(ji,jj) - zzeT )
#endif
            zeG(ji,jj) = zzeG   ;   zeB(ji,jj) = zzeB          ! Green ; Blue store at jk+1 w-level
            zeU(ji,jj) = zzeU   ;   zeT(ji,jj) = zzeT          !    UV ; total  -        -      -
         END_2D
      END DO
      !
      DO jk = nkG+1, nkU                  !* down to UV extinction *!   (< ~1000 m : B + UV, IR+RG removed from calculation)
         ik = MIN( jk , ipk )
         DO_2D( 0, 0, 0, 0 )
            !                                      !- inverse of RGB attenuation lengths
            zlogc     = zc(ji,jj,ik,0)
            zCb       = 0.768 + zlogc * ( 0.087 - zlogc * ( 0.179 + zlogc * 0.025 ) )
            zCmax     = 0.299 - zlogc * ( 0.289 - zlogc * 0.579 )
            zpsimax   = 0.6   - zlogc * ( 0.640 - zlogc * ( 0.021 + zlogc * 0.115 ) )
            ! zdelpsi = 0.710 + zlogc * ( 0.159 + zlogc * 0.021 )
            zCze   = zc(ji,jj,ik,1)
            zrdpsi = zc(ji,jj,ik,2)                               ! 1/zdelpsi
            zpsi   = zc(ji,jj,ik,3) * gdepw(ji,jj,jk+1,Kmm)         ! gdepw/zze
!!st05            zpsi   = zc(ji,jj,ik,3) * gdepw(ji,jj,jk,Kmm)         ! gdepw/zze
            !                                                  ! make sure zchl value is such that: 0.03 < zchl < 10. 
            zchl = MAX(  0.03_wp , MIN( zCze * ( zCb + zCmax * EXP( -( (zpsi - zpsimax) * zrdpsi )**2 ) ) , 10._wp )  )
            !                                                  ! Convert chlorophyll value to attenuation coefficient
            itab = NINT( 41 + 20.*LOG10(zchl) + 1.e-15 )       ! look-up table index
            r1_LB = rktab(1,itab)    ;   r1_LU = rktab(4,itab) ! Blue ; UV
            !
            !                                      !- fluxes at jk+1 w-level
            ze3t = e3t(ji,jj,jk,Kmm)
            zzeB = zeB(ji,jj) * EXP( - ze3t * r1_LB )          ! Blue
            zzeU = zeU(ji,jj) * EXP( - ze3t * r1_LU )          !   UV 
            zzeT = ( zzeB + zzeU ) * wmask(ji,jj,jk+1)         ! Total             -      -
#if defined key_RK3
            !                                      !- RK3 : temperature trend at jk t-level
            pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs) + r1_rho0_rcp * ( zeT(ji,jj) - zzeT ) / ze3t
#else
            !                                      !- MLF : heat content trend due to Qsr flux (qsr_hc)
            qsr_hc(ji,jj,jk) = r1_rho0_rcp * ( zeT(ji,jj) - zzeT )
#endif
            zeB(ji,jj) = zzeB  ;  zeU(ji,jj) = zzeU            ! Blue ; UV store at jk+1 w-level
            zeT(ji,jj) = zzeT                                  ! total       -        -      -
         END_2D
      END DO
      !
      DO jk = nkU+1, nkB                  !* down to Blue extinction *!   (< ~1300 m : B, IR+RG+UV removed from calculation)
         ik = MIN( jk , ipk )
         DO_2D( 0, 0, 0, 0 )
            !                                      !- inverse of RGB attenuation lengths
            zlogc     = zc(ji,jj,ik,0)
            zCb       = 0.768 + zlogc * ( 0.087 - zlogc * ( 0.179 + zlogc * 0.025 ) )
            zCmax     = 0.299 - zlogc * ( 0.289 - zlogc * 0.579 )
            zpsimax   = 0.6   - zlogc * ( 0.640 - zlogc * ( 0.021 + zlogc * 0.115 ) )
            ! zdelpsi = 0.710 + zlogc * ( 0.159 + zlogc * 0.021 )
            zCze   = zc(ji,jj,ik,1)
            zrdpsi = zc(ji,jj,ik,2)                               ! 1/zdelpsi
            zpsi   = zc(ji,jj,ik,3) * gdepw(ji,jj,jk+1,Kmm)         ! gdepw/zze
!!st05            zpsi   = zc(ji,jj,ik,3) * gdepw(ji,jj,jk,Kmm)         ! gdepw/zze
            !                                                  ! make sure zchl value is such that: 0.03 < zchl < 10. 
            zchl = MAX(  0.03_wp , MIN( zCze * ( zCb + zCmax * EXP( -( (zpsi - zpsimax) * zrdpsi )**2 ) ) , 10._wp )  )
            !                                                  ! Convert chlorophyll value to attenuation coefficient
            itab = NINT( 41 + 20.*LOG10(zchl) + 1.e-15 )       ! look-up table index
            r1_LB = rktab(1,itab)                              ! Blue
            !
            !                                      !- fluxes at jk+1 w-level
            ze3t = e3t(ji,jj,jk,Kmm)
            zzeB = zeB(ji,jj) * EXP( - ze3t * r1_LB )          ! Blue
            zzeT = ( zzeB ) * wmask(ji,jj,jk+1)         ! Total             -      -
#if defined key_RK3
            !                                      !- RK3 : temperature trend at jk t-level
            pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs) + r1_rho0_rcp * ( zeT(ji,jj) - zzeT ) / ze3t
#else
            !                                      !- MLF : heat content trend due to Qsr flux (qsr_hc)
            qsr_hc(ji,jj,jk) = r1_rho0_rcp * ( zeT(ji,jj) - zzeT )
#endif
            zeB(ji,jj) = zzeB                                 ! Blue ; UV store at jk+1 w-level
            zeT(ji,jj) = zzeT                                 ! total       -        -      -
         END_2D
      END DO

      !
      DEALLOCATE( zc )
      !
   END SUBROUTINE qsr_5BDc


   SUBROUTINE qsr_5BD( kt, Kmm, pts, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE qsr_5BD  ***
      !!
      !! ** Purpose :   IR-Red-Green-Blue-UV solar radiation using chlorophyll data
      !!
      !! ** Method  : The profile of the solar radiation within the ocean is defined
      !!      through 2 wavebands (rn_si0,rn_si1) or 3 wavebands (RGB) and a ratio rn_abs
      !!      Considering the 2 wavebands case:
      !!         I(k) = Qsr*( rn_abs*EXP(z(k)/rn_si0) + (1.-rn_abs)*EXP(z(k)/rn_si1) )
      !!         The temperature trend associated with the solar radiation penetration
      !!         is given by : zta = 1/e3t dk[ I ] / (rho0*Cp)
      !!         At the bottom, boudary condition for the radiation is no flux :
      !!      all heat which has not been absorbed in the above levels is put
      !!      in the last ocean level.
      !!         The computation is only done down to the level where
      !!      I(k) < 1.e-15 W/m2 (i.e. over the top nk levels) .
      !!
      !! ** Action  : - update ta with the penetrative solar radiation trend
      !!              - send  trend for further diagnostics (l_trdtra=T)
      !!
      !! Reference  : Lengaigne et al. 2007, Clim. Dyn., V28, 5, 503-516.
      !!              Morel, A. et Berthon, JF, 1989, Limnol Oceanogr 34(8), 1545-1562
      !!----------------------------------------------------------------------
      INTEGER,                                   INTENT(in   ) ::   kt, Kmm, Krhs   ! ocean time-step and time-level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts,jpt), INTENT(inout) ::   pts         ! active tracers and RHS of tracer equation
      !!
      INTEGER  ::   ji, jj, jk               ! dummy loop indices
      REAL(wp) ::   zze0, zzeR, zzeG, zzeB, zzeU, zzeT        !    -         -
      REAL(wp) ::   zz0 , zz1 , zz2 , ze3t                    !    -         -
      REAL(wp), DIMENSION(T2D(0))   ::   ze0, zeR, zeG, zeB, zeU, zeT
      !!----------------------------------------------------------------------
      !
      !
      !                       !===================================================!
      !                       !==  IR-R-G-B-UV fluxes with constant chlorophyll ==!
      !                       !============================********===============!
      !
      !                             !=  surface light  =!
      !
      zz0 =           rn_abs              ! Infrared absorption
      zz1 =         ( rn_par ) / 3._wp    ! R-G-B equi-partition
      zz2 =    ( 1._wp - rn_abs - rn_par) ! uv partition
      !
      DO_2D( 0, 0, 0, 0 )                 ! surface light
         ze0(ji,jj) = zz0 * qsr(ji,jj)   ;   zeR(ji,jj) = zz1 * qsr(ji,jj)    ! IR    ; Red
         zeG(ji,jj) = zz1 * qsr(ji,jj)   ;   zeB(ji,jj) = zz1 * qsr(ji,jj)    ! Green ; Blue
         zeU(ji,jj) = zz2 * qsr(ji,jj)   ;   zeT(ji,jj) =       qsr(ji,jj)    ! UV    ; Total
      END_2D
      !              
      !                             !=  interior light  =!
      !
      DO jk = 1, nk0                      !* near surface layers *!   (< ~12 meters : IR + RGB + UV)
         DO_2D( 0, 0, 0, 0 )
            ze3t = e3t(ji,jj,jk,Kmm)
            zze0 = ze0(ji,jj) * EXP( - ze3t*r1_si0 )   ;   zzeR = zeR(ji,jj) * EXP( - ze3t*r1_LR )   ! IR    ; Red  at jk+1 w-level
            zzeG = zeG(ji,jj) * EXP( - ze3t*r1_LG  )   ;   zzeB = zeB(ji,jj) * EXP( - ze3t*r1_LB )   ! Green ; Blue      -      -
            zzeU = zeU(ji,jj) * EXP( - ze3t*r1_LU )                                                  !    UV             -      -
            zzeT = ( zze0 + zzeU + zzeB + zzeG + zzeR ) * wmask(ji,jj,jk+1)                          ! Total             -      -
            !
#if defined key_RK3
            !                                      !- RK3 : temperature trend at jk t-level
            pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs) + r1_rho0_rcp * ( zeT(ji,jj) - zzeT ) / ze3t
#else
            !                                      !- MLF : heat content trend due to Qsr flux (qsr_hc)
            qsr_hc(ji,jj,jk) = r1_rho0_rcp * ( zeT(ji,jj) - zzeT )
#endif
            ze0(ji,jj) = zze0   ;   zeR(ji,jj) = zzeR           ! IR    ; Red  store at jk+1 w-level
            zeG(ji,jj) = zzeG   ;   zeB(ji,jj) = zzeB           ! Green ; Blue   -        -      -
            zeU(ji,jj) = zzeU   ;   zeT(ji,jj) = zzeT           ! UV    ; total  -        -      -
         END_2D
         !
      END DO
      !
      DO jk = nk0+1, nkR                  !* down to Red extinction *!   (< ~71 meters : RGB + UV, IR removed from calculation)
         DO_2D( 0, 0, 0, 0 )
            ze3t = e3t(ji,jj,jk,Kmm)
            zzeR = zeR(ji,jj) * EXP( - ze3t*r1_LR )   ;   zzeU = zeU(ji,jj) * EXP( - ze3t*r1_LU )   ! Red   ;   UV at jk+1 w-level
            zzeG = zeG(ji,jj) * EXP( - ze3t*r1_LG )   ;   zzeB = zeB(ji,jj) * EXP( - ze3t*r1_LB )   ! Green ; Blue      -      -
            zzeT = ( zzeR + zzeG + zzeB + zzeU ) * wmask(ji,jj,jk+1)                                ! Total             -      -
            !
#if defined key_RK3
            !                                      !- RK3 : temperature trend at jk t-level
            pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs) + r1_rho0_rcp * ( zeT(ji,jj) - zzeT ) / ze3t
#else
            !                                      !- MLF : heat content trend due to Qsr flux (qsr_hc)
            qsr_hc(ji,jj,jk) = r1_rho0_rcp * ( zeT(ji,jj) - zzeT )
#endif
            zeR(ji,jj) = zzeR   ;   zeU(ji,jj) = zzeU          ! Red   ;   UV store at jk+1 w-level
            zeG(ji,jj) = zzeG   ;   zeB(ji,jj) = zzeB          ! Green ; Blue   -        -      -
            zeT(ji,jj) = zzeT                                  ! total          -        -      -
         END_2D
      END DO
      !
      DO jk = nkR+1, nkG                  !* down to Green extinction *!   (< ~350 m : GB + UV , IR+R removed from calculation)
         DO_2D( 0, 0, 0, 0 )
            ze3t = e3t(ji,jj,jk,Kmm)
            zzeG = zeG(ji,jj) * EXP( - ze3t * r1_LG )   ;   zzeB = zeB(ji,jj) * EXP( - ze3t * r1_LB ) ! Green ; Blue
            zzeU = zeU(ji,jj) * EXP( - ze3t * r1_LU )                                                 !    UV
            zzeT = ( zzeG + zzeB + zzeU ) * wmask(ji,jj,jk+1)                                         ! Total             -      -
#if defined key_RK3
            !                                      !- RK3 : temperature trend at jk t-level
            pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs) + r1_rho0_rcp * ( zeT(ji,jj) - zzeT ) / ze3t
#else
            !                                      !- MLF : heat content trend due to Qsr flux (qsr_hc)
            qsr_hc(ji,jj,jk) = r1_rho0_rcp * ( zeT(ji,jj) - zzeT )
#endif
            zeG(ji,jj) = zzeG   ;   zeB(ji,jj) = zzeB          ! Green ; Blue store at jk+1 w-level
            zeU(ji,jj) = zzeU   ;   zeT(ji,jj) = zzeT          !    UV ; total  -        -      -
         END_2D
      END DO
      !
      DO jk = nkG+1, nkU                  !* down to UV extinction *!   (< ~1000 m : B + UV, IR+RG removed from calculation)
         DO_2D( 0, 0, 0, 0 )
            ze3t = e3t(ji,jj,jk,Kmm)
            zzeB = zeB(ji,jj) * EXP( - ze3t * r1_LB )          ! Blue
            zzeU = zeU(ji,jj) * EXP( - ze3t * r1_LU )          !   UV 
            zzeT = ( zzeB + zzeU ) * wmask(ji,jj,jk+1)         ! Total             -      -
#if defined key_RK3
            !                                      !- RK3 : temperature trend at jk t-level
            pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs) + r1_rho0_rcp * ( zeT(ji,jj) - zzeT ) / ze3t
#else
            !                                      !- MLF : heat content trend due to Qsr flux (qsr_hc)
            qsr_hc(ji,jj,jk) = r1_rho0_rcp * ( zeT(ji,jj) - zzeT )
#endif
            zeB(ji,jj) = zzeB  ;  zeU(ji,jj) = zzeU            ! Blue ; UV store at jk+1 w-level
            zeT(ji,jj) = zzeT                                  ! total       -        -      -
         END_2D
      END DO
      !
      DO jk = nkU+1, nkB                  !* down to Blue extinction *!   (< ~1300 m : B, IR+RG+UV removed from calculation)
         DO_2D( 0, 0, 0, 0 )
            ze3t = e3t(ji,jj,jk,Kmm)
            zzeB = zeB(ji,jj) * EXP( - ze3t * r1_LB )          ! Blue
            zzeT = ( zzeB ) * wmask(ji,jj,jk+1)         ! Total             -      -
#if defined key_RK3
            !                                      !- RK3 : temperature trend at jk t-level
            pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs) + r1_rho0_rcp * ( zeT(ji,jj) - zzeT ) / ze3t
#else
            !                                      !- MLF : heat content trend due to Qsr flux (qsr_hc)
            qsr_hc(ji,jj,jk) = r1_rho0_rcp * ( zeT(ji,jj) - zzeT )
#endif
            zeB(ji,jj) = zzeB                                 ! Blue ; UV store at jk+1 w-level
            zeT(ji,jj) = zzeT                                 ! total       -        -      -
         END_2D
      END DO
      !
   END SUBROUTINE qsr_5BD


   FUNCTION qsr_ext_lev( pL, pfr ) RESULT( klev )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE trc_oce_ext_lev  ***
      !!       
      !! ** Purpose :   compute the maximum level of light penetration
      !!          
      !! ** Method  :   the function provides the level at which irradiance, I, 
      !!              has a negligible effect on temperature.
      !!                T(n+1)-T(n) = ∆t dk[I] / ( rho0 Cp e3t_k )  
      !!              I(k) has a negligible effect on temperature at level k if:
      !!                ∆t I(k) / ( rho0*Cp*e3t_k ) <= 1.e-15 °C
      !!              with I(z) = Qsr*pfr*EXP(-z/L), therefore :
      !!                z >= L * LOG( 1.e-15 * rho0*Cp*e3t_k / ( ∆t*Qsr*pfr ) )
      !!              with Qsr being the maximum normal surface irradiance at sea 
      !!              level (~1000 W/m2).
      !!                # pL is the longest depth of extinction:
      !!                   - pL = 23.00 m (2 bands case)
      !!                   - pL = 48.24 m (3 bands case: blue waveband & 0.03 mg/m2 for the chlorophyll)
      !!                # pfr is the fraction of solar radiation which penetrates,
      !!                considering Qsr=1000 W/m2 and rn_abs = 0.58:
      !!                   - Qsr*pfr0 = Qsr *    rn_abs    = 580 W/m2   (top absorbtion)
      !!                   - Qsr*pfr1 = Qsr * (1-rn_abs)   = 420 W/m2 (2 bands case)
      !!                   - Qsr*pfr1 = Qsr * (1-rn_abs)/3 = 140 W/m2 (3 bands case & equi-partition)
      !!
      !!----------------------------------------------------------------------
      INTEGER              ::   klev   ! result: maximum level of light penetration
      REAL(wp), INTENT(in) ::   pL     ! depth of extinction
      REAL(wp), INTENT(in) ::   pfr    ! frac. solar radiation which penetrates 
      !
      INTEGER  ::   jk                 ! dummy loop index
      REAL(wp) ::   zcoef              ! local scalar
      REAL(wp) ::   zhext              ! deepest depth until which light penetrates
      REAL(wp) ::   ze3t , zdw         ! max( e3t_k ) and min( w-depth_k+1 )
      REAL(wp) ::   zprec = 10.e-15_wp ! required precision 
      REAL(wp) ::   zQmax= 1000._wp    ! maximum normal surface irradiance at sea level (W/m2)
      !!----------------------------------------------------------------------
      !
      zcoef    =  zprec * rho0_rcp / ( rDt * zQmax * pfr)
      !
#if defined key_vco_1d || defined key_vco_1d3d
      ! z- or zps coordinate (use 1D ref vertcial coordinate)
         klev = jpkm1                              ! Level of light extinction zco / zps
         DO jk = jpkm1, 1, -1
            zdw  = gdepw_1d(jk+1)                  ! max w-depth at jk+1 level
            ze3t =   e3t_1d(jk  )                  ! minimum e3t at jk   level
            zhext =  - pL * LOG( zcoef * ze3t )    ! extinction depth
            IF( zdw >= zhext )   klev = jk         ! last T-level reached by Qsr
         END DO
#else
      ! s- or s-z- coordinate (use 3D vertical coordinate)
         klev = jpkm1                              ! Level of light extinction 
         DO jk = jpkm1, 1, -1    ! 
            IF( SUM( tmask(:,:,jk) ) > 0 ) THEN    ! ocean point at that level
               zdw  = MAXVAL( gdepw_0(:,:,jk+1) *       wmask(:,:,jk)         )    ! max w-depth at jk+1 level
               ze3t = MINVAL(   e3t_0(:,:,jk  ) , mask=(wmask(:,:,jk+1)==1)   )    ! minimum e3t at jk   level
               zhext =  - pL * LOG( zcoef * ze3t )                                 ! extinction depth
               IF( zdw >= zhext )   klev = jk                                      ! last T-level reached by Qsr
            ELSE                                   ! only land point at level jk
               klev = jk                                                           ! local domain sea-bed level 
            ENDIF
         END DO
         CALL mpp_max('tra_qsr', klev)             ! needed for reproducibility   !!st may be modified to avoid this comm.
         !                                                                        !!st use ssmask to remove the comm ?
#endif
      !
!!st      IF(lwp) WRITE(numout,*) '                level of e3t light extinction = ', klev, ' ref depth = ', gdepw_1d(klev+1), ' m'
   END FUNCTION qsr_ext_lev


   SUBROUTINE tra_qsr_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_qsr_init  ***
      !!
      !! ** Purpose :   Initialization for the penetrative solar radiation
      !!
      !! ** Method  :   The profile of solar radiation within the ocean is set
      !!      from two length scale of penetration (rn_si0,rn_si1) and a ratio
      !!      (rn_abs). These parameters are read in the namtra_qsr namelist. The
      !!      default values correspond to clear water (type I in Jerlov'
      !!      (1968) classification.
      !!         called by tra_qsr at the first timestep (nit000)
      !!
      !! ** Action  : - initialize rn_si0, rn_si1 and rn_abs
      !!
      !! Reference : Jerlov, N. G., 1968 Optical Oceanography, Elsevier, 194pp.
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk                  ! dummy loop indices
      INTEGER  ::   ios, ierror, ioptio   ! local integer
      REAL(wp) ::   zLU, zLB, zLG, zLR    ! local scalar
      REAL(wp) ::   zVlp, zVlu, zchl      !   -      -
      !
      CHARACTER(len=100) ::   cn_dir   ! Root directory for location of ssr files
      TYPE(FLD_N)        ::   sn_chl   ! informations about the chlorofyl field to be read
      !!
      NAMELIST/namtra_qsr/  sn_chl, cn_dir, ln_qsr_rgb, ln_qsr_2bd, ln_qsr_5bd, ln_qsr_bio,  &
         &                  nn_chldta, rn_abs, rn_par, rn_si0, rn_si1
      !!----------------------------------------------------------------------
      !
      READ_NML_REF(numnam,namtra_qsr)
      READ_NML_CFG(numnam,namtra_qsr)
      IF(lwm) WRITE ( numond, namtra_qsr )
      !
      IF(lwp) THEN            !**  control print  **!
         WRITE(numout,*)
         WRITE(numout,*) 'tra_qsr_init : penetration of the surface solar radiation'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtra_qsr : set the parameter of penetration'
         WRITE(numout,*) '      RGB (Red-Green-Blue) light penetration           ln_qsr_rgb = ', ln_qsr_rgb
         WRITE(numout,*) '      2 band               light penetration           ln_qsr_2bd = ', ln_qsr_2bd
         WRITE(numout,*) '      UV-R-G-B-IR          light penetration           ln_qsr_5bd = ', ln_qsr_5bd
         WRITE(numout,*) '      bio-model            light penetration           ln_qsr_bio = ', ln_qsr_bio
         WRITE(numout,*) '      RGB : 3D/Surface Chl data or Cst value (2,1,0)   nn_chldta  = ', nn_chldta
         WRITE(numout,*) '      RGB & 2 bands: fraction of light (rn_si1)        rn_abs     = ', rn_abs
         WRITE(numout,*) '      RGB & 2 bands: shortess attenuation depth        rn_si0     = ', rn_si0
         WRITE(numout,*) '      2 bands: longest attenuation depth               rn_si1     = ', rn_si1
         WRITE(numout,*) '      5 bands: fraction of PAR                         rn_par     = ', rn_par
         WRITE(numout,*)
      ENDIF
      !
      ioptio = 0              !**  Parameter control  **!
      IF( ln_qsr_rgb  )   ioptio = ioptio + 1
      IF( ln_qsr_2bd  )   ioptio = ioptio + 1
      IF( ln_qsr_5bd  )   ioptio = ioptio + 1
      IF( ln_qsr_bio  )   ioptio = ioptio + 1
      !
      IF( ioptio /= 1 )   CALL ctl_stop( 'Choose ONE type of light penetration in namelist namtra_qsr',  &
         &                               ' 2 bands, 3 RGB bands, 5 UV-R-G-B-IR bands or bio-model light penetration' )
      !
      IF( ln_qsr_rgb .AND.   nn_chldta == 0                       )   nqsr = np_RGB
      IF( ln_qsr_rgb .AND. ( nn_chldta == 1 .OR. nn_chldta == 2 ) )   nqsr = np_RGBc
      IF( ln_qsr_2bd                                              )   nqsr = np_2BD
      IF( ln_qsr_5bd .AND.   nn_chldta == 0                       )   nqsr = np_5BD
      IF( ln_qsr_5bd .AND. ( nn_chldta == 1 .OR. nn_chldta == 2 ) )   nqsr = np_5BDc
      IF( ln_qsr_bio                                              )   nqsr = np_BIO
      !
      !                       !**  Initialisation  **!
      !
      !                                !==  Infrared attenuation  ==!   (all schemes)
      !                                !============================!
      !
      r1_si0 = 1._wp / rn_si0                ! inverse of infrared attenuation length
      !
      nk0 = qsr_ext_lev( rn_si0, rn_abs )    ! level of light extinction
      !
      IF(lwp) WRITE(numout,*) '   ==>>>   Infrared light attenuation'
      IF(lwp) WRITE(numout,*) '              level of infrared extinction = ', nk0, ' ref depth = ', gdepw_1d(nk0+1), ' m'
      IF(lwp) WRITE(numout,*)
      !
      SELECT CASE( nqsr )
      !
      CASE( np_RGBc, np_RGB )          !==  Red-Green-Blue light attenuation  ==!   (Chl data or constant)
         !                             !========================================!
         !
         IF( nqsr == np_RGB ) THEN   ;   zchl   = 0.05        ! constant Chl value
         ELSE                        ;   zchl   = 0.03        ! minimum  Chl value
         ENDIF
         zchl   = MAX( 0.03_wp , MIN( zchl , 10._wp) )     ! NB. make sure that chosen value verifies: 0.03 < zchl < 10
         nc_tab = NINT( 41 + 20.*LOG10(zchl) + 1.e-15 )    ! Convert Chl value to attenuation coefficient look-up table index
         !
         CALL trc_oce_tab( rktab )                 ! tabulated attenuation coef.
         !
         zVlp =  ( 1._wp - rn_abs ) / 3._wp        ! visible light equi-partition
         !
         !     1 / length          !   attenuation  length   !         attenuation level
         r1_LR = rktab(3,nc_tab)   ;   zLR = 1._wp / r1_LR   ;   nkR = qsr_ext_lev( zLR, zVlp )   ! Red   
         r1_LG = rktab(2,nc_tab)   ;   zLG = 1._wp / r1_LG   ;   nkG = qsr_ext_lev( zLG, zVlp )   ! Green
         r1_LB = rktab(1,nc_tab)   ;   zLB = 1._wp / r1_LB   ;   nkB = qsr_ext_lev( zLB, zVlp )   ! Blue
         !
         nkV = nkB                                 ! maximum level of light penetration
         !
         IF( nqsr == np_RGB ) THEN
            IF(lwp) WRITE(numout,*) '   ==>>>   RGB:  light attenuation with a constant Chlorophyll = ', zchl
         ELSE
            IF(lwp) WRITE(numout,*) '   ==>>>   RGB:  light attenuation using Chlorophyll data with min(Chl) = ', zchl
         ENDIF            
         IF(lwp) WRITE(numout,*) '                 level of Red   extinction = ', nkR, ' ref depth = ', gdepw_1d(nkR+1), ' m'
         IF(lwp) WRITE(numout,*) '                 level of Green extinction = ', nkG, ' ref depth = ', gdepw_1d(nkG+1), ' m'
         IF(lwp) WRITE(numout,*) '                 level of Blue  extinction = ', nkB, ' ref depth = ', gdepw_1d(nkB+1), ' m'
         IF(lwp) WRITE(numout,*)
         !
         IF( nqsr == np_RGBc ) THEN                ! Chl data : set sf_chl structure
            IF(lwp) WRITE(numout,*) '   ==>>>   Chlorophyll read in a file'
            ALLOCATE( sf_chl(1), STAT=ierror )
            IF( ierror > 0 ) THEN
               CALL ctl_stop( 'tra_qsr_init: unable to allocate sf_chl structure' )   ;   RETURN
            ENDIF
            IF    ( nn_chldta == 1 ) THEN
               ALLOCATE( sf_chl(1)%fnow(jpi,jpj,1) )
               IF( sn_chl%ln_tint )   ALLOCATE( sf_chl(1)%fdta(jpi,jpj,1,2) )
            ELSEIF( nn_chldta == 2 ) THEN
               ALLOCATE( sf_chl(1)%fnow(jpi,jpj,jpk) )
               IF( sn_chl%ln_tint )   ALLOCATE( sf_chl(1)%fdta(jpi,jpj,jpk,2) )
            ENDIF
            !                                        ! fill sf_chl with sn_chl and control print
            CALL fld_fill( sf_chl, (/ sn_chl /), cn_dir, 'tra_qsr_init',                             &
               &           'Solar penetration function of read chlorophyll', 'namtra_qsr' , no_print )
         ENDIF
         !
      CASE( np_2BD )                   !==  2 bands light attenuation (IR+ visible light) ==!
         !
         !
         IF( lk_top ) CALL trc_oce_tab( rktab )      ! tabulated attenuation coef.
         !
         r1_si1 = 1._wp / rn_si1                     ! inverse of visible light attenuation
         zVlp =  ( 1._wp - rn_abs )                  ! visible light partition
         nkV  = qsr_ext_lev( rn_si1, zVlp )          ! level of visible light extinction
         !
         IF(lwp) WRITE(numout,*) '   ==>>>   2 bands attenuation (Infrared + Visible light) '
         IF(lwp) WRITE(numout,*) '                level of visible light extinction = ', nkV, ' ref depth = ', gdepw_1d(nkV+1), ' m'
         IF(lwp) WRITE(numout,*)
         !
      CASE( np_5BDc, np_5BD )          !==  UV-Red-Green-Blue-IR light attenuation  ==!   (Chl data or constant)
         !                             !==============================================!
         !
         IF( nqsr == np_5BD ) THEN   ;   zchl   = 0.05        ! constant Chl value
         ELSE                        ;   zchl   = 0.03        ! minimum  Chl value
         ENDIF
         zchl   = MAX( 0.03_wp , MIN( zchl , 10._wp) )     ! NB. make sure that chosen value verifies: 0.03 < zchl < 10
         nc_tab = NINT( 41 + 20.*LOG10(zchl) + 1.e-15 )    ! Convert Chl value to attenuation coefficient look-up table index
         !
         CALL trc_oce_tab( rktab )                 ! tabulated attenuation coef.
         !
         zVlp =  ( rn_par ) / 3._wp         ! visible light equi-partition
         zVlu =  ( 1._wp - rn_abs - rn_par) !      uv light      partition
         !
         !     1 / length          !   attenuation  length   !         attenuation level
         r1_LR = rktab(3,nc_tab)   ;   zLR = 1._wp / r1_LR   ;   nkR = qsr_ext_lev( zLR, zVlp )   ! Red   
         r1_LG = rktab(2,nc_tab)   ;   zLG = 1._wp / r1_LG   ;   nkG = qsr_ext_lev( zLG, zVlp )   ! Green
         r1_LB = rktab(1,nc_tab)   ;   zLB = 1._wp / r1_LB   ;   nkB = qsr_ext_lev( zLB, zVlp )   ! Blue
         r1_LU = rktab(4,nc_tab)   ;   zLU = 1._wp / r1_LU   ;   nkU = qsr_ext_lev( zLU, zVlu )   ! UV
         !
         nkV = nkB                                 ! maximum level of light penetration
         !
         IF( nqsr == np_5BD ) THEN
            IF(lwp) WRITE(numout,*) '   ==>>>   5BD:  light attenuation with a constant Chlorophyll = ', zchl
         ELSE
            IF(lwp) WRITE(numout,*) '   ==>>>   5BD:  light attenuation using Chlorophyll data with min(Chl) = ', zchl
         ENDIF
         IF(lwp) WRITE(numout,*) '                 level of Red   extinction = ', nkR, ' ref depth = ', gdepw_1d(nkR+1), ' m'
         IF(lwp) WRITE(numout,*) '                 level of Green extinction = ', nkG, ' ref depth = ', gdepw_1d(nkG+1), ' m'
         IF(lwp) WRITE(numout,*) '                 level of Blue  extinction = ', nkB, ' ref depth = ', gdepw_1d(nkB+1), ' m'
         IF(lwp) WRITE(numout,*) '                 level of   UV  extinction = ', nkU, ' ref depth = ', gdepw_1d(nkU+1), ' m'
         IF(lwp) WRITE(numout,*)
         !
         IF( nqsr == np_5BDc ) THEN                ! Chl data : set sf_chl structure
            IF(lwp) WRITE(numout,*) '   ==>>>   Chlorophyll read in a file'
            ALLOCATE( sf_chl(1), STAT=ierror )
            IF( ierror > 0 ) THEN
               CALL ctl_stop( 'tra_qsr_init: unable to allocate sf_chl structure' )   ;   RETURN
            ENDIF
            IF    ( nn_chldta == 1 ) THEN
               ALLOCATE( sf_chl(1)%fnow(jpi,jpj,1) )
               IF( sn_chl%ln_tint )   ALLOCATE( sf_chl(1)%fdta(jpi,jpj,1,2) )
            ELSEIF( nn_chldta == 2 ) THEN
               ALLOCATE( sf_chl(1)%fnow(jpi,jpj,jpk) )
               IF( sn_chl%ln_tint )   ALLOCATE( sf_chl(1)%fdta(jpi,jpj,jpk,2) )
            ENDIF
            !                                        ! fill sf_chl with sn_chl and control print
            CALL fld_fill( sf_chl, (/ sn_chl /), cn_dir, 'tra_qsr_init',                             &
               &           'Solar penetration function of read chlorophyll', 'namtra_qsr' , no_print )
         ENDIF
         !
      CASE( np_BIO )                   !==  BIO light penetration  ==!
         !
         IF(lwp) WRITE(numout,*) '   ==>>>   bio-model light penetration'
         IF( .NOT.lk_top )   CALL ctl_stop( 'No bio model : ln_qsr_bio = true impossible ' )
         !
         CALL trc_oce_tab( rktab )                 ! tabulated attenuation coef.
         !
         nkV = trc_oce_ext_lev( r_si2, 33._wp )    ! maximum level of light extinction
         !
         IF(lwp) WRITE(numout,*) '        level of light extinction = ', nkV, ' ref depth = ', gdepw_1d(nkV+1), ' m'
         !
      END SELECT
      !jth this needs to be the the maximum across all procs
      call mpp_max('traqsr',nksr)
      !
      nksr = nkV       ! name of max level of light extinction used in traatf(_qco).F90
      !
#if ! defined key_RK3
      qsr_hc(:,:,:) = 0._wp      ! MLF : now qsr heat content set to zero where it will not be computed
#endif
      !
      !                          ! Sea-ice :   1st ocean level attenuation coefficient (used in sbcssm)
      IF( iom_varid( numror, 'fraqsr_1lev', ldstop = .FALSE. ) > 0 ) THEN
         CALL iom_get( numror, jpdom_auto, 'fraqsr_1lev'  , fraqsr_1lev  )
      ELSE
         fraqsr_1lev(:,:) = 1._wp   ! default : no penetration
      ENDIF
      !
   END SUBROUTINE tra_qsr_init

   !!======================================================================
END MODULE traqsr

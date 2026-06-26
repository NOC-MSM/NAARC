MODULE icefrm
   !!======================================================================
   !!                       ***  MODULE  icefrm  ***
   !! Sea-ice form drag param from CICE
   !!======================================================================
   !! History :  4.2        ! 2023  (D. Schroeder)       Form drag from CICE5
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------

   USE ice                  ! sea-ice variables
   USE in_out_manager       ! I/O manager
   USE iom                  ! for iom_put
   USE prtctl               ! print control
   USE lbclnk               ! lateral boundary conditions (or MPP link)
   USE timing               ! timing
   USE phycst               ! physical constants
   USE eosbn2               ! equation of state
   USE sbcblk , ONLY : nn_frm, csw, csa, cra, crw, cfa, cfw, rn_Cd_i, rn_frm_ht0 ! FORM drag param

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_frm        ! called by icestp.F90

   ! Sail/keel parameters    
   REAL(wp), PARAMETER :: zazlpha   = 0._wp     ! weight functions for area of ridged ice 
   REAL(wp), PARAMETER :: zbeta2   = 0.75_wp       
   REAL(wp), PARAMETER :: zbeta    = 0.5_wp     ! power exponent appearing in astar and 
                                                ! L=Lmin(A*/(A*-A))**beta [0,1]
   REAL(wp), PARAMETER :: tanar    = 0.4_wp     ! sail slope 
   REAL(wp), PARAMETER :: tanak    = 0.4_wp     ! keel slope
   REAL(wp), PARAMETER :: hkoverhr = 4._wp      ! hkeel/hridge ratio
   REAL(wp), PARAMETER :: dkoverdr = 1._wp      ! dkeel/distrdg ratio
   REAL(wp), PARAMETER :: zphir    = 0.8_wp     ! porosity of ridges
   REAL(wp), PARAMETER :: zphik    = 0.8_wp     ! porosity of keels

   ! Floe parameters   
   REAL(wp), PARAMETER :: zLmin    = 8._wp      ! min length of floe (m) [5,100]
   REAL(wp), PARAMETER :: zLmax    = 300._wp    ! max length of floe (m) [30,3000]
   REAL(wp), PARAMETER :: zLmoy    = 300._wp    ! average length of floe (m) [30,1000]

   ! Melt pond parameters   
   REAL(wp), PARAMETER :: zlpmin    = 2.26_wp   ! min pond length (m) see Eq. 17 [1,10]
   REAL(wp), PARAMETER :: zlpmax    = 24.63_wp  ! max pond length (m) see Eq. 17 [10,100]
   REAL(wp), PARAMETER :: cpa       = 0.2_wp    ! ratio of local form drag over 
                                                ! geometrical parameter [0,1]

  ! Limits for stability
   REAL(wp), PARAMETER :: camax    = 0.0028_wp  ! 0.0112_wp  ! Maximum for atmospheric drag
   REAL(wp), PARAMETER :: cwmax    = 0.01_wp    ! 0.02_wp    ! Maximum for ocean drag
   REAL(wp), PARAMETER :: camin    = 0.0007_wp  ! Minimum for atmospheric drag
   REAL(wp), PARAMETER :: cwmin    = 0.0025_wp  ! Minimum for ocean drag

   ! Other parameters   
   REAL(wp), PARAMETER :: sHGB     = 0.18_wp    ! attenuation parameter for the sheltering function
   REAL(wp), PARAMETER :: mrdg     = 20._wp     ! screening effect see Lu2011 [5,50]
   REAL(wp), PARAMETER :: mrdgo    = 10._wp     ! screening effect see Lu2011 [5,50]
   REAL(wp), PARAMETER :: zref     = 10._wp     ! reference height
   REAL(wp), PARAMETER :: sl       = 22._wp     ! Sheltering parameter Lupkes2012 [10,30]
   REAL(wp), PARAMETER :: ocnruf   = 0.000327_wp! ocean surface roughness (m)
   REAL(wp), PARAMETER :: iceruf   = 0.0005_wp  ! ice surface roughness (m)
   REAL(wp), PARAMETER :: ocnrufi  = 1._wp/ocnruf  ! inverse ocean roughness
   REAL(wp), PARAMETER :: icerufi  = 1._wp/iceruf  ! inverse ice roughness

   ! Local variables
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zardg_drag    ! ridged ice concentration
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zvrdg_drag    ! ridged ice thickness
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zhfreebd     ! freeboard (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zhdraft      ! draft of ice + snow column (Stoessel1993)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zhridge      ! ridge height
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zdistrdg     ! distance between ridges
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zhkeel       ! keel depth
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zdkeel       ! distance between keels
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zlfloe       ! floe length
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zdfloe       ! distance between floes
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zdrag_ia_skin ! neutral skin drag coefficient
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zdrag_ia_floe ! neutral floe edge drag coefficient
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zdrag_ia_pond ! neutral pond edge drag coefficient
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zdrag_ia_rdg  ! neutral ridge drag coefficient
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zdrag_io_skin ! skin drag coefficient
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zdrag_io_floe ! floe edge drag coefficient
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zdrag_io_keel ! keel drag coefficient

   REAL(wp) :: zai    ! ice area 
   REAL(wp) :: zaii   ! ice area inverse
   REAL(wp) :: ztmp1  ! temporary value for ridges and keels conditions
   REAL(wp) :: zsca   ! wind attenuation function
   REAL(wp) :: zscw   ! ocean attenuation function
   REAL(wp) :: ztecar ! temporary value (ridge)
   REAL(wp) :: ztecwk ! temporary value (keel)
   REAL(wp) :: ztecaf ! temporary value (floe atmosphere)
   REAL(wp) :: ztecwf ! temporary value (floe ocean)
   REAL(wp) :: zastar ! new constant for form drag
   REAL(wp) :: zlp    ! pond length (m)
   REAL(wp) :: zapond !  melt pond fraction of a grid cell


   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2018)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_frm( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ice_frm  ***
      !!                   
      !! ** Purpose :   Neutral drag coefficients for turbulent momentum
      !!                exchange between ocean - sea-ice and atmosphere sea-ice 
      !!                calculated from ridge height, distance, floe size and pond
      !!                fraction based upon Tsamados et al. (2014), 
      !!                JPO, DOI: 10.1175/JPO-D-13-0215.1.
      !!
      !! ** Action  :   Calculates drag_io and drag_ia 
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER  ::   ji , jj                                        ! dummy loop indices
      INTEGER  ::   ij                                             ! combined i and j indices
      INTEGER  ::   icells                                         ! number of cells where the drag coefficients will be computed
      INTEGER, ALLOCATABLE, SAVE, DIMENSION(:)  ::   indxi , indxj ! compressed i and j indices
      !
      REAL(wp), PARAMETER :: ardgi    = -0.004216_wp ! parameter used to estimate ardg value (intercept)
      REAL(wp), PARAMETER :: ardgi1   =  0.157815_wp ! parameter used to estimate ardg value (fist degree)
      REAL(wp), PARAMETER :: ardgi2   = -0.005002_wp ! parameter used to estimate ardg value (second degree)
      REAL(wp), PARAMETER :: vrdgi    = -0.059959_wp ! parameter used to estimate vrdg value (intercept)
      REAL(wp), PARAMETER :: vrdgi1   =  0.501150_wp ! parameter used to estimate vrdg value (first degree)
      REAL(wp), PARAMETER :: vrdgi2   =  0.077504_wp ! parameter used to estimate vrdg value (second degree)

      REAL(wp), DIMENSION(jpi,jpj) ::   zmsk00

      !
      IF( ln_timing )   CALL timing_start('ice_frm')
      !
      ALLOCATE ( indxi(jpj*jpi), indxj(jpj*jpi) )
      !
      IF( kt == nit000 ) THEN   ! at first time-step
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'ice_frm: sea-ice form drag param'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
         IF(lwp) WRITE(numout,*)
      ENDIF
      ! 
      
      zmsk00(:,:) = MERGE( 1._wp, 0._wp, at_i(:,:) >= epsi06  )

      ! David >
      ALLOCATE( zardg_drag  (jpi,jpj), zvrdg_drag   (jpi,jpj), zhfreebd(jpi,jpj) ,   &
                zhdraft    (jpi,jpj) , zhridge     (jpi,jpj) , zdistrdg(jpi,jpj) ,   &
                zhkeel     (jpi,jpj) , zdkeel      (jpi,jpj) , zlfloe(jpi,jpj) ,   &
                zdfloe     (jpi,jpj) , zdrag_ia_skin(jpi,jpj), zdrag_ia_floe(jpi,jpj) ,   &
                zdrag_ia_rdg(jpi,jpj), zdrag_ia_pond(jpi,jpj), zdrag_io_skin(jpi,jpj) ,   &
                zdrag_io_floe(jpi,jpj),zdrag_io_keel(jpi,jpj) )


      !----------------------------------------------------------!
      !   Evaluation of ridged ice concentration and thickness   !
      !----------------------------------------------------------!
      ! Note: Currently there is no tracer for level or ridged-ice in SI3.
      ! Therefore, we preliminary use this approximation from Félicien Tournay
      ! and Antoine Barthélemy in 2017, at UCLouvain. They used model outputs
      ! from CICE to fit polynomials giving the volume and concentration of
      ! deformed ice as functions of sea ice properties. 

       DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            zardg_drag(ji,jj) = ardgi2 * ( vt_i(ji,jj) ** 2._wp ) + ardgi1 * vt_i(ji,jj) + ardgi
            zvrdg_drag(ji,jj) = vrdgi2 * ( vt_i(ji,jj) ** 2._wp ) + vrdgi1 * vt_i(ji,jj) + vrdgi
            zardg_drag(ji,jj) = min( max( 0._wp , zardg_drag(ji,jj) ) , 1._wp          )
            zvrdg_drag(ji,jj) = min( max( 0._wp , zvrdg_drag(ji,jj) ) , vt_i(ji,jj) )

            ! For small values of ice volume, vrdg_drag can be zero and ardg_drag very small.
            ! To avoid problems with skin drags, ardg_drag is set to zero when vrdg_drag is.
            IF (zvrdg_drag(ji,jj) < epsi10) THEN
               zardg_drag(ji,jj) = 0._wp
            ENDIF
       END_2D

      !-----------------------------------!
      !   Initialize                      !
      !-----------------------------------!

      zastar = 1._wp / ( 1._wp - ( zLmin / zLmax ) ** ( 1._wp / zbeta ) )

       DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            drag_ia(:,:) = rn_Cd_i 
            drag_io(:,:) = rn_cio
            zhfreebd(:,:) = 0._wp
            zhdraft (:,:) = 0._wp
            zhridge (:,:) = 0._wp
            zdistrdg(:,:) = 0._wp
            zhkeel  (:,:) = 0._wp
            zdkeel  (:,:) = 0._wp
            zlfloe  (:,:) = 0._wp
            zdfloe  (:,:) = 0._wp
            zdrag_ia_skin(:,:) = 0._wp
            zdrag_ia_floe(:,:) = 0._wp
            zdrag_ia_pond(:,:) = 0._wp
            zdrag_ia_rdg (:,:) = 0._wp
            zdrag_io_skin(:,:) = 0._wp
            zdrag_io_floe(:,:) = 0._wp
            zdrag_io_keel(:,:) = 0._wp
      END_2D

      !------------------------------------------!
      !   Identify cells with nonzero ice area   !
      !------------------------------------------!

      icells = 0
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
          IF (at_i(ji,jj) > 0.01_wp) THEN
            icells = icells + 1
            indxi(icells) = ji
            indxj(icells) = jj
          ENDIF
      END_2D

      DO ij = 1, icells
         ji = indxi(ij)
         jj = indxj(ij)

      !--------------------------------!
      !   Compute average quantities   !
      !--------------------------------!

         zai  = at_i(ji,jj)
         zaii = 1._wp / zai

         !!! melt ponds !!
         zapond = SUM(a_ip(ji,jj,:))

         !!! draft and freeboard !!!
         zhdraft(ji,jj) = ( rhoi * vt_i(ji,jj) + rhos * vt_s(ji,jj) ) * zaii / rho0 ! without ponds
         zhfreebd(ji,jj) = (vt_i(ji,jj) + vt_s(ji,jj) ) * zaii - zhdraft(ji,jj)

         ! Do not allow draft larger than ice thickness (see Eq. 28)  
         IF (zhdraft(ji,jj) >= vt_i(ji,jj) * zaii) THEN
            ! replace excess snow with ice so hi (sea ice volume) ~= zhdraft
            zhfreebd(ji,jj) = ( zhdraft(ji,jj) * zai * ( 1._wp - rhoi / rho0 ) &
              &      + (vt_s(ji,jj) - ( vt_i(ji,jj) - zhdraft(ji,jj) * zai ) * rhoi / rhos ) &
              &      * ( 1._wp - rhos / rho0 ) ) * zaii ! Stoessel1993  
         ENDIF

         !!! floes !!!
         ! floe size parameterization see Eq. 13
         zlfloe(ji,jj) = zLmin * ( zastar / ( zastar - zai ) ) ** zbeta

         ! distance between floes parameterization see Eq. 14
         zdfloe(ji,jj) = zlfloe(ji,jj) * ( 1._wp / sqrt( zai ) - 1._wp )

         !!! ridges and keels !!!
         ! hridge, hkeel, distrdg and dkeel estimates for 
         ! simple triangular geometry as in the CICE code

         IF (zardg_drag(ji,jj) > 0.01_wp) THEN

            zhridge(ji,jj) = zvrdg_drag(ji,jj) / zardg_drag(ji,jj) * 2._wp                    &
              &              * ( zazlpha + zbeta2 * hkoverhr / dkoverdr * tanar / tanak )    &
              &              / ( zphir + zphik * tanar / tanak * hkoverhr ** 2._wp / dkoverdr )
            zdistrdg(ji,jj) = 2._wp * zhridge(ji,jj) * zai / zardg_drag(ji,jj)               &
              &              * ( zazlpha / tanar + zbeta2 / tanak * hkoverhr / dkoverdr )
            zhkeel(ji,jj) = hkoverhr * zhridge(ji,jj)
            zdkeel(ji,jj) = dkoverdr * zdistrdg(ji,jj)

            ! Use the height of ridges relative to the mean freeboard of
            ! the pack.  

            ztmp1 = max( 0._wp , zhridge(ji,jj) - zhfreebd(ji,jj) )

      !----------------------!
      !   Skin drag (atmo)   !
      !----------------------!

            zdrag_ia_skin(ji,jj) = csa * ( 1._wp - mrdg * ( ztmp1 / zdistrdg(ji,jj) ) )
            zdrag_ia_skin(ji,jj) = max( min( zdrag_ia_skin(ji,jj) , camax ) , 0._wp )

      !-------------------------!
      !   Ridge effect (atmo)   !
      !-------------------------!

            IF (ztmp1 > epsi10) THEN
               zsca = 1._wp - exp( -sHGB * zdistrdg(ji,jj) / ztmp1 )
               ztecar = cra * 0.5_wp
               zdrag_ia_rdg(ji,jj) = ztecar * ztmp1 / zdistrdg(ji,jj) * zsca                &
                  &      *  ( log( ztmp1 * icerufi ) / log( zref * icerufi ) ) ** 2._wp
               zdrag_ia_rdg(ji,jj) = min( zdrag_ia_rdg(ji,jj) , camax )
            ENDIF

            ztmp1 = max( 0._wp , zhkeel(ji,jj) - zhdraft(ji,jj) )

      !----------------------------------!
      !   Skin drag bottom ice (ocean)   !
      !----------------------------------!

            zdrag_io_skin(ji,jj) = csw * (1._wp - mrdgo * ( ztmp1 / zdkeel(ji,jj) ) )
            zdrag_io_skin(ji,jj) = max( min( zdrag_io_skin(ji,jj) , cwmax) , 0._wp)

      !-------------------------!
      !   Keel effect (ocean)   !
      !-------------------------!

            IF (ztmp1 > epsi10) THEN
               zscw = 1._wp - exp( -sHGB * zdkeel(ji,jj) / ztmp1 )
               ztecwk = crw * 0.5_wp
               zdrag_io_keel(ji,jj) = ztecwk * ztmp1 / zdkeel(ji,jj) * zscw                &
                 &        *  ( log( ztmp1 * icerufi ) / log( zref * icerufi ) ) ** 2._wp
               zdrag_io_keel(ji,jj) = max( min( zdrag_io_keel(ji,jj) , cwmax ) , 0._wp )
            ENDIF

         ELSE                ! ardg_drag > 0.01_wp

            !!! skin drag in case of ardg_drag < 0.01 but at_i > 0.01 !!!
            ! this is needed when the ice concentration is high but the ice volume is small
            ! otherwise skin drag is set to zero even though there is undeformed ice

            zdrag_ia_skin(ji,jj) = csa
            zdrag_io_skin(ji,jj) = csw

         ENDIF               ! ardg_drag > 0.01_wp

      !------------------------------!
      ! Floe edge drag effect (atmo) !
      !------------------------------!

         IF (zhfreebd(ji,jj) > epsi10) THEN
            zsca = 1._wp - exp( -sl * zbeta * ( 1._wp - zai ) )
            ztecaf = cfa * 0.5_wp                                                 &
               &      * ( log( zhfreebd(ji,jj) * ocnrufi ) / log( zref * ocnrufi ) ) ** 2._wp * zsca
            zdrag_ia_floe(ji,jj) = ztecaf * zhfreebd(ji,jj) / zlfloe(ji,jj)
            zdrag_ia_floe(ji,jj) = max( min( zdrag_ia_floe(ji,jj) , camax ) , 0._wp )
         ENDIF

      !-------------------------!
      ! Pond edge effect (atmo) !
      !-------------------------!

        IF (zhfreebd(ji,jj) > epsi10) THEN
           zsca = ( zapond ) ** ( 1._wp / ( zref * zbeta ) )
           zlp  = zlpmin * ( 1._wp - zapond ) + zlpmax * zapond
           zdrag_ia_pond(ji,jj) = cpa * 0.5_wp * zsca * zapond * zhfreebd(ji,jj) / zlp &
             &          * ( log ( zhfreebd(ji,jj) * ocnrufi ) / log( zref * ocnrufi ) ) ** 2._wp
           zdrag_ia_pond(ji,jj) = max( min( zdrag_ia_pond(ji,jj) , camax ) , 0._wp )
        ENDIF

      !-------------------------------!
      ! Floe edge drag effect (ocean) !
      !-------------------------------!

         IF (zhdraft(ji,jj) > epsi10) THEN
            zscw = 1._wp - exp( -sl * zbeta * ( 1._wp - zai ) )
            ztecwf = cfw * 0.5_wp & 
              &   * ( log( zhdraft(ji,jj) * ocnrufi ) / log( zref * ocnrufi ) ) ** 2._wp  * zscw
            zdrag_io_floe(ji,jj) = ztecwf * zhdraft(ji,jj) / zlfloe(ji,jj)
            zdrag_io_floe(ji,jj) = max( min( zdrag_io_floe(ji,jj) , cwmax) , 0._wp )
         ENDIF

      !-------------------------------!
      ! Total drag coefficient (atmo) !
      !-------------------------------!

         drag_ia(ji,jj) = zdrag_ia_skin(ji,jj) + zdrag_ia_floe(ji,jj) & 
            &              + zdrag_ia_rdg(ji,jj)  + zdrag_ia_pond(ji,jj)
         drag_ia(ji,jj) = min( drag_ia(ji,jj) , camax )
         drag_ia(ji,jj) = max( drag_ia(ji,jj) , camin )

      !--------------------------------!
      ! Total drag coefficient (ocean) !
      !--------------------------------!
  
         drag_io(ji,jj) = zdrag_io_skin(ji,jj) + zdrag_io_floe(ji,jj) + zdrag_io_keel(ji,jj)
         drag_io(ji,jj) = min( drag_io(ji,jj) , cwmax )
         drag_io(ji,jj) = max( drag_io(ji,jj) , cwmin )
      
      !--------------------------------!
      ! JDHA adjustment close to coast !
      !--------------------------------!
         
         IF (ht_0(ji,jj) <= rn_frm_ht0) drag_io(ji,jj) = rn_cio

      ENDDO               ! ij


      CALL iom_put( 'drag_io', drag_io * zmsk00 )
      CALL iom_put( 'drag_io_skin', zdrag_io_skin * zmsk00 )
      CALL iom_put( 'drag_io_keel', zdrag_io_keel * zmsk00 )
      CALL iom_put( 'drag_io_floe', zdrag_io_floe * zmsk00 )

      CALL iom_put( 'drag_ia', drag_ia * zmsk00 )
      CALL iom_put( 'drag_ia_skin', zdrag_ia_skin * zmsk00 )
      CALL iom_put( 'drag_ia_rdg', zdrag_ia_rdg * zmsk00 )
      CALL iom_put( 'drag_ia_floe', zdrag_ia_floe * zmsk00 )
      CALL iom_put( 'drag_ia_pond', zdrag_ia_pond * zmsk00 )

      CALL iom_put( 'ardg_drag', zardg_drag * zmsk00 )
      CALL iom_put( 'vrdg_drag', zvrdg_drag * zmsk00 )

      !
      DEALLOCATE(indxi, indxj)
      DEALLOCATE(zardg_drag, zvrdg_drag, zhfreebd, zhdraft, zhridge, zdistrdg, zhkeel, &
       &         zdkeel, zlfloe, zdfloe, zdrag_ia_skin, zdrag_ia_floe, zdrag_ia_rdg, &
       &         zdrag_ia_pond, zdrag_io_skin, zdrag_io_floe, zdrag_io_keel)
      !
      IF( ln_timing )   CALL timing_stop('ice_frm')
      !
   END SUBROUTINE ice_frm

#else
    !!----------------------------------------------------------------------
    !!   Default option         Empty Module                No sea-ice model
    !!----------------------------------------------------------------------
CONTAINS
    SUBROUTINE ice_frm( kt )         ! Empty routine
       WRITE(*,*) 'ice_frm: You should not have seen this print! error?', kt
    END SUBROUTINE ice_frm
#endif


END MODULE icefrm

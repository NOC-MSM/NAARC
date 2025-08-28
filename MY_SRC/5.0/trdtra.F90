MODULE trdtra
   !!======================================================================
   !!                       ***  MODULE  trdtra  ***
   !! Ocean diagnostics:  ocean tracers trends pre-processing
   !!=====================================================================
   !! History :  3.3  !  2010-06  (C. Ethe) creation for the TRA/TRC merge
   !!            3.5  !  2012-02  (G. Madec) update the comments 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   trd_tra       : pre-process the tracer trends
   !!   trd_tra_adv   : transform a div(U.T) trend into a U.grad(T) trend
   !!   trd_tra_mng   : tracer trend manager: dispatch to the diagnostic modules
   !!   trd_tra_iom   : output 3D tracer trends using IOM
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean domain
   USE domutl, ONLY : lbnd_ij
   USE sbc_oce        ! surface boundary condition: ocean
   USE zdf_oce        ! ocean vertical physics
   USE trd_oce        ! trends: ocean variables
   USE trdtrc         ! ocean passive mixed layer tracers trends 
   USE trdglo         ! trends: global domain averaged
   USE trdpen         ! trends: Potential ENergy
   USE trdmxl         ! ocean active mixed layer tracers trends 
   USE ldftra         ! ocean active tracers lateral physics
   USE ldfslp
   USE zdfddm         ! vertical physics: double diffusion
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trd_tra   ! called by all tra_... modules

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   trdtx, trdty, trdt   ! use to store the temperature trends
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   avt_evd  ! store avt_evd to calculate EVD trend

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION trd_tra_alloc()
      !!---------------------------------------------------------------------
      !!                  ***  FUNCTION trd_tra_alloc  ***
      !!---------------------------------------------------------------------
      ALLOCATE( trdtx(A2D(0),jpk), trdty(A2D(0),jpk), trdt(A2D(0),jpk), avt_evd(A2D(0),jpk), STAT=trd_tra_alloc )
      !
      CALL mpp_sum ( 'trdtra', trd_tra_alloc )
      IF( trd_tra_alloc /= 0 )   CALL ctl_stop( 'STOP', 'trd_tra_alloc: failed to allocate arrays' )
   END FUNCTION trd_tra_alloc


   SUBROUTINE trd_tra( kt, Kmm, Krhs, ctype, ktra, ktrd, ptrd, pu, ptra )
      !!
      INTEGER                   , INTENT(in)           ::   kt      ! time step
      CHARACTER(len=3)          , INTENT(in)           ::   ctype   ! tracers trends type 'TRA'/'TRC'
      INTEGER                   , INTENT(in)           ::   ktra    ! tracer index
      INTEGER                   , INTENT(in)           ::   ktrd    ! tracer trend index
      INTEGER                   , INTENT(in)           ::   Kmm, Krhs ! time level indices
      REAL(wp), DIMENSION(:,:,:), INTENT(in)           ::   ptrd    ! tracer trend  or flux
      REAL(wp), DIMENSION(:,:,:), INTENT(in), OPTIONAL ::   pu      ! now velocity
      REAL(wp), DIMENSION(:,:,:), INTENT(in), OPTIONAL ::   ptra    ! now tracer variable
      !
      INTEGER, DIMENSION(2) :: itu, ittra
      !!
      IF( PRESENT(pu)   ) THEN ; itu(:)   = lbnd_ij(pu)   ; ELSE ; itu(:)   = 1 ; ENDIF
      IF( PRESENT(ptra) ) THEN ; ittra(:) = lbnd_ij(ptra) ; ELSE ; ittra(:) = 1 ; ENDIF

      CALL trd_tra_t( kt, Kmm, Krhs, ctype, ktra, ktrd, ptrd, lbnd_ij(ptrd), pu, itu(:), ptra, ittra(:) )
   END SUBROUTINE trd_tra


   SUBROUTINE trd_tra_t( kt, Kmm, Krhs, ctype, ktra, ktrd, ptrd, kttrd, pu, ktu, ptra, kttra )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_tra  ***
      !! 
      !! ** Purpose : pre-process tracer trends
      !!
      !! ** Method  : - mask the trend
      !!              - advection (ptra present) converte the incoming flux (U.T) 
      !!              into trend (U.T => -U.grat(T)=div(U.T)-T.div(U)) through a 
      !!              call to trd_tra_adv
      !!              - 'TRA' case : regroup T & S trends
      !!              - send the trends to trd_tra_mng (trdtrc) for further processing
      !!----------------------------------------------------------------------
      INTEGER,  DIMENSION(2)              , INTENT(in)           ::   kttrd, ktu, kttra
      INTEGER                             , INTENT(in)           ::   kt         ! time step
      CHARACTER(len=3)                    , INTENT(in)           ::   ctype      ! tracers trends type 'TRA'/'TRC'
      INTEGER                             , INTENT(in)           ::   ktra       ! tracer index
      INTEGER                             , INTENT(in)           ::   ktrd       ! tracer trend index
      INTEGER                             , INTENT(in)           ::   Kmm, Krhs  ! time level indices
      REAL(wp), DIMENSION(AB2D(kttrd),JPK), INTENT(in)           ::   ptrd       ! tracer trend  or flux
      REAL(wp), DIMENSION(AB2D(ktu),  JPK), INTENT(in), OPTIONAL ::   pu         ! now velocity
      REAL(wp), DIMENSION(AB2D(kttra),JPK), INTENT(in), OPTIONAL ::   ptra       ! now tracer variable
      !
      INTEGER  ::   ji,jj,jk    ! loop indices
      INTEGER  ::   i01   ! 0 or 1
      REAL(wp) ::  z1e3w
      REAL(wp),              DIMENSION(T2D(0),jpk) ::   ztrds             ! 3D workspace
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   zwt, zws, ztrdt   ! 3D workspace
      !!----------------------------------------------------------------------
      !      
      IF( .NOT. ALLOCATED( trdtx ) ) THEN      ! allocate trdtra arrays
         IF( trd_tra_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trd_tra : unable to allocate arrays' )
         avt_evd(T2D(0),:) = 0._wp
      ENDIF
      !
      i01 = COUNT( (/ PRESENT(pu) .OR. ( ktrd /= jptra_xad .AND. ktrd /= jptra_yad .AND. ktrd /= jptra_zad ) /) )
      !
      IF( ctype == 'TRA' .AND. ktra == jp_tem ) THEN   !==  Temperature trend  ==!
         !
         SELECT CASE( ktrd*i01 )
         !                            ! advection: transform the advective flux into a trend
         CASE( jptra_xad )   ;   CALL trd_tra_adv( ptrd, pu, ptra, 'X', trdtx, Kmm ) 
         CASE( jptra_yad )   ;   CALL trd_tra_adv( ptrd, pu, ptra, 'Y', trdty, Kmm ) 
         CASE( jptra_zad )   ;   CALL trd_tra_adv( ptrd, pu, ptra, 'Z', trdt, Kmm )
         CASE( jptra_bbc,    &        ! qsr, bbc: on temperature only, send to trd_tra_mng
            &  jptra_qsr )   ;   trdt(T2D(0),:) = ptrd(T2D(0),:) * tmask(T2D(0),:)
                                 ztrds(:,:,:) = 0._wp
                                 CALL trd_tra_mng( trdt, ztrds, ktrd, kt, Kmm, Krhs )
 !!gm Gurvan, verify the jptra_evd trend please !
         CASE( jptra_evd )   ;   avt_evd(T2D(0),:) = ptrd(T2D(0),:) * tmask(T2D(0),:)
         CASE DEFAULT                 ! other trends: masked trends
            trdt(T2D(0),:) = ptrd(T2D(0),:) * tmask(T2D(0),:)              ! mask & store
         END SELECT
         !
      ENDIF

      IF( ctype == 'TRA' .AND. ktra == jp_sal ) THEN      !==  Salinity trends  ==!
         !
         SELECT CASE( ktrd*i01 )
         !                            ! advection: transform the advective flux into a trend
         !                            !            and send T & S trends to trd_tra_mng
         CASE( jptra_xad  )   ;   CALL trd_tra_adv( ptrd , pu  , ptra, 'X'  , ztrds, Kmm, Krhs ) 
                                  CALL trd_tra_mng( trdtx, ztrds, ktrd, kt, Kmm   )
         CASE( jptra_yad  )   ;   CALL trd_tra_adv( ptrd , pu  , ptra, 'Y'  , ztrds, Kmm, Krhs ) 
                                  CALL trd_tra_mng( trdty, ztrds, ktrd, kt, Kmm   )
         CASE( jptra_zad  )   ;   CALL trd_tra_adv( ptrd , pu  , ptra, 'Z'  , ztrds, Kmm, Krhs ) 
                                  CALL trd_tra_mng( trdt , ztrds, ktrd, kt, Kmm   )
         CASE( jptra_zdfp )           ! diagnose the "PURE" Kz trend (here: just before the swap)
            !                         ! iso-neutral diffusion case otherwise jptra_zdf is "PURE"
            ALLOCATE( zwt(T2D(0),jpk), zws(T2D(0),jpk), ztrdt(T2D(0),jpk) )
            !
            zwt(:,:, 1 ) = 0._wp   ;   zws(:,:, 1 ) = 0._wp            ! vertical diffusive fluxes
            zwt(:,:,jpk) = 0._wp   ;   zws(:,:,jpk) = 0._wp
            DO_3D( 0, 0, 0, 0, 2, jpk )
               z1e3w = ( 1._wp / e3w(ji,jj,jk,Kmm) ) * tmask(ji,jj,jk)
               zwt(ji,jj,jk) = avt(ji,jj,jk) * ( ts(ji,jj,jk-1,jp_tem,Krhs) - ts(ji,jj,jk,jp_tem,Krhs) ) * z1e3w
               zws(ji,jj,jk) = avs(ji,jj,jk) * ( ts(ji,jj,jk-1,jp_sal,Krhs) - ts(ji,jj,jk,jp_sal,Krhs) ) * z1e3w
            END_3D
            !
            ztrdt(:,:,jpk) = 0._wp   ;   ztrds(:,:,jpk) = 0._wp
            DO_3D( 0, 0, 0, 0, 1, jpkm1 )
               z1e3w = 1._wp / e3w(ji,jj,jk,Kmm) 
               ztrdt(ji,jj,jk) = ( zwt(ji,jj,jk) - zwt(ji,jj,jk+1) ) * z1e3w
               ztrds(ji,jj,jk) = ( zws(ji,jj,jk) - zws(ji,jj,jk+1) ) * z1e3w
            END_3D
            CALL trd_tra_mng( ztrdt, ztrds, jptra_zdfp, kt, Kmm, Krhs )  
            !
            !                         ! Also calculate EVD trend at this point. 
            DO_3D( 0, 0, 0, 0, 2, jpk )
               z1e3w = ( 1._wp / e3w(ji,jj,jk,Kmm) ) * tmask(ji,jj,jk)
               zwt(ji,jj,jk) = avt_evd(ji,jj,jk) * ( ts(ji,jj,jk-1,jp_tem,Krhs) - ts(ji,jj,jk,jp_tem,Krhs) ) * z1e3w
               zws(ji,jj,jk) = avt_evd(ji,jj,jk) * ( ts(ji,jj,jk-1,jp_sal,Krhs) - ts(ji,jj,jk,jp_sal,Krhs) ) * z1e3w
            END_3D
            !
            DO_3D( 0, 0, 0, 0, 1, jpkm1 )
               z1e3w = 1._wp / e3w(ji,jj,jk,Kmm) 
               ztrdt(ji,jj,jk) = ( zwt(ji,jj,jk) - zwt(ji,jj,jk+1) ) * z1e3w
               ztrds(ji,jj,jk) = ( zws(ji,jj,jk) - zws(ji,jj,jk+1) ) * z1e3w
            END_3D
            CALL trd_tra_mng( ztrdt, ztrds, jptra_evd, kt, Kmm, Krhs )  
            !
            DEALLOCATE( zwt, zws, ztrdt )
            !
         CASE DEFAULT                 ! other trends: mask and send T & S trends to trd_tra_mng
            ztrds(:,:,:) = ptrd(T2D(0),:) * tmask(T2D(0),:)
            CALL trd_tra_mng( trdt, ztrds, ktrd, kt, Kmm, Krhs )  
         END SELECT
      ENDIF

      IF( ctype == 'TRC' ) THEN                           !==  passive tracer trend  ==!
         !
         SELECT CASE( ktrd*i01 )
         !                            ! advection: transform the advective flux into a masked trend
         CASE( jptra_xad )   ;   CALL trd_tra_adv( ptrd , pu , ptra, 'X', ztrds, Kmm ) 
         CASE( jptra_yad )   ;   CALL trd_tra_adv( ptrd , pu , ptra, 'Y', ztrds, Kmm ) 
         CASE( jptra_zad )   ;   CALL trd_tra_adv( ptrd , pu , ptra, 'Z', ztrds, Kmm ) 
         CASE DEFAULT                 ! other trends: just masked 
                                 ztrds(:,:,:) = ptrd(T2D(0),:) * tmask(T2D(0),:)
         END SELECT
         !                            ! send trend to trd_trc
         CALL trd_trc( ztrds, ktra, ktrd, kt, Kmm ) 
         !
      ENDIF
      !
   END SUBROUTINE trd_tra_t


   SUBROUTINE trd_tra_adv( pf, pu, pt, cdir, ptrd, Kmm )
      !!
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pf      ! advective flux in one direction
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pu      ! now velocity   in one direction
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pt      ! now or before tracer
      CHARACTER(len=1)          , INTENT(in   ) ::   cdir    ! X/Y/Z direction
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   ptrd    ! advective trend in one direction
      INTEGER,  INTENT(in)                      ::   Kmm     ! time level index
      !!
      CALL trd_tra_adv_t( pf, lbnd_ij(pf), pu, lbnd_ij(pu), pt, lbnd_ij(pt), cdir, ptrd, lbnd_ij(ptrd), Kmm )
   END SUBROUTINE trd_tra_adv


   SUBROUTINE trd_tra_adv_t( pf, ktf, pu, ktu, pt, ktt, cdir, ptrd, kttrd, Kmm )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_tra_adv  ***
      !! 
      !! ** Purpose :   transformed a advective flux into a masked advective trends
      !!
      !! ** Method  :   use the following transformation: -div(U.T) = - U grad(T) + T.div(U)
      !!       i-advective trends = -un. di-1[T] = -( di-1[fi] - tn di-1[un] )
      !!       j-advective trends = -un. di-1[T] = -( dj-1[fi] - tn dj-1[un] )
      !!       k-advective trends = -un. di+1[T] = -( dk+1[fi] - tn dk+1[un] )
      !!                where fi is the incoming advective flux.
      !!----------------------------------------------------------------------
      INTEGER,  DIMENSION(2)              , INTENT(in   ) ::   ktf, ktu, ktt, kttrd
      REAL(wp), DIMENSION(AB2D(ktf),  JPK), INTENT(in   ) ::   pf      ! advective flux in one direction
      REAL(wp), DIMENSION(AB2D(ktu),  JPK), INTENT(in   ) ::   pu      ! now velocity   in one direction
      REAL(wp), DIMENSION(AB2D(ktt),  JPK), INTENT(in   ) ::   pt      ! now or before tracer
      CHARACTER(len=1)                    , INTENT(in   ) ::   cdir    ! X/Y/Z direction
      REAL(wp), DIMENSION(AB2D(kttrd),JPK), INTENT(  out) ::   ptrd    ! advective trend in one direction
      INTEGER                             , INTENT(in   ) ::   Kmm     ! time level index
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ii, ij, ik   ! index shift as function of the direction
      !!----------------------------------------------------------------------
      !
      SELECT CASE( cdir )             ! shift depending on the direction
      CASE( 'X' )   ;   ii = 1   ;   ij = 0   ;   ik = 0      ! i-trend
      CASE( 'Y' )   ;   ii = 0   ;   ij = 1   ;   ik = 0      ! j-trend
      CASE( 'Z' )   ;   ii = 0   ;   ij = 0   ;   ik =-1      ! k-trend
      END SELECT
      !
      !                               ! set to zero uncomputed values
      ptrd(:,:,:) = 0._wp 
      !
      DO_3D( 0, 0, 0, 0, 1, jpkm1 )   ! advective trend
         ptrd(ji,jj,jk) = - (     pf (ji,jj,jk) - pf (ji-ii,jj-ij,jk-ik)                        &
           &                  - ( pu(ji,jj,jk) - pu(ji-ii,jj-ij,jk-ik) ) * pt(ji,jj,jk)  )   &
           &              * r1_e1e2t(ji,jj) / e3t(ji,jj,jk,Kmm) * tmask(ji,jj,jk)
      END_3D
      !
   END SUBROUTINE trd_tra_adv_t


   SUBROUTINE trd_tra_mng( ptrdx, ptrdy, ktrd, kt, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_tra_mng  ***
      !! 
      !! ** Purpose :   Dispatch all tracer trends computation, e.g. 3D output,
      !!                integral constraints, potential energy, and/or 
      !!                mixed layer budget.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(T2D(0),jpk), INTENT(inout) ::   ptrdx   ! Temperature or U trend
      REAL(wp), DIMENSION(T2D(0),jpk), INTENT(inout) ::   ptrdy   ! Salinity    or V trend
      INTEGER                        , INTENT(in   ) ::   ktrd    ! tracer trend index
      INTEGER                        , INTENT(in   ) ::   kt      ! time step
      INTEGER                        , INTENT(in   ) ::   Kmm, Krhs   ! time level index
      !!----------------------------------------------------------------------

      !                   ! 3D output of tracers trends using IOM interface
      IF( ln_tra_trd )   CALL trd_tra_iom ( ptrdx, ptrdy, ktrd, kt, Kmm )

      !                   ! Integral Constraints Properties for tracers trends           
      IF( ln_glo_trd )   CALL trd_glo( ptrdx, ptrdy, ktrd, 'TRA', kt, Kmm )

      !                   ! Potential ENergy trends
      IF( ln_PE_trd  )   CALL trd_pen( ptrdx, ptrdy, ktrd, kt, rDt, Kmm )

      !                   ! Mixed layer trends for active tracers
      IF( ln_tra_mxl )   THEN   
         !-----------------------------------------------------------------------------------------------
         ! W.A.R.N.I.N.G :
         ! jptra_ldf : called by traldf.F90
         !                 at this stage we store:
         !                  - the lateral geopotential diffusion (here, lateral = horizontal)
         !                  - and the iso-neutral diffusion if activated 
         ! jptra_zdf : called by trazdf.F90
         !                 * in case of iso-neutral diffusion we store the vertical diffusion component in the 
         !                   lateral trend including the K_z contrib, which will be removed later (see trd_mxl)
         !-----------------------------------------------------------------------------------------------

         SELECT CASE ( ktrd )
         CASE ( jptra_xad )        ;   CALL trd_tra_mxl( ptrdx, ptrdy, jpmxl_xad, kt, Kmm )   ! zonal    advection
         CASE ( jptra_yad )        ;   CALL trd_tra_mxl( ptrdx, ptrdy, jpmxl_yad, kt, Kmm )   ! merid.   advection
         CASE ( jptra_zad )        ;   CALL trd_tra_mxl( ptrdx, ptrdy, jpmxl_zad, kt, Kmm )   ! vertical advection
         CASE ( jptra_ldf )        ;   CALL trd_tra_mxl( ptrdx, ptrdy, jpmxl_ldf, kt, Kmm )   ! lateral  diffusion
         CASE ( jptra_bbl )        ;   CALL trd_tra_mxl( ptrdx, ptrdy, jpmxl_bbl, kt, Kmm )   ! bottom boundary layer
         CASE ( jptra_zdf )
            IF( ln_traldf_iso ) THEN ; CALL trd_tra_mxl( ptrdx, ptrdy, jpmxl_ldf, kt, Kmm )   ! lateral  diffusion (K_z)
            ELSE                   ;   CALL trd_tra_mxl( ptrdx, ptrdy, jpmxl_zdf, kt, Kmm )   ! vertical diffusion (K_z)
            ENDIF
         CASE ( jptra_dmp )        ;   CALL trd_tra_mxl( ptrdx, ptrdy, jpmxl_dmp, kt, Kmm )   ! internal 3D restoring (tradmp)
         CASE ( jptra_qsr )        ;   CALL trd_tra_mxl( ptrdx, ptrdy, jpmxl_for, kt, Kmm )   ! air-sea : penetrative sol radiat
         CASE ( jptra_nsr )        ;   ptrdx(:,:,2:jpk) = 0._wp   ;   ptrdy(:,:,2:jpk) = 0._wp
                                       CALL trd_tra_mxl( ptrdx, ptrdy, jpmxl_for, kt, Kmm )   ! air-sea : non penetr sol radiation
         CASE ( jptra_bbc )        ;   CALL trd_tra_mxl( ptrdx, ptrdy, jpmxl_bbc, kt, Kmm )   ! bottom bound cond (geoth flux)
         CASE ( jptra_npc )        ;   CALL trd_tra_mxl( ptrdx, ptrdy, jpmxl_npc, kt, Kmm )   ! non penetr convect adjustment
         CASE ( jptra_atf )        ;   CALL trd_tra_mxl( ptrdx, ptrdy, jpmxl_atf, kt, Kmm )   ! asselin time filter (last trend)
                                   !
                                       CALL trd_mxl( kt, rDt, Kmm, Krhs )                     ! trends: Mixed-layer (output) 
         END SELECT
         !
      ENDIF
      !
   END SUBROUTINE trd_tra_mng


   SUBROUTINE trd_tra_iom( ptrdx, ptrdy, ktrd, kt, Kmm )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_tra_iom  ***
      !! 
      !! ** Purpose :   output 3D tracer trends using IOM
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(T2D(0),jpk), INTENT(inout) ::   ptrdx   ! Temperature or U trend
      REAL(wp), DIMENSION(T2D(0),jpk), INTENT(inout) ::   ptrdy   ! Salinity    or V trend
      INTEGER                        , INTENT(in   ) ::   ktrd    ! tracer trend index
      INTEGER                        , INTENT(in   ) ::   kt      ! time step
      INTEGER                        , INTENT(in   ) ::   Kmm     ! time level index
      !!
      INTEGER  :: ji, jj
      REAL(wp) :: z1e3
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   z2dx, z2dy   ! 2D workspace 
      !!----------------------------------------------------------------------
      !
!!gm Rq: mask the trends already masked in trd_tra, but lbc_lnk should probably be added
      !
      ! Trends evaluated every time step that could go to the standard T file and can be output every ts into a 1ts file if 1ts output is selected
      SELECT CASE( ktrd )
      ! This total trend is done every time step
      CASE( jptra_tot  )   ;   CALL iom_put( "ttrd_tot" , ptrdx )           ! model total trend
                               CALL iom_put( "strd_tot" , ptrdy )
      END SELECT
      !
      ! These trends are done every second time step. When 1ts output is selected must go different (2ts) file from standard T-file
      IF( MOD( kt, 2 ) == 0 ) THEN
         SELECT CASE( ktrd )
         CASE( jptra_xad  )   ;   CALL iom_put( "ttrd_xad"  , ptrdx )        ! x- horizontal advection
                                  CALL iom_put( "strd_xad"  , ptrdy )
         CASE( jptra_yad  )   ;   CALL iom_put( "ttrd_yad"  , ptrdx )        ! y- horizontal advection
                                  CALL iom_put( "strd_yad"  , ptrdy )
         CASE( jptra_zad  )   ;   CALL iom_put( "ttrd_zad"  , ptrdx )        ! z- vertical   advection
                                  CALL iom_put( "strd_zad"  , ptrdy )
                                  IF( lk_linssh ) THEN                   ! cst volume : adv flux through z=0 surface
                                     ALLOCATE( z2dx(T2D(0)), z2dy(T2D(0)) )
                                     DO_2D( 0, 0, 0, 0 )
                                        z1e3 = ww(ji,jj,1) / e3t(ji,jj,1,Kmm)
                                        z2dx(ji,jj) = ts(ji,jj,1,jp_tem,Kmm) * z1e3
                                        z2dy(ji,jj) = ts(ji,jj,1,jp_sal,Kmm) * z1e3
                                     END_2D
                                     CALL iom_put( "ttrd_sad", z2dx )
                                     CALL iom_put( "strd_sad", z2dy )
                                     DEALLOCATE( z2dx, z2dy )
                                  ENDIF
         CASE( jptra_totad  ) ;   CALL iom_put( "ttrd_totad", ptrdx )        ! total   advection
                                  CALL iom_put( "strd_totad", ptrdy )
         CASE( jptra_ldf  )   ;   CALL iom_put( "ttrd_ldf"  , ptrdx )        ! lateral diffusion
                                  CALL iom_put( "strd_ldf"  , ptrdy )
         CASE( jptra_zdf  )   ;   CALL iom_put( "ttrd_zdf"  , ptrdx )        ! vertical diffusion (including Kz contribution)
                                  CALL iom_put( "strd_zdf"  , ptrdy )
         CASE( jptra_zdfp )   ;   CALL iom_put( "ttrd_zdfp" , ptrdx )        ! PURE vertical diffusion (no isoneutral contribution)
                                  CALL iom_put( "strd_zdfp" , ptrdy )
         CASE( jptra_evd )    ;   CALL iom_put( "ttrd_evd"  , ptrdx )        ! EVD trend (convection)
                                  CALL iom_put( "strd_evd"  , ptrdy )
         CASE( jptra_dmp  )   ;   CALL iom_put( "ttrd_dmp"  , ptrdx )        ! internal restoring (damping)
                                  CALL iom_put( "strd_dmp"  , ptrdy )
         CASE( jptra_bbl  )   ;   CALL iom_put( "ttrd_bbl"  , ptrdx )        ! bottom boundary layer
                                  CALL iom_put( "strd_bbl"  , ptrdy )
         CASE( jptra_npc  )   ;   CALL iom_put( "ttrd_npc"  , ptrdx )        ! static instability mixing
                                  CALL iom_put( "strd_npc"  , ptrdy )
         CASE( jptra_bbc  )   ;   CALL iom_put( "ttrd_bbc"  , ptrdx )        ! geothermal heating   (only on temperature)
         CASE( jptra_nsr  )   ;   CALL iom_put( "ttrd_qns"  , ptrdx(:,:,1) ) ! surface forcing + runoff (ln_rnf=T)
                                  CALL iom_put( "strd_cdt"  , ptrdy(:,:,1) )        ! output as 2D surface fields
         CASE( jptra_qsr  )   ;   CALL iom_put( "ttrd_qsr"  , ptrdx )        ! penetrative solar radiat. (only on temperature)
         END SELECT
         ! the Asselin filter trend  is also every other time step but needs to be lagged one time step
         ! Even when 1ts output is selected can go to the same (2ts) file as the trends plotted every even time step.
      ELSE IF( MOD( kt, 2 ) == 1 ) THEN
         SELECT CASE( ktrd )
         CASE( jptra_atf  )   ;   CALL iom_put( "ttrd_atf" , ptrdx )        ! asselin time Filter
                                  CALL iom_put( "strd_atf" , ptrdy )
         END SELECT
      END IF
      !
   END SUBROUTINE trd_tra_iom

   !!======================================================================
END MODULE trdtra

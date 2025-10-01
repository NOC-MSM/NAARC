MODULE traldf
   !!======================================================================
   !!                       ***  MODULE  traldf  ***
   !! Ocean Active tracers : lateral diffusive trends
   !!=====================================================================
   !! History :  9.0  ! 2005-11  (G. Madec)  Original code
   !!  NEMO      3.0  ! 2008-01  (C. Ethe, G. Madec)  merge TRC-TRA
   !!            3.7  ! 2013-12  (G. Madec) remove the optional computation from T & S anomaly profiles and traldf_bilapg
   !!             -   ! 2013-12  (F. Lemarie, G. Madec)  triad operator (Griffies) + Method of Stabilizing Correction
   !!             -   ! 2014-01  (G. Madec, S. Masson)  restructuration/simplification of lateral diffusive operators
   !!            4.5  ! 2022-08  (G, Madec, S. Techene)  refactorization to reduce local memory usage
   !!                 !                                + add tra_ldf_iso_a33 routine and traldf_iso_h/z.h90 files
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_ldf       : update the tracer trend with the lateral diffusion trend
   !!   tra_ldf_init  : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE ldftra         ! lateral diffusion: eddy diffusivity & EIV coeff.
   USE ldfslp         ! lateral diffusion: iso-neutral slope
   USE traldf_iso     ! lateral diffusion: laplacian iso-neutral standard operator  (traldf_iso_lap/_blp  routines)
   USE traldf_lev     ! lateral diffusion: laplacian iso-level            operator  (traldf_lap/_blp      routines)
   USE traldf_triad   ! lateral diffusion: laplacian iso-neutral triad    operator  (tra_ldf_triad(_blp)  routines)
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! ocean active tracers trends
   !
   USE prtctl         ! Print control
   USE in_out_manager ! I/O manager
   USE iom            ! I/O library
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_ldf        ! called by step.F90
   PUBLIC   tra_ldf_init   ! called by nemogcm.F90

   LOGICAL  ::   l_ptr   ! flag to compute the diffusive part of poleward heat & salt transport
   LOGICAL  ::   l_hst   ! flag to compute the diffusive part of heat and salt transport

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_ldf( kt, Kbb, Kmm, pts, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf  ***
      !!
      !! ** Purpose :   compute the lateral ocean tracer physics.
      !!----------------------------------------------------------------------
      INTEGER,                                   INTENT(in   ) :: kt              ! ocean time-step index
      INTEGER,                                   INTENT(in   ) :: Kbb, Kmm, Krhs  ! ocean time level indices
      INTEGER  ::  ji, jj, jk   ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts,jpt), INTENT(inout) :: pts             ! active tracers and RHS of tracer equation
      !!
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   zTtrd, zStrd
      !!----------------------------------------------------------------------
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 .AND. lwp )  THEN
            WRITE(numout,*)
            SELECT CASE ( nldf_tra )   !* compute lateral mixing trend and add it to the general trend
               CASE ( np_lap    )   ;   WRITE(numout,*) 'traldf_lev_lap   : iso-level laplacian diffusive operator'
               CASE ( np_lap_i  )   ;   WRITE(numout,*) 'traldf_iso_lap   : iso-neutral standard laplacian diffusive operator'
               CASE ( np_lap_it )   ;   WRITE(numout,*) 'traldf_triad_lap : iso-neutral triad laplacian diffusive operator'
               CASE ( np_blp    )   ;   WRITE(numout,*) 'traldf_lev_blp   : iso-level bilaplacian diffusive operator'
               CASE ( np_blp_i  )   ;   WRITE(numout,*) 'traldf_iso_blp   : iso-neutral bilaplacian diffusive operator'
               CASE ( np_blp_it )   ;   WRITE(numout,*) 'traldf_triad_blp : iso-neutral triad bilaplacian diffusive operator'
            END SELECT
            WRITE(numout,*) '~~~~~~~~~~~~~~~~ '
         ENDIF
      ENDIF
      !
      IF( ln_timing )   CALL timing_start('tra_ldf')
      !
      IF( l_trdtra )   THEN                    !* Save ta and sa trends
         ALLOCATE( zTtrd(T2D(0),jpk) , zStrd(T2D(0),jpk) )
         zTtrd(:,:,:) = pts(T2D(0),:,jp_tem,Krhs)
         zStrd(:,:,:) = pts(T2D(0),:,jp_sal,Krhs)
      ENDIF
      !
      SELECT CASE ( nldf_tra )   !* compute lateral mixing trend and add it to the general trend
      !                                !-  laplacian  - !
      CASE ( np_lap    )                     ! level operator
         CALL traldf_lev_lap  ( kt, Kbb, Kmm, pts, Krhs, l_ptr, l_hst )
      CASE ( np_lap_i  )                     ! standard iso-neutral operator
         CALL traldf_iso_lap  ( kt, Kbb, Kmm, pts, Krhs, l_ptr, l_hst )
      CASE ( np_lap_it )                     ! laplacian: triad iso-neutral operator
         CALL traldf_triad_lap( kt, Kmm, nit000,'TRA', ahtu, ahtv, pts(:,:,:,:,Kbb),    &
            &                                    pts(:,:,:,:,Kbb), pts(:,:,:,:,Krhs), jpts,  1 )
      !                                !- bilaplacian - !
      CASE ( np_blp    )                     ! iso-level operators
         CALL traldf_lev_blp   ( kt, Kbb, Kmm, pts, Krhs, l_ptr, l_hst )
      CASE ( np_blp_i  )                     ! standard iso-neutral operator
         CALL traldf_iso_blp   ( kt, Kbb, Kmm, pts, Krhs, l_ptr, l_hst )
      CASE ( np_blp_it )                     ! bilaplacian: iso-level & iso-neutral operators
         CALL traldf_triad_blp( kt, Kmm, nit000,'TRA', ahtu, ahtv, pts(:,:,:,:,Kbb),    &
            &                                     pts(:,:,:,:,Krhs),             jpts )
      END SELECT
      !
      IF( l_trdtra )   THEN                    !* save the horizontal diffusive trends for further diagnostics
         zTtrd(:,:,:) = pts(T2D(0),:,jp_tem,Krhs) - zTtrd(:,:,:)
         zStrd(:,:,:) = pts(T2D(0),:,jp_sal,Krhs) - zStrd(:,:,:)
         CALL trd_tra( kt, Kmm, Krhs, 'TRA', jp_tem, jptra_ldf, zTtrd )
         CALL trd_tra( kt, Kmm, Krhs, 'TRA', jp_sal, jptra_ldf, zStrd )
         DEALLOCATE( zTtrd, zStrd )
      ENDIF
      
      !DO_2D( nn_hls, nn_hls, nn_hls, nn_hls)
      !   IF( mjg0(jj) == 3605 ) THEN
      !      DO jk = 1, mbkt(ji,jj)
      !         pts(ji,jj,jk,jp_sal,Krhs) = MIN(35._wp, pts(ji,jj,jk,jp_sal,Krhs))
      !      END DO
      !   ENDIF
      !END_2D
      
      !                                        !* print mean trends (used for debugging)
      IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab3d_1=pts(:,:,:,jp_tem,Krhs), clinfo1=' ldf  - Ta: ', mask1=tmask, &
         &                                  tab3d_2=pts(:,:,:,jp_sal,Krhs), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      !
      IF( ln_timing )   CALL timing_stop('tra_ldf')
      !
   END SUBROUTINE tra_ldf


   SUBROUTINE tra_ldf_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_init  ***
      !!
      !! ** Purpose :   Choice of the operator for the lateral tracer diffusion
      !!
      !! ** Method  :   set nldf_tra from the namtra_ldf logicals
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio, ierr   ! temporary integers
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                     !==  Namelist print  ==!
         WRITE(numout,*)
         WRITE(numout,*) 'tra_ldf_init : lateral tracer diffusive operator'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtra_ldf: already read in ldftra module'
         WRITE(numout,*) '      see ldf_tra_init report for lateral mixing parameters'
         WRITE(numout,*)
         !
         SELECT CASE( nldf_tra )             ! print the choice of operator
         CASE( np_no_ldf )   ;   WRITE(numout,*) '   ==>>>   NO lateral diffusion'
         CASE( np_lap    )   ;   WRITE(numout,*) '   ==>>>   laplacian iso-level operator'
         CASE( np_lap_i  )   ;   WRITE(numout,*) '   ==>>>   Rotated laplacian operator (standard)'
         CASE( np_lap_it )   ;   WRITE(numout,*) '   ==>>>   Rotated laplacian operator (triad)'
         CASE( np_blp    )   ;   WRITE(numout,*) '   ==>>>   bilaplacian iso-level operator'
         CASE( np_blp_i  )   ;   WRITE(numout,*) '   ==>>>   Rotated bilaplacian operator (standard)'
         CASE( np_blp_it )   ;   WRITE(numout,*) '   ==>>>   Rotated bilaplacian operator (triad)'
         END SELECT
      ENDIF
      !
      l_ptr = .FALSE.                        ! set flag for heat & salt diffusive diagnostics
      l_hst = .FALSE.
      IF(   ( iom_use( 'sophtldf'  ) .OR. iom_use(  'sopstldf' )  )   )   l_ptr = .TRUE.   ! diffusive poleward transport
      IF(   ( iom_use("uadv_heattr") .OR. iom_use("vadv_heattr") .OR. &                 
         &    iom_use("uadv_salttr") .OR. iom_use("vadv_salttr")  )   )   l_hst = .TRUE.   ! vertically cumulated diffusive fluxes
      !
   END SUBROUTINE tra_ldf_init

   !!======================================================================
END MODULE traldf

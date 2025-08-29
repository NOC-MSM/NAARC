MODULE bdydyn
   !!======================================================================
   !!                       ***  MODULE  bdydyn  ***
   !! Unstructured Open Boundary Cond. :   Apply boundary conditions to velocities
   !!======================================================================
   !! History :  1.0  !  2005-02  (J. Chanut, A. Sellar)  Original code
   !!             -   !  2007-07  (D. Storkey) Move Flather implementation to separate routine.
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.2  !  2008-04  (R. Benshila) consider velocity instead of transport 
   !!            3.3  !  2010-09  (E.O'Dea) modifications for Shelf configurations 
   !!            3.3  !  2010-09  (D.Storkey) add ice boundary conditions
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!----------------------------------------------------------------------
   !!   bdy_dyn        : split velocities into barotropic and baroclinic parts
   !!                    and call bdy_dyn2d and bdy_dyn3d to apply boundary
   !!                    conditions
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE bdy_oce         ! ocean open boundary conditions
   USE bdydyn2d        ! open boundary conditions for barotropic solution
   USE bdydyn3d        ! open boundary conditions for baroclinic velocities
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  !

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_dyn    ! routine called in dyn_nxt

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_dyn( kt, Kbb, puu, pvv, Kaa, dyn3d_only )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn  ***
      !!
      !! ** Purpose : - Wrapper routine for bdy_dyn2d and bdy_dyn3d.
      !!
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT(in)    ::   kt           ! Main time step counter
      INTEGER                             , INTENT(in)    ::   Kbb, Kaa     ! Ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::   puu, pvv     ! Ocean velocities (to be updated at open boundaries)
      LOGICAL, OPTIONAL                   , INTENT(in)    ::   dyn3d_only   ! T => only update baroclinic velocities
      !
      INTEGER ::   jk, ii, ij, ib_bdy, ib, igrd     ! Loop counter
      LOGICAL ::   ll_dyn2d, ll_dyn3d, ll_orlanski
      REAL(wp), DIMENSION(jpi,jpj) :: zua2d, zva2d     ! after barotropic velocities
      !!----------------------------------------------------------------------
      !
      ll_dyn2d = .true.
      ll_dyn3d = .true.
      !
      IF( PRESENT(dyn3d_only) ) THEN
         IF( dyn3d_only )   ll_dyn2d = .false.
      ENDIF
      !
      ll_orlanski = .false.
      DO ib_bdy = 1, nb_bdy
         IF ( cn_dyn2d(ib_bdy) == 'orlanski' .OR. cn_dyn2d(ib_bdy) == 'orlanski_npo' &
     &   .OR. cn_dyn3d(ib_bdy) == 'orlanski' .OR. cn_dyn3d(ib_bdy) == 'orlanski_npo')   ll_orlanski = .true.
      END DO

      !-------------------------------------------------------
      ! Split velocities into barotropic and baroclinic parts
      !-------------------------------------------------------

      !                          ! "After" velocities: 
      zua2d(:,:) = 0._wp
      zva2d(:,:) = 0._wp     
      DO jk = 1, jpkm1
         zua2d(:,:) = zua2d(:,:) + e3u(:,:,jk,Kaa) * puu(:,:,jk,Kaa) * umask(:,:,jk)
         zva2d(:,:) = zva2d(:,:) + e3v(:,:,jk,Kaa) * pvv(:,:,jk,Kaa) * vmask(:,:,jk)
      END DO
      zua2d(:,:) = zua2d(:,:) * r1_hu(:,:,Kaa)
      zva2d(:,:) = zva2d(:,:) * r1_hv(:,:,Kaa)

      DO jk = 1 , jpkm1
         puu(:,:,jk,Kaa) = ( puu(:,:,jk,Kaa) - zua2d(:,:) ) * umask(:,:,jk)
         pvv(:,:,jk,Kaa) = ( pvv(:,:,jk,Kaa) - zva2d(:,:) ) * vmask(:,:,jk)
      END DO


      IF( ll_orlanski ) THEN     ! "Before" velocities (Orlanski condition only) 
         DO jk = 1 , jpkm1
            puu(:,:,jk,Kbb) = ( puu(:,:,jk,Kbb) - uu_b(:,:,Kbb) ) * umask(:,:,jk)
            pvv(:,:,jk,Kbb) = ( pvv(:,:,jk,Kbb) - vv_b(:,:,Kbb) ) * vmask(:,:,jk)
         END DO
      ENDIF

      !-------------------------------------------------------
      ! Apply boundary conditions to barotropic and baroclinic
      ! parts separately
      !-------------------------------------------------------

      IF( ll_dyn2d )   CALL bdy_dyn2d( kt, zua2d, zva2d, uu_b(:,:,Kbb), vv_b(:,:,Kbb), r1_hu(:,:,Kaa), r1_hv(:,:,Kaa), ssh(:,:,Kaa) )

      IF( ll_dyn3d )   CALL bdy_dyn3d( kt, Kbb, puu, pvv, Kaa )

      !-------------------------------------------------------
      ! Recombine velocities
      !-------------------------------------------------------
      !
      DO jk = 1 , jpkm1
         puu(:,:,jk,Kaa) = ( puu(:,:,jk,Kaa) + zua2d(:,:) ) * umask(:,:,jk)
         pvv(:,:,jk,Kaa) = ( pvv(:,:,jk,Kaa) + zva2d(:,:) ) * vmask(:,:,jk)
      END DO
      !
      IF ( ll_orlanski ) THEN
         DO jk = 1 , jpkm1
            puu(:,:,jk,Kbb) = ( puu(:,:,jk,Kbb) + uu_b(:,:,Kbb) ) * umask(:,:,jk)
            pvv(:,:,jk,Kbb) = ( pvv(:,:,jk,Kbb) + vv_b(:,:,Kbb) ) * vmask(:,:,jk)
         END DO
      END IF
      !
   END SUBROUTINE bdy_dyn

   !!======================================================================
END MODULE bdydyn

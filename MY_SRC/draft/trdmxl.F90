MODULE trdmxl
   !!======================================================================
   !!                       ***  MODULE  trdmxl  ***
   !! Ocean diagnostics:  mixed layer T-S trends 
   !!======================================================================
   !! History :  OPA  !  1995-04  (J. Vialard)    Original code
   !!                 !  1997-02  (E. Guilyardi)  Adaptation global + base cmo
   !!                 !  1999-09  (E. Guilyardi)  Re-writing + netCDF output
   !!   NEMO     1.0  !  2002-06  (G. Madec)      F90: Free form and module
   !!             -   !  2004-08  (C. Talandier)  New trends organization
   !!            2.0  !  2005-05  (C. Deltel)     Diagnose trends of time averaged ML T & S
   !!            3.5  !  2012-03  (G. Madec)      complete reorganisation + change in the time averaging
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   trd_mxl          : T and S cumulated trends averaged over the mixed layer
   !!   trd_mxl_zint     : T and S trends vertical integration
   !!   trd_mxl_init     : initialization step
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE trd_oce         ! trends: ocean variables
   USE trdmxl_oce      ! ocean variables trends
   USE ldftra          ! lateral diffusion: eddy diffusivity & EIV coeff.
   USE zdf_oce         ! ocean vertical physics
   USE phycst          ! Define parameters for the routines
   USE dianam          ! build the name of file (routine)
   USE ldfslp          ! iso-neutral slopes 
   USE zdfmxl          ! mixed layer depth
   USE zdfddm          ! ocean vertical physics: double diffusion
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE trdmxl_rst      ! restart for diagnosing the ML trends
   !
   USE in_out_manager  ! I/O manager
   USE ioipsl          ! NetCDF library
   USE prtctl          ! Print control
   USE restart         ! for lrst_oce
   USE lib_mpp         ! MPP library
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trd_mxl        ! routine called by step.F90
   PUBLIC   trd_mxl_init   ! routine called by opa.F90
   PUBLIC   trd_tra_mxl    ! routine called by tracers routines

   INTEGER ::   nkstp      ! current time step 
   
   CHARACTER (LEN=40) ::  clhstnam         ! name of the trends NetCDF file
   INTEGER ::   nh_t, nmoymltrd
   INTEGER ::   nidtrd
   INTEGER ::   ndimtrd1                        
   INTEGER ::   ionce, icount                   

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE trd_tra_mxl( ptrdx, ptrdy, ktrd, kt, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_tra_mng  ***
      !! 
      !! ** Purpose :   Dispatch all trends computation, e.g. 3D output, integral
      !!                constraints, barotropic vorticity, kinetic enrgy, 
      !!                potential energy, and/or mixed layer budget.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   ptrdx   ! Temperature or U trend 
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   ptrdy   ! Salinity    or V trend
      INTEGER                   , INTENT(in   ) ::   ktrd    ! tracer trend index
      INTEGER                   , INTENT(in   ) ::   kt      ! time step index
      INTEGER                   , INTENT(in   ) ::   Kmm     ! time level index
      !
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      !!----------------------------------------------------------------------

      !                         !==============================!
      IF ( kt /= nkstp ) THEN   !=  1st call at kt time step  =!
         !                      !==============================!
         !

         !!gm BUG?
         !!gm CAUTION:  double check the definition of nmln  it is the nb of w-level, not t-level I guess

         ! ... Set nmxl(ji,jj) = index of first T point below control surf. or outside mixed-layer
         IF( nn_ctls == 0 ) THEN       ! * control surface = mixed-layer with density criterion
            nmxl(:,:) = nmln(:,:)    ! array nmln computed in zdfmxl.F90
         ELSEIF( nn_ctls == 1 ) THEN   ! * control surface = read index from file
            nmxl(:,:) = nbol(:,:)
         ELSEIF( nn_ctls >= 2 ) THEN   ! * control surface = model level
            nn_ctls = MIN( nn_ctls, jpktrd - 1 )
            nmxl(:,:) = nn_ctls + 1
         ENDIF

         nkstp = kt
         !                          !==  reset trend arrays to zero  ==!
         DO_3D( 0, 0, 0, 0, 1, jpktrd )
            tmltrd(ji,jj,jk) = 0._wp
            smltrd(ji,jj,jk) = 0._wp
         END_3D
         !
         !                          !==  now ML weights for vertical averaging  ==!
         DO_3D( 0, 0, 0, 0, 1, jpktrd )  ! initialize wkx with vertical scale factor in mixed-layer
            IF( jk - nmxl(ji,jj) < 0 )   THEN
               wkx(ji,jj,jk) = e3t(ji,jj,jk,Kmm) * tmask(ji,jj,jk)
            ELSE
               wkx(ji,jj,jk) = 0._wp
            ENDIF
         END_3D
         !
         hmxl(:,:) = 0._wp               ! NOW mixed-layer depth
         DO_3D( 0, 0, 0, 0, 1, jpktrd )
            hmxl(ji,jj)    = hmxl(ji,jj) + wkx(ji,jj,jk)
            wkx (ji,jj,jk) = wkx (ji,jj,jk) / MAX( 1.e-20_wp, hmxl(ji,jj) ) * tmask(ji,jj,1)
         END_3D
         !
         !                          !==  Vertically averaged T and S  ==!
         tml(:,:) = 0._wp   ;   sml(:,:) = 0._wp
         DO_3D( 0, 0, 0, 0, 1, jpktrd )
            tml(ji,jj) = tml(ji,jj) + wkx(ji,jj,jk) * ts(ji,jj,jk,jp_tem,Kmm)
            sml(ji,jj) = sml(ji,jj) + wkx(ji,jj,jk) * ts(ji,jj,jk,jp_sal,Kmm)
         END_3D
         !
      ENDIF

      ! mean now trends over the now ML 
      DO_2D( 0, 0, 0, 0 )
         tmltrd(ji,jj,ktrd) = tmltrd(ji,jj,ktrd) + ptrdx(ji,jj,jk) * wkx(ji,jj,jk)   ! temperature
         smltrd(ji,jj,ktrd) = smltrd(ji,jj,ktrd) + ptrdy(ji,jj,jk) * wkx(ji,jj,jk)   ! salinity
      END_2D

 
!!gm to be put juste before the output !
!      ! Lateral boundary conditions
!      CALL lbc_lnk( 'trdmxl', tmltrd(:,:,jl), 'T', 1.0_wp , smltrd(:,:,jl), 'T', 1.0_wp )
!!gm end



         SELECT CASE( ktrd )
         CASE( jptra_npc  )               ! non-penetrative convection: regrouped with zdf
!!gm : to be completed ! 
!		   IF( ....
!!gm end
         CASE( jptra_zdfp )               ! iso-neutral diffusion: "pure" vertical diffusion
!                                   ! regroup iso-neutral diffusion in one term
         tmltrd(:,:,jpmxl_ldf) = tmltrd(:,:,jpmxl_ldf) + ( tmltrd(:,:,jpmxl_zdf) - tmltrd(:,:,jpmxl_zdfp) )
         smltrd(:,:,jpmxl_ldf) = smltrd(:,:,jpmxl_ldf) + ( smltrd(:,:,jpmxl_zdf) - smltrd(:,:,jpmxl_zdfp) )
         !                                   ! put in zdf the dia-neutral diffusion
         tmltrd(:,:,jpmxl_zdf) = tmltrd(:,:,jpmxl_zdfp)
         smltrd(:,:,jpmxl_zdf) = smltrd(:,:,jpmxl_zdfp)
         IF( ln_zdfnpc ) THEN
            tmltrd(:,:,jpmxl_zdf) = tmltrd(:,:,jpmxl_zdf) + tmltrd(:,:,jpmxl_npc)
            smltrd(:,:,jpmxl_zdf) = smltrd(:,:,jpmxl_zdf) + smltrd(:,:,jpmxl_npc)
         ENDIF
         !
      CASE( jptra_atf  )               ! last trends of the current time step: perform the time averaging & output
         !
         ! after ML           :   zhmla                      NB will be swaped to provide hmln and hmlb
         !
         ! entrainement ent_1 :   tb_mln - tb_mlb        ==>> use previous timestep ztn_mla = tb_mln
         !                                                    "     "         "     tn_mln = tb_mlb  (unfiltered tb!)
         !                                                   NB: tn_mln itself comes from the 2 time step before (ta_mla)
         !
         ! atf trend          :   ztbf_mln - tb_mln      ==>> use previous timestep tn_mla = tb_mln
         !                                                   need to compute tbf_mln, using the current tb
         !                                                   which is the before fitered tracer
         !
         ! entrainement ent_2 :   zta_mla - zta_mln      ==>> need to compute zta_mla and zta_mln
         !
         ! time averaging     :   mean: CALL trd_mean( kt, ptrd, ptrdm )
         !                              and out put the starting mean value and the total trends
         !                              (i.e. difference between starting and ending values)
         !                        hat : CALL trd_hat ( kt, ptrd, ptrdm )
         !                              and output the starting hat value and the total hat trends
         !
         ! swaps              :   hmlb   <==   hmln   <== zhmla
         !                        tb_mlb <==  tn_mln  <== zta_mla
         !                        tb_mln <== ztn_mla     ==>> now T over after h, need to be computed here 
         !                                                    to be used at next time step (unfiltered before)
         !
      END SELECT
      !
   END SUBROUTINE trd_tra_mxl


   SUBROUTINE trd_mean( kt, ptrd, ptrdm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_mean  ***
      !! 
      !! ** Purpose :   
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   ptrd    ! trend at kt
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   ptrdm   ! cumulative trends at kt
      INTEGER                   , INTENT(in   ) ::   kt      ! time step index
      !!----------------------------------------------------------------------
      !
      IF ( kt == nn_it000 )   ptrdm(:,:,:) = 0._wp
      !
      ptrdm(:,:,:) = ptrdm(:,:,:) + ptrd(:,:,:)
      !
      IF ( MOD( kt - nn_it000 + 1, nn_trd ) == 0 ) THEN
         !
         ! call iom put????  avec en argument le tableau de nom des trends?
         !
      ENDIF
      !
   END SUBROUTINE trd_mean


   SUBROUTINE trd_mxl( kt, p2dt, Kmm, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_mxl  ***
      !! 
      !! ** Purpose :  Compute and cumulate the mixed layer trends over an analysis
      !!               period, and write NetCDF outputs.
      !!
      !! ** Method/usage :
      !!          The stored trends can be chosen twofold (according to the ln_trdmxl_instant 
      !!          logical namelist variable) :
      !!          1) to explain the difference between initial and final 
      !!             mixed-layer T & S (where initial and final relate to the
      !!             current analysis window, defined by nn_trd in the namelist)
      !!          2) to explain the difference between the current and previous 
      !!             TIME-AVERAGED mixed-layer T & S (where time-averaging is
      !!             performed over each analysis window).
      !!
      !! ** Consistency check : 
      !!        If the control surface is fixed ( nn_ctls > 1 ), the residual term (dh/dt
      !!        entrainment) should be zero, at machine accuracy. Note that in the case
      !!        of time-averaged mixed-layer fields, this residual WILL NOT BE ZERO
      !!        over the first two analysis windows (except if restart).
      !!        N.B. For ORCA2_ICE, use e.g. nn_trd=5, rn_ucf=1., nn_ctls=8
      !!             for checking residuals.
      !!             On a NEC-SX5 computer, this typically leads to:
      !!                   O(1.e-20) temp. residuals (tml_res) when ln_trdmxl_instant=.false.
      !!                   O(1.e-21) temp. residuals (tml_res) when ln_trdmxl_instant=.true.
      !!
      !! ** Action :
      !!       At each time step, mixed-layer averaged trends are stored in the 
      !!       tmltrd(:,:,jpmxl_xxx) array (see trdmxl_oce.F90 for definitions of jpmxl_xxx).
      !!       This array is known when trd_mxl is called, at the end of the stp subroutine, 
      !!       except for the purely vertical K_z diffusion term, which is embedded in the
      !!       lateral diffusion trend.
      !!
      !!       In I), this K_z term is diagnosed and stored, thus its contribution is removed
      !!       from the lateral diffusion trend.
      !!       In II), the instantaneous mixed-layer T & S are computed, and misc. cumulative
      !!       arrays are updated.
      !!       In III), called only once per analysis window, we compute the total trends,
      !!       along with the residuals and the Asselin correction terms.
      !!       In IV), the appropriate trends are written in the trends NetCDF file.
      !!
      !! References :  Vialard et al.,2001, JPO.
      !!----------------------------------------------------------------------
      INTEGER , INTENT(in   ) ::   kt     ! ocean time-step index
      REAL(wp), INTENT(in   ) ::   p2dt   ! time step  [s]
      INTEGER , INTENT(in   ) ::   Kmm, Krhs ! time level index
      !
      INTEGER :: ji, jj, jk, jl, ik, it, itmod
      LOGICAL :: lldebug = .FALSE.
      REAL(wp) :: zavt, zfn, zfn2
      !                                              ! z(ts)mltot : dT/dt over the anlysis window (including Asselin)
      !                                              ! z(ts)mlres : residual = dh/dt entrainment term
      REAL(wp), DIMENSION(T2D(0)  )   ::  ztmltot , zsmltot , ztmlres , zsmlres , ztmlatf , zsmlatf
      REAL(wp), DIMENSION(T2D(0)  )   ::  ztmltot2, zsmltot2, ztmlres2, zsmlres2, ztmlatf2, zsmlatf2, ztmltrdm2, zsmltrdm2
      REAL(wp), DIMENSION(T2D(0),jpk) ::  ztmltrd2, zsmltrd2   ! only needed for mean diagnostics
      !!----------------------------------------------------------------------

      ! ======================================================================
      ! I. Separate and remove the K_z trend from the iso-neutral diffusion term (if any)
      ! ======================================================================
      IF( ln_traldf_iso ) THEN
         !
         ! Diagnose the purely vertical (K_z) diffusion trend
         ! ... These terms can be estimated by flux computation at the lower boundary of the ML 
         !     (we compute (-1/h) * K_z * d_z( T ) and (-1/h) * K_z * d_z( S ))

         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            ik=nmln(ji,jj)
            tmltrd(ji,jj,jpmxl_zdf) = - avt(ji,jj,ik) * ( ts(ji,jj,ik-1,jp_tem,Krhs) - ts(ji,jj,ik,jp_tem,Krhs) )   &
               &        / e3w(ji,jj,ik,Kmm) * tmask(ji,jj,ik) / MAX( 1., hmxl(ji,jj) )
            smltrd(ji,jj,jpmxl_zdf) = - avs(ji,jj,ik) * ( ts(ji,jj,ik-1,jp_sal,Krhs) - ts(ji,jj,ik,jp_sal,Krhs) )   &
               &        / e3w(ji,jj,ik,Kmm) * tmask(ji,jj,ik) / MAX( 1., hmxl(ji,jj) )
         END_2D
         !
         ! ... Remove this K_z trend from the iso-neutral diffusion term (if any)
         tmltrd(:,:,jpmxl_ldf) = tmltrd(:,:,jpmxl_ldf) - tmltrd(:,:,jpmxl_zdf)
         smltrd(:,:,jpmxl_ldf) = smltrd(:,:,jpmxl_ldf) - smltrd(:,:,jpmxl_zdf)
      END IF  

      ! ======================================================================
      ! II. Cumulate the trends over the analysis window
      ! ======================================================================

      ztmltrd2(:,:,:) = 0.e0   ;    zsmltrd2(:,:,:) = 0.e0  ! <<< reset arrays to zero
      ztmltot2(:,:)   = 0.e0   ;    zsmltot2(:,:)   = 0.e0
      ztmlres2(:,:)   = 0.e0   ;    zsmlres2(:,:)   = 0.e0
      ztmlatf2(:,:)   = 0.e0   ;    zsmlatf2(:,:)   = 0.e0

      ! II.1 Set before values of vertically average T and S 
      ! ----------------------------------------------------
      IF( kt > nit000 ) THEN
         !   ... temperature ...                    ... salinity ...
         tmlb   (:,:) = tml   (:,:)             ;   smlb   (:,:) = sml   (:,:)
         tmlatfn(:,:) = tmltrd(:,:,jpmxl_atf)   ;   smlatfn(:,:) = smltrd(:,:,jpmxl_atf)
      END IF


      ! II.3 Initialize mixed-layer "before" arrays for the 1rst analysis window    
      ! ------------------------------------------------------------------------
      IF( kt == 2 ) THEN  !  i.e. ( .NOT. ln_rstart ).AND.( kt == nit000 + 1)
         !
         !   ... temperature ...                ... salinity ...
         tmlbb  (:,:) = tmlb   (:,:)   ;   smlbb  (:,:) = smlb   (:,:)
         tmlbn  (:,:) = tml    (:,:)   ;   smlbn  (:,:) = sml    (:,:)
         tmlatfb(:,:) = tmlatfn(:,:)   ;   smlatfb(:,:) = smlatfn(:,:)
         
         tmltrd_csum_ub (:,:,:) = 0.e0  ;   smltrd_csum_ub (:,:,:) = 0.e0
         tmltrd_atf_sumb(:,:)   = 0.e0  ;   smltrd_atf_sumb(:,:)   = 0.e0

         hmxlbn(:,:) = hmxl(:,:)

         IF( sn_cfctl%l_prtctl ) THEN
            WRITE(numout,*) '             we reach kt == nit000 + 1 = ', nit000+1
            CALL prt_ctl(tab2d_1=tmlbb   , clinfo1=' tmlbb   -   : ', mask1=tmask)
            CALL prt_ctl(tab2d_1=tmlbn   , clinfo1=' tmlbn   -   : ', mask1=tmask)
            CALL prt_ctl(tab2d_1=tmlatfb , clinfo1=' tmlatfb -   : ', mask1=tmask)
         END IF
         !
      END IF

      IF( ( ln_rstart ) .AND. ( kt == nit000 ) .AND. sn_cfctl%l_prtctl ) THEN
         IF( ln_trdmxl_instant ) THEN
            WRITE(numout,*) '             restart from kt == nit000 = ', nit000
            CALL prt_ctl(tab2d_1=tmlbb   , clinfo1=' tmlbb   -   : ', mask1=tmask)
            CALL prt_ctl(tab2d_1=tmlbn   , clinfo1=' tmlbn   -   : ', mask1=tmask)
            CALL prt_ctl(tab2d_1=tmlatfb , clinfo1=' tmlatfb -   : ', mask1=tmask)
         ELSE
            WRITE(numout,*) '             restart from kt == nit000 = ', nit000
            CALL prt_ctl(tab2d_1=tmlbn          , clinfo1=' tmlbn           -  : ', mask1=tmask)
            CALL prt_ctl(tab2d_1=hmxlbn         , clinfo1=' hmxlbn          -  : ', mask1=tmask)
            CALL prt_ctl(tab2d_1=tml_sumb       , clinfo1=' tml_sumb        -  : ', mask1=tmask)
            CALL prt_ctl(tab2d_1=tmltrd_atf_sumb, clinfo1=' tmltrd_atf_sumb -  : ', mask1=tmask)
            CALL prt_ctl(tab3d_1=tmltrd_csum_ub , clinfo1=' tmltrd_csum_ub  -  : ', mask1=tmask, kdim=1)
         END IF
      END IF

      ! II.4 Cumulated trends over the analysis period
      ! ----------------------------------------------
      !
      !         [  1rst analysis window ] [     2nd analysis window     ]                       
      !
      !     o---[--o-----o-----o-----o--]-[--o-----o-----o-----o-----o--]---o-----o--> time steps
      !                          nn_trd                           2*nn_trd       etc.
      !     1      2     3     4    =5 e.g.                          =10
      !
      IF( ( kt >= 2 ).OR.( ln_rstart ) ) THEN
         !
         nmoymltrd = nmoymltrd + 1
         
         ! ... Cumulate over BOTH physical contributions AND over time steps
         DO jl = 1, jpltrd
            tmltrdm(:,:) = tmltrdm(:,:) + tmltrd(:,:,jl)
            smltrdm(:,:) = smltrdm(:,:) + smltrd(:,:,jl)
         END DO

         ! ... Special handling of the Asselin trend 
         tmlatfm(:,:) = tmlatfm(:,:) + tmlatfn(:,:)
         smlatfm(:,:) = smlatfm(:,:) + smlatfn(:,:)

         ! ... Trends associated with the time mean of the ML T/S
         tmltrd_sum    (:,:,:) = tmltrd_sum    (:,:,:) + tmltrd    (:,:,:) ! tem
         tmltrd_csum_ln(:,:,:) = tmltrd_csum_ln(:,:,:) + tmltrd_sum(:,:,:)
         tml_sum       (:,:)   = tml_sum       (:,:)   + tml       (:,:)
         smltrd_sum    (:,:,:) = smltrd_sum    (:,:,:) + smltrd    (:,:,:) ! sal
         smltrd_csum_ln(:,:,:) = smltrd_csum_ln(:,:,:) + smltrd_sum(:,:,:)
         sml_sum       (:,:)   = sml_sum       (:,:)   + sml       (:,:)
         hmxl_sum      (:,:)   = hmxl_sum      (:,:)   + hmxl      (:,:)   ! rmxl
         !
      END IF

      ! ======================================================================
      ! III. Prepare fields for output (get here ONCE PER ANALYSIS PERIOD)
      ! ======================================================================

      ! Convert to appropriate physical units
      ! N.B. It may be useful to check IOIPSL time averaging with :
      !      tmltrd (:,:,:) = 1. ; smltrd (:,:,:) = 1.
      tmltrd(:,:,:) = tmltrd(:,:,:) * rn_ucf   ! (actually needed for 1:jpltrd-1, but trdmxl(:,:,jpltrd)
      smltrd(:,:,:) = smltrd(:,:,:) * rn_ucf   !  is no longer used, and is reset to 0. at next time step)
      
      ! define time axis
      it = kt
      itmod = kt - nit000 + 1

      MODULO_NTRD : IF( MOD( itmod, nn_trd ) == 0 ) THEN        ! nitend MUST be multiple of nn_trd
         !      
         zfn  = REAL( MAX(1,nmoymltrd), wp )   ;   zfn2 = zfn * zfn
         
         ! III.1 Prepare fields for output ("instantaneous" diagnostics) 
         ! -------------------------------------------------------------
         
         !-- Compute total trends
         ztmltot(:,:) = ( tml(:,:) - tmlbn(:,:) + tmlb(:,:) - tmlbb(:,:) ) / p2dt
         zsmltot(:,:) = ( sml(:,:) - smlbn(:,:) + smlb(:,:) - smlbb(:,:) ) / p2dt
         
         !-- Compute residuals
         ztmlres(:,:) = ztmltot(:,:) - ( tmltrdm(:,:) - tmlatfn(:,:) + tmlatfb(:,:) )
         zsmlres(:,:) = zsmltot(:,:) - ( smltrdm(:,:) - smlatfn(:,:) + smlatfb(:,:) )
      
         !-- Diagnose Asselin trend over the analysis window 
         ztmlatf(:,:) = tmlatfm(:,:) - tmlatfn(:,:) + tmlatfb(:,:)
         zsmlatf(:,:) = smlatfm(:,:) - smlatfn(:,:) + smlatfb(:,:)
         
         ! III.2 Prepare fields for output ("mean" diagnostics) 
         ! ----------------------------------------------------
         
         !-- Update the ML depth time sum (to build the Leap-Frog time mean)
         hmxl_sum(:,:) = hmxlbn(:,:) + 2._wp * ( hmxl_sum(:,:) - hmxl(:,:) ) + hmxl(:,:)

         !-- Compute temperature total trends
         tml_sum (:,:) = tmlbn(:,:) + 2._wp * ( tml_sum(:,:) - tml(:,:) ) + tml(:,:)
         ztmltot2(:,:) = ( tml_sum(:,:) - tml_sumb(:,:) ) / p2dt    ! now in degC/s
         
         !-- Compute salinity total trends
         sml_sum (:,:) = smlbn(:,:) + 2._wp * ( sml_sum(:,:) - sml(:,:) ) + sml(:,:)
         zsmltot2(:,:) = ( sml_sum(:,:) - sml_sumb(:,:) ) / p2dt    ! now in psu/s
         
         !-- Compute temperature residuals
         DO jl = 1, jpltrd
            ztmltrd2(:,:,jl) = tmltrd_csum_ub(:,:,jl) + tmltrd_csum_ln(:,:,jl)
         END DO

         ztmltrdm2(:,:) = 0._wp
         DO jl = 1, jpltrd
            ztmltrdm2(:,:) = ztmltrdm2(:,:) + ztmltrd2(:,:,jl)
         END DO

         ztmlres2(:,:) =  ztmltot2(:,:)  -       &
              ( ztmltrdm2(:,:) - tmltrd_sum(:,:,jpmxl_atf) + tmltrd_atf_sumb(:,:) )
         
         !-- Compute salinity residuals
         DO jl = 1, jpltrd
            zsmltrd2(:,:,jl) = smltrd_csum_ub(:,:,jl) + smltrd_csum_ln(:,:,jl)
         END DO

         zsmltrdm2(:,:) = 0._wp
         DO jl = 1, jpltrd
            zsmltrdm2(:,:) = zsmltrdm2(:,:) + zsmltrd2(:,:,jl)
         END DO

         zsmlres2(:,:) =  zsmltot2(:,:)  -       &
              ( zsmltrdm2(:,:) - smltrd_sum(:,:,jpmxl_atf) + smltrd_atf_sumb(:,:) )
         
         !-- Diagnose Asselin trend over the analysis window
         ztmlatf2(:,:) = ztmltrd2(:,:,jpmxl_atf) - tmltrd_sum(:,:,jpmxl_atf) + tmltrd_atf_sumb(:,:)
         zsmlatf2(:,:) = zsmltrd2(:,:,jpmxl_atf) - smltrd_sum(:,:,jpmxl_atf) + smltrd_atf_sumb(:,:)
         
         ! III.3 Time evolution array swap
         ! -------------------------------
         
         ! For T/S instantaneous diagnostics 
         !   ... temperature ...               ... salinity ...
         tmlbb  (:,:) = tmlb   (:,:)  ;   smlbb  (:,:) = smlb   (:,:)
         tmlbn  (:,:) = tml    (:,:)  ;   smlbn  (:,:) = sml    (:,:)
         tmlatfb(:,:) = tmlatfn(:,:)  ;   smlatfb(:,:) = smlatfn(:,:)

         ! For T mean diagnostics 
         tmltrd_csum_ub (:,:,:) = zfn * tmltrd_sum(:,:,:) - tmltrd_csum_ln(:,:,:)
         tml_sumb       (:,:)   = tml_sum(:,:)
         tmltrd_atf_sumb(:,:)   = tmltrd_sum(:,:,jpmxl_atf)
         
         ! For S mean diagnostics 
         smltrd_csum_ub (:,:,:) = zfn * smltrd_sum(:,:,:) - smltrd_csum_ln(:,:,:)
         sml_sumb       (:,:)   = sml_sum(:,:)
         smltrd_atf_sumb(:,:)   = smltrd_sum(:,:,jpmxl_atf)
         
         ! ML depth
         hmxlbn         (:,:)   = hmxl    (:,:)
         
         IF( sn_cfctl%l_prtctl ) THEN
            IF( ln_trdmxl_instant ) THEN
               CALL prt_ctl(tab2d_1=tmlbb   , clinfo1=' tmlbb   -   : ', mask1=tmask)
               CALL prt_ctl(tab2d_1=tmlbn   , clinfo1=' tmlbn   -   : ', mask1=tmask)
               CALL prt_ctl(tab2d_1=tmlatfb , clinfo1=' tmlatfb -   : ', mask1=tmask)
            ELSE
               CALL prt_ctl(tab2d_1=tmlbn          , clinfo1=' tmlbn           -  : ', mask1=tmask)
               CALL prt_ctl(tab2d_1=hmxlbn         , clinfo1=' hmxlbn          -  : ', mask1=tmask)
               CALL prt_ctl(tab2d_1=tml_sumb       , clinfo1=' tml_sumb        -  : ', mask1=tmask)
               CALL prt_ctl(tab2d_1=tmltrd_atf_sumb, clinfo1=' tmltrd_atf_sumb -  : ', mask1=tmask)
               CALL prt_ctl(tab3d_1=tmltrd_csum_ub , clinfo1=' tmltrd_csum_ub  -  : ', mask1=tmask, kdim=1)
            END IF
         END IF

         ! III.4 Convert to appropriate physical units
         ! -------------------------------------------

         !    ... temperature ...                         ... salinity ...
         ztmltot (:,:)   = ztmltot(:,:)   * rn_ucf/zfn  ; zsmltot (:,:)   = zsmltot(:,:)   * rn_ucf/zfn
         ztmlres (:,:)   = ztmlres(:,:)   * rn_ucf/zfn  ; zsmlres (:,:)   = zsmlres(:,:)   * rn_ucf/zfn
         ztmlatf (:,:)   = ztmlatf(:,:)   * rn_ucf/zfn  ; zsmlatf (:,:)   = zsmlatf(:,:)   * rn_ucf/zfn

         tml_sum (:,:)   = tml_sum (:,:)  /  (2*zfn) ; sml_sum (:,:)   = sml_sum (:,:)  /  (2*zfn)
         ztmltot2(:,:)   = ztmltot2(:,:)  * rn_ucf/zfn2 ; zsmltot2(:,:)   = zsmltot2(:,:)  * rn_ucf/zfn2
         ztmltrd2(:,:,:) = ztmltrd2(:,:,:)* rn_ucf/zfn2 ; zsmltrd2(:,:,:) = zsmltrd2(:,:,:)* rn_ucf/zfn2
         ztmlatf2(:,:)   = ztmlatf2(:,:)  * rn_ucf/zfn2 ; zsmlatf2(:,:)   = zsmlatf2(:,:)  * rn_ucf/zfn2
         ztmlres2(:,:)   = ztmlres2(:,:)  * rn_ucf/zfn2 ; zsmlres2(:,:)   = zsmlres2(:,:)  * rn_ucf/zfn2

         hmxl_sum(:,:)   = hmxl_sum(:,:)  /  (2*zfn)  ! similar to tml_sum and sml_sum

         ! * Debugging information *
         IF( lldebug .AND. lwp ) THEN
            !
            WRITE(numout,*)
            WRITE(numout,*) 'trd_mxl : write trends in the Mixed Layer for debugging process:'
            WRITE(numout,*) '~~~~~~~  '
            WRITE(numout,*) '          TRA kt = ', kt, 'nmoymltrd = ', nmoymltrd
            WRITE(numout,*)
            WRITE(numout,*) '          >>>>>>>>>>>>>>>>>>  TRA TEMPERATURE <<<<<<<<<<<<<<<<<<'
            WRITE(numout,*) '          TRA ztmlres    : ', SUM(ztmlres(:,:))
            WRITE(numout,*) '          TRA ztmltot    : ', SUM(ztmltot(:,:))
            WRITE(numout,*) '          TRA tmltrdm    : ', SUM(tmltrdm(:,:))
            WRITE(numout,*) '          TRA tmlatfb    : ', SUM(tmlatfb(:,:))
            WRITE(numout,*) '          TRA tmlatfn    : ', SUM(tmlatfn(:,:))
            DO jl = 1, jpltrd
               WRITE(numout,*) '          * TRA TREND INDEX jpmxl_xxx = jl = ', jl, &
                    & ' tmltrd : ', SUM(tmltrd(:,:,jl))
            END DO
            WRITE(numout,*) '          TRA ztmlres (jpi/2,jpj/2) : ', ztmlres (jpi/2,jpj/2)
            WRITE(numout,*) '          TRA ztmlres2(jpi/2,jpj/2) : ', ztmlres2(jpi/2,jpj/2)
            WRITE(numout,*)
            WRITE(numout,*) '          >>>>>>>>>>>>>>>>>>  TRA SALINITY <<<<<<<<<<<<<<<<<<'
            WRITE(numout,*) '          TRA zsmlres    : ', SUM(zsmlres(:,:))
            WRITE(numout,*) '          TRA zsmltot    : ', SUM(zsmltot(:,:))
            WRITE(numout,*) '          TRA smltrdm    : ', SUM(smltrdm(:,:))
            WRITE(numout,*) '          TRA smlatfb    : ', SUM(smlatfb(:,:))
            WRITE(numout,*) '          TRA smlatfn    : ', SUM(smlatfn(:,:))
            DO jl = 1, jpltrd
               WRITE(numout,*) '          * TRA TREND INDEX jpmxl_xxx = jl = ', jl, &
                    & ' smltrd : ', SUM(smltrd(:,:,jl))
            END DO
            WRITE(numout,*) '          TRA zsmlres (jpi/2,jpj/2) : ', zsmlres (jpi/2,jpj/2)
            WRITE(numout,*) '          TRA zsmlres2(jpi/2,jpj/2) : ', zsmlres2(jpi/2,jpj/2)
            !
         END IF
         !
      END IF MODULO_NTRD

      ! ======================================================================
      ! IV. Write trends in the NetCDF file
      ! ======================================================================

      !-- Write the trends for T/S instantaneous diagnostics 

      ! clem => these fields do not exist in field_def
      
      IF( ln_trdmxl_instant ) THEN           

         CALL iom_put( "mxl_depth", hmxl(:,:) )
         
         !................................. ( ML temperature ) ...................................
         
         !-- Output the fields
         CALL iom_put( "tml"     , tml    (:,:) ) 
         CALL iom_put( "tml_tot" , ztmltot(:,:) ) 
         CALL iom_put( "tml_res" , ztmlres(:,:) ) 
         
         DO jl = 1, jpltrd - 1
            CALL iom_put( trim("tml"//ctrd(jl,2)), tmltrd (:,:,jl) )
         END DO
         
         CALL iom_put( trim("tml"//ctrd(jpmxl_atf,2)), ztmlatf(:,:) )
         
         !.................................. ( ML salinity ) .....................................
         
         !-- Output the fields
         CALL iom_put( "sml"    , sml    (:,:) )
         CALL iom_put( "sml_tot", zsmltot(:,:) ) 
         CALL iom_put( "sml_res", zsmlres(:,:) ) 
         
         DO jl = 1, jpltrd - 1
            CALL iom_put( trim("sml"//ctrd(jl,2)), smltrd(:,:,jl) )
         END DO
         
         CALL iom_put( trim("sml"//ctrd(jpmxl_atf,2)), zsmlatf(:,:) )


         
      ELSE      !-- Write the trends for T/S mean diagnostics 

         CALL iom_put( "mxl_depth", hmxl_sum(:,:) )
         
         !................................. ( ML temperature ) ...................................
         
         !-- Output the fields
         CALL iom_put( "tml"     , tml_sum (:,:) ) 
         CALL iom_put( "tml_tot" , ztmltot2(:,:) ) 
         CALL iom_put( "tml_res" , ztmlres2(:,:) ) 

         DO jl = 1, jpltrd - 1
            CALL iom_put( trim("tml"//ctrd(jl,2)), ztmltrd2(:,:,jl) )
         END DO
         
         CALL iom_put( trim("tml"//ctrd(jpmxl_atf,2)), ztmlatf2(:,:) )
         
         !.................................. ( ML salinity ) .....................................
                     
         !-- Output the fields
         CALL iom_put( "sml"    , sml_sum (:,:) )
         CALL iom_put( "sml_tot", zsmltot2(:,:) ) 
         CALL iom_put( "sml_res", zsmlres2(:,:) ) 

         DO jl = 1, jpltrd - 1
            CALL iom_put( trim("sml"//ctrd(jl,2)), zsmltrd2(:,:,jl) )
         END DO
         
         CALL iom_put( trim("sml"//ctrd(jpmxl_atf,2)), zsmlatf2(:,:) )
         !
      END IF
      !

      IF( MOD( itmod, nn_trd ) == 0 ) THEN
         !
         ! III.5 Reset cumulative arrays to zero
         ! -------------------------------------
         nmoymltrd = 0
         
         !   ... temperature ...               ... salinity ...
         tmltrdm        (:,:)   = 0.e0   ;   smltrdm        (:,:)   = 0.e0
         tmlatfm        (:,:)   = 0.e0   ;   smlatfm        (:,:)   = 0.e0
         tml_sum        (:,:)   = 0.e0   ;   sml_sum        (:,:)   = 0.e0
         tmltrd_csum_ln (:,:,:) = 0.e0   ;   smltrd_csum_ln (:,:,:) = 0.e0
         tmltrd_sum     (:,:,:) = 0.e0   ;   smltrd_sum     (:,:,:) = 0.e0

         hmxl_sum       (:,:)   = 0.e0
         !
      END IF

      ! ======================================================================
      ! V. Write restart file
      ! ======================================================================

      IF( lrst_oce )   CALL trd_mxl_rst_write( kt ) 

      !
   END SUBROUTINE trd_mxl


   SUBROUTINE trd_mxl_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_mxl_init  ***
      !! 
      !! ** Purpose :   computation of vertically integrated T and S budgets
      !!      from ocean surface down to control surface (NetCDF output)
      !!----------------------------------------------------------------------
      INTEGER  ::   jl     ! dummy loop indices
      INTEGER  ::   inum   ! logical unit
      INTEGER  ::   ios    ! local integer
      REAL(wp) ::   zjulian, zsto, zout
      CHARACTER (LEN=40) ::   clop
      CHARACTER (LEN=12) ::   clmxl, cltu, clsu
      !!
      NAMELIST/namtrd_mxl/ nn_trd , cn_trdrst_in , ln_trdmxl_restart,       &
         &                 nn_ctls, cn_trdrst_out, ln_trdmxl_instant, rn_ucf, rn_rho_c
      !!----------------------------------------------------------------------
      !
      READ_NML_REF(numnam,namtrd_mxl)
      READ_NML_CFG(numnam,namtrd_mxl)
      IF(lwm) WRITE( numond, namtrd_mxl )
      !
      IF(lwp) THEN                      ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' trd_mxl_init : Mixed-layer trends'
         WRITE(numout,*) ' ~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtrd : set trends parameters'
         WRITE(numout,*) '      frequency of trends diagnostics (glo)      nn_trd             = ', nn_trd
         WRITE(numout,*) '      density criteria used to defined the MLD   rn_rho_c           = ', rn_rho_c
         WRITE(numout,*) '      control surface type            (mld)      nn_ctls            = ', nn_ctls
         WRITE(numout,*) '      restart for ML diagnostics                 ln_trdmxl_restart  = ', ln_trdmxl_restart
         WRITE(numout,*) '      instantaneous or mean ML T/S               ln_trdmxl_instant  = ', ln_trdmxl_instant
         WRITE(numout,*) '      unit conversion factor                     rn_ucf             = ', rn_ucf
         WRITE(numout,*) '      criteria to compute the MLD                rn_rho_c           = ', rn_rho_c
      ENDIF



      ! I.1 Check consistency of user defined preferences
      ! -------------------------------------------------
      
      IF ( rn_rho_c /= rho_c )   CALL ctl_warn( 'Unless you have good reason to do so, you should use the value ',    &
         &                                      'defined in zdfmxl.F90 module to calculate the mixed layer depth' )

      IF( MOD( nitend, nn_trd ) /= 0 ) THEN
         WRITE(ctmp1,*) '                Your nitend parameter, nitend = ', nitend
         WRITE(ctmp2,*) '                is no multiple of the trends diagnostics frequency        '
         WRITE(ctmp3,*) '                          you defined, nn_trd   = ', nn_trd
         WRITE(ctmp4,*) '                This will not allow you to restart from this simulation.  '
         WRITE(ctmp5,*) '                You should reconsider this choice.                        ' 
         WRITE(ctmp6,*) 
         WRITE(ctmp7,*) '                N.B. the nitend parameter is also constrained to be a     '
         WRITE(ctmp8,*) '                     multiple of the nn_fsbc parameter '
         CALL ctl_stop( ctmp1, ctmp2, ctmp3, ctmp4, ctmp5, ctmp6, ctmp7, ctmp8 )
      END IF

      !                                   ! allocate trdmxl arrays
      IF( trdmxl_oce_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trd_mxl_init : unable to allocate trdmxl_oce arrays' )



      nkstp     = nit000 - 1              ! current time step indicator initialization




      ! I.2 Initialize arrays to zero or read a restart file
      ! ----------------------------------------------------

      nmoymltrd = 0

      !     ... temperature ...                  ... salinity ...
      tml            (:,:)   = 0.e0    ;    sml            (:,:)   = 0.e0     ! inst.
      tmltrdm        (:,:)   = 0.e0    ;    smltrdm        (:,:)   = 0.e0
      tmlatfm        (:,:)   = 0.e0    ;    smlatfm        (:,:)   = 0.e0
      tml_sum        (:,:)   = 0.e0    ;    sml_sum        (:,:)   = 0.e0     ! mean
      tmltrd_sum     (:,:,:) = 0.e0    ;    smltrd_sum     (:,:,:) = 0.e0
      tmltrd_csum_ln (:,:,:) = 0.e0    ;    smltrd_csum_ln (:,:,:) = 0.e0

      hmxl           (:,:)   = 0.e0            
      hmxl_sum       (:,:)   = 0.e0

      IF( ln_rstart .AND. ln_trdmxl_restart ) THEN
         CALL trd_mxl_rst_read
      ELSE
         !     ... temperature ...                  ... salinity ...
         tmlb           (:,:)   = 0.e0    ;    smlb           (:,:)   = 0.e0  ! inst.
         tmlbb          (:,:)   = 0.e0    ;    smlbb          (:,:)   = 0.e0  
         tmlbn          (:,:)   = 0.e0    ;    smlbn          (:,:)   = 0.e0 
         tmlatfb        (:,:)   = 0.e0    ;    smlatfb        (:,:)   = 0.e0
         tml_sumb       (:,:)   = 0.e0    ;    sml_sumb       (:,:)   = 0.e0  ! mean
         tmltrd_csum_ub (:,:,:) = 0.e0    ;    smltrd_csum_ub (:,:,:) = 0.e0
         tmltrd_atf_sumb(:,:)   = 0.e0    ;    smltrd_atf_sumb(:,:)   = 0.e0  
      END IF

      icount = 1   ;   ionce  = 1                            ! open specifier

      ! I.3 Read control surface from file ctlsurf_idx
      ! ----------------------------------------------
 
      IF( nn_ctls == 1 ) THEN
         CALL ctl_opn( inum, 'ctlsurf_idx', 'OLD', 'UNFORMATTED', 'SEQUENTIAL', -1, numout, lwp )
         READ ( inum, * ) nbol
         CLOSE( inum )
      END IF

      ! ======================================================================
      ! II. netCDF output initialization
      ! ======================================================================

      ! clmxl = legend root for netCDF output
      IF( nn_ctls == 0 ) THEN      ! control surface = mixed-layer with density criterion
         clmxl = 'Mixed Layer '  !                   (array nmln computed in zdfmxl.F90)
      ELSE IF( nn_ctls == 1 ) THEN ! control surface = read index from file 
         clmxl = '      Bowl '
      ELSE IF( nn_ctls >= 2 ) THEN ! control surface = model level
         WRITE(clmxl,'(A10,I2,1X)') 'Levels 1 -', nn_ctls
      END IF



      ! II.3 Define the T grid trend file (nidtrd)
      ! ------------------------------------------
      !-- Define long and short names for the NetCDF output variables
      !       ==> choose them according to trdmxl_oce.F90 <==

      ctrd(jpmxl_xad,1) = " Zonal advection"                  ;   ctrd(jpmxl_xad,2) = "_xad"
      ctrd(jpmxl_yad,1) = " Meridional advection"             ;   ctrd(jpmxl_yad,2) = "_yad"
      ctrd(jpmxl_zad,1) = " Vertical advection"               ;   ctrd(jpmxl_zad,2) = "_zad"
      ctrd(jpmxl_ldf,1) = " Lateral diffusion"                ;   ctrd(jpmxl_ldf,2) = "_ldf"
      ctrd(jpmxl_for,1) = " Forcing"                          ;   ctrd(jpmxl_for,2) = "_for"
      ctrd(jpmxl_zdf,1) = " Vertical diff. (Kz)"              ;   ctrd(jpmxl_zdf,2) = "_zdf"
      ctrd(jpmxl_bbc,1) = " Geothermal flux"                  ;   ctrd(jpmxl_bbc,2) = "_bbc"
      ctrd(jpmxl_bbl,1) = " Adv/diff. Bottom boundary layer"  ;   ctrd(jpmxl_bbl,2) = "_bbl"
      ctrd(jpmxl_dmp,1) = " Tracer damping"                   ;   ctrd(jpmxl_dmp,2) = "_dmp"
      ctrd(jpmxl_npc,1) = " Non penetrative convec. adjust."  ;   ctrd(jpmxl_npc,2) = "_npc"
      ctrd(jpmxl_zdfp,1) = " iso-neut.: pure vert. diffusion" ;  ctrd(jpmxl_zdfp,2) = "_zdfp"
      ctrd(jpmxl_atf,1) = " Asselin time filter"              ;   ctrd(jpmxl_atf,2) = "_atf"
                                                                  

      !-- Define physical units
      IF     ( rn_ucf == 1.       ) THEN   ;   cltu = "degC/s"     ;   clsu = "p.s.u./s"
      ELSEIF ( rn_ucf == 3600.*24.) THEN   ;   cltu = "degC/day"   ;   clsu = "p.s.u./day"
      ELSE                                 ;   cltu = "unknown?"   ;   clsu = "unknown?"
      END IF
      !
   END SUBROUTINE trd_mxl_init

   !!======================================================================
END MODULE trdmxl

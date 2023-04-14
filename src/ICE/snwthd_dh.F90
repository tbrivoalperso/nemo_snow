MODULE snwthd_dh
   !!======================================================================
   !!                       ***  MODULE snwthd_dh ***
   !!   seaice : snow thermodynamic growth and melt
   !!======================================================================
   !! History :       !  2003-05  (M. Vancoppenolle) Original code in 1D
   !!                 !  2005-06  (M. Vancoppenolle) 3D version
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   snw_thd_dh        : vertical sea-ice growth and melt
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE ice            ! sea-ice: variables
   USE ice1D          ! sea-ice: thermodynamics variables
   USE icevar         ! for CALL ice_var_snwblow
   USE snwent         ! snow enthalpy remapping
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   snw_thd_dh        ! called by snw_thd

   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: snwthd_dh.F90 Theo $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE snw_thd_dh( zq_rema, zevap_rema , zh_s, ze_s)
      !!------------------------------------------------------------------
      !!                ***  ROUTINE snw_thd_dh  ***
      !!
      !! ** Purpose :   Snow thickness changes due to growth/melting
      !!
      !! ** Method  :   Snow surface melting arises from imbalance in surface fluxes
      !!                Snow thickness can increase by precipitation and decrease by sublimation
      !!
      !!                - Compute available flux of heat for surface ablation
      !!                - Compute snow and sea ice enthalpies
      !!                - Surface ablation and sublimation
      !!                - Returns the remaining heat and mass fluxes after snow melting 
      !!                  / sublimation, and the snow enthalpy and thicknesses profiles
      !!                - Enthalpy is NOT REMAPPED in this routine, but later in icethd_dh
      !!                - ice / snow conversion is not computed here, but later in icethd_dh too
      !! 
      !! ** Notes     : - h=max(0,h+dh) are often used to ensure positivity of h.
      !!                very small negative values can occur otherwise (e.g. -1.e-20)
      !!                - The routine is simply an extraction of snow in the previous 
      !!                  icethd_dh.F90 (from v4.2-stable) routine
      !! References : Bitz and Lipscomb, 1999, J. Geophys. Res.
      !!              Fichefet T. and M. Maqueda 1997, J. Geophys. Res., 102(C6), 12609-12646
      !!              Vancoppenolle, Fichefet and Bitz, 2005, Geophys. Res. Let.
      !!              Vancoppenolle et al.,2009, Ocean Modelling
      !!------------------------------------------------------------------
      REAL(wp), DIMENSION(jpij), INTENT(out) ::   zq_rema     ! remaining heat flux from snow melting       (J.m-2)
      REAL(wp), DIMENSION(jpij), INTENT(out) ::   zevap_rema  ! remaining mass flux from snow sublimation   (kg.m-2)
      REAL(wp), DIMENSION(jpij,0:nlay_s  ), INTENT(out) ::   zh_s      ! snw layer thickness (m) 
      REAL(wp), DIMENSION(jpij,0:nlay_s  ), INTENT(out) ::   ze_s      ! snw layer enthalpy (J.m-3)

!
      INTEGER  ::   ji, jk       ! dummy loop indices
      INTEGER  ::   iter         ! local integer

      REAL(wp) ::   zdum

      REAL(wp), DIMENSION(jpij) ::   zq_top      ! heat for surface ablation                   (J.m-2)
      REAL(wp), DIMENSION(jpij) ::   zdeltah
      REAL(wp), DIMENSION(jpij) ::   zsnw        ! distribution of snow after wind blowing

      !!------------------------------------------------------------------
      ! Initialise remaining heat and mass fluxes after melt and sublimation
      zq_rema(1:npti)    = 0._wp
      zevap_rema(1:npti) = 0._wp

      !
      ! initialize snw layer thicknesses and enthalpies
      zh_s(1:npti,0) = 0._wp
      ze_s(1:npti,0) = 0._wp
      DO jk = 1, nlay_s
         DO ji = 1, npti
            zh_s(ji,jk) = h_s_1d(ji) * r1_nlay_s
            ze_s(ji,jk) = e_s_1d(ji,jk)
         END DO
      END DO
      !
      !                       ! ============================================== !
      !                       ! Available heat for surface ablation !
      !                       ! ============================================== !
      !
      IF( ln_cndflx .AND. .NOT.ln_cndemulate ) THEN
         !
         DO ji = 1, npti
            zq_top(ji)     = MAX( 0._wp, qml_ice_1d(ji) * rDt_ice )
         END DO
         !
      ELSE
         !
         DO ji = 1, npti
            zdum           = qns_ice_1d(ji) + qsr_ice_1d(ji) - qtr_ice_top_1d(ji) - qcn_ice_top_1d(ji)
            qml_ice_1d(ji) = zdum * MAX( 0._wp , SIGN( 1._wp, t_su_1d(ji) - rt0 ) )
            zq_top(ji)     = MAX( 0._wp, qml_ice_1d(ji) * rDt_ice )
         END DO
         !
      ENDIF
      !
      !                       ! ============ !
      !                       !     Snow     !
      !                       ! ============ !
      !
      ! Internal melting
      ! ----------------
      ! IF snow temperature is above freezing point, THEN snow melts (should not happen but sometimes it does)
      DO jk = 1, nlay_s
         DO ji = 1, npti
            IF( t_s_1d(ji,jk) > rt0 ) THEN
               hfx_res_1d    (ji) = hfx_res_1d    (ji) - ze_s(ji,jk) * zh_s(ji,jk) * a_i_1d(ji) * r1_Dt_ice   ! heat flux to the ocean [W.m-2], < 0
               wfx_snw_sum_1d(ji) = wfx_snw_sum_1d(ji) + rhos        * zh_s(ji,jk) * a_i_1d(ji) * r1_Dt_ice   ! mass flux
               ! updates
               dh_s_mlt(ji)    =             dh_s_mlt(ji) - zh_s(ji,jk)
               h_s_1d  (ji)    = MAX( 0._wp, h_s_1d  (ji) - zh_s(ji,jk) )
               zh_s    (ji,jk) = 0._wp
               ze_s    (ji,jk) = 0._wp
            END IF
         END DO
      END DO

     IF( .NOT. ln_snwext ) THEN ! Nb: this part of the code is the same as in the snwthd_snwfl routine. We keep it here to be consistent with 4.2.stable version
         ! Snow precipitation
         !-------------------
         CALL ice_var_snwblow( 1._wp - at_i_1d(1:npti), zsnw(1:npti) )   ! snow distribution over ice after wind blowing

         DO ji = 1, npti
            IF( sprecip_1d(ji) > 0._wp ) THEN
               zh_s(ji,0) = zsnw(ji) * sprecip_1d(ji) * rDt_ice * r1_rhos / at_i_1d(ji)   ! thickness of precip
               ze_s(ji,0) = MAX( 0._wp, - qprec_ice_1d(ji) )                              ! enthalpy of the precip (>0, J.m-3)
               !
               hfx_spr_1d(ji) = hfx_spr_1d(ji) + ze_s(ji,0) * zh_s(ji,0) * a_i_1d(ji) * r1_Dt_ice   ! heat flux from snow precip (>0, W.m-2)
               wfx_spr_1d(ji) = wfx_spr_1d(ji) - rhos       * zh_s(ji,0) * a_i_1d(ji) * r1_Dt_ice   ! mass flux, <0
               !
               ! update thickness
               h_s_1d(ji) = h_s_1d(ji) + zh_s(ji,0)
            ENDIF
         END DO
      ENDIF
      ! Snow melting
      ! ------------
      ! If heat still available (zq_top > 0)
      ! then all snw precip has been melted and we need to melt more snow
      DO jk = 0, nlay_s
         DO ji = 1, npti
            IF( zh_s(ji,jk) > 0._wp .AND. zq_top(ji) > 0._wp ) THEN
               !
               rswitch = MAX( 0._wp , SIGN( 1._wp , ze_s(ji,jk) - epsi20 ) )
               zdum    = - rswitch * zq_top(ji) / MAX( ze_s(ji,jk), epsi20 )   ! thickness change
               zdum    = MAX( zdum , - zh_s(ji,jk) )                           ! bound melting

               hfx_snw_1d    (ji) = hfx_snw_1d    (ji) - ze_s(ji,jk) * zdum * a_i_1d(ji) * r1_Dt_ice   ! heat used to melt snow(W.m-2, >0)
               wfx_snw_sum_1d(ji) = wfx_snw_sum_1d(ji) - rhos        * zdum * a_i_1d(ji) * r1_Dt_ice   ! snow melting only = water into the ocean

               ! updates available heat + thickness
               dh_s_mlt(ji)    =              dh_s_mlt(ji)    + zdum
               zq_top  (ji)    = MAX( 0._wp , zq_top  (ji)    + zdum * ze_s(ji,jk) )
               h_s_1d  (ji)    = MAX( 0._wp , h_s_1d  (ji)    + zdum )
               zh_s    (ji,jk) = MAX( 0._wp , zh_s    (ji,jk) + zdum )
!!$               IF( zh_s(ji,jk) == 0._wp )   ze_s(ji,jk) = 0._wp

               !
            ENDIF
         END DO
      END DO
      DO ji = 1, npti
               zq_rema (ji) = zq_top (ji) ! remaining heat at the end of the routine in J.m-2 (used to melt ice later on)
      END DO

      ! Snow sublimation
      !-----------------
      ! qla_ice is always >=0 (upwards), heat goes to the atmosphere, therefore snow sublimates
      !    comment: not counted in mass/heat exchange in iceupdate.F90 since this is an exchange with atm. (not ocean)
      zdeltah   (1:npti) = 0._wp ! total snow thickness that sublimates, < 0
      zevap_rema(1:npti) = 0._wp
      DO ji = 1, npti
         zdeltah   (ji) = MAX( - evap_ice_1d(ji) * r1_rhos * rDt_ice, - h_s_1d(ji) )   ! amount of snw that sublimates, < 0
         zevap_rema(ji) = evap_ice_1d(ji) * rDt_ice + zdeltah(ji) * rhos               ! remaining evap in kg.m-2 (used for ice sublimation later on)
      END DO

      DO jk = 0, nlay_s
         DO ji = 1, npti
            zdum = MAX( -zh_s(ji,jk), zdeltah(ji) ) ! snow layer thickness that sublimates, < 0
            !
            hfx_sub_1d    (ji) = hfx_sub_1d    (ji) + ze_s(ji,jk) * zdum * a_i_1d(ji) * r1_Dt_ice  ! Heat flux of snw that sublimates [W.m-2], < 0
            wfx_snw_sub_1d(ji) = wfx_snw_sub_1d(ji) - rhos        * zdum * a_i_1d(ji) * r1_Dt_ice  ! Mass flux by sublimation

            ! update thickness
            h_s_1d(ji)    = MAX( 0._wp , h_s_1d(ji)    + zdum )
            zh_s  (ji,jk) = MAX( 0._wp , zh_s  (ji,jk) + zdum )
!!$            IF( zh_s(ji,jk) == 0._wp )   ze_s(ji,jk) = 0._wp

            ! update sublimation left
            zdeltah(ji) = MIN( zdeltah(ji) - zdum, 0._wp )
         END DO
      END DO

     IF( ln_snwext ) THEN
   
        ! Remapping of snw enthalpy on a regular grid
        !--------------------------------------------
         CALL snw_ent( zh_s, ze_s, e_s_1d)
   
         ! recalculate t_s_1d from e_s_1d
         DO jk = 1, nlay_s
            DO ji = 1,npti
               IF( h_s_1d(ji) > 0._wp ) THEN
                  t_s_1d(ji,jk) = rt0 + ( - e_s_1d(ji,jk) * r1_rhos * r1_rcpi + rLfus * r1_rcpi )
               ELSE
                  t_s_1d(ji,jk) = rt0
               ENDIF
            END DO
         END DO
   
      ENDIF

      !
      !
   END SUBROUTINE snw_thd_dh

#else
   !!----------------------------------------------------------------------
   !!   Default option                                NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE snwthd_dh

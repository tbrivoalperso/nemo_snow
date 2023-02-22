MODULE snwthd_dh
   !!======================================================================
   !!                       ***  MODULE snwthd_dh ***
   !!   seaice : thermodynamic growth and melt
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
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   snw_thd_dh        ! called by snw_thd

   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: snwthd_dh.F90 14686 2021-04-08 15:36:01Z clem $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE snw_thd_dh( zq_rema, zevap_rema )
      !!------------------------------------------------------------------
      !!                ***  ROUTINE snw_thd_dh  ***
      !!
      !! ** Purpose :   compute ice and snow thickness changes due to growth/melting
      !!
      !! ** Method  :   Ice/Snow surface melting arises from imbalance in surface fluxes
      !!                Bottom accretion/ablation arises from flux budget
      !!                Snow thickness can increase by precipitation and decrease by sublimation
      !!                If snow load excesses Archmiede limit, snow-ice is formed by
      !!                the flooding of sea-water in the snow
      !!
      !!                - Compute available flux of heat for surface ablation
      !!                - Compute snow and sea ice enthalpies
      !!                - Surface ablation and sublimation
      !!                - Bottom accretion/ablation
      !!                - Snow ice formation
      !!
      !! ** Note     :  h=max(0,h+dh) are often used to ensure positivity of h.
      !!                very small negative values can occur otherwise (e.g. -1.e-20)
      !!
      !! References : Bitz and Lipscomb, 1999, J. Geophys. Res.
      !!              Fichefet T. and M. Maqueda 1997, J. Geophys. Res., 102(C6), 12609-12646
      !!              Vancoppenolle, Fichefet and Bitz, 2005, Geophys. Res. Let.
      !!              Vancoppenolle et al.,2009, Ocean Modelling
      !!------------------------------------------------------------------
      REAL(wp), DIMENSION(jpij), INTENT(out) ::   zq_rema     ! remaining heat flux from snow melting       (J.m-2)
      REAL(wp), DIMENSION(jpij), INTENT(out) ::   zevap_rema  ! remaining mass flux from snow sublimation   (kg.m-2)
      !
      INTEGER  ::   ji, jk       ! dummy loop indices
      INTEGER  ::   iter         ! local integer

      REAL(wp) ::   zdum

      REAL(wp), DIMENSION(jpij) ::   zq_top      ! heat for surface ablation                   (J.m-2)
      REAL(wp), DIMENSION(jpij) ::   zdeltah
      REAL(wp), DIMENSION(jpij) ::   zsnw        ! distribution of snow after wind blowing

      REAL(wp), DIMENSION(jpij,0:nlay_s  ) ::   zh_s      ! snw layer thickness (m)
      REAL(wp), DIMENSION(jpij,0:nlay_s  ) ::   ze_s      ! snw layer enthalpy (J.m-3)

      !!------------------------------------------------------------------

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

               zq_rema (ji) = zq_top (ji) ! remaining heat at the end of the routine in J.m-2 (used to melt ice later on)
               !
            ENDIF
         END DO
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

      !
      !
      ! Remapping of snw enthalpy on a regular grid
      !--------------------------------------------
      CALL snw_ent( zh_s, ze_s, e_s_1d )

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

   END SUBROUTINE snw_thd_dh

   SUBROUTINE snw_ent( ph_old, pe_old, pe_new )
      !!-------------------------------------------------------------------
      !!               ***   ROUTINE snw_ent  ***
      !!
      !! ** Purpose :
      !!           This routine computes new vertical grids in the snow,
      !!           and consistently redistributes temperatures.
      !!           Redistribution is made so as to ensure to energy conservation
      !!
      !!
      !! ** Method  : linear conservative remapping
      !!
      !! ** Steps : 1) cumulative integrals of old enthalpies/thicknesses
      !!            2) linear remapping on the new layers
      !!
      !! ------------ cum0(0)                        ------------- cum1(0)
      !!                                    NEW      -------------
      !! ------------ cum0(1)               ==>      -------------
      !!     ...                                     -------------
      !! ------------                                -------------
      !! ------------ cum0(nlay_s+1)                 ------------- cum1(nlay_s)
      !!
      !!
      !! References : Bitz & Lipscomb, JGR 99; Vancoppenolle et al., GRL, 2005
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(jpij,0:nlay_s), INTENT(in   ) ::   ph_old             ! old thicknesses (m)
      REAL(wp), DIMENSION(jpij,0:nlay_s), INTENT(in   ) ::   pe_old             ! old enthlapies (J.m-3)
      REAL(wp), DIMENSION(jpij,1:nlay_s), INTENT(inout) ::   pe_new             ! new enthlapies (J.m-3, remapped)
      !
      INTEGER  :: ji         !  dummy loop indices
      INTEGER  :: jk0, jk1   !  old/new layer indices
      !
      REAL(wp), DIMENSION(jpij,0:nlay_s+1) ::   zeh_cum0, zh_cum0   ! old cumulative enthlapies and layers interfaces
      REAL(wp), DIMENSION(jpij,0:nlay_s)   ::   zeh_cum1, zh_cum1   ! new cumulative enthlapies and layers interfaces
      REAL(wp), DIMENSION(jpij)            ::   zhnew               ! new layers thicknesses
      !!-------------------------------------------------------------------

      !--------------------------------------------------------------------------
      !  1) Cumulative integral of old enthalpy * thickness and layers interfaces
      !--------------------------------------------------------------------------
      zeh_cum0(1:npti,0) = 0._wp
      zh_cum0 (1:npti,0) = 0._wp
      DO jk0 = 1, nlay_s+1
         DO ji = 1, npti
            zeh_cum0(ji,jk0) = zeh_cum0(ji,jk0-1) + pe_old(ji,jk0-1) * ph_old(ji,jk0-1)
            zh_cum0 (ji,jk0) = zh_cum0 (ji,jk0-1) + ph_old(ji,jk0-1)
         END DO
      END DO

      !------------------------------------
      !  2) Interpolation on the new layers
      !------------------------------------
      ! new layer thickesses
      DO ji = 1, npti
         zhnew(ji) = SUM( ph_old(ji,0:nlay_s) ) * r1_nlay_s
      END DO

      ! new layers interfaces
      zh_cum1(1:npti,0) = 0._wp
      DO jk1 = 1, nlay_s
         DO ji = 1, npti
            zh_cum1(ji,jk1) = zh_cum1(ji,jk1-1) + zhnew(ji)
         END DO
      END DO

      zeh_cum1(1:npti,0:nlay_s) = 0._wp
      ! new cumulative q*h => linear interpolation
      DO jk0 = 1, nlay_s+1
         DO jk1 = 1, nlay_s-1
            DO ji = 1, npti
               IF( zh_cum1(ji,jk1) <= zh_cum0(ji,jk0) .AND. zh_cum1(ji,jk1) > zh_cum0(ji,jk0-1) ) THEN
                  zeh_cum1(ji,jk1) = ( zeh_cum0(ji,jk0-1) * ( zh_cum0(ji,jk0) - zh_cum1(ji,jk1  ) ) +  &
                     &                 zeh_cum0(ji,jk0  ) * ( zh_cum1(ji,jk1) - zh_cum0(ji,jk0-1) ) )  &
                     &             / ( zh_cum0(ji,jk0) - zh_cum0(ji,jk0-1) )
               ENDIF
            END DO
         END DO
      END DO
      ! to ensure that total heat content is strictly conserved, set:
      zeh_cum1(1:npti,nlay_s) = zeh_cum0(1:npti,nlay_s+1)

      ! new enthalpies
      DO jk1 = 1, nlay_s
         DO ji = 1, npti
            rswitch      = MAX( 0._wp , SIGN( 1._wp , zhnew(ji) - epsi20 ) )
            pe_new(ji,jk1) = rswitch * ( zeh_cum1(ji,jk1) - zeh_cum1(ji,jk1-1) ) / MAX( zhnew(ji), epsi20 )
         END DO
      END DO

   END SUBROUTINE snw_ent


#else
   !!----------------------------------------------------------------------
   !!   Default option                                NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE snwthd_dh

MODULE icethd_dh
   !!======================================================================
   !!                       ***  MODULE icethd_dh ***
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
   !!   ice_thd_dh        : vertical sea-ice growth and melt
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE ice            ! sea-ice: variables
   USE ice1D          ! sea-ice: thermodynamics variables
   USE icethd_sal     ! sea-ice: salinity profiles
   USE icevar         ! for CALL ice_var_snwblow
   USE icectl         ! sea-ice: control print
   USE snwthd_dh       ! Changes in height due to snow melt (& snowfall if ln_snwext=F)
   USE snwent         ! snow enthalpy remapping

   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)

   USE MODE_SNOW3L    ! For isbaes

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_thd_dh        ! called by ice_thd

   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: icethd_dh.F90 14686 2021-04-08 15:36:01Z clem $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_thd_dh( isnow, zq_rema, zevap_rema, zh_s, ze_s )
      !!------------------------------------------------------------------
      !!                ***  ROUTINE ice_thd_dh  ***
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
      REAL(wp), DIMENSION(jpij), INTENT(in)    ::   isnow       ! presence of snow or not
      REAL(wp), DIMENSION(jpij), INTENT(inout) ::   zq_rema     ! remaining heat flux from snow melting       (J.m-2)
      REAL(wp), DIMENSION(jpij), INTENT(inout) ::   zevap_rema  ! remaining mass flux from snow sublimation   (kg.m-2)
      REAL(wp), DIMENSION(jpij,0:nlay_s  ), INTENT(inout) ::   zh_s      ! snw layer thickness (m) 
      REAL(wp), DIMENSION(jpij,0:nlay_s  ), INTENT(inout) ::   ze_s      ! snw layer enthalpy (J.m-3)

      INTEGER  ::   ji, jk       ! dummy loop indices
      INTEGER  ::   iter         ! local integer

      REAL(wp) ::   ztmelts      ! local scalar
      REAL(wp) ::   zdum
      REAL(wp) ::   zfracs       ! fractionation coefficient for bottom salt entrapment
      REAL(wp) ::   zswi1        ! switch for computation of bottom salinity
      REAL(wp) ::   zswi12       ! switch for computation of bottom salinity
      REAL(wp) ::   zswi2        ! switch for computation of bottom salinity
      REAL(wp) ::   zgrr         ! bottom growth rate
      REAL(wp) ::   zt_i_new     ! bottom formation temperature
      REAL(wp) ::   z1_rho       ! 1/(rhos+rho0-rhoi)

      REAL(wp) ::   zQm          ! enthalpy exchanged with the ocean (J/m2), >0 towards the ocean
      REAL(wp) ::   zEi          ! specific enthalpy of sea ice (J/kg)
      REAL(wp) ::   zEw          ! specific enthalpy of exchanged water (J/kg)
      REAL(wp) ::   zdE          ! specific enthalpy difference (J/kg)
      REAL(wp) ::   zfmdt        ! exchange mass flux x time step (J/m2), >0 towards the ocean

      REAL(wp), DIMENSION(jpij) ::   zq_top      ! heat for surface ablation                   (J.m-2)
      REAL(wp), DIMENSION(jpij) ::   zq_bot      ! heat for bottom ablation                    (J.m-2)
      REAL(wp), DIMENSION(jpij) ::   zf_tt       ! Heat budget to determine melting or freezing(W.m-2)
      REAL(wp), DIMENSION(jpij) ::   zdeltah
      REAL(wp), DIMENSION(jpij) ::   zsnw        ! distribution of snow after wind blowing

      INTEGER , DIMENSION(jpij,nlay_i)     ::   icount    ! number of layers vanishing by melting
      REAL(wp), DIMENSION(jpij,0:nlay_i+1) ::   zh_i      ! ice layer thickness (m)


      REAL(wp), DIMENSION(jpij, nlay_s) ::   rho_s_bef    ! For isbaes
      REAL(wp), DIMENSION(jpij, nlay_s) ::   eh_s_1d      ! 
      REAL(wp), DIMENSION(jpij, nlay_s) ::   oh_s_1d      ! 

      REAL(wp) ::   zswitch_sal

      REAL(wp), DIMENSION(jpij) ::   z1_rho_isbaes       ! 1/(rhos+rho0-rhoi) (isbaes formalism)
      REAL(wp), DIMENSION(jpij) ::   dh_snowice_tmp       ! Temporary snow ice conversion height (saved for isbaes) 
      INTEGER  ::   num_iter_max      ! Heat conservation
      !!------------------------------------------------------------------

      ! Discriminate between time varying salinity and constant
      SELECT CASE( nn_icesal )                  ! varying salinity or not
         CASE( 1 , 3 )   ;   zswitch_sal = 0._wp   ! prescribed salinity profile
         CASE( 2 )       ;   zswitch_sal = 1._wp   ! varying salinity profile
      END SELECT

      ! Theo : snowfall / melt is now computed in snw_thd_dh 
      !
      !                       ! ============================================== !
      !                       !               Snowfall / melt                  !
      !                       ! ============================================== !
      !

      IF( .NOT. (ln_snwext) )  CALL snw_thd_dh(isnow, zq_rema, zevap_rema, zh_s, ze_s)
      DO ji = 1, npti
         zq_top(ji) = zq_rema(ji)
      END DO
      !
      ! initialize ice layer thicknesses and enthalpies
      eh_i_old(1:npti,0:nlay_i+1) = 0._wp
      h_i_old (1:npti,0:nlay_i+1) = 0._wp
      zh_i    (1:npti,0:nlay_i+1) = 0._wp
      DO jk = 1, nlay_i
         DO ji = 1, npti
            eh_i_old(ji,jk) = h_i_1d(ji) * r1_nlay_i * e_i_1d(ji,jk)
            h_i_old (ji,jk) = h_i_1d(ji) * r1_nlay_i
            zh_i    (ji,jk) = h_i_1d(ji) * r1_nlay_i
         END DO
      END DO
      !
      IF( ln_snwext .OR. ln_isbaes) THEN
              
         ! initialize snw layer thicknesses and enthalpies
         zh_s(1:npti,0) = 0._wp
         ze_s(1:npti,0) = 0._wp
         DO jk = 1, nlay_s
            DO ji = 1, npti
               IF(ln_isbaes) THEN
                  zh_s(ji,jk) = dh_s_1d(ji,jk) 
                  ze_s(ji,jk) = e_s_1d(ji,jk)
               ELSE        
                  zh_s(ji,jk) = h_s_1d(ji) * r1_nlay_s
                  ze_s(ji,jk) = e_s_1d(ji,jk)
               ENDIF
            END DO
         END DO
      ENDIF
      PRINT*,'H_s_1d before sni', h_s_1d(1)
      !
      !                       ! ============================================== !
      !                       ! Available heat for surface and bottom ablation !
      !                       ! ============================================== !
      !
      !
      DO ji = 1, npti
         zf_tt(ji)         = qcn_ice_bot_1d(ji) + qsb_ice_bot_1d(ji) + fhld_1d(ji) + qtr_ice_bot_1d(ji) * frq_m_1d(ji)
         zq_bot(ji)        = MAX( 0._wp, zf_tt(ji) * rDt_ice )
      END DO
      !
      !                       ! ============ !
      !                       !     Ice      !
      !                       ! ============ !

      ! Surface ice melting
      !--------------------
      DO jk = 1, nlay_i
         DO ji = 1, npti
            ztmelts = - rTmlt * sz_i_1d(ji,jk)   ! Melting point of layer k [C]

            IF( t_i_1d(ji,jk) >= (ztmelts+rt0) ) THEN   !-- Internal melting

               zEi            = - e_i_1d(ji,jk) * r1_rhoi             ! Specific enthalpy of layer k [J/kg, <0]
               zdE            =   0._wp                               ! Specific enthalpy difference (J/kg, <0)
               !                                                          set up at 0 since no energy is needed to melt water...(it is already melted)
               zdum           = MIN( 0._wp , - zh_i(ji,jk) )          ! internal melting occurs when the internal temperature is above freezing
               !                                                          this should normally not happen, but sometimes, heat diffusion leads to this
               zfmdt          = - zdum * rhoi                         ! Recompute mass flux [kg/m2, >0]
               !
               dh_i_itm(ji)   = dh_i_itm(ji) + zdum                   ! Cumulate internal melting
               !
               hfx_res_1d(ji) = hfx_res_1d(ji) + zEi  * zfmdt             * a_i_1d(ji) * r1_Dt_ice    ! Heat flux to the ocean [W.m-2], <0
               !                                                                                          ice enthalpy zEi is "sent" to the ocean
               wfx_res_1d(ji) = wfx_res_1d(ji) - rhoi * zdum              * a_i_1d(ji) * r1_Dt_ice    ! Mass flux
               sfx_res_1d(ji) = sfx_res_1d(ji) - rhoi * zdum * s_i_1d(ji) * a_i_1d(ji) * r1_Dt_ice    ! Salt flux
               !                                                                                          using s_i_1d and not sz_i_1d(jk) is ok
            ELSE                                        !-- Surface melting

               zEi            = - e_i_1d(ji,jk) * r1_rhoi             ! Specific enthalpy of layer k [J/kg, <0]
               zEw            =    rcp * ztmelts                      ! Specific enthalpy of resulting meltwater [J/kg, <0]
               zdE            =    zEi - zEw                          ! Specific enthalpy difference < 0

               zfmdt          = - zq_top(ji) / zdE                    ! Mass flux to the ocean [kg/m2, >0]

               zdum           = - zfmdt * r1_rhoi                     ! Melt of layer jk [m, <0]

               zdum           = MIN( 0._wp , MAX( zdum , - zh_i(ji,jk) ) )    ! Melt of layer jk cannot exceed the layer thickness [m, <0]

               zq_top(ji)     = MAX( 0._wp , zq_top(ji) - zdum * rhoi * zdE ) ! update available heat

               dh_i_sum(ji)   = dh_i_sum(ji) + zdum                   ! Cumulate surface melt

               zfmdt          = - rhoi * zdum                         ! Recompute mass flux [kg/m2, >0]

               zQm            = zfmdt * zEw                           ! Energy of the melt water sent to the ocean [J/m2, <0]

               hfx_thd_1d(ji) = hfx_thd_1d(ji) + zEw  * zfmdt             * a_i_1d(ji) * r1_Dt_ice    ! Heat flux [W.m-2], < 0
               hfx_sum_1d(ji) = hfx_sum_1d(ji) - zdE  * zfmdt             * a_i_1d(ji) * r1_Dt_ice    ! Heat flux used in this process [W.m-2], > 0
               wfx_sum_1d(ji) = wfx_sum_1d(ji) - rhoi * zdum              * a_i_1d(ji) * r1_Dt_ice    ! Mass flux
               sfx_sum_1d(ji) = sfx_sum_1d(ji) - rhoi * zdum * s_i_1d(ji) * a_i_1d(ji) * r1_Dt_ice    ! Salt flux >0
               !                                                                                          using s_i_1d and not sz_i_1d(jk) is ok)
            END IF
            ! update thickness
            zh_i(ji,jk) = MAX( 0._wp, zh_i(ji,jk) + zdum )
            h_i_1d(ji)  = MAX( 0._wp, h_i_1d(ji)  + zdum )
            !
            ! update heat content (J.m-2) and layer thickness
            eh_i_old(ji,jk) = eh_i_old(ji,jk) + zdum * e_i_1d(ji,jk)
            h_i_old (ji,jk) = h_i_old (ji,jk) + zdum
            !
            !
            ! Ice sublimation
            ! ---------------
            zdum               = MAX( - zh_i(ji,jk) , - zevap_rema(ji) * r1_rhoi )
            !
            hfx_sub_1d(ji)     = hfx_sub_1d(ji)     + e_i_1d(ji,jk) * zdum              * a_i_1d(ji) * r1_Dt_ice ! Heat flux [W.m-2], < 0
            wfx_ice_sub_1d(ji) = wfx_ice_sub_1d(ji) - rhoi          * zdum              * a_i_1d(ji) * r1_Dt_ice ! Mass flux > 0
            sfx_sub_1d(ji)     = sfx_sub_1d(ji)     - rhoi          * zdum * s_i_1d(ji) * a_i_1d(ji) * r1_Dt_ice ! Salt flux >0
            !                                                                                                      clem: flux is sent to the ocean for simplicity
            !                                                                                                            but salt should remain in the ice except
            !                                                                                                            if all ice is melted. => must be corrected
            ! update remaining mass flux and thickness
            zevap_rema(ji) = zevap_rema(ji) + zdum * rhoi
            zh_i(ji,jk)    = MAX( 0._wp, zh_i(ji,jk) + zdum )
            h_i_1d(ji)     = MAX( 0._wp, h_i_1d(ji)  + zdum )
            dh_i_sub(ji)   = dh_i_sub(ji) + zdum

            ! update heat content (J.m-2) and layer thickness
            eh_i_old(ji,jk) = eh_i_old(ji,jk) + zdum * e_i_1d(ji,jk)
            h_i_old (ji,jk) = h_i_old (ji,jk) + zdum

            ! record which layers have disappeared (for bottom melting)
            !    => icount=0 : no layer has vanished
            !    => icount=5 : 5 layers have vanished
            rswitch       = MAX( 0._wp , SIGN( 1._wp , - zh_i(ji,jk) ) )
            icount(ji,jk) = NINT( rswitch )

         END DO
      END DO
      PRINT*,'zevap rema',zevap_rema
      ! remaining "potential" evap is sent to ocean
      DO ji = 1, npti
         wfx_err_sub_1d(ji) = wfx_err_sub_1d(ji) - zevap_rema(ji) * a_i_1d(ji) * r1_Dt_ice  ! <=0 (net evap for the ocean in kg.m-2.s-1)
      END DO


      ! Ice Basal growth
      !------------------
      ! Basal growth is driven by heat imbalance at the ice-ocean interface,
      ! between the inner conductive flux  (qcn_ice_bot), from the open water heat flux
      ! (fhld) and the sensible ice-ocean flux (qsb_ice_bot).
      ! qcn_ice_bot is positive downwards. qsb_ice_bot and fhld are positive to the ice

      ! If salinity varies in time, an iterative procedure is required, because
      ! the involved quantities are inter-dependent.
      ! Basal growth (dh_i_bog) depends upon new ice specific enthalpy (zEi),
      ! which depends on forming ice salinity (s_i_new), which depends on dh/dt (dh_i_bog)
      ! -> need for an iterative procedure, which converges quickly

      num_iter_max = 1
      IF( nn_icesal == 2 )   num_iter_max = 5  ! salinity varying in time

      DO ji = 1, npti
         IF(  zf_tt(ji) < 0._wp  ) THEN
            DO iter = 1, num_iter_max   ! iterations

               ! New bottom ice salinity (Cox & Weeks, JGR88 )
               !--- zswi1  if dh/dt < 2.0e-8
               !--- zswi12 if 2.0e-8 < dh/dt < 3.6e-7
               !--- zswi2  if dh/dt > 3.6e-7
               zgrr     = MIN( 1.0e-3, MAX ( dh_i_bog(ji) * r1_Dt_ice , epsi10 ) )
               zswi2    = MAX( 0._wp , SIGN( 1._wp , zgrr - 3.6e-7 ) )
               zswi12   = MAX( 0._wp , SIGN( 1._wp , zgrr - 2.0e-8 ) ) * ( 1.0 - zswi2 )
               zswi1    = 1. - zswi2 * zswi12
               zfracs   = MIN( zswi1  * 0.12 + zswi12 * ( 0.8925 + 0.0568 * LOG( 100.0 * zgrr ) )   &
                  &          + zswi2  * 0.26 / ( 0.26 + 0.74 * EXP ( - 724300.0 * zgrr ) )  , 0.5 )

               s_i_new(ji)    = zswitch_sal * zfracs * sss_1d(ji) + ( 1. - zswitch_sal ) * s_i_1d(ji)  ! New ice salinity

               ztmelts        = - rTmlt * s_i_new(ji)                                                  ! New ice melting point (C)

               zt_i_new       = zswitch_sal * t_bo_1d(ji) + ( 1. - zswitch_sal) * t_i_1d(ji, nlay_i)

               zEi            = rcpi * ( zt_i_new - (ztmelts+rt0) ) &                                  ! Specific enthalpy of forming ice (J/kg, <0)
                  &             - rLfus * ( 1.0 - ztmelts / ( MIN( zt_i_new - rt0, -epsi10 ) ) ) + rcp * ztmelts

               zEw            = rcp  * ( t_bo_1d(ji) - rt0 )                                           ! Specific enthalpy of seawater (J/kg, < 0)

               zdE            = zEi - zEw                                                              ! Specific enthalpy difference (J/kg, <0)

               dh_i_bog(ji)   = rDt_ice * MAX( 0._wp , zf_tt(ji) / ( zdE * rhoi ) )

            END DO
            ! Contribution to Energy and Salt Fluxes
            zfmdt = - rhoi * dh_i_bog(ji)                                                              ! Mass flux x time step (kg/m2, < 0)

            hfx_thd_1d(ji) = hfx_thd_1d(ji) + zEw  * zfmdt                      * a_i_1d(ji) * r1_Dt_ice   ! Heat flux to the ocean [W.m-2], >0
            hfx_bog_1d(ji) = hfx_bog_1d(ji) - zdE  * zfmdt                      * a_i_1d(ji) * r1_Dt_ice   ! Heat flux used in this process [W.m-2], <0
            wfx_bog_1d(ji) = wfx_bog_1d(ji) - rhoi * dh_i_bog(ji)               * a_i_1d(ji) * r1_Dt_ice   ! Mass flux, <0
            sfx_bog_1d(ji) = sfx_bog_1d(ji) - rhoi * dh_i_bog(ji) * s_i_new(ji) * a_i_1d(ji) * r1_Dt_ice   ! Salt flux, <0

            ! update thickness
            zh_i(ji,nlay_i+1) = zh_i(ji,nlay_i+1) + dh_i_bog(ji)
            h_i_1d(ji)        = h_i_1d(ji)        + dh_i_bog(ji)

            ! update heat content (J.m-2) and layer thickness
            eh_i_old(ji,nlay_i+1) = eh_i_old(ji,nlay_i+1) + dh_i_bog(ji) * (-zEi * rhoi)
            h_i_old (ji,nlay_i+1) = h_i_old (ji,nlay_i+1) + dh_i_bog(ji)

         ENDIF

      END DO

      ! Ice Basal melt
      !---------------
      DO jk = nlay_i, 1, -1
         DO ji = 1, npti
            IF(  zf_tt(ji)  >  0._wp  .AND. jk > icount(ji,jk) ) THEN   ! do not calculate where layer has already disappeared by surface melting

               ztmelts = - rTmlt * sz_i_1d(ji,jk)  ! Melting point of layer jk (C)

               IF( t_i_1d(ji,jk) >= (ztmelts+rt0) ) THEN   !-- Internal melting

                  zEi            = - e_i_1d(ji,jk) * r1_rhoi     ! Specific enthalpy of melting ice (J/kg, <0)
                  zdE            = 0._wp                         ! Specific enthalpy difference   (J/kg, <0)
                  !                                                  set up at 0 since no energy is needed to melt water...(it is already melted)
                  zdum           = MIN( 0._wp , - zh_i(ji,jk) )  ! internal melting occurs when the internal temperature is above freezing
                  !                                                  this should normally not happen, but sometimes, heat diffusion leads to this
                  dh_i_itm (ji)  = dh_i_itm(ji) + zdum
                  !
                  zfmdt          = - zdum * rhoi                 ! Mass flux x time step > 0
                  !
                  hfx_res_1d(ji) = hfx_res_1d(ji) + zEi  * zfmdt             * a_i_1d(ji) * r1_Dt_ice   ! Heat flux to the ocean [W.m-2], <0
                  !                                                                                         ice enthalpy zEi is "sent" to the ocean
                  wfx_res_1d(ji) = wfx_res_1d(ji) - rhoi * zdum              * a_i_1d(ji) * r1_Dt_ice   ! Mass flux
                  sfx_res_1d(ji) = sfx_res_1d(ji) - rhoi * zdum * s_i_1d(ji) * a_i_1d(ji) * r1_Dt_ice   ! Salt flux
                  !                                                                                         using s_i_1d and not sz_i_1d(jk) is ok
               ELSE                                        !-- Basal melting

                  zEi            = - e_i_1d(ji,jk) * r1_rhoi                       ! Specific enthalpy of melting ice (J/kg, <0)
                  zEw            = rcp * ztmelts                                   ! Specific enthalpy of meltwater (J/kg, <0)
                  zdE            = zEi - zEw                                       ! Specific enthalpy difference   (J/kg, <0)

                  zfmdt          = - zq_bot(ji) / zdE                              ! Mass flux x time step (kg/m2, >0)

                  zdum           = - zfmdt * r1_rhoi                               ! Gross thickness change

                  zdum           = MIN( 0._wp , MAX( zdum, - zh_i(ji,jk) ) )       ! bound thickness change

                  zq_bot(ji)     = MAX( 0._wp , zq_bot(ji) - zdum * rhoi * zdE )   ! update available heat. MAX is necessary for roundup errors

                  dh_i_bom(ji)   = dh_i_bom(ji) + zdum                             ! Update basal melt

                  zfmdt          = - zdum * rhoi                                   ! Mass flux x time step > 0

                  zQm            = zfmdt * zEw                                     ! Heat exchanged with ocean

                  hfx_thd_1d(ji) = hfx_thd_1d(ji) + zEw  * zfmdt             * a_i_1d(ji) * r1_Dt_ice   ! Heat flux to the ocean [W.m-2], <0
                  hfx_bom_1d(ji) = hfx_bom_1d(ji) - zdE  * zfmdt             * a_i_1d(ji) * r1_Dt_ice   ! Heat used in this process [W.m-2], >0
                  wfx_bom_1d(ji) = wfx_bom_1d(ji) - rhoi * zdum              * a_i_1d(ji) * r1_Dt_ice   ! Mass flux
                  sfx_bom_1d(ji) = sfx_bom_1d(ji) - rhoi * zdum * s_i_1d(ji) * a_i_1d(ji) * r1_Dt_ice   ! Salt flux
                  !                                                                                         using s_i_1d and not sz_i_1d(jk) is ok
               ENDIF
               ! update thickness
               zh_i(ji,jk) = MAX( 0._wp, zh_i(ji,jk) + zdum )
               h_i_1d(ji)  = MAX( 0._wp, h_i_1d(ji)  + zdum )
               !
               ! update heat content (J.m-2) and layer thickness
               eh_i_old(ji,jk) = eh_i_old(ji,jk) + zdum * e_i_1d(ji,jk)
               h_i_old (ji,jk) = h_i_old (ji,jk) + zdum
            ENDIF
         END DO
      END DO
      PRINT*,'ICE HEIGHT',h_i_1d(1)
      IF(ln_isbaes) THEN
         ! With isbaes, snow grid height are not constant along the vertical
         ! Thus, changes in snow height due to snow to ice conversion may impact cells with different size on the vertical
         ! Therefore, enthalpies and snow age are multiplied by the grid thicknesses and we use these quantities instead      
         DO ji = 1, npti
            eh_s_1d(ji,:) = e_s_1D(ji,:) * dh_s_1d(ji,:) !* a_i_1d(ji)
            oh_s_1d(ji,:) = o_s_1d(ji,:) * dh_s_1d(ji,:) !* a_i_1d(ji)
         END DO
      ENDIF

      ! Remove snow if ice has melted entirely
      ! --------------------------------------
      DO jk = 0, nlay_s
         DO ji = 1,npti
            IF( h_i_1d(ji) == 0._wp ) THEN
               ! mass & energy loss to the ocean
               hfx_res_1d(ji) = hfx_res_1d(ji) - ze_s(ji,jk) * zh_s(ji,jk) * a_i_1d(ji) * r1_Dt_ice  ! heat flux to the ocean [W.m-2], < 0
               IF( ln_isbaes) THEN
                   ! Mass flux is computed from 3D density arrays instead of constant density
                   wfx_res_1d(ji) = wfx_res_1d(ji) + rho_s_1d(ji,jk)        * zh_s(ji,jk) * a_i_1d(ji) * r1_Dt_ice  ! mass flux

                   dh_s_1d(ji,jk) = 0._wp
                   swe_s_1d(ji,jk) = 0._wp
               ELSE
                    wfx_res_1d(ji) = wfx_res_1d(ji) + rhos        * zh_s(ji,jk) * a_i_1d(ji) * r1_Dt_ice  ! mass flux
               ENDIF

               wfx_res_1d(ji) = wfx_res_1d(ji) + rhos        * zh_s(ji,jk) * a_i_1d(ji) * r1_Dt_ice  ! mass flux

               ! update thickness and energy
               h_s_1d(ji)    = 0._wp
               ze_s  (ji,jk) = 0._wp
               zh_s  (ji,jk) = 0._wp
            ENDIF
         END DO
      END DO

      ! Snow load on ice
      ! -----------------
      ! When snow load exceeds Archimede's limit and sst is positive,
      ! snow-ice formation (next bloc) can lead to negative ice enthalpy.
      ! Therefore we consider here that this excess of snow falls into the ocean
      IF(ln_isbaes) THEN
         DO ji = 1, npti
             IF(h_s_1d(ji) > 0.00000000000000001) THEN
                zdeltah(ji) = h_s_1d(ji) + h_i_1d(ji) * (rhoi-rho0) * (1._wp / (SUM(rho_s_1d(ji,:)*dh_s_1d(ji,:))/h_s_1d(ji)))
             ELSE
                zdeltah(ji) = 0._wp
             ENDIF
         END DO
      ELSE
         zdeltah(1:npti) = h_s_1d(1:npti) + h_i_1d(1:npti) * (rhoi-rho0) * r1_rhos
      ENDIF
      
      DO jk = 0, nlay_s
         DO ji = 1, npti
            IF( zdeltah(ji) > 0._wp .AND. sst_1d(ji) > 0._wp ) THEN
               ! snow layer thickness that falls into the ocean
               zdum = MIN( zdeltah(ji) , zh_s(ji,jk) )
               ! mass & energy loss to the ocean
               hfx_res_1d(ji) = hfx_res_1d(ji) - ze_s(ji,jk) * zdum * a_i_1d(ji) * r1_Dt_ice  ! heat flux to the ocean [W.m-2], < 0
               IF(ln_isbaes) THEN
                  wfx_res_1d(ji) = wfx_res_1d(ji) + rho_s_1d(ji,jk)        * zdum * a_i_1d(ji) * r1_Dt_ice  ! mass flux
                  ! update thickness and energy
                  eh_s_1d(ji,jk)= MAX( 0._wp, eh_s_1d(ji,jk) - zdum * e_s_1d(ji,jk))
               ELSE
                  wfx_res_1d(ji) = wfx_res_1d(ji) + rhos        * zdum * a_i_1d(ji) * r1_Dt_ice  ! mass flux
               ENDIF 
               h_s_1d(ji)    = MAX( 0._wp, h_s_1d(ji)  - zdum )
               zh_s  (ji,jk) = MAX( 0._wp, zh_s(ji,jk) - zdum )
               ! update snow thickness that still has to fall
               zdeltah(ji)   = MAX( 0._wp, zdeltah(ji) - zdum )
            ENDIF
         END DO
      END DO

      ! Snow-Ice formation
      ! ------------------
      ! When snow load exceeds Archimede's limit, snow-ice interface goes down under sea-level,
      ! flooding of seawater transforms snow into ice. Thickness that is transformed is dh_snowice (positive for the ice)
      z1_rho = 1._wp / ( rhos+rho0-rhoi )
      IF(ln_isbaes) THEN
         DO ji = 1, npti
           IF(h_s_1d(ji) > 0.00000000000000001) THEN
               z1_rho_isbaes(ji) = 1._wp / ( SUM(rho_s_1d(ji,:)*dh_s_1d(ji,:))/h_s_1d(ji) + rho0-rhoi )
           ELSE
               z1_rho_isbaes(ji) = 0._wp ! Not used anyway, so we can assign a random value
           ENDIF
         END DO
      ENDIF

      zdeltah(1:npti) = 0._wp
      DO ji = 1, npti
         !
         IF(h_s_1d(ji) > 0.00000000000000001) THEN
            IF(ln_isbaes) THEN
               ! dh_snowice is computed from a variying density
               dh_snowice(ji) = MAX( 0._wp , ( SUM(rho_s_1d(ji,:) * dh_s_1d(ji,:)) + (rhoi-rho0) * h_i_1d(ji) ) * z1_rho_isbaes(ji) )
            ELSE     
               dh_snowice(ji) = MAX( 0._wp , ( rhos * h_s_1d(ji) + (rhoi-rho0) * h_i_1d(ji) ) * z1_rho )
            ENDIF
            h_i_1d(ji)    = h_i_1d(ji) + dh_snowice(ji)
            h_s_1d(ji)    = h_s_1d(ji) - dh_snowice(ji)
            ! Mass flux: All snow is thrown in the ocean, and seawater is taken to replace the volume
            wfx_sni_1d    (ji) = wfx_sni_1d    (ji) - dh_snowice(ji) * rhoi * a_i_1d(ji) * r1_Dt_ice

            IF(ln_isbaes) THEN
               ! We compute heat & salt fluxes from a varying density
               zdum = SUM(rho_s_1d(ji,:)*dh_s_1d(ji,:))/h_s_1d(ji)
               zfmdt          = ( zdum - rhoi ) * dh_snowice(ji)    ! <0
               zEw            = rcp * sst_1d(ji)
               zQm            = zfmdt * zEw

               hfx_thd_1d(ji) = hfx_thd_1d(ji)  + zEw        * zfmdt * a_i_1d(ji) * r1_Dt_ice ! Heat flux
               sfx_sni_1d(ji) = sfx_sni_1d(ji) + sss_1d(ji) * zfmdt * a_i_1d(ji) * r1_Dt_ice ! Salt flux

                  ! Case constant salinity in time: virtual salt flux to keep salinity constant
                  IF( nn_icesal /= 2 )  THEN
                     sfx_bri_1d(ji) = sfx_bri_1d(ji) - sss_1d(ji) * zfmdt                 * a_i_1d(ji) * r1_Dt_ice  &  ! put back sss_m     into the ocean
                        &                            - s_i_1d(ji) * dh_snowice(ji) * rhoi * a_i_1d(ji) * r1_Dt_ice     ! and get  rn_icesal from the ocean
                  ENDIF
            ELSE
               zfmdt          = ( rhos - rhoi ) * dh_snowice(ji)    ! <0
               zEw            = rcp * sst_1d(ji)
               zQm            = zfmdt * zEw

               hfx_thd_1d(ji) = hfx_thd_1d(ji)  + zEw        * zfmdt * a_i_1d(ji) * r1_Dt_ice ! Heat flux
               sfx_sni_1d(ji) = sfx_sni_1d(ji) + sss_1d(ji) * zfmdt * a_i_1d(ji) * r1_Dt_ice ! Salt flux

                  ! Case constant salinity in time: virtual salt flux to keep salinity constant
                  IF( nn_icesal /= 2 )  THEN
                     sfx_bri_1d(ji) = sfx_bri_1d(ji) - sss_1d(ji) * zfmdt                 * a_i_1d(ji) * r1_Dt_ice  &  ! put back sss_m     into the ocean
                        &                            - s_i_1d(ji) * dh_snowice(ji) * rhoi * a_i_1d(ji) * r1_Dt_ice     ! and get  rn_icesal from the ocean
                  ENDIF

                  wfx_snw_sni_1d(ji) = wfx_snw_sni_1d(ji) + dh_snowice(ji) * rhos * a_i_1d(ji) * r1_Dt_ice
            ENDIF        
         
            ! update thickness
            zh_i(ji,0)  = zh_i(ji,0) + dh_snowice(ji)
            zdeltah(ji) =              dh_snowice(ji)

            ! update heat content (J.m-2) and layer thickness
            h_i_old (ji,0) = h_i_old (ji,0) + dh_snowice(ji)
            eh_i_old(ji,0) = eh_i_old(ji,0) + zfmdt * zEw           ! 1st part (sea water enthalpy)

         ELSE ! IF(h_s_1d(ji) > 0.00000000000000001
            ! we set fluxes to 0. (this is just to avoid nans)    
            zdeltah(ji) = 0._wp
            dh_snowice(ji) = 0._wp
            wfx_sni_1d(ji) = 0._wp
            wfx_snw_sni_1d(ji) = 0._wp
         ENDIF 
      END DO
      !
      DO jk = nlay_s, 1, -1   ! flooding of snow starts from the base
         DO ji = 1, npti
            zdum           = MIN( zdeltah(ji), zh_s(ji,jk) )     ! amount of snw that floods, > 0
            zh_s(ji,jk)    = MAX( 0._wp, zh_s(ji,jk) - zdum )    ! remove some snow thickness
            eh_i_old(ji,0) = eh_i_old(ji,0) + zdum * ze_s(ji,jk) ! 2nd part (snow enthalpy)

            IF(ln_isbaes) THEN
               ! When ln_isbaes=T, remove enthalpy per snow layer
               ! => This is necessary as thicknesses are constant 
               ! 
               eh_s_1d(ji,jk) = eh_s_1d(ji,jk) - zdum * ze_s(ji,jk)

               ! Mass flux is computed from 3D densities
               wfx_snw_sni_1d(ji) = wfx_snw_sni_1d(ji) + zdum * rho_s_1d(ji,jk) *a_i_1d(ji) * r1_Dt_ice
            ENDIF
            ! update dh_snowice
            zdeltah(ji)    = MAX( 0._wp, zdeltah(ji) - zdum )
         END DO
      END DO

      IF(ln_isbaes) THEN
           ! If ln_isbaes:
           ! - We compute a new grid according to the new snow height after snow / ice conversion
           ! - Then, we regrid snow quantities on this new grid
           ! This is necessary, since those quantities will be advected later on

           DO ji = 1, npti
              h_s_1d(ji) = SUM(zh_s(ji,:)) ! Compute new total snow height
           END DO
           rho_s_bef(:,:) = rho_s_1d(:,:)
           
           ! We recompute the isbaes grid from the new snow height            
           CALL SNOW3LGRID(dh_s_1d(:,1:nlay_s),h_s_1d(:),PSNOWDZ_OLD=zh_s(:,1:nlay_s))

           !
           ! Mass/Heat redistribution:
           !
           ! Regrid quantities on the new grid 
           CALL SNOW3LTRANSF(h_s_1d(:),zh_s(:,1:nlay_s), dh_s_1d(:,:),rho_s_1d(:,:),eh_s_1d(:,:),oh_s_1d(:,:))
           !

           DO ji = 1, npti
              IF(h_s_1d(ji) > 0.00000000000000001) THEN
                 ! Recompute snow enthalpies and age 
                 dh_s_1d(ji,:) = zh_s(ji,1:nlay_s)
   
                 e_s_1d(ji,:) = eh_s_1d(ji,:) / (dh_s_1d(ji,:)) 
                 o_s_1d(ji,:) = oh_s_1d(ji,:) / (dh_s_1d(ji,:)) 
              ELSE
                 ! We reset rho to previous value when snow height = 0 just to avoid having nans
                 rho_s_1d(:,:) = rho_s_bef(:,:) 
              ENDIF

           END DO
       ENDIF
           ! recalculate t_s_1d from e_s_1d
           DO jk = 1, nlay_s
              DO ji = 1,npti
                 IF( h_s_1d(ji) > 0._wp ) THEN
                    IF(ln_isbaes) THEN
                       t_s_1d(ji,jk) = rt0 + ( - e_s_1d(ji,jk) * (1 / rho_s_1d(ji,jk)) * r1_rcpi + rLfus * r1_rcpi ) 
                    ELSE        
                       t_s_1d(ji,jk) = rt0 + ( - e_s_1d(ji,jk) * r1_rhos * r1_rcpi + rLfus * r1_rcpi )
                    ENDIF
                 ELSE
                    t_s_1d(ji,jk) = rt0
                 ENDIF
              END DO
           END DO
      PRINT*,'H_s_1d after sni', h_s_1d(1)

!      !
!      !
!!!!$      ! --- Update snow diags --- !
!!!!$      !!clem: this is wrong. dh_s_tot is not used anyway
!!!!$      DO ji = 1, npti
!!!!$         dh_s_tot(ji) = dh_s_tot(ji) + dh_s_mlt(ji) + zdeltah(ji) + zdh_s_sub(ji) - dh_snowice(ji)
!!!!$      END DO
!!      !
!!      !
!      ! Remapping of snw enthalpy on a regular grid
!      !--------------------------------------------
!      CALL snw_ent( zh_s, ze_s, e_s_1d)
!
!      ! recalculate t_s_1d from e_s_1d
!      DO jk = 1, nlay_s
!         DO ji = 1,npti
!            IF( h_s_1d(ji) > 0._wp ) THEN
!               IF(ln_isbaes) THEN
!                  t_s_1d(ji,jk) = rt0 + ( - e_s_1d(ji,jk) * (1 / rho_s_1d(ji,jk)) * r1_rcpi + rLfus * r1_rcpi ) 
!               ELSE        
!                  t_s_1d(ji,jk) = rt0 + ( - e_s_1d(ji,jk) * r1_rhos * r1_rcpi + rLfus * r1_rcpi )
!               ENDIF
!            ELSE
!               t_s_1d(ji,jk) = rt0
!            ENDIF
!         END DO
!      END DO
!
      ! Note: remapping of ice enthalpy is done in icethd.F90

      ! --- ensure that a_i = 0 & h_s = 0 where h_i = 0 ---
      WHERE( h_i_1d(1:npti) == 0._wp )
         a_i_1d (1:npti) = 0._wp
         h_s_1d (1:npti) = 0._wp
         t_su_1d(1:npti) = rt0
      END WHERE

   END SUBROUTINE ice_thd_dh

#else
   !!----------------------------------------------------------------------
   !!   Default option                                NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icethd_dh

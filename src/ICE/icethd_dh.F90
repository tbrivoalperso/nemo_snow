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
   USE snwthd_iceconv
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)

#if defined key_isbaes   
   USE MODE_SNOW3L    ! For isbaes
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_thd_dh        ! called by ice_thd

   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: icethd_dh.F90 14686 2021-04-08 15:36:01Z clem $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_thd_dh( isnow, zq_rema, zevap_rema,thickness_si, mass_si, enthalpy_si, zh_s, ze_s )
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
      REAL(wp), DIMENSION(jpij), INTENT(inout) ::   thickness_si 
      REAL(wp), DIMENSION(jpij), INTENT(inout) ::   mass_si
      REAL(wp), DIMENSION(jpij), INTENT(inout) ::   enthalpy_si

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

      REAL(wp) ::   zswitch_sal

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

      ! Remove snow if ice has melted entirely
      ! --------------------------------------
      DO jk = 0, nlay_s
         DO ji = 1,npti
            IF( h_i_1d(ji) == 0._wp ) THEN
               ! mass & energy loss to the ocean
               hfx_res_1d(ji) = hfx_res_1d(ji) - e_s_1d(ji,jk) * r1_Dt_ice  ! heat flux to the ocean [W.m-2], < 0
#if defined key_isbaes
               IF( ln_isbaes) THEN
                   ! Mass flux is computed from 3D density arrays instead of constant density
                   wfx_res_1d(ji) = wfx_res_1d(ji) + rho_s_1d(ji,jk)        * zh_s(ji,jk) * a_i_1d(ji) * r1_Dt_ice  ! mass flux

                   dh_s_1d(ji,jk) = 0._wp
                   swe_s_1d(ji,jk) = 0._wp
               ELSE
                    wfx_res_1d(ji) = wfx_res_1d(ji) + rhos        * zh_s(ji,jk) * a_i_1d(ji) * r1_Dt_ice  ! mass flux
               ENDIF
               h_s_1d(ji)    = 0._wp
               e_s_1d(ji,jk)    = 0._wp
               ze_s  (ji,jk) = 0._wp
               zh_s  (ji,jk) = 0._wp
               dh_s_1d(ji,jk) = 0._wp
               rhov_s_1d(ji,jk) = 0._wp

#else
               wfx_res_1d(ji) = wfx_res_1d(ji) + rhos        * zh_s(ji,jk) * a_i_1d(ji) * r1_Dt_ice  ! mass flux

               ! update thickness and energy
               h_s_1d(ji)    = 0._wp
               ze_s  (ji,jk) = 0._wp
               zh_s  (ji,jk) = 0._wp
#endif
            ENDIF
         END DO
      END DO


!#if defined key_isbaes
!      IF(.NOT. ln_isbaes) THEN
!#endif
         PRINT*,'MASS 1', SUM(rho_s_1d(1,:) * dh_s_1d(1,:) * a_i_1d(1))
         PRINT*,'enthalpy before snowice', SUM(e_s_1d(1,:))
         DO ji = 1, npti
            ! We initialise thickness_si, mass_si, enthalpy_si here so that a ln_snwiceconv could be easily added
            mass_si(ji) = 0; enthalpy_si(ji) = 0; thickness_si(ji) = 0; ! fields sent to sea ice model
         END DO

         CALL snw_thd_iceconv( isnow, zh_s, ze_s, thickness_si, mass_si, enthalpy_si )
         PRINT*,'enthalpy after snowice', SUM(e_s_1d(1,:))

         PRINT*,'MASS 2', SUM(rho_s_1d(1,:) * dh_s_1d(1,:) * a_i_1d(1))
         PRINT*,'thickness_si, mass_si, enthalpy_si', thickness_si(1), mass_si(1), enthalpy_si(1)
         ! Snow-Ice formation
         ! ------------------
         ! When snow load exceeds Archimede's limit, snow-ice interface goes down under sea-level,
         ! flooding of seawater transforms snow into ice. Thickness that is transformed is dh_snowice (positive for the ice)
         z1_rho = 1._wp / ( rhos+rho0-rhoi )
         zdeltah(1:npti) = 0._wp
         DO ji = 1, npti
            !
            dh_snowice(ji) = thickness_si(ji) !MAX( 0._wp , ( rhos * h_s_1d(ji) + (rhoi-rho0) * h_i_1d(ji) ) * z1_rho )

            h_i_1d(ji)    = h_i_1d(ji) + dh_snowice(ji)
            !h_s_1d(ji)    = h_s_1d(ji) - dh_snowice(ji)

            ! Contribution to energy flux to the ocean [J/m2], >0 (if sst<0)
            zfmdt          =  mass_si(ji) - rhoi * dh_snowice(ji)    ! <0
            zEw            = rcp * sst_1d(ji)
            zQm            = zfmdt * zEw

            hfx_thd_1d(ji) = hfx_thd_1d(ji) + zEw        * zfmdt * a_i_1d(ji) * r1_Dt_ice ! Heat flux
            sfx_sni_1d(ji) = sfx_sni_1d(ji) + sss_1d(ji) * zfmdt * a_i_1d(ji) * r1_Dt_ice ! Salt flux

            ! Case constant salinity in time: virtual salt flux to keep salinity constant
            IF( nn_icesal /= 2 )  THEN
               sfx_bri_1d(ji) = sfx_bri_1d(ji) - sss_1d(ji) * zfmdt                 * a_i_1d(ji) * r1_Dt_ice  &  ! put back sss_m     into the ocean
                  &                            - s_i_1d(ji) * dh_snowice(ji) * rhoi * a_i_1d(ji) * r1_Dt_ice     ! and get  rn_icesal from the ocean
            ENDIF

            ! Mass flux: All snow is thrown in the ocean, and seawater is taken to replace the volume
            wfx_sni_1d    (ji) = wfx_sni_1d    (ji) - dh_snowice(ji) * rhoi * a_i_1d(ji) * r1_Dt_ice
            wfx_snw_sni_1d(ji) = wfx_snw_sni_1d(ji) + mass_si(ji) * a_i_1d(ji) * r1_Dt_ice
            ! update thickness
            zh_i(ji,0)  = zh_i(ji,0) + dh_snowice(ji)

            ! update heat content (J.m-2) and layer thickness
            h_i_old (ji,0) = h_i_old (ji,0) + dh_snowice(ji)
            eh_i_old(ji,0) = eh_i_old(ji,0) + zfmdt * zEw           ! 1st part (sea water enthalpy)
         END DO
         
         ! Now, the loop over the snow layers is done in snwthd_iceconv
         DO ji = 1, npti
            eh_i_old(ji,0) = eh_i_old(ji,0) + enthalpy_si(ji) 
         END DO
!



!         
!!!!$         ! --- Update snow diags --- !
!!!!$         !!clem: this is wrong. dh_s_tot is not used anyway
!!!!$         DO ji = 1, npti
!!!!$            dh_s_tot(ji) = dh_s_tot(ji) + dh_s_mlt(ji) + zdeltah(ji) + zdh_s_sub(ji) - dh_snowice(ji)
!!!!$         END DO
!!         !
!!         !

!#if defined key_isbaes
!      END IF
!#endif
         ! Note: remapping of ice enthalpy is done in icethd.F90

      ! --- ensure that a_i = 0 & h_s = 0 where h_i = 0 ---
      WHERE( h_i_1d(1:npti) == 0._wp )
         a_i_1d (1:npti) = 0._wp
         h_s_1d (1:npti) = 0._wp
         t_su_1d(1:npti) = rt0
      END WHERE
#if defined key_isbaes
      IF(ln_isbaes) THEN
         DO jk = 1, nlay_s
            WHERE( h_i_1d(1:npti) == 0._wp )
               dh_s_1d (1:npti,jk) = 0._wp ! ln_isbaes
            END WHERE
         END DO
      ENDIF
#endif
   END SUBROUTINE ice_thd_dh

#else
   !!----------------------------------------------------------------------
   !!   Default option                                NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icethd_dh

MODULE snwthd_iceconv
   !!======================================================================
   !!                       ***  MODULE snwthd_iceconv ***
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
   !!   snw_thd_iceconv        : vertical sea-ice growth and melt
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

#if defined key_isbaes   
   USE MODE_SNOW3L    ! For isbaes
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   snw_thd_iceconv        ! called by ice_thd

   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: snwthd_iceconv.F90 14686 2021-04-08 15:36:01Z clem $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE snw_thd_iceconv( isnow, zh_s, ze_s, thickness_si, mass_si, enthalpy_si )
      !!------------------------------------------------------------------
      !!                ***  ROUTINE snw_thd_iceconv  ***
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
      REAL(wp), DIMENSION(jpij,0:nlay_s  ), INTENT(inout) ::   zh_s      ! snw layer thickness (m) 
      REAL(wp), DIMENSION(jpij,0:nlay_s  ), INTENT(inout) ::   ze_s      ! snw layer enthalpy (J.m-3)

      REAL(wp), DIMENSION(jpij), INTENT(inout) ::   thickness_si  ! Thickness removed to snow by snow to ice conversion (m) 
      REAL(wp), DIMENSION(jpij), INTENT(inout) ::   mass_si       ! Mass removed to snow by snow to ice conversion (kg)      
      REAL(wp), DIMENSION(jpij), INTENT(inout) ::   enthalpy_si   ! Enthalpy removed to snow by snow to ice conversion (J/m2) 


      
      INTEGER  ::   ji, jk       ! dummy loop indices

      REAL(wp) ::   zfrac, zden, zdh, zdum      ! local scalar

      REAL(wp), DIMENSION(jpij)          ::   zhi, zdeltah
      REAL(wp), DIMENSION(jpij,nlay_s)   ::   zm_s, mass_snow, zrhos
      REAL(wp), DIMENSION(jpij,nlay_s)   ::   dh_sni
      !!------------------------------------------------------------------

! Snow to ice conversion      
! Initialize fields

      
#if defined key_isbaes      
     IF( ln_snwext .OR. ln_isbaes) THEN
#else
     IF( ln_snwext) THEN
#endif
         ! initialize snw layer thicknesses and enthalpies
         zh_s(1:npti,0) = 0._wp
         ze_s(1:npti,0) = 0._wp
         DO jk = 1, nlay_s
            DO ji = 1, npti
#if defined key_isbaes            
               zh_s(ji,jk) = dh_s_1d(ji,jk)
               ze_s(ji,jk) = e_s_1d(ji,jk)
#else
               zh_s(ji,jk) = h_s_1d(ji) * r1_nlay_s
               ze_s(ji,jk) = e_s_1d(ji,jk)
#endif
            END DO
         END DO
      ENDIF

      DO ji = 1, npti      
         DO jk = nlay_s, 1, -1
#if defined key_isbaes
            mass_snow(ji,jk) = rho_s_1d(ji,jk) * dh_s_1d(ji,jk)
            zrhos(ji,jk) = rho_s_1d(ji,jk)
#else
            mass_snow(ji,jk) = rhos * h_s_1d(ji) * r1_nlay_s !* a_i_1d(ji)
            zrhos(ji,jk) = rhos
#endif
            dh_sni(ji,jk) = 0._wp
         END DO
      END DO

      DO ji = 1, npti
      !   mass_si(ji) = 0; enthalpy_si(ji) = 0; thickness_si(ji) = 0; ! fields sent to sea ice model
         zhi(ji) = h_i_1d(ji); ! pseudo ice thickness,
      END DO

      DO ji = 1, npti
         zdh = 0._wp
         IF((SUM(zh_s(ji,:)) .ne. 0._wp) .OR. (h_i_1d(ji) .ne. 0._wp)) THEN     
            ! Loop over snow layers
            DO jk = nlay_s, 1, -1
               IF((zh_s(ji,jk)) >  0._wp) THEN
                   
                  !!! Calculate snow ice formation for layer k
                  zden = zrhos(ji,jk) + rho0 - rhoi ! denominator
                  zdh  = 0. !MAX( SUM(mass_snow(ji,:)) + ( rhoi - rho0 ) * zhi(ji) , 0._wp ) / zden ! snow ice for layer jk
                  dh_sni(ji,jk) = MIN( zdh, zh_s(ji,jk) ) ! snow ice cannot exceed available thickness
                  zfrac = dh_sni(ji,jk) / zh_s(ji,jk)       ! fraction of layer k consumed by snow ice (between 0 and 1)

                  !!! Cumulate over snow layers the mass, thickness and enthalpy of snow contributing to snow ice
                  thickness_si(ji) = thickness_si(ji) + dh_sni(ji,jk)

                  mass_si(ji)      = mass_si(ji)      + zfrac * mass_snow(ji,jk)
                  enthalpy_si(ji)  = enthalpy_si(ji)  + dh_sni(ji,jk) * e_s_1d(ji,jk)    

                  !!! Update thickness, mass and enthalpy of layer
                  zh_s(ji,jk)     = zh_s(ji,jk) - dh_sni(ji,jk)
!                  h_s_1d(ji) = h_s_1d(ji) - dh_sni(ji,jk) 
                ENDIF
            END DO
         END IF
      END DO

      DO ji = 1, npti
#if defined key_isbaes      
         dh_s_1d(ji,1:nlay_s) = zh_s(ji,1:nlay_s)
         h_s_1d(ji) = SUM(dh_s_1d(ji,1:nlay_s))
         !e_s_1d(ji,1:nlay_s) = ze_s(ji,1:nlay_s)
         !rhov_s_1d(ji,1:nlay_s) = mass_snow(ji,1:nlay_s)
#else            
         h_s_1d(ji) = SUM(zh_s(ji,:))
         dh_s_1d(ji,1:nlay_s) = h_s_1d(ji) * r1_nlay_s
#endif         
      END DO

#if defined key_isbaes      

#else
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
#endif      
!      IF(h_s_1d(1) .ne. h_s_1d(9) ) STOP

   END SUBROUTINE snw_thd_iceconv

#else
   !!----------------------------------------------------------------------
   !!   Default option                                NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE snwthd_iceconv

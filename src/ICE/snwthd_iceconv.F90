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

      REAL(wp), DIMENSION(jpij), INTENT(out) ::   thickness_si  ! Thickness removed to snow by snow to ice conversion (m) 
      REAL(wp), DIMENSION(jpij), INTENT(out) ::   mass_si       ! Mass removed to snow by snow to ice conversion (kg)      
      REAL(wp), DIMENSION(jpij), INTENT(out) ::   enthalpy_si   ! Enthalpy removed to snow by snow to ice conversion (J/m2) 


      
      INTEGER  ::   ji, jk       ! dummy loop indices

      REAL(wp) ::   zfrac, zden, zdh, zdum      ! local scalar

      REAL(wp), DIMENSION(jpij)          ::   zhi, zdeltah
      REAL(wp), DIMENSION(jpij,nlay_s)   ::   zm_s, mass_snow, zrhos
      REAL(wp), DIMENSION(jpij,nlay_s)   ::   dh_snowice
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
               IF(ln_isbaes) THEN
                  zh_s(ji,jk) = dh_s_1d(ji,jk)
                  ze_s(ji,jk) = e_s_1d(ji,jk)
               ELSE
                  zh_s(ji,jk) = h_s_1d(ji) * r1_nlay_s
                  ze_s(ji,jk) = e_s_1d(ji,jk)
               ENDIF
#else
               zh_s(ji,jk) = h_s_1d(ji) * r1_nlay_s
               ze_s(ji,jk) = e_s_1d(ji,jk)
#endif
            END DO
         END DO
      ENDIF

      DO ji = 1, npti      
#if defined key_isbaes
         DO jk = nlay_s, 1, -1
            IF(ln_isbaes) THEN
               mass_snow(ji,jk) = rho_s_1d(ji,jk) * dh_s_1d(ji,jk)
               zrhos(ji,jk) = rho_s_1d(ji,jk)
            ELSE
               mass_snow(ji,jk) = rhos * h_s_1d(ji) * r1_nlay_s !* a_i_1d(ji)
               zrhos(ji,jk) = rhos
            ENDIF
               zm_s(ji,jk) = mass_snow(ji,jk)                ! pseudo snow mass
         END DO
#else
         DO jk = nlay_s, 1, -1
            mass_snow(ji,jk) = rhos * h_s_1d(ji) * r1_nlay_s !* a_i_1d(ji)
            zrhos(ji,jk) = rhos
            zm_s(ji,jk) = mass_snow(ji,jk)                ! pseudo snow mass
         END DO
#endif
      END DO

#if defined key_isbaes
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
#else
      zdeltah(1:npti) = h_s_1d(1:npti) + h_i_1d(1:npti) * (rhoi-rho0) * r1_rhos
#endif
      DO jk = 0, nlay_s
         DO ji = 1, npti
            IF( zdeltah(ji) > 0._wp .AND. sst_1d(ji) > 0._wp ) THEN
               ! snow layer thickness that falls into the ocean
               zdum = MIN( zdeltah(ji) , zh_s(ji,jk) )
               ! mass & energy loss to the ocean
               hfx_res_1d(ji) = hfx_res_1d(ji) - ze_s(ji,jk) * zdum * a_i_1d(ji) * r1_Dt_ice  ! heat flux to the ocean [W.m-2], < 0

#if defined key_isbaes
               IF(ln_isbaes) THEN
                  wfx_res_1d(ji) = wfx_res_1d(ji) + rho_s_1d(ji,jk)        * zdum * a_i_1d(ji) * r1_Dt_ice  ! mass flux
                  ! update thickness and energy
                  e_s_1d(ji,jk)= MAX( 0._wp, e_s_1d(ji,jk) - (zdum / dh_s_1d(ji,jk))  * e_s_1d(ji,jk))
               ELSE
                  wfx_res_1d(ji) = wfx_res_1d(ji) + rhos        * zdum * a_i_1d(ji) * r1_Dt_ice  ! mass flux
               ENDIF
#else
               wfx_res_1d(ji) = wfx_res_1d(ji) + rhos        * zdum * a_i_1d(ji) * r1_Dt_ice  ! mass flux
#endif
               h_s_1d(ji)    = MAX( 0._wp, h_s_1d(ji)  - zdum )
               zh_s  (ji,jk) = MAX( 0._wp, zh_s(ji,jk) - zdum )
               ! update snow thickness that still has to fall
               zdeltah(ji)   = MAX( 0._wp, zdeltah(ji) - zdum )
            ENDIF
         END DO
      END DO

      DO ji = 1, npti
         mass_si(ji) = 0; enthalpy_si(ji) = 0; thickness_si(ji) = 0; ! fields sent to sea ice model
         zhi = h_i_1d(ji); ! pseudo ice thickness, 
      END DO

      DO ji = 1, npti
         zdh = 0._wp
         IF((SUM(zh_s(ji,:)) .ne. 0._wp) .OR. (h_i_1d(ji) .ne. 0._wp)) THEN     
            ! Loop over snow layers
            DO jk = nlay_s, 1, -1
               IF((zh_s(ji,jk)) >  0._wp) THEN
                   
                  !!! Calculate snow ice formation for layer k
                  zden = zrhos(ji,jk) + rho0 - rhoi ! denominator
                  zdh  = MAX( SUM(mass_snow(ji,:)) + ( rhoi - rho0 ) * zhi(ji) , 0._wp ) / zden  ! snow ice for layer jk
                  dh_snowice(ji,jk) = MIN( zdh, zh_s(ji,jk) ) ! snow ice cannot exceed available thickness
                  zfrac = dh_snowice(ji,jk) / zh_s(ji,jk)       ! fraction of layer k consumed by snow ice (between 0 and 1)

                  !!! Cumulate over snow layers the mass, thickness and enthalpy of snow contributing to snow ice
                  thickness_si(ji) = thickness_si(ji) + dh_snowice(ji,jk)

                  mass_si(ji)      = mass_si(ji)      + zfrac * mass_snow(ji,jk)
#if defined key_isbaes                  
                  IF(ln_isbaes) THEN
                      enthalpy_si(ji)  = enthalpy_si(ji)  + dh_snowice(ji,jk) * e_s_1d(ji,jk) 
                      !ze_s(ji,jk) = ze_s(ji,jk) - dh_snowice(ji,jk) * e_s_1d(ji,jk) 
                      mass_snow(ji,jk)     = mass_snow(ji,jk)     * ( 1.0 - zfrac ) 
                  ELSE    
                      enthalpy_si(ji)  = enthalpy_si(ji)  + dh_snowice(ji,jk) * e_s_1d(ji,jk)    
                  ENDIF
#endif

                  !!! Update pseudo thickness and snow mass for layer jk (needed to converge to zero-freeboard)
                  zhi(ji) = zhi(ji) + dh_snowice(ji,jk)
                  zm_s(ji,jk) = zm_s(ji,jk) - zfrac * mass_snow(ji,jk)

                  !!! Update thickness, mass and enthalpy of layer
                  zh_s(ji,jk)     = zh_s(ji,jk) - dh_snowice(ji,jk)
                  !mass_snow(ji,jk)     = mass_snow(ji,jk)     * ( 1.0 - zfrac )
                ENDIF
            END DO
         END IF
      END DO

      DO ji = 1, npti
#if defined key_isbaes      
         IF(ln_isbaes) THEN
            dh_s_1d(ji,1:nlay_s) = zh_s(ji,1:nlay_s)
            h_s_1d(ji) = SUM(dh_s_1d(ji,:))
            e_s_1d(ji,1:nlay_s) = ze_s(ji,1:nlay_s)
            !rhov_s_1d(ji,1:nlay_s) = mass_snow(ji,1:nlay_s)
         ELSE
#endif            
         h_s_1d(ji) = SUM(zh_s(ji,:))
         dh_s_1d(ji,1:nlay_s) = h_s_1d(ji) * r1_nlay_s
#if defined key_isbaes      
         ENDIF
#endif         
      END DO

            
      IF(.NOT. ln_isbaes) THEN 
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
      
!      IF(h_s_1d(1) .ne. h_s_1d(9) ) STOP

   END SUBROUTINE snw_thd_iceconv

#else
   !!----------------------------------------------------------------------
   !!   Default option                                NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE snwthd_iceconv

MODULE snwthd_snwfl
   !!======================================================================
   !!                       ***  MODULE snwthd_snwfl ***
   !!   seaice : thicknesses changes due to snowfall 
   !!======================================================================
   !! History :       !  2003-05  (M. Vancoppenolle) Original code in 1D
   !!                 !  2005-06  (M. Vancoppenolle) 3D version
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   snw_thd_snwfl        : vertical sea-ice growth and melt
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

   PUBLIC   snw_thd_snwfl        ! called by snw_thd

   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: snwthd_snwfl.F90 Theo $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE snw_thd_snwfl
      !!------------------------------------------------------------------
      !!                ***  ROUTINE snw_thd_snwfl  ***
      !!
      !! ** Purpose :   Snow thickness changes due to snowfall
      !!
      !! ** Method  :   Snow thickness increase by precipitation 
      !!
      !!                - Compute change in snow height due to snowfall 
      !!                - Snow enthalpy is remapped at the end of the routine 
      !! ** Notes   :   This is an extraction of the part of the code related to
      !                 snowfall in the icethd_dh routine of NEMO4.2-stable.
      !! References : Bitz and Lipscomb, 1999, J. Geophys. Res.
      !!              Fichefet T. and M. Maqueda 1997, J. Geophys. Res., 102(C6), 12609-12646
      !!              Vancoppenolle, Fichefet and Bitz, 2005, Geophys. Res. Let.
      !!              Vancoppenolle et al.,2009, Ocean Modelling
      !!------------------------------------------------------------------
      REAL(wp), DIMENSION(jpij,0:nlay_s  ) ::   zh_s      ! snw layer thickness (m) 
      REAL(wp), DIMENSION(jpij,0:nlay_s  ) ::   ze_s      ! snw layer enthalpy (J.m-3)

!
      INTEGER  ::   ji, jk       ! dummy loop indices
      INTEGER  ::   iter         ! local integer

      REAL(wp) ::   zdum

      REAL(wp), DIMENSION(jpij) ::   zsnw        ! distribution of snow after wind blowing

      !!------------------------------------------------------------------
      ! Initialise remaining heat and mass fluxes after melt and sublimation
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

! THEO Nb: remapping of the snow enthalpy is done later on in icethd_dh      
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
   END SUBROUTINE snw_thd_snwfl

#else
   !!----------------------------------------------------------------------
   !!   Default option                                NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE snwthd_snwfl

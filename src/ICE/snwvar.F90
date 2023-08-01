MODULE snwvar
   !!======================================================================
   !!                       ***  MODULE snwvar ***
   !!   sea-ice:  series of functions to transform or compute ice variables
   !!======================================================================
   !! History :   -   !  2006-01  (M. Vancoppenolle) Original code
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!
   !!                 There are three sets of variables
   !!                 VGLO : global variables of the model
   !!                        - v_i (jpi,jpj,jpl)
   !!                        - v_s (jpi,jpj,jpl)
   !!                        - a_i (jpi,jpj,jpl)
   !!                        - t_s (jpi,jpj,jpl)
   !!                        - e_i (jpi,jpj,nlay_i,jpl)
   !!                        - e_s (jpi,jpj,nlay_s,jpl)
   !!                        - sv_i(jpi,jpj,jpl)
   !!                        - oa_i(jpi,jpj,jpl)
   !!                 VEQV : equivalent variables sometimes used in the model
   !!                        - h_i(jpi,jpj,jpl)
   !!                        - h_s(jpi,jpj,jpl)
   !!                        - t_i(jpi,jpj,nlay_i,jpl)
   !!                        ...
   !!                 VAGG : aggregate variables, averaged/summed over all
   !!                        thickness categories
   !!                        - vt_i(jpi,jpj)
   !!                        - vt_s(jpi,jpj)
   !!                        - at_i(jpi,jpj)
   !!                        - st_i(jpi,jpj)
   !!                        - et_s(jpi,jpj)  total snow heat content
   !!                        - et_i(jpi,jpj)  total ice thermal content
   !!                        - sm_i(jpi,jpj)  mean ice salinity
   !!                        - tm_i(jpi,jpj)  mean ice temperature
   !!                        - tm_s(jpi,jpj)  mean snw temperature
   !!----------------------------------------------------------------------
   !!   snw_var_enthalpy  : compute ice and snow enthalpies from temperature
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants (ocean directory)
   USE sbc_oce , ONLY : sss_m, ln_ice_embd, nn_fsbc
   USE ice            ! sea-ice: variables
   USE ice1D          ! sea-ice: thermodynamics variables
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   snw_var_enthalpy

   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: snwvar.F90 15385 2021-10-15 13:52:48Z clem $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE snw_var_enthalpy
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE snw_var_enthalpy ***
      !!
      !! ** Purpose :   Computes sea ice energy of melting q_i (J.m-3) from temperature
      !!
      !! ** Method  :   Formula (Bitz and Lipscomb, 1999)
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jk   ! dummy loop indices
      REAL(wp) ::   ztmelts  ! local scalar
      !!-------------------------------------------------------------------
      !
      DO jk = 1, nlay_s             ! Snow energy of melting
         DO ji = 1, npti
            IF(ln_isbaes) THEN
               e_s_1d(ji,jk) = rho_s_1d(ji,jk) * ( rcpi * ( rt0 - t_s_1d(ji,jk) ) + rLfus )
            ELSE        
               e_s_1d(ji,jk) = rhos * ( rcpi * ( rt0 - t_s_1d(ji,jk) ) + rLfus )
            ENDIF
         END DO
      END DO
      !
   END SUBROUTINE snw_var_enthalpy

#else
   !!----------------------------------------------------------------------
   !!   Default option         Dummy module           NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE snwvar

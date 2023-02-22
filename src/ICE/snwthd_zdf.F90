MODULE snwthd_zdf
   !!======================================================================
   !!                  ***  MODULE icethd   ***
   !!   sea-ice : master routine for thermodynamics
   !!======================================================================
   !! History :  1.0  !  2000-01  (M.A. Morales Maqueda, H. Goosse, T. Fichefet) original code 1D
   !!            4.0  !  2018     (many people)       SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_thd       : thermodynamics of sea ice
   !!   ice_thd_init  : initialisation of sea-ice thermodynamics
   !!----------------------------------------------------------------------
   USE phycst         ! physical constants
   USE dom_oce        ! ocean space and time domain variables
   USE ice            ! sea-ice: variables
   USE ice1D          ! sea-ice: thermodynamics variables
   USE icevar         ! sea-ice: operations
   USE icectl         ! sea-ice: control print

   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   snw_thd_zdf         ! called by snwthd module

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: snwthd.F90 15388 2023-02-17 11:33:47Z Theo $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE snw_thd_zdf( zradtr_s, zradab_s, za_s_fra )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE snw_thd  ***
      !!
      !! ** Purpose : computes the time evolution of snow temperature profiles,
      !!              using the original Bitz and Lipscomb (1999) algorithm
      !!
      !! ** Method : solves the heat equation diffusion in the snow with a Neumann
      !!             boundary condition. The numerical scheme is an iterative Crank-
      !!             Nicolson on a non-uniform multilayer grid in the snow system. 
      !!             
      !! ** Action : Not finished yet            
      !!             
      !!            
      !!             
      !!             
      !!             
      !!             
      !!-------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(jpij,0:nlay_s), INTENT(out) ::   zradtr_s  ! Radiation transmited through the snow
      REAL(wp), DIMENSION(jpij,0:nlay_s), INTENT(out) ::   zradab_s  ! Radiation absorbed in the snow 
      REAL(wp), DIMENSION(jpij), INTENT(out) ::   za_s_fra    ! ice fraction covered by snow
      !
      INTEGER  :: ji, jj, jk, jl   ! dummy loop indices
      !
      REAL(wp) ::   zhs_ssl   =  0.03_wp      ! surface scattering layer in the snow
      REAL(wp) ::   zh_min    =  1.e-3_wp     ! minimum ice/snow thickness for conduction
      !
      REAL(wp), DIMENSION(jpij) ::   zraext_s     ! extinction coefficient of radiation in the snow
      REAL(wp), DIMENSION(jpij)          ::   isnow       ! snow presence (1) or not (0)
      REAL(wp), DIMENSION(jpij)          ::   isnow_comb  ! snow presence for met-office
      REAL(wp), DIMENSION(jpij) ::   zh_s, z1_h_s ! snow layer thickness
      REAL(wp), DIMENSION(jpij,nlay_s)   ::   ztsold      ! Old temperature in the snow

      ! calculate ice fraction covered by snow for radiation
      CALL ice_var_snwfra( h_s_1d(1:npti), za_s_fra(1:npti) )

      !------------------
      ! 1) Initialization
      !------------------
      !
      ! extinction radiation in the snow
      IF    ( nn_qtrice == 0 ) THEN   ! constant
         zraext_s(1:npti) = rn_kappa_s
      ELSEIF( nn_qtrice == 1 ) THEN   ! depends on melting/freezing conditions
         WHERE( t_su_1d(1:npti) < rt0 )   ;   zraext_s(1:npti) = rn_kappa_sdry   ! no surface melting
         ELSEWHERE                        ;   zraext_s(1:npti) = rn_kappa_smlt   !    surface melting
         END WHERE
      ENDIF

      ! thicknesses
      DO ji = 1, npti
         ! ice thickness
         ! snow thickness
         IF( h_s_1d(ji) > 0._wp ) THEN
            zh_s  (ji) = MAX( zh_min , h_s_1d(ji) ) * r1_nlay_s ! set a minimum thickness for conduction
            z1_h_s(ji) = 1._wp / zh_s(ji)                       !       it must be very small
            isnow (ji) = 1._wp
         ELSE
            zh_s  (ji) = 0._wp
            z1_h_s(ji) = 0._wp
            isnow (ji) = 0._wp
         ENDIF
         ! for Met-Office
         IF( h_s_1d(ji) < zh_min ) THEN
            isnow_comb(ji) = h_s_1d(ji) / zh_min
         ELSE
            isnow_comb(ji) = 1._wp
         ENDIF
      END DO

      ztsold (1:npti,:) = t_s_1d(1:npti,:)   ! Old snow temperature

      !-------------
      ! 2) Radiation
      !-------------
      ! --- Transmission/absorption of solar radiation in the ice --- !
      zradtr_s(1:npti,0) = qtr_ice_top_1d(1:npti)
      DO jk = 1, nlay_s
         DO ji = 1, npti
            !                             ! radiation transmitted below the layer-th snow layer
            zradtr_s(ji,jk) = zradtr_s(ji,0) * EXP( - zraext_s(ji) * MAX( 0._wp, zh_s(ji) * REAL(jk) - zhs_ssl ) )
            !                             ! radiation absorbed by the layer-th snow layer
            zradab_s(ji,jk) = zradtr_s(ji,jk-1) - zradtr_s(ji,jk)
         END DO
      END DO

      !
   END SUBROUTINE snw_thd_zdf
#else
   !!----------------------------------------------------------------------
   !!   Default option         Dummy module          NO  SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE snwthd_zdf

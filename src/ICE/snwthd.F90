MODULE snwthd
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
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants (ocean directory) 
   USE ice             ! sea-ice: variables
   USE ice1D          ! sea-ice: thermodynamics variables
   USE icectl         ! sea-ice: control print

   USE snwthd_zdf      ! snow: vertical diffusion 
   USE snwthd_dh       ! snow: growing / melting 
   !
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE lib_fortran     ! fortran utilities (glob_sum + no signed zero)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   snw_thd         ! called by ice_thd module
   PUBLIC   snw_thd_init    ! called by ice_thd_init

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: snwthd.F90 15388 2023-02-17 11:33:47Z Theo $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE snw_thd( zradtr_s, zradab_s, za_s_fra, zq_rema, zevap_rema ) 
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE snw_thd  ***
      !!
      !! ** Purpose : This routine manages snow thermodynamics in detached mode
      !!
      !! ** Action : 
      !!             
      !!             
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
      REAL(wp), DIMENSION(jpij), INTENT(out) ::   zq_rema     ! remaining heat flux from snow melting       (J.m-2)
      REAL(wp), DIMENSION(jpij), INTENT(out) ::   zevap_rema  ! remaining mass flux from snow sublimation   (kg.m-2)
      !
      !!-------------------------------------------------------------------
      ! controls
      IF( ln_timing    )   CALL timing_start('snwthd')                                                             ! timing
      IF( ln_icediachk )   CALL ice_cons_hsm(0, 'snwthd', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      IF( ln_icediachk )   CALL ice_cons2D  (0, 'snwthd',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft) ! conservation

      e_s_1d_old(:,:) = e_s_1d(:,:)    
      !------------------
      ! 1) Thermodynamics 
      !------------------

      CALL snw_thd_zdf( zradtr_s, zradab_s, za_s_fra )

      !------------------
      ! 2) Snowfall / melt 
      !------------------

      IF( ln_icedH )   CALL snw_thd_dh( zq_rema, zevap_rema )
      !
      !
      IF( ln_icediachk )   CALL ice_cons_hsm(1, 'snwthd', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
      IF( ln_icediachk )   CALL ice_cons2D  (1, 'snwthd',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft)
      !
      IF( ln_timing )   CALL timing_stop('snwthd')                                        ! timing
      !
   END SUBROUTINE snw_thd

   SUBROUTINE snw_thd_init
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE ice_thd_init ***
      !!
      !! ** Purpose :   Physical constants and parameters associated with
      !!                snow thermodynamics (detached mode)
      !!
      !! ** Method  :   Read the namthd namelist and check the parameters
      !!                called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namthd
      !!-------------------------------------------------------------------
      INTEGER  ::   ios   ! Local integer output status for namelist read
      !!
      !
      NAMELIST/namthd_snw/ rn_cnd_s, rn_kappa_s, rn_kappa_smlt, rn_kappa_sdry
      !!-------------------------------------------------------------------
      !
      READ  ( numnam_ice_ref, namthd_snw, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namthd_snw in reference namelist' )
      READ  ( numnam_ice_cfg, namthd_snw, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namthd_snw in configuration namelist' )
      IF(lwm) WRITE( numoni, namthd_snw )
      !
      IF(lwp) THEN                          ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'snw_thd: Snow parameters'
         WRITE(numout,*) '~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namthd_snw:'
         WRITE(numout,*) '      thermal conductivity in the snow                          rn_cnd_s      = ', rn_cnd_s
         WRITE(numout,*) '      extinction radiation parameter in snw      (nn_qtrice=0)  rn_kappa_s    = ', rn_kappa_s
         WRITE(numout,*) '      extinction radiation parameter in melt snw (nn_qtrice=1)  rn_kappa_smlt = ', rn_kappa_smlt
         WRITE(numout,*) '      extinction radiation parameter in dry  snw (nn_qtrice=1)  rn_kappa_sdry = ', rn_kappa_sdry
      ENDIF


      !
      !
   END SUBROUTINE snw_thd_init

#else
   !!----------------------------------------------------------------------
   !!   Default option         Dummy module          NO  SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE snwthd

MODULE snwthd
   !!======================================================================
   !!                  ***  MODULE icethd   ***
   !!   sea-ice : master routine for snow thermodynamics
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
   USE snwthd_dh       ! snow: melting 
   USE snwthd_snwfl      ! snowfall

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

   SUBROUTINE snw_thd( zradtr_s, zradab_s, za_s_fra, qcn_snw_bot_1d, isnow, & 
                       zq_rema, zevap_rema, zh_s, ze_s ) 
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE snw_thd  ***
      !!
      !! ** Purpose : This routine manages snow thermodynamics in detached mode
      !!
      !! ** Action : - call snw_thd_zdf ( vertical diffusion of heat in snow)
      !!             - call snw_thd_dh ( changes in heat and snow thickness with 
      !!               snowfall / melt / sublimation 
      !! ** Returns  - The radiative fluxes transmitted and absorbed through snow
      !!               (zradtr_s & zradab_s)
      !!             - The ice fraction covered by snow (za_s_fra)
      !!             - The conduction flux at snow / ice interface (qcn_snw_bot_1d)
      !!             - The snow presence (1) or not (0) at the beggining of the 
      !!               procedure (isnow)
      !!             - the remaining heat and mass fluxes after changes in snow 
      !!               height due to snowfall/melt/sublimation (zq_rema & zevap_rema)
      !!             - New thicknesses and enthalpy profiles to be remapped later on
      !!               in icethd_dh (zh_s & ze_s)
      !!                          
      !! ** Global variables altered after this routine (but not passed as arguments):
      !!             - e_s_1d, t_s_1d, h_s_1d, t_su_1d (surface T°), t_si_1d (snw / SI interfacial T°)     
      !!             - Heat & mass fluxes:  qml_ice_1d (heat remaining for melting the surface), 
      !!               hfx_res_1d (HF to the ocean after snow internal melting), wfx_snw_sum_1d (mass
      !!               flux to the ocean after melt), hfx_spr_1d (HF from precip), wfx_spr_1d (MF from
      !!               precip), hfx_snw_1d (HF from surface melting), wfx_snw_sum_1d (MF from surface 
      !!               melting), hfx_sub_1d (HF from snw sublimation), wfx_snw_sub_1d (MF from sub.)
      !!                             
      !!             
      !!             
      !!             
      !! ** Note:    This routine outputs zh_s, which is the thicknesses of the snow layer, 
      !!              which are all equal to the total snow thickness / nlay_s. We keep it 
      !!              that way so that it is consistent with ISBA-ES coupling.     
      !!             
      !!-------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(jpij,0:nlay_s), INTENT(out) ::   zradtr_s  ! Radiation transmited through the snow
      REAL(wp), DIMENSION(jpij,0:nlay_s), INTENT(out) ::   zradab_s  ! Radiation absorbed in the snow 
      REAL(wp), DIMENSION(jpij), INTENT(out) ::   za_s_fra    ! ice fraction covered by snow
      REAL(wp), DIMENSION(jpij), INTENT(out) ::   qcn_snw_bot_1d    ! Conduction flux at snow / ice interface
      REAL(wp), DIMENSION(jpij), INTENT(out) ::   isnow       ! snow presence (1) or not (0)
      REAL(wp), DIMENSION(jpij), INTENT(out) ::   zq_rema     ! remaining heat flux from snow melting       (J.m-2)
      REAL(wp), DIMENSION(jpij), INTENT(out) ::   zevap_rema  ! remaining mass flux from snow sublimation   (kg.m-2)
      REAL(wp), DIMENSION(jpij,0:nlay_s), INTENT(out) ::   zh_s  ! Thicknesses of the snow layers (m)
      REAL(wp), DIMENSION(jpij,0:nlay_s), INTENT(out) ::   ze_s  ! Snow enthalpy per unit volume of the snow layers  

      !
      !!-------------------------------------------------------------------
      ! controls
      IF( ln_timing    )   CALL timing_start('snwthd')                                                             ! timing

      !------------------
      ! 1) Snowfall 
      !------------------

      ! We compute the changes in snow thickness due to snowfall first. This
      ! avoid unconsistencies between the surface temperature used to compute
      ! the available heat for surface melt (in snwthd_dh). This way the
      ! available heat for surface melting (computed in snwthd_dh) is always
      ! computed from the surface ice or snow T° at time=t+1, regardless of the
      ! presence of snow or not at time=t.
      ! 
      IF( ln_icedH )   CALL snw_thd_snwfl ! Change in snow thickness and enthalpy due to snowfall 

      !------------------
      ! 2) Thermodynamics 
      !------------------

      IF( .NOT.ln_cndflx ) THEN                           ! No conduction flux ==> default option
         CALL snw_thd_zdf( np_cnd_OFF, zradtr_s, zradab_s, za_s_fra, qcn_snw_bot_1d, isnow)
      ELSEIF( ln_cndflx .AND. .NOT.ln_cndemulate ) THEN   ! Conduction flux as surface boundary condition ==> Met Office default option
         CALL snw_thd_zdf( np_cnd_ON, zradtr_s, zradab_s, za_s_fra, qcn_snw_bot_1d, isnow )
      ELSEIF( ln_cndflx .AND.      ln_cndemulate ) THEN   ! Conduction flux is emulated 
         CALL snw_thd_zdf( np_cnd_EMU, zradtr_s, zradab_s, za_s_fra, qcn_snw_bot_1d, isnow )
         CALL snw_thd_zdf( np_cnd_ON, zradtr_s, zradab_s, za_s_fra, qcn_snw_bot_1d, isnow )
      ENDIF

      !------------------
      ! 3) Snow melt & sublimation 
      !------------------

      IF( ln_icedH )   CALL snw_thd_dh( zq_rema, zevap_rema, zh_s, ze_s)
      !
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

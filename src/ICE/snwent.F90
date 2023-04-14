MODULE snwent
   !!======================================================================
   !!                  ***  MODULE snwent   ***
   !!   sea-ice : remapping of snow enthalpy 
   !!======================================================================
   !! History :  1.0  !  2000-01  (M.A. Morales Maqueda, H. Goosse, T. Fichefet) original code 1D
   !!            4.0  !  2018     (many people)       SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   snw_ent       : snow enthalpy remapping 
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants (ocean directory) 
   USE ice             ! sea-ice: variables
   USE ice1D          ! sea-ice: thermodynamics variables

  ! !
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE lib_fortran     ! fortran utilities (glob_sum + no signed zero)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   snw_ent         ! called by ice_thd_dh, snwent_dh & snwent_snwfl module

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: snwent.F90 15388 2023-02-17 11:33:47Z Theo $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

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
       REAL(wp), DIMENSION(jpij,0:nlay_s), INTENT(in   ) ::   ph_old ! old thicknesses (m)
       REAL(wp), DIMENSION(jpij,0:nlay_s), INTENT(in   ) ::   pe_old ! old enthlapies (J.m-3)
       REAL(wp), DIMENSION(jpij,1:nlay_s), INTENT(inout) ::   pe_new ! new enthlapies (J.m-3, remapped)
       !
       INTEGER  :: ji         !  dummy loop indices
       INTEGER  :: jk0, jk1   !  old/new layer indices
       !
       REAL(wp), DIMENSION(jpij,0:nlay_s+1) ::   zeh_cum0, zh_cum0   ! old cumulative enthlapies and layers interfaces
       REAL(wp), DIMENSION(jpij,0:nlay_s)   ::   zeh_cum1, zh_cum1   ! new cumulative enthlapies and layers interfaces
       REAL(wp), DIMENSION(jpij)            ::   zhnew               ! new layers thicknesses
       !!-------------------------------------------------------------------
 
       !--------------------------------------------------------------------------
       !  1) Cumulative integral of old enthalpy * thickness and layers
       !  interfaces
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
   !!   Default option         Dummy module          NO  SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE snwent

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
   USE snwvar         ! sea-ice: operations
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

   SUBROUTINE snw_thd_zdf(k_cnd, zradtr_s, zradab_s, za_s_fra, qcn_snw_bot_1d, isnow)
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE snw_thd  ***
      !!
      !! ** Purpose : computes the time evolution of snow temperature profiles,
      !!              using the original Bitz and Lipscomb (1999) algorithm
      !!
      !! ** Method : - solves the heat equation diffusion in the snow with a Neumann
      !!             boundary condition at the top. The numerical scheme is an iterative Crank-
      !!             Nicolson on a non-uniform multilayer grid in the snow system. 
      !!             - the method is similar as the one used in icethd_zdf_bl99.F90 when
      !!             ln_snwext=false. The only difference is that the heat diffusion is computed
      !!             using the T° of the first ice level at time=t instead of t+1. 
      !!             
      !! ** Action : - uses a similar method as in icethd_zdf_bl99.F90 to resolve the heat eq.
      !!             in the snow
      !!             - returns the radiation transmitted and absorbed through the snow, the
      !!             ice fraction covered by snow and the conduction flux at snow / ice 
      !!             interface: Fc,si = - Kappa_int * (T_i1 - T_s) (see eq. 3.61 in LIM3 book 
      !!             as an example)
      !!            
      !!-------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   k_cnd     ! conduction flux (off, on, emulated)
      REAL(wp), DIMENSION(jpij,0:nlay_s), INTENT(out) ::   zradtr_s  ! Radiation transmited through the snow
      REAL(wp), DIMENSION(jpij,0:nlay_s), INTENT(out) ::   zradab_s  ! Radiation absorbed in the snow 
      REAL(wp), DIMENSION(jpij), INTENT(out) ::   za_s_fra    ! ice fraction covered by snow
      REAL(wp), DIMENSION(jpij), INTENT(out) ::   qcn_snw_bot_1d ! Conduction flux at snow / ice interface
      REAL(wp), DIMENSION(jpij), INTENT(out) ::   isnow       ! snow presence (1) or not (0)
      !
      REAL(wp), DIMENSION(jpij,0:nlay_i) ::   zradtr_i    ! Radiation transmitted through the ice
      REAL(wp), DIMENSION(jpij,0:nlay_i) ::   zradab_i    ! Radiation absorbed in the ice

      INTEGER ::   ji, jk         ! spatial loop index
      INTEGER ::   jm             ! current reference number of equation
      INTEGER ::   jm_mint, jm_maxt
      INTEGER ::   iconv          ! number of iterations in iterative procedure
      INTEGER ::   iconv_max = 50 ! max number of iterations in iterative procedure
      !
      INTEGER, DIMENSION(jpij) ::   jm_min    ! reference number of top equation
      INTEGER, DIMENSION(jpij) ::   jm_max    ! reference number of bottom equation

      LOGICAL, DIMENSION(jpij) ::   l_T_converged   ! true when T converges (per grid point)

      !
      REAL(wp) ::   zg1s      =  2._wp        ! for the tridiagonal system
      REAL(wp) ::   zg1      =  2._wp        ! for the tridiagonal system

      REAL(wp) ::   zgamma    =  18009._wp    ! for specific heat
      REAL(wp) ::   zbeta     =  0.117_wp     ! for thermal conductivity (could be 0.13)
      REAL(wp) ::   zkimin    =  0.10_wp      ! minimum ice thermal conductivity
      REAL(wp) ::   ztsu_err  =  1.e-5_wp     ! range around which t_su is considered at 0C

      REAL(wp) ::   zdti_bnd  =  1.e-4_wp     ! maximal authorized error on temperature

      REAL(wp) ::   zhs_ssl   =  0.03_wp      ! surface scattering layer in the snow
      REAL(wp) ::   zh_min    =  1.e-3_wp     ! minimum ice/snow thickness for conduction

      REAL(wp) ::   zdti_max                  ! current maximal error on temperature
      REAL(wp) ::   zcpi                      ! Ice specific heat
      REAL(wp) ::   zhfx_err, zdq             ! diag errors on heat

      !
      REAL(wp), DIMENSION(jpij) ::   zraext_s     ! extinction coefficient of radiation in the snow
      REAL(wp), DIMENSION(jpij) ::   ztsub        ! surface temperature at previous iteration
      REAL(wp), DIMENSION(jpij) ::   zh_s, z1_h_s ! snow layer thickness
      REAL(wp), DIMENSION(jpij) ::   zqns_ice_b   ! solar radiation absorbed at the surface
      REAL(wp), DIMENSION(jpij) ::   zfnet        ! surface flux function
      REAL(wp), DIMENSION(jpij) ::   zdqns_ice_b  ! derivative of the surface flux function
      !
      REAL(wp), DIMENSION(jpij       )   ::   ztsuold     ! Old surface temperature in the ice

      REAL(wp), DIMENSION(jpij,nlay_s)   ::   ztsold      ! Old temperature in the snow
      REAL(wp), DIMENSION(jpij,nlay_s)   ::   ztsb        ! Temporary temperature in the snow to check the convergence
      REAL(wp), DIMENSION(jpij,0:nlay_s) ::   zkappa_s    ! Kappa factor in the snow
      REAL(wp), DIMENSION(jpij,0:nlay_i) ::   zkappa_i    ! Kappa factor in the ice

      REAL(wp), DIMENSION(jpij,0:nlay_s) ::   zeta_s      ! Eta factor in the snow
      REAL(wp), DIMENSION(jpij) ::   zkappa_si    ! Kappa factor at snow / ice interface

      REAL(wp), DIMENSION(jpij)          ::   zkappa_comb ! Combined snow and ice surface conductivity
      REAL(wp), DIMENSION(jpij)          ::   zq_ini      ! diag errors on heat
      REAL(wp), DIMENSION(jpij)          ::   zghe        ! G(he), th. conduct enhancement factor, mono-cat
      
      ! Ice variables (1st layer only)
      REAL(wp), DIMENSION(jpij) ::   zh_i, z1_h_i ! ice layer thickness
      REAL(wp), DIMENSION(jpij,0:nlay_i) ::   ztcond_i    ! Ice thermal conductivity  
      REAL(wp), DIMENSION(jpij,0:nlay_i) ::   ztcond_i_cp ! copy
      ! Nb: Theo : we keep ztcond_i ztcond_i_cp 2D here bc it is needed if virtual_itd = T
      
      REAL(wp), DIMENSION(jpij)          ::   isnow_comb  ! snow presence for met-office
      REAL(wp), DIMENSION(jpij,nlay_i+nlay_s+1)   ::   zindterm    ! 'Ind'ependent term
      REAL(wp), DIMENSION(jpij,nlay_i+nlay_s+1)   ::   zindtbis    ! Temporary 'ind'ependent term
      REAL(wp), DIMENSION(jpij,nlay_i+nlay_s+1)   ::   zdiagbis    ! Temporary 'dia'gonal term
      REAL(wp), DIMENSION(jpij,nlay_i+nlay_s+1,3) ::   ztrid       ! Tridiagonal system terms

      ! Mono-category
      REAL(wp) ::   zepsilon   ! determines thres. above which computation of G(h) is done
      REAL(wp) ::   zhe        ! dummy factor
      REAL(wp) ::   zcnd_i     ! mean sea ice thermal conductivity
      !!------------------------------------------------------------------

      ! --- diag error on heat diffusion - PART 1 --- !
      DO ji = 1, npti
         zq_ini(ji) = ( SUM( e_i_1d(ji,1:nlay_i) ) * h_i_1d(ji) * r1_nlay_i  +  &
            &           SUM( e_s_1d(ji,1:nlay_s) ) * h_s_1d(ji) * r1_nlay_s )
      END DO

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
         IF( h_i_1d(ji) > 0._wp ) THEN
            zh_i  (ji) = MAX( zh_min , h_i_1d(ji) ) * r1_nlay_i ! set a minimum thickness for conduction
            z1_h_i(ji) = 1._wp / zh_i(ji)                       !       it must be very small
         ELSE
            zh_i  (ji) = 0._wp
            z1_h_i(ji) = 0._wp
         ENDIF

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

      ! Store initial temperatures and non solar heat fluxes
      IF( k_cnd == np_cnd_OFF .OR. k_cnd == np_cnd_EMU ) THEN
         ztsub      (1:npti) = t_su_1d(1:npti)                          ! surface temperature at iteration n-1
         ztsuold    (1:npti) = t_su_1d(1:npti)                          ! surface temperature initial value
       !  t_su_1d    (1:npti) = MIN( t_su_1d(1:npti), rt0 - ztsu_err )   ! required to leave the choice between melting or not
         WHERE(isnow(1:npti) == 1._wp) t_su_1d    (1:npti) = MIN( t_su_1d(1:npti), rt0 - ztsu_err )   ! required to leave the choice between melting or not
         zdqns_ice_b(1:npti) = dqns_ice_1d(1:npti)                      ! derivative of incoming nonsolar flux
         zqns_ice_b (1:npti) = qns_ice_1d(1:npti)                       ! store previous qns_ice_1d value
         !
      ENDIF
      
      ztsold (1:npti,:) = t_s_1d(1:npti,:)   ! Old snow temperature

      !-------------
      ! 2) Radiation
      !-------------
      ! --- Transmission/absorption of solar radiation in the snw --- !
      zradtr_s(1:npti,0) = qtr_ice_top_1d(1:npti)
      DO jk = 1, nlay_s
         DO ji = 1, npti
            !                             ! radiation transmitted below the layer-th snow layer
            zradtr_s(ji,jk) = zradtr_s(ji,0) * EXP( - zraext_s(ji) * MAX( 0._wp, zh_s(ji) * REAL(jk) - zhs_ssl ) )
            !                             ! radiation absorbed by the layer-th snow layer
            zradab_s(ji,jk) = zradtr_s(ji,jk-1) - zradtr_s(ji,jk)
         END DO
      END DO

      iconv = 0          ! number of iterations
      !
      l_T_converged(:) = .FALSE.
      ! Convergence calculated until all sub-domain grid points have converged
      ! Calculations keep going for all grid points until sub-domain convergence (vectorisation optimisation)
      ! but values are not taken into account (results independant of MPI partitioning)
      !
      !                                                                            !============================!
      DO WHILE ( ( .NOT. ALL (l_T_converged(1:npti)) ) .AND. iconv < iconv_max )   ! Iterative procedure begins !
         !                                                                         !============================!
         iconv = iconv + 1
         !
         ztsb(1:npti,:) = t_s_1d(1:npti,:)
         !
         !--------------------------------
         ! 3) Sea ice thermal conductivity of the first ice layer
         !--------------------------------
         IF( ln_cndi_U64 ) THEN         !-- Untersteiner (1964) formula: k = k0 + beta.S/T
            !
            DO ji = 1, npti
               ztcond_i_cp(ji,0)      = rcnd_i + zbeta * sz_i_1d(ji,1)      / MIN( -epsi10, t_i_1d(ji,1) - rt0 )
               ztcond_i_cp(ji,nlay_i) = rcnd_i + zbeta * sz_i_1d(ji,nlay_i) / MIN( -epsi10, t_bo_1d(ji)  - rt0 )
            END DO
            DO jk = 1, nlay_i-1
               DO ji = 1, npti
                  ztcond_i_cp(ji,jk) = rcnd_i + zbeta * 0.5_wp * ( sz_i_1d(ji,jk) + sz_i_1d(ji,jk+1) ) /  &
                     &                    MIN( -epsi10, 0.5_wp * (  t_i_1d(ji,jk) +  t_i_1d(ji,jk+1) ) - rt0 )
               END DO
            END DO
            !
         ELSEIF( ln_cndi_P07 ) THEN     !-- Pringle et al formula: k = k0 + beta1.S/T - beta2.T
            !
            DO ji = 1, npti
               ztcond_i_cp(ji,0)      = rcnd_i + 0.09_wp  *  sz_i_1d(ji,1)      / MIN( -epsi10, t_i_1d(ji,1) - rt0 )  &
                  &                            - 0.011_wp * ( t_i_1d(ji,1) - rt0 )
               ztcond_i_cp(ji,nlay_i) = rcnd_i + 0.09_wp  *  sz_i_1d(ji,nlay_i) / MIN( -epsi10, t_bo_1d(ji)  - rt0 )  &
                  &                            - 0.011_wp * ( t_bo_1d(ji) - rt0 )
            END DO
            DO jk = 1, nlay_i-1
               DO ji = 1, npti
                  ztcond_i_cp(ji,jk) = rcnd_i + 0.09_wp  *   0.5_wp * ( sz_i_1d(ji,jk) + sz_i_1d(ji,jk+1) ) /       &
                     &                         MIN( -epsi10, 0.5_wp * (  t_i_1d(ji,jk) +  t_i_1d(ji,jk+1) ) - rt0 ) &
                     &                        - 0.011_wp * ( 0.5_wp * (  t_i_1d(ji,jk) +  t_i_1d(ji,jk+1) ) - rt0 )
               END DO
            END DO
            !
         ENDIF

         ! Variable used after iterations
         ! Value must be frozen after convergence for MPP independance reason
         DO ji = 1, npti
            IF ( .NOT. l_T_converged(ji) ) &
               ztcond_i(ji,:) = MAX( zkimin, ztcond_i_cp(ji,:) )
         END DO

         !
         !--- G(he) : enhancement of thermal conductivity in mono-category case
         ! Computation of effective thermal conductivity G(h)
         ! Used in mono-category case only to simulate an ITD implicitly
         ! Fichefet and Morales Maqueda, JGR 1997
         zghe(1:npti) = 1._wp
         !
         IF( ln_virtual_itd ) THEN
            !
            zepsilon = 0.1_wp
            DO ji = 1, npti
               zcnd_i = SUM( ztcond_i(ji,:) ) / REAL( nlay_i+1, wp )                                ! Mean sea ice thermal conductivity
               zhe = ( rn_cnd_s * h_i_1d(ji) + zcnd_i * h_s_1d(ji) ) / ( rn_cnd_s + zcnd_i )        ! Effective thickness he (zhe)
               IF( zhe >=  zepsilon * 0.5_wp * EXP(1._wp) )  &
                  &   zghe(ji) = MIN( 2._wp, 0.5_wp * ( 1._wp + LOG( 2._wp * zhe / zepsilon ) ) )   ! G(he)
            END DO
            !
         ENDIF
         !
         !-----------------
         ! 4) kappa factors
         !-----------------
         !--- Snow
         ! Variable used after iterations
         ! Value must be frozen after convergence for MPP independance reason
         DO jk = 0, nlay_s-1
            DO ji = 1, npti
               IF ( .NOT. l_T_converged(ji) ) &
                  zkappa_s(ji,jk) = zghe(ji) * rn_cnd_s * z1_h_s(ji)
            END DO
         END DO
         DO ji = 1, npti   ! Snow-ice interface
            IF ( .NOT. l_T_converged(ji) ) &
               zkappa_s(ji,nlay_s) = isnow(ji) * zghe(ji) * rn_cnd_s * ztcond_i(ji,0) &
                  &                            / ( 0.5_wp * ( ztcond_i(ji,0) * zh_s(ji) + rn_cnd_s * zh_i(ji) ) )
         END DO
        !--- Ice
         ! Variable used after iterations
         ! Value must be frozen after convergence for MPP independance reason
         DO jk = 0, nlay_i
            DO ji = 1, npti
               IF ( .NOT. l_T_converged(ji) ) &
                  zkappa_i(ji,jk) = zghe(ji) * ztcond_i(ji,jk) * z1_h_i(ji)
            END DO
         END DO
         DO ji = 1, npti   ! Snow-ice interface
            IF ( .NOT. l_T_converged(ji) ) THEN
               ! If there is snow then use the same snow-ice interface conductivity for the top layer of ice
               IF( h_s_1d(ji) > 0._wp )   zkappa_i(ji,0) = zkappa_s(ji,nlay_s)
           ENDIF
         END DO

         !
         !--------------------------------------
         ! 5) Sea ice specific heat, eta factors
         !--------------------------------------

         DO jk = 1, nlay_s
            DO ji = 1, npti
               zeta_s(ji,jk) = rDt_ice * r1_rhos * r1_rcpi * z1_h_s(ji)
            END DO
         END DO
         !  
         !
         !----------------------------------------!
         !                                        !
         !   Conduction flux is off or emulated   !
         !                                        !
         !----------------------------------------!
         !
         IF( k_cnd == np_cnd_OFF .OR. k_cnd == np_cnd_EMU ) THEN
            !
            ! ==> The original BL99 temperature computation is used
            !       (with qsr_ice, qns_ice and dqns_ice as inputs)
            !
            !----------------------------
            ! 6) surface flux computation
            !----------------------------
            ! update of the non solar flux according to the update in T_su
            DO ji = 1, npti
               ! Variable used after iterations
               ! Value must be frozen after convergence for MPP independance reason
               IF ( ( .NOT. l_T_converged(ji) ) .AND. (isnow(ji) == 1 ) ) &
                  qns_ice_1d(ji) = qns_ice_1d(ji) + dqns_ice_1d(ji) * ( t_su_1d(ji) - ztsub(ji) )
            END DO

            DO ji = 1, npti
               zfnet(ji) = qsr_ice_1d(ji) - qtr_ice_top_1d(ji) + qns_ice_1d(ji) ! net heat flux = net - transmitted solar + non solar
            END DO
            !
            !----------------------------
            ! 7) tridiagonal system terms
            !----------------------------
            ! layer denotes the number of the layer in the snow or in the ice
            ! jm denotes the reference number of the equation in the tridiagonal
            ! system, terms of tridiagonal system are indexed as following :
            ! 1 is subdiagonal term, 2 is diagonal and 3 is superdiagonal one

            ! ice interior terms (top equation has the same form as the others)
            ztrid   (1:npti,:,:) = 0._wp
            zindterm(1:npti,:)   = 0._wp
            zindtbis(1:npti,:)   = 0._wp
            zdiagbis(1:npti,:)   = 0._wp

            DO ji = 1, npti
               ! snow interior terms (bottom equation has the same form as the others)
                                               !---------------------!
               IF( h_s_1d(ji) > 0._wp ) THEN   !  snow-covered cells !
                  !                            !---------------------!
                  DO jm = 3, nlay_s + 1
                     jk = jm - 1
                     ztrid   (ji,jm,1) =       - zeta_s(ji,jk) *   zkappa_s(ji,jk-1)
                     ztrid   (ji,jm,2) = 1._wp + zeta_s(ji,jk) * ( zkappa_s(ji,jk-1) + zkappa_s(ji,jk) )
                     ztrid   (ji,jm,3) =       - zeta_s(ji,jk) *                       zkappa_s(ji,jk)
                     zindterm(ji,jm)   = ztsold(ji,jk) + zeta_s(ji,jk) * zradab_s(ji,jk)
                  END DO

                  IF( t_su_1d(ji) < rt0 ) THEN   !--  case 1 : no surface melting

                     jm_min(ji) = 1
                     jm_max(ji) = nlay_s + 1

                     ! surface equation
                     ztrid   (ji,1,1) = 0._wp
                     ztrid   (ji,1,2) = zdqns_ice_b(ji) - zg1s * zkappa_s(ji,0)
                     ztrid   (ji,1,3) =                   zg1s * zkappa_s(ji,0)
                     zindterm(ji,1)   = zdqns_ice_b(ji) * t_su_1d(ji) - zfnet(ji)

                     ! first layer of snow equation
                     ztrid   (ji,2,1) =       - zeta_s(ji,1) *                    zkappa_s(ji,0) * zg1s
                     ztrid   (ji,2,2) = 1._wp + zeta_s(ji,1) * ( zkappa_s(ji,1) + zkappa_s(ji,0) * zg1s )
                     ztrid   (ji,2,3) =       - zeta_s(ji,1) *   zkappa_s(ji,1)
                     zindterm(ji,2)   = ztsold(ji,1) + zeta_s(ji,1) * zradab_s(ji,1)

                  ELSE                            !--  case 2 : surface is melting
                     !
                     jm_min(ji) = 2
                     jm_max(ji) = nlay_s + 1

                     ! first layer of snow equation
                     ztrid   (ji,2,1) = 0._wp
                     ztrid   (ji,2,2) = 1._wp + zeta_s(ji,1) * ( zkappa_s(ji,1) + zkappa_s(ji,0) * zg1s )
                     ztrid   (ji,2,3) =       - zeta_s(ji,1) *   zkappa_s(ji,1)
                     zindterm(ji,2)   = ztsold(ji,1) + zeta_s(ji,1) * ( zradab_s(ji,1) + zkappa_s(ji,0) * zg1s * t_su_1d(ji) )
                  ENDIF

                  !
                  zindtbis(ji,jm_min(ji)) = zindterm(ji,jm_min(ji))
                  zdiagbis(ji,jm_min(ji)) = ztrid   (ji,jm_min(ji),2)
                  !
               ENDIF
            END DO
            !
            !------------------------------
            ! 8) tridiagonal system solving
            !------------------------------
            ! Solve the tridiagonal system with Gauss elimination method.
            ! Thomas algorithm, from Computational fluid Dynamics, J.D. ANDERSON, McGraw-Hill 1984
               DO ji = 1, npti
                  IF( h_s_1d(ji) > 0._wp ) THEN
                     jm_mint = jm_min(ji)
                     jm_maxt = jm_max(ji)
                     DO jm = jm_mint+1, jm_maxt
                        zdiagbis(ji,jm) = ztrid   (ji,jm,2) - ztrid(ji,jm,1) * ztrid   (ji,jm-1,3) / zdiagbis(ji,jm-1)
                        zindtbis(ji,jm) = zindterm(ji,jm  ) - ztrid(ji,jm,1) * zindtbis(ji,jm-1  ) / zdiagbis(ji,jm-1)
                     END DO
                  ENDIF
               END DO

               ! snow temperatures
               DO ji = 1, npti
                  ! Variables used after iterations
                  ! Value must be frozen after convergence for MPP independance reason
                  IF ( .NOT. l_T_converged(ji) .AND. h_s_1d(ji) > 0._wp ) &
                     &   t_s_1d(ji,nlay_s) = ( zindtbis(ji,nlay_s+1) - ztrid(ji,nlay_s+1,3) * t_i_1d(ji,1) ) / zdiagbis(ji,nlay_s+1)

             END DO

               !!clem SNWLAY
               DO jm = nlay_s, 2, -1
                  DO ji = 1, npti
                     jk = jm - 1
                     IF ( .NOT. l_T_converged(ji) .AND. h_s_1d(ji) > 0._wp ) &
                        &   t_s_1d(ji,jk) = ( zindtbis(ji,jm) - ztrid(ji,jm,3) * t_s_1d(ji,jk+1) ) / zdiagbis(ji,jm)

                  END DO
               END DO

               ! surface temperature
               DO ji = 1, npti
                  IF( .NOT. l_T_converged(ji) ) THEN
                     ztsub(ji) = t_su_1d(ji)
                     IF( t_su_1d(ji) < rt0 ) THEN
                        ! recompute t_su_1d only if isnow = 1. Without the presence of snow, t_su_1d is computed in
                        !  ice_thd_zdf_bl99_snwext routine
                        IF( h_s_1d(ji) > 0._wp ) t_su_1d(ji) = ( zindtbis(ji,jm_min(ji)) - ztrid(ji,jm_min(ji),3) *  &
                           &          ( isnow(ji) * t_s_1d(ji,1) ) ) / zdiagbis(ji,jm_min(ji))
                     ENDIF
                  ENDIF
               END DO
               !
               !--------------------------------------------------------------
               ! 9) Has the scheme converged?, end of the iterative procedure
               !--------------------------------------------------------------
               ! check that nowhere it has started to melt
               ! zdti_max is a measure of error, it has to be under zdti_bnd

               DO ji = 1, npti

                  zdti_max = 0._wp

                  IF ( .NOT. l_T_converged(ji) ) THEN
                     ! Check convergence on surface T° only in the presence of snow. 

                     IF( h_s_1d(ji) > 0._wp ) THEN
                        t_su_1d(ji) = MAX( MIN( t_su_1d(ji) , rt0 ) , rt0 - 100._wp ) 
                        zdti_max    = MAX( zdti_max, ABS( t_su_1d(ji) - ztsub(ji) ) )    
                        DO jk = 1, nlay_s
                           t_s_1d(ji,jk) = MAX( MIN( t_s_1d(ji,jk), rt0 ), rt0 - 100._wp )
                           zdti_max      = MAX ( zdti_max , ABS( t_s_1d(ji,jk) - ztsb(ji,jk) ) )
                        END DO
                     ENDIF

                     ! convergence test
                     IF( ln_zdf_chkcvg ) THEN
                        tice_cvgerr_1d(ji) = zdti_max
                        tice_cvgstp_1d(ji) = REAL(iconv)
                     ENDIF
                     
                     ! if there is no snow, l_T_converged will automatically be T
                     IF( zdti_max < zdti_bnd )   l_T_converged(ji) = .TRUE.

                  ENDIF
            END DO

         ENDIF ! k_cnd
      END DO  ! End of the do while iterative procedure

      !
      !-----------------------------
      ! 10) Fluxes at the interfaces
      !-----------------------------

      !
      ! --- calculate conduction fluxes at the snow / ice interface 
      !     
      DO ji = 1, npti
         qcn_snw_bot_1d(ji) = 0._wp
         ! The kappa factor is obtained as in eq. 3.61 in LIM3 book by equalling the lowermost heat
         ! conductive flux in the snow and the uppermost one in the sea ice
         zkappa_si(ji) = isnow(ji) * rn_cnd_s * ztcond_i(ji,0) &
                  &                            / ( 0.5_wp * ( ztcond_i(ji,0) * zh_s(ji) + rn_cnd_s * zh_i(ji) ) )

         ! Conduction flux at snow ice - interval : Fc,si = - Kappa_int * (T_i1 - T_s)
         qcn_snw_bot_1d(ji) = - isnow(ji) * zkappa_si(ji) * ( t_i_1d(ji, 1) - t_s_1d (ji,nlay_s) ) 
         !IF(isnow(ji) == 1) qcn_ice_top_1d(ji) = qcn_snw_bot_1d(ji)  
      END DO

      ! We compute the conduction flux at the snow / ice surface here because it
      ! is used later on in snwthd_dh even when isnow=0.
      IF( k_cnd == np_cnd_OFF .OR. k_cnd == np_cnd_EMU ) THEN
         !
         DO ji = 1, npti
            IF(isnow(ji) == 1._wp) qcn_ice_top_1d(ji) = -           isnow(ji)   * zkappa_s(ji,0) * zg1s * ( t_s_1d(ji,1) - t_su_1d(ji) ) &
               &                 - ( 1._wp - isnow(ji) ) * zkappa_i(ji,0) * zg1 * ( t_i_1d(ji,1) - t_su_1d(ji) )
         END DO
         !
      ELSEIF( k_cnd == np_cnd_ON ) THEN
         !
         DO ji = 1, npti
            qcn_ice_top_1d(ji) = qcn_ice_1d(ji)
         END DO
         !
      ENDIF
      !

      ! --- Diagnose the heat loss due to changing non-solar / conduction flux --- !
      !
      IF( k_cnd == np_cnd_OFF .OR. k_cnd == np_cnd_EMU ) THEN
         !
         DO ji = 1, npti
            hfx_err_difs_1d(ji) = isnow(ji) * (hfx_err_difs_1d(ji) - ( qns_ice_1d(ji) - zqns_ice_b(ji) ) * a_i_1d(ji))
         END DO
         !
      ENDIF
      !
      ! --- Diagnose the heat loss due to non-fully converged temperature solution (should not be above 10-4 W-m2) --- !
      !
      IF( k_cnd == np_cnd_OFF .OR. k_cnd == np_cnd_ON ) THEN
         CALL snw_var_enthalpy
      
        ! zhfx_err = correction on the diagnosed heat flux due to non-convergence of the algorithm used to solve heat equation
        DO ji = 1, npti
           zdq = - zq_ini(ji) + ( SUM( e_i_1d(ji,1:nlay_i) ) * h_i_1d(ji) * r1_nlay_i +  &
               &                   SUM( e_s_1d(ji,1:nlay_s) ) * h_s_1d(ji) * r1_nlay_s )

           IF( isnow(ji) == 1._wp ) THEN

              ! Nb: this part does not pass in debug mode (floating overflow)
              ! Diagnostics are only computed in snow here
              IF( k_cnd == np_cnd_OFF ) THEN

                 IF( t_su_1d(ji) < rt0 ) THEN  ! case T_su < 0degC
                    zhfx_err = ( qns_ice_1d(ji)     + qsr_ice_1d(ji)     - zradtr_s(ji,nlay_s) - qcn_snw_bot_1d(ji)  &
                       &       + zdq * r1_Dt_ice ) * a_i_1d(ji)
                 ELSE                          ! case T_su = 0degC
                    zhfx_err = ( qcn_ice_top_1d(ji) + qtr_ice_top_1d(ji) - zradtr_s(ji,nlay_s) - qcn_snw_bot_1d(ji)  &
                       &       + zdq * r1_Dt_ice ) * a_i_1d(ji)
                 ENDIF

              ELSEIF( k_cnd == np_cnd_ON ) THEN

                 zhfx_err    = ( qcn_ice_top_1d(ji) + qtr_ice_top_1d(ji) - zradtr_s(ji,nlay_s) - qcn_snw_bot_1d(ji)  &
                    &          + zdq * r1_Dt_ice ) * a_i_1d(ji)

              ENDIF
              !
              ! total heat sink to be sent to the ocean
              hfx_err_difs_1d(ji) = isnow(ji) * (hfx_err_difs_1d(ji) + zhfx_err)

              !
              ! hfx_difs = Heat flux diagnostic of sensible heat used to warm/cool ice in W.m-2
              hfx_difs_1d(ji) = hfx_difs_1d(ji) - zdq * r1_Dt_ice * a_i_1d(ji)
              !
           ENDIF
        END DO
         
      ENDIF

      !
      !--------------------------------------------------------------------
      ! 11) reset inner snow temperature
      !--------------------------------------------------------------------
      !
      IF( k_cnd == np_cnd_EMU ) THEN
         ! Restore temperatures to their initial values
         t_s_1d    (1:npti,:) = ztsold        (1:npti,:)
      ENDIF
      !
      ! --- SIMIP diagnostics
      !
      DO ji = 1, npti
         !--- Snow-ice interfacial temperature (diagnostic SIMIP)
         IF( h_s_1d(ji) >= zhs_ssl ) THEN
            t_si_1d(ji) = (   rn_cnd_s       * h_i_1d(ji) * r1_nlay_i * t_s_1d(ji,nlay_s)   &
               &            + ztcond_i(ji,1) * h_s_1d(ji) * r1_nlay_s * t_i_1d(ji,1)      ) &
               &          / ( rn_cnd_s       * h_i_1d(ji) * r1_nlay_i &
               &            + ztcond_i(ji,1) * h_s_1d(ji) * r1_nlay_s )
         ELSE
            t_si_1d(ji) = t_su_1d(ji)
         ENDIF
      END DO
      !
      !
   END SUBROUTINE snw_thd_zdf
#else
   !!----------------------------------------------------------------------
   !!   Default option         Dummy module          NO  SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE snwthd_zdf

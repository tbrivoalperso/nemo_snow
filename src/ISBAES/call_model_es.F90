!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
!================================================================
SUBROUTINE CALL_MODEL(JI,KSIZE2,PTSTEP, ZP_A_S_FRA, ZP_SNOWBLOW, ZP_PA, ZP_RADXS, ZP_Q_REMA,  ZP_EVAP_REMA, ZP_BDG)
!
USE MODD_CSTS !,       ONLY : XLMTT, XLSTT
USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_MEB_PAR
USE MODD_SNOW_PAR
USE ice1D ! SI3 1D variables
USE ice   ! SI3 variables 
USE sbc_phy        ! Catalog of functions for physical/meteorological parameters in the marine boundary layer
USE par_kind 
USE snwvar
USE icevar
USE MODI_SNOW3L
USE MODE_SNOW3L
USE sbcblk
USE in_out_manager ! I/O manager

IMPLICIT NONE
!
INTEGER, INTENT(IN) :: JI 
INTEGER, INTENT(IN) :: KSIZE2 ! Number of snow layers

REAL(wp), INTENT(IN) :: PTSTEP ! Time step of the integration 
REAL(wp), INTENT(IN) ::   ZP_A_S_FRA    ! ice fraction covered by snow
REAL(wp), INTENT(IN) ::   ZP_SNOWBLOW   ! snowfall fraction blown by wind
REAL(wp), INTENT(IN) ::   ZP_PA   ! Pressure at atmospheric level

REAL(wp),DIMENSION(1), INTENT(out) ::   ZP_RADXS    ! Radiation transmited through the snow
REAL(wp),DIMENSION(1), INTENT(out) ::   ZP_Q_REMA     ! remaining heat flux from snow melting       (J.m-2)
REAL(wp),DIMENSION(1), INTENT(out) ::   ZP_EVAP_REMA  ! remaining mass flux from snow sublimation   (kg.m-2)
REAL(wp),DIMENSION(1), INTENT(out) ::   ZP_BDG        ! Heat budget diagnostic 

!
REAL(wp), DIMENSION(1,KSIZE2) :: ZP_SNOWSWE
REAL(wp), DIMENSION(1,KSIZE2) :: ZP_SNOWDZ
REAL(wp), DIMENSION(1,KSIZE2) :: ZP_SNOWRHO
REAL(wp), DIMENSION(1,KSIZE2) :: ZP_SNOWHEAT
REAL(wp), DIMENSION(1,KSIZE2) :: ZP_SNOWTEMP
REAL(wp), DIMENSION(1,KSIZE2) :: ZP_SNOWLIQ
REAL(wp), DIMENSION(1,KSIZE2) :: ZP_SNOWGRAN1
REAL(wp), DIMENSION(1,KSIZE2) :: ZP_SNOWGRAN2
REAL(wp), DIMENSION(1,KSIZE2) :: ZP_SNOWHIST
REAL(wp), DIMENSION(1,KSIZE2) :: ZP_SNOWAGE
REAL(wp),DIMENSION(1)         :: ZP_SNOWALB
REAL(wp),DIMENSION(1)         :: ZP_SWNETSNOW
REAL(wp),DIMENSION(1)         :: ZP_SWNETSNOWS
REAL(wp),DIMENSION(1)         :: ZP_LWNETSNOW
REAL(wp),DIMENSION(1)         :: ZP_PS
REAL(wp),DIMENSION(1)         :: ZP_SRSNOW
REAL(wp),DIMENSION(1)         :: ZP_RRSNOW
REAL(wp),DIMENSION(1)         :: ZP_PSN3L
REAL(wp),DIMENSION(1)         :: ZP_TA
REAL(wp),DIMENSION(1)         :: ZP_CT
REAL(wp),DIMENSION(1)         :: ZP_TG
REAL(wp),DIMENSION(1)         :: ZP_D_G
REAL(wp),DIMENSION(1)         :: ZP_SOILD
REAL(wp),DIMENSION(1)         :: ZP_DELHEATG
REAL(wp),DIMENSION(1)         :: ZP_DELHEATG_SFC
REAL(wp),DIMENSION(1)         :: ZP_SW_RAD
REAL(wp),DIMENSION(1)         :: ZP_QA
REAL(wp),DIMENSION(1)         :: ZP_LVTT
REAL(wp),DIMENSION(1)         :: ZP_LSTT
REAL(wp),DIMENSION(1)         :: ZP_VMOD
REAL(wp),DIMENSION(1)         :: ZP_LW_RAD
REAL(wp),DIMENSION(1)         :: ZP_RHOA
REAL(wp),DIMENSION(1)         :: ZP_UREF
REAL(wp),DIMENSION(1)         :: ZP_EXNS
REAL(wp),DIMENSION(1)         :: ZP_EXNA
REAL(wp),DIMENSION(1)         :: ZP_DIRCOSZW
REAL(wp),DIMENSION(1)         :: ZP_ZREF
REAL(wp),DIMENSION(1)         :: ZP_Z0NAT
REAL(wp),DIMENSION(1)         :: ZP_Z0HNAT
REAL(wp),DIMENSION(1)         :: ZP_Z0EFF
REAL(wp),DIMENSION(1)         :: ZP_ALB
REAL(wp),DIMENSION(1)         :: ZP_SOILCOND
REAL(wp),DIMENSION(1)         :: ZP_THRUFAL
REAL(wp),DIMENSION(1)         :: ZP_GRNDFLUX
! REAL(wp)        :: ZP_FLSN_COR   ! Not used 
REAL(wp),DIMENSION(1)         :: ZP_RESTOREN
REAL(wp),DIMENSION(1)         :: ZP_EVAPCOR
REAL(wp),DIMENSION(1)         :: ZP_SOILCOR
REAL(wp),DIMENSION(1)         :: ZP_GFLXCOR
REAL(wp),DIMENSION(1)         :: ZP_RNSNOW
REAL(wp),DIMENSION(1)         :: ZP_HSNOW
REAL(wp),DIMENSION(1)         :: ZP_GFLUXSNOW
REAL(wp),DIMENSION(1)         :: ZP_DELHEATN
REAL(wp),DIMENSION(1)         :: ZP_DELHEATN_SFC
REAL(wp),DIMENSION(1)         :: ZP_DELPHASEN
REAL(wp),DIMENSION(1)         :: ZP_DELPHASEN_SFC
REAL(wp),DIMENSION(1)         :: ZP_MELTSTOT
REAL(wp),DIMENSION(1)         :: ZP_SNREFREEZ

REAL(wp),DIMENSION(1)         :: ZP_SNOWSFCH
REAL(wp),DIMENSION(1)         :: ZP_HPSNOW
REAL(wp),DIMENSION(1)         :: ZP_LES3L
REAL(wp),DIMENSION(1)         :: ZP_LEL3L
REAL(wp),DIMENSION(1)         :: ZP_EVAP
REAL(wp),DIMENSION(1)         :: ZP_SNDRIFT
REAL(wp),DIMENSION(1)         :: ZP_RI
REAL(wp),DIMENSION(1)         :: ZP_QS
REAL(wp),DIMENSION(1)         :: ZP_EMISNOW
REAL(wp),DIMENSION(1)         :: ZP_CDSNOW
REAL(wp),DIMENSION(1)         :: ZP_USTARSNOW
REAL(wp),DIMENSION(1)         :: ZP_CHSNOW
REAL(wp),DIMENSION(1)         :: ZP_SNOWHMASS
REAL(wp),DIMENSION(1)         :: ZP_VEGTYPE
REAL(wp),DIMENSION(1)         :: ZP_FOREST
REAL(wp),DIMENSION(1)         :: ZP_PEW_A_COEF
REAL(wp),DIMENSION(1)         :: ZP_PEW_B_COEF
REAL(wp),DIMENSION(1)         :: ZP_PET_A_COEF
REAL(wp),DIMENSION(1)         :: ZP_PET_B_COEF
REAL(wp),DIMENSION(1)         :: ZP_PEQ_A_COEF
REAL(wp),DIMENSION(1)         :: ZP_PEQ_B_COEF
REAL(wp),DIMENSION(1)         :: ZP_ZENITH
REAL(wp),DIMENSION(1)         :: ZP_LAT,ZP_LON
REAL(wp),DIMENSION(1)         :: ZP_PSN_INV
REAL(wp),DIMENSION(1)         :: ZP_PSN
REAL(wp),DIMENSION(1)         :: ZP_PSN_GFLXCOR
REAL(wp),DIMENSION(1)         :: ZP_WORK
REAL(wp),DIMENSION(1) :: ZP_DELHEAT_SNWFL, ZP_DELHEAT_SUB, ZP_DELHEAT_MLT, ZP_DELHEAT_DIF 

!

REAL(wp) :: h_s_bef
REAL(wp), DIMENSION(1,KSIZE2) :: dh_s_bef


REAL(wp), PARAMETER :: ZDEPTHABS = 0.60 ! m
REAL(wp) ::  ZSCAP
REAL(wp)          ::   zq_ini      ! diag errors on heat
REAL(wp)          ::   zm_ini      ! diag errors on mass
REAL(wp) ::   zdq             ! diag errors on heat
REAL(wp) ::   zdm             ! diag errors on mass

!
INTEGER :: JWRK 


! some options
 CHARACTER(LEN=100)       :: HSNOWRES
!                                      HSNOWRES  = ISBA-SNOW3L turbulant exchange option
!                                      'DEF' = Default: Louis (ISBA: Noilhan and Mahfouf 1996)
!                                      'RIL' = Limit Richarson number under very stable
!                                              conditions (currently testing)
CHARACTER(LEN=100)        :: HIMPLICIT_WIND   ! wind implicitation option
!                                           ! 'OLD' = direct
!                                           ! 'NEW' = Taylor serie, order 1
LOGICAL                 :: OMEB       ! True = coupled to MEB. This means surface fluxes ae IMPOSED
!                                     ! as an upper boundary condition to the explicit snow schemes. 
!                                     ! If = False, then energy
!                                     ! budget and fluxes are computed herein.

LOGICAL                 :: OSI3       ! Coupled with SI3 
!                                     
!                                     
!                                     
CHARACTER(4)            :: OSNOWDRIFT
LOGICAL                 :: OSNOWDRIFT_SUBLIM ! activate snowdrift, sublimation during drift

TYPE(DATE_TIME)         :: TPTIME      ! current date and time


! Define options
HSNOWRES = 'DEF'
HIMPLICIT_WIND = 'OLD'
OMEB = .false.
OSI3 = .true.
!OMEB = OSI3
OSNOWDRIFT = 'NONE' !'DFLT' 
OSNOWDRIFT_SUBLIM = .false.


! Initialise constants (in modd_csts)
XPI         = 2.*ASIN(1.)
XDAY   = 86400.
XSIYEA = 365.25*XDAY*2.*XPI/ 6.283076
XSIDAY = XDAY/(1.+XDAY/XSIYEA)

XKARMAN     = 0.4
XBOLTZ      = 1.380658E-23
XLIGHTSPEED = 299792458.
XPLANCK     = 6.6260755E-34
XAVOGADRO   = 6.0221367E+23

XRADIUS = 6371229.
XG      = 9.80665
XOMEGA = 2.*XPI/XSIDAY

XP00 = 1.E5

XSTEFAN = ( 2.* XPI**5 / 15. ) * ( (XBOLTZ / XPLANCK)* XBOLTZ ) * (XBOLTZ/(XLIGHTSPEED*XPLANCK))**2
XI0     = 1370.
XMD    = 28.9644E-3
XMV    = 18.0153E-3
XRD    = XAVOGADRO * XBOLTZ / XMD
XRV    = XAVOGADRO * XBOLTZ / XMV
XCPD   = 7.* XRD /2.
XCPV   = 4.* XRV
XRHOLW = 1000.
XCL    = 4.218E+3
XCI    = 2.106E+3
XTT    = rt0 ! We use the SI3 parameter for the 0°C 
XTTSI  = XTT - 1.8
XICEC  = 0.5

XTTS   = XTT*(1-XICEC) + XTTSI*XICEC
XLVTT  = 2.5008E+6
XLSTT  = 2.8345E+6
XLMTT  = XLSTT - XLVTT
XESTT  = 611.14

XGAMW  = (XCL - XCPV) / XRV
XBETAW = (XLVTT/XRV) + (XGAMW * XTT)
XALPW  = LOG(XESTT) + (XBETAW /XTT) + (XGAMW *LOG(XTT))
XGAMI  = (XCI - XCPV) / XRV
XBETAI = (XLSTT/XRV) + (XGAMI * XTT)
XALPI  = LOG(XESTT) + (XBETAI /XTT) + (XGAMI *LOG(XTT))

XTH00 = 300.
XRHOLI = 917.
XCONDI = 2.22
NDAYSEC = 24*3600 ! Number of seconds in a day
XSURF_TINY    = 1.0e-80
XSURF_TINY_12 = SQRT    (XSURF_TINY    )
XSURF_EPSILON = EPSILON (XSURF_EPSILON ) * 10.0

! Constants in MEB_PAR

XTAU_LW=0.5
XRAGNC_FACTOR=200.
XKDELTA_WR=0.25
!
ZP_UREF         = 10. ! atm. level for wind => zu => No flux computation no need 
ZP_DIRCOSZW     = 1. ! Cosine of the angle between the normal to the surface and the vertical => = 1 (bertrand) 
ZP_ZREF         = 2. ! atm. level for temp. and humidity 

!
! Initialize:
!
ZP_PSN_GFLXCOR  = 0.
ZP_WORK         = 0.
ZP_SOILD        = 0.
!
! pack the variables
!
! --- diag error on heat diffusion - PART 1 --- !
zq_ini = SUM( e_s_1d(JI,1:nlay_s))!  * dh_s_1d(JI,1:nlay_s) )!* r1_nlay_s 
zm_ini = SUM(rho_s_1d(JI,1:nlay_s) * dh_s_1d(JI,1:nlay_s)) * a_i_1d(JI) 

!WHERE (e_s_1d(:,:) .eq. XUNDEF) e_s_1d(:,:) = 0._wp
!WHERE (e_s_1d(:,:) .eq. XUNDEF) dh_s_1d(:,:) = 0._wp

!DO JWRK=1,KSIZE2
!  IF(dh_s_1d(JI,JWRK) .eq. 0.) THEN
!     dh_s_1d(
!     swe_s_1d(JI,JWRK) = XUNDEF
!     rho_s_1d(JI,JWRK) = XUNDEF
!     e_s_1d(JI,JWRK) = XUNDEF
!     o_s_1d(JI,JWRK) = XUNDEF
!     t_s_1d(JI,JWRK)   = XUNDEF
!     lwc_s_1d(JI,JWRK)   = XUNDEF
!  ENDIF
!ENDDO

DO JWRK=1,KSIZE2
     IF (e_s_1d(JI,JWRK) .eq. 0.) dh_s_1d(JI,JWRK) = 0._wp
     IF(dh_s_1d(JI,JWRK) .eq. 0.) THEN
        ZP_SNOWSWE (1,JWRK) = 0. 
        ZP_SNOWRHO (1,JWRK) = 400. 
        ZP_SNOWTEMP(1,JWRK) = 273.15 
        ZP_SNOWAGE (1,JWRK) = 10.
        ZP_SNOWLIQ (1,JWRK) = 0.
        ZP_SNOWDZ  (1,JWRK) = 0. 
        ZP_SNOWHEAT(1,JWRK) = 0.
     ELSE
 
        ZP_SNOWSWE (1,JWRK) = rho_s_1d(JI,JWRK) * dh_s_1d(JI,JWRK) !swe_s_1d(JI,JWRK) ! Snow layer(s) liquid Water Equivalent (SWE:kg m-2) 
        ZP_SNOWRHO (1,JWRK) = rho_s_1d(JI,JWRK) ! Snow layer(s) averaged density (kg/m3) 
        ZSCAP     = SNOW3LSCAP(ZP_SNOWRHO(1,JWRK))
        ZP_SNOWTEMP(1,JWRK) = t_s_1d(JI,JWRK)   ! Snow temperature => °C ou K ?
        ZP_SNOWAGE (1,JWRK) = o_s_1d (JI,JWRK)  ! Snow age (verifier si c'est x area ou pas)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        ZP_SNOWLIQ (1,JWRK) = lwc_s_1d(JI,JWRK) ! Diagnostique => Pas besoin d'advecter 
        ZP_SNOWDZ  (1,JWRK) = dh_s_1d(JI,JWRK)       ! Snow layer(s) thickness (m) (per layer 
!        ! Compute the snow heat from SI3 T°
!         ZP_SNOWHEAT(JI,JWRK) = ZP_SNOWDZ(JI,JWRK)*( !         ZSCAP*(ZP_SNOWTEMP(JI,JWRK)-XTT)        &
!                      - XLMTT*ZP_SNOWRHO(JI,JWRK) ) +
!                      XLMTT*XRHOLW*ZP_SNOWLIQ(JI,JWRK)
        ZP_SNOWHEAT(1,JWRK) = - e_s_1d(JI,JWRK) !* dh_s_1d(JI,JWRK) ! * a_i_1d(JI) !- e_s_1d(JI,JWRK) * dh_s_1d(JI,JWRK)  * a_i_1d(JI)  ! Snow layer(s) heat content (J/m2) (=> verifier correspondance unités) 
     ENDIF
ENDDO
!
DO JWRK=1,KSIZE2
   ZP_SNOWGRAN1(1,JWRK) = XUNDEF ! Not used
   ZP_SNOWGRAN2(1,JWRK) = XUNDEF ! Not used
   ZP_SNOWHIST (1,JWRK) = XUNDEF ! Not used
ENDDO
!  
ZP_D_G          = h_i_1d(JI) * r1_nlay_i ! Assumed first soil layer thickness (m) 
!
h_s_bef            = SUM(ZP_SNOWDZ  (1,:)) ! Save height for later
ZP_LVTT    (1)     = XLVTT  ! Fourni par modd_csts 
ZP_LSTT    (1)     = XLSTT  ! Fourni par modd_csts 
ZP_EMISNOW (1)     = 0.99 ! Snow Emissivity 
ZP_SNOWALB (1)     = albs_isbaes_1d (JI) ! Snow albedo
ZP_PSN3L   (1)     = 1. !ZP_A_S_FRA      (JI) * a_i_1d(JI) ! Total Snow fraction : ZP_A_S_FRA => à mettre en intent=IN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
ZP_Z0NAT   (1)     = 1e-03 ! Values from surfex
ZP_Z0HNAT  (1)     = 1e-04 ! Values from surfex
ZP_Z0EFF   (1)     = ZP_Z0NAT(1)  ! effective roughness length for momentum => equal to z0  
ZP_RNSNOW  (1)     = (1. - albs_isbaes_1d(JI)) * qsr_ice_isbaes_1d(JI) + qlw_ice_isbaes_1d(JI) ! net radiative flux from snow => qsr_ice ?? => verifier avec VIRGINIE !!!!!!!!!!!!!!!!!!!
ZP_HSNOW   (1)     = qsb_ice_isbaes_1d(JI) 
ZP_HPSNOW  (1)     = qprec_ice_1d(JI) ! heat release from rainfall 

ZP_PS          (1) = slp_isbaes_1d(JI)      ! pressure at the surface => slp !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
ZP_SRSNOW      (1) = snow_isbaes_1d(JI) * ZP_SNOWBLOW / at_i_1d(JI)      ! Snow rate => sprecip_1d (Kg/m2/s) => 
ZP_CT          (1) = 1. / (rcpi * SUM(rho_s_1d(JI,:)) ) !inverse of the product of snow heat capacity and layer thickness [(m2 K)/J] 
ZP_DELHEATG    (1) = 0. ! ground heat content change (diagnostic) (W/m2) JUST NEED TO DECLARE POINTER  
ZP_DELHEATG_SFC(1) = 0. ! ground heat content change in sfc only (diagnostic) (W/m2) JUST NEED TO DECLARE POINTER
ZP_SW_RAD      (1) = qsr_ice_isbaes_1d(JI) !* a_i_1d(JI) ! qtr_ice_top_1d(JI) !/ (1 - albs_isbaes_1d (JI)) !qsr_ice_isbaes_1d(JI) ! Incoming solar radiation 
ZP_QA          (1) = qair_isbaes_1d(JI) ! air humidity at atm. level (kg/kg) => pqair (specific Q in SI3, whatabout ISBA?)
ZP_VMOD        (1) = wndm_isbaes_1d(JI) ! module of the horizontal wind => wndm_ice 
ZP_LW_RAD      (1) = qlwdwn_ice_isbaes_1d(JI) ! NOT USED IF FLUX ARE PRESCRIBED qlw_ice_isbaes_1d(JI) 
ZP_RHOA        (1) = rho_air_isbaes_1d(JI) ! => rhoa => No flux computation no need
ZP_EXNS        (1) = (ZP_PS(1)/XP00)**(XRD/XCPD) ! Exner function at sea surface => No flux computation no need
ZP_EXNA        (1) = (ZP_PA/XP00)**(XRD/XCPD) ! Exner function at atm level => No flux computation no need
ZP_TA          (1) = tair_isbaes_1d(JI) * ZP_EXNA(1)      ! DOIT ETRE LA TEMPERATURE ABSOLUE = multiplier tpot par fonction exner air temperature at atm. level 
ZP_TG          (1) = t_i_1d(JI,1)       * ZP_EXNS(1) ! Ground T° => 1st ice level 
ZP_ALB         (1) = albi_isbaes_1d(JI) ! green areas albedo => snow free albedo   
ZP_RRSNOW      (1) = rain_isbaes_1d(JI) !* a_i_1d(JI) !* ZP_A_S_FRA(JI) ! rain rate over snow [kg/(m2 s)] !!!!!!! MULTIPLIER PAR LA FRACTION DE NEIGE 
ZP_SOILCOND    (1) = cnd_i_isbaes_1d(JI) ! Temporary heat conductivity of litter + soil => conductivity of 1st ice layer 
ZP_PEW_A_COEF  (1) = 0. ! XUNDEF !0. !ZP_VMOD(JI) ! 0. ! Coeff flux => No flux computation no need ! B COEFF 0 et A COEFF 1
ZP_PEW_B_COEF  (1) = ZP_VMOD(1) ! => No flux computation no need
ZP_PET_A_COEF  (1) =  0. 
ZP_PET_B_COEF  (1) =  ZP_TA(1) / (ZP_PA/XP00)**(XRD/XCPD) 
ZP_PEQ_A_COEF  (1) =  0. 
ZP_PEQ_B_COEF  (1) =  ZP_QA(1)
ZP_LAT         (1) = gphit_1d(JI) 
ZP_LON         (1) = glamt_1d(JI) ! => glamt_1d
ZP_ZENITH      (1) = 0. ! solar zenith angle
ZP_GRNDFLUX    (1) = 0. !qcn_snw_bot_1d(JI) ! snow-ground flux before correction (W m-2) 
ZP_DELHEATN    (1) = 0. ! total snow heat content change in the surface layer (W m-2)
ZP_DELHEATN_SFC(1) = 0. ! total snow heat content change during the timestep (W m-2)
ZP_SNOWSFCH    (1) = 0. ! snow surface layer pseudo-heating term owing to changes in grid thickness (W m-2) ! DIAG => PAS BESOIN
ZP_LES3L       (1) = 0. ! qla_ice_isbaes_1d(JI) ! hfx_sub_1d(JI) ! ATTENTION INVERSE AVEC SUBLIMATION !evaporation heat flux from snow (W/m2) !!!!!!!!!!!!!!!! EAU LIQUIDE QUI SEVAPORE ????????????
ZP_LEL3L       (1) = 0. ! ATTENTION INVERSE !! hfx_sub_1d(JI) ! sublimation (W/m2) ! A
ZP_EVAP        (1) = 0. !ZP_LES3L(JI) + ZP_LEL3L(JI) ! total evaporative flux (kg/m2/s) LES + LEL ?????
ZP_SWNETSNOW   (1) = qsr_ice_isbaes_1d(JI) * (1. - albs_isbaes_1d(JI)) ! net shortwave radiation entering top of snowpack (W m-2)
ZP_SWNETSNOWS  (1) = qsr_ice_isbaes_1d(JI) * (1. - albs_isbaes_1d(JI)) ! net shortwave radiation in uppermost layer of snowpack
!/!\ /!\/!\/!\/!\/!\/!\/!\/!\ SWNETSNOWS must take into account solar penetration => for now it is qsr * (1 -albedo)
ZP_LWNETSNOW    (1)= qlw_ice_isbaes_1d(JI) ! net longwave radiation entering top of snowpack
ZP_MELTSTOT     (1)= 0.0
ZP_SNREFREEZ    (1)= 0.0
ZP_DELPHASEN    (1)= 0.0 
ZP_DELPHASEN_SFC(1)= 0.0
!
ZP_VEGTYPE(1)  = 1. ! fraction of permanet snow/ice => ???? ! Mettre 1 => bertrand 
ZP_FOREST(1)   = 0. ! => ?????
!
!
! ===============================================================
! conversion of snow heat from J/m3 into J/m2
!WHERE(ZP_SNOWSWE(:,:)>0.) &
!  ZP_SNOWHEAT(:,:) = ZP_SNOWHEAT(:,:) / ZP_SNOWRHO (:,:) * ZP_SNOWSWE (:,:) 
! ===============================================================
!
ZP_PSN_INV       = 0.
ZP_PSN           = ZP_PSN3L 
!
!
CALL SNOW3L(JI, HSNOWRES, TPTIME, OMEB, OSI3, HIMPLICIT_WIND,                   &
           ZP_PEW_A_COEF, ZP_PEW_B_COEF,                                 &
           ZP_PET_A_COEF, ZP_PEQ_A_COEF,ZP_PET_B_COEF, ZP_PEQ_B_COEF,    &
           ZP_SNOWSWE, ZP_SNOWRHO, ZP_SNOWHEAT, ZP_SNOWALB,              &
           ZP_SNOWGRAN1, ZP_SNOWGRAN2, ZP_SNOWHIST, ZP_SNOWAGE, PTSTEP,  &
           ZP_PS, ZP_SRSNOW, ZP_RRSNOW, ZP_PSN3L, ZP_TA, ZP_TG,     &
           ZP_SW_RAD, ZP_QA, ZP_VMOD, ZP_LW_RAD, ZP_RHOA, ZP_UREF,       &
           ZP_EXNS, ZP_EXNA, ZP_DIRCOSZW, ZP_ZREF, ZP_Z0NAT, ZP_Z0EFF,   &
           ZP_Z0HNAT, ZP_ALB, ZP_SOILCOND, ZP_D_G,                  &
           ZP_LVTT, ZP_LSTT, ZP_SNOWLIQ,                                 &
           ZP_SNOWTEMP, ZP_SNOWDZ, ZP_THRUFAL, ZP_MELTSTOT, ZP_SNREFREEZ,&
           ZP_GRNDFLUX, ZP_EVAPCOR, ZP_SOILCOR, ZP_GFLXCOR, ZP_SNOWSFCH, &
           ZP_DELHEATN, ZP_DELHEATN_SFC, ZP_DELPHASEN, ZP_DELPHASEN_SFC, &
           ZP_SWNETSNOW, ZP_SWNETSNOWS, ZP_LWNETSNOW, ZP_RESTOREN,       &
           ZP_RNSNOW, ZP_HSNOW, ZP_GFLUXSNOW, ZP_HPSNOW, ZP_LES3L,       &
           ZP_LEL3L, ZP_EVAP, ZP_SNDRIFT, ZP_RI,                         &
           ZP_EMISNOW, ZP_CDSNOW, ZP_USTARSNOW,                          &
           ZP_CHSNOW, ZP_SNOWHMASS, ZP_QS, ZP_RADXS, ZP_VEGTYPE,  ZP_FOREST,       &
           ZP_ZENITH, ZP_LAT, ZP_LON, OSNOWDRIFT,OSNOWDRIFT_SUBLIM,      &
           ZP_DELHEAT_SNWFL, ZP_DELHEAT_SUB, ZP_DELHEAT_MLT, ZP_DELHEAT_DIF)

! unpack variables
!
!h_s_1d(JI) = 0.
h_s_1d(JI) = 0.
DO JWRK=1,KSIZE2

!     IF(ZP_SNOWSWE(1,JWRK) .eq. XUNDEF) swe_s_1d(JI,JWRK)   = 0.
!     IF(ZP_SNOWRHO(1,JWRK) .eq. XUNDEF) rho_s_1d(JI,JWRK)   = 330.
!     IF(ZP_SNOWHEAT(1,JWRK).eq. XUNDEF) e_s_1d(JI,JWRK)     = 0.
!     IF(ZP_SNOWTEMP(1,JWRK).eq. XUNDEF) t_s_1d(JI,JWRK)     = 273.15
!     IF(ZP_SNOWLIQ(1,JWRK) .eq. XUNDEF) lwc_s_1d(JI,JWRK)   = 0.
!     IF(ZP_SNOWDZ(1,JWRK).eq. XUNDEF) dh_s_1d(JI,JWRK)   = 0.

     swe_s_1d(JI,JWRK) = ZP_SNOWSWE  (1,JWRK) ! Should we advect it ???
     rho_s_1d(JI,JWRK) = ZP_SNOWRHO  (1,JWRK)
     e_s_1d(JI,JWRK)   =  - ZP_SNOWHEAT(1,JWRK) !/ (ZP_SNOWDZ   (JI,JWRK) )!* a_i_1d(JI))
     o_s_1d(JI,JWRK)   = ZP_SNOWAGE  (1,JWRK)
     t_s_1d(JI,JWRK)   = ZP_SNOWTEMP (1,JWRK)
     lwc_s_1d(JI,JWRK) = ZP_SNOWLIQ  (1,JWRK) ! No need because it is a diagnostic
     dh_s_1d(JI,JWRK)  = ZP_SNOWSWE(1,JWRK)/ZP_SNOWRHO(1,JWRK) !ZP_SNOWDZ   (JI,JWRK) 

ENDDO
h_s_1d(JI) = SUM(dh_s_1d(JI,:))

!Ensure that if no snow, all variables equals to 0

!DO JWRK=1,KSIZE2
!  IF(dh_s_1d(JI,JWRK) .eq. 0.) THEN
!
!     swe_s_1d(JI,JWRK) = 0.
!     rho_s_1d(JI,JWRK) = 0.
!     e_s_1d(JI,JWRK) = 0.
!     o_s_1d(JI,JWRK) = 0.
!     t_s_1d(JI,JWRK)   = 0.
!     lwc_s_1d(JI,JWRK)   = 0.
!  ENDIF
!ENDDO

zdq = SUM( e_s_1d(JI,1:nlay_s)) - zq_ini  !  * dh_s_1d(JI,1:nlay_s) )
zdm = - zm_ini + SUM(rho_s_1d(JI,1:nlay_s) * dh_s_1d(JI,1:nlay_s)) * a_i_1d(JI) 


t_su_1d(JI) = t_s_1d(JI,1)
albs_isbaes_1d(JI)   = ZP_SNOWALB(1)     
!   PEK%TSNOW%EMIS(JI)  = ZP_EMISNOW     (JI) ! Constante
!  DMK%XCDSNOW   (JI)  = ZP_CDSNOW      (JI) ! OUT
!  DMK%XUSTARSNOW(JI)  = ZP_USTARSNOW   (JI) ! OUT
!  DMK%XCHSNOW   (JI)  = ZP_CHSNOW      (JI) ! OUT
!  DMK%XSNOWHMASS(JI)  = ZP_SNOWHMASS   (JI) 
!  DMK%XRNSNOW   (JI)  = ZP_RNSNOW      (JI)
!  DMK%XHSNOW    (JI)  = ZP_HSNOW       (JI)
!  DMK%XHPSNOW  (JI)   = ZP_HPSNOW      (JI)
!  DMK%XGFLUXSNOW(JI)  = ZP_GFLUXSNOW   (JI)
  !
!  PDELHEATG    (JI)   = ZP_DELHEATG    (JI)
!  PDELHEATG_SFC(JI)   = ZP_DELHEATG_SFC(JI)

!  PRI          (JI)   = ZP_RI          (JI)
!  PQS          (JI)   = ZP_QS          (JI)
qcn_snw_bot_1d(JI)  = (ZP_GRNDFLUX(1) ) !    + ZP_GFLXCOR(1))     ! Somme des flux radiatifs et convectif ??
  
  !  PFLSN_COR     (JI)  = ZP_FLSN_COR    (JI) ! Not used
!  PDELHEATN    (JI)   = ZP_DELHEATN    (JI)
!  PDELHEATN_SFC(JI)   = ZP_DELHEATN_SFC(JI)
!  PSNOWSFCH    (JI)   = ZP_SNOWSFCH    (JI) ! Check avec bertrand
!  PGSFCSNOW    (JI)   = ZP_RESTOREN    (JI)  !heat flux between the surface and sub-surface 
!                                                  snow layers (W/m2)
!  PLES3L       (JI)   = ZP_LES3L       (JI)
!  PLEL3L       (JI)   = ZP_LEL3L       (JI)
!  PEVAP        (JI)   = ZP_EVAP        (JI)
!  ZSOILCOR     (JI)   = ZP_SOILCOR     (JI)  ! PSOILCOR => C'est quoi la diff avec EVAPCOR => Faire somme des deux pour zevap_ram ??
  !
!  qsr_ice_1d(JI) = ZP_SWNETSNOW   (JI)
   !qns_ice_1d(JI) = ZP_LWNETSNOW   (JI)
!  ZSWNET_NS  (JI) = ZP_SWNETSNOWS  (JI)
!  ZLWNET_N   (JI) = ZP_LWNETSNOW   (JI)
   
! Heat fluxes for budget diagnostics
  ! We put 0 to all heat fluxes except from hfx_snw that we now consider as the total heat content change in the snow
hfx_sub_1d(JI) = hfx_sub_1d(JI)   - ZP_DELHEAT_SUB(1) * r1_Dt_ice ! Minus factor following SI3 conventions
hfx_difs_1d(JI) = hfx_difs_1d(JI) + ZP_DELHEAT_DIF(1) * r1_Dt_ice 
hfx_spr_1d(JI)  = hfx_spr_1d(JI)  - ZP_DELHEAT_SNWFL(1) * r1_Dt_ice ! Minus factor following SI3 conventions
hfx_snw_1d(JI)  = hfx_snw_1d(JI)  + ZP_DELHEAT_MLT(1) * r1_Dt_ice
hfx_res_1d(JI)  = hfx_res_1d(JI)  + zdq * r1_Dt_ice + (ZP_DELHEAT_SUB(1) + ZP_DELHEAT_DIF(1) + ZP_DELHEAT_SNWFL(1) + ZP_DELHEAT_MLT(1) )* r1_Dt_ice ! Minus factor because enthalpy is negative per convention in ISBAES 
!PRINT*,'zdq + ZP_DELHEAT_SUB(1) + ZP_DELHEAT_DIF(1) + ZP_DELHEAT_SNWFL(1) + ZP_DELHEAT_MLT(1)',(zdq + ZP_DELHEAT_SUB(1) + &
!& ZP_DELHEAT_DIF(1) + ZP_DELHEAT_SNWFL(1) + ZP_DELHEAT_MLT(1))  * r1_Dt_ice
!PRINT*,'ZP_GFLXCOR(:)',ZP_GFLXCOR(:)
!PRINT*,'ZP_SOILCOR(:)',ZP_SOILCOR(:)
!PRINT*,'ZP_EVAPCOR(:)',ZP_EVAPCOR(:)
!PRINT*,'ZP_SNOWSFCH(:)',ZP_SNOWSFCH(:)
!PRINT*,'ZP_GRNDFLUX(1)',ZP_GRNDFLUX(1)
!PRINT*,'BUDGET ', ZP_GFLUXSNOW(1) + ZP_SNOWHMASS(1)* r1_Dt_ice - ZP_GRNDFLUX(1)
!PRINT*,'ZP_DELHEAT_SNWFL(1)',ZP_DELHEAT_SNWFL(1),'ZP_SNOWHMASS(1)',ZP_SNOWHMASS
!PRINT*,'ZP_DELHEAT_SUB(1)',ZP_DELHEAT_SUB(1), 'ZP_DELHEAT_MLT(1)', ZP_DELHEAT_MLT(1), 'ZP_DELHEAT_MLT(1) + ZP_DELHEAT_SUB(1)', &
!        & ZP_DELHEAT_MLT(1) + ZP_DELHEAT_SUB(1),'ZP_DELPHASEN(1)',ZP_DELPHASEN(1) 
!PRINT*,'ZP_DELHEAT_DIF(1)',ZP_DELHEAT_DIF(1),'ZP_DELHEATN(1)',ZP_DELHEATN(1), 'ZP_GFLUXSNOW(:)', ZP_GFLUXSNOW(:)
!PRINT*,'ZP_DELHEATN(:)+ZP_DELPHASEN(:)-ZP_GFLUXSNOW(:)+ZP_GRNDFLUX(:)+ZP_GFLXCOR(:)',ZP_DELHEATN(:)+ZP_DELPHASEN(:)-ZP_GFLUXSNOW(:)+ZP_GRNDFLUX(:)+ZP_GFLXCOR(:)
!PRINT*,'ZP_DELHEATN_SFC(:)-ZP_SNOWSFCH(:)+ZP_SWNETSNOW(:)-ZP_SWNETSNOWS(:)-ZP_GFLUXSNOW(:)+ZP_RESTOREN(:)', ZP_DELHEATN_SFC(:)-ZP_SNOWSFCH(:)+ZP_SWNETSNOW(:)-ZP_SWNETSNOWS(:)-ZP_GFLUXSNOW(:)+ZP_RESTOREN(:)
wfx_spr_1d    (JI)   = wfx_spr_1d    (JI) -(ZP_SRSNOW(1) + ZP_RRSNOW(1) ) * a_i_1d(JI) ! METTRE EVAP A SUB 
wfx_snw_sub_1d   (JI)   =  wfx_snw_sub_1d   (JI) + ((ZP_PSN3L(1)*ZP_LES3L(1)/XLSTT) - (ZP_EVAPCOR(1) + ZP_SOILCOR(1))) * a_i_1d(JI)
!wfx_snw_sub_1d   (JI)   =  wfx_snw_sub_1d   (JI) + (ZP_PSN3L(1)*ZP_LES3L(1)/XLSTT) * a_i_1d(JI)
wfx_snw_sum_1d(JI)   = wfx_snw_sum_1d(JI) + ZP_THRUFAL(1) * a_i_1d(JI) ! rate that liquid water leaves snow pack (kg/(m2 s)): 
wfx_res_1d(JI)   = wfx_res_1d(JI) - zdm * r1_Dt_ice - (wfx_spr_1d(JI) + wfx_snw_sub_1d(JI) + wfx_snw_sum_1d(JI))

ZP_Q_REMA   =  (ZP_GFLXCOR(1)      + ZP_RADXS(1) ) * rDt_ice ! En J / m2
!ZP_EVAP_REMA      = (ZP_EVAP(1) - wfx_snw_sub_1d(JI)) * rDt_ice + (ZP_EVAPCOR(1) + ZP_SOILCOR(1)) * rDt_ice
ZP_EVAP_REMA      =  (ZP_EVAPCOR(1) + ZP_SOILCOR(1)) * rDt_ice

ZP_BDG = - zdq * r1_Dt_ice - (ZP_GFLUXSNOW(1) + ZP_SNOWHMASS(1)* r1_Dt_ice - ZP_GRNDFLUX(1) - ZP_GFLXCOR(1) - ZP_RADXS(1))
PRINT*,ZP_BDG
!- ZP_GFLUXSNOW(1) + - ZP_SNOWHMASS(1)* r1_Dt_ice + ZP_GRNDFLUX(1) + ZP_GFLXCOR(1) + ZP_RADXS(1)  
!Surface total flux => Qtot = qns_ice_1d(ji) + qsr_ice_1d(ji) - qtr_ice_top_1d(ji) - qcn_ice_top_1d(ji)                                            
qla_ice_isbaes_1d(JI) = ZP_LES3L(1) + ZP_LEL3L(1)
qsb_ice_isbaes_1d(JI) = ZP_HSNOW(1)
qlw_ice_isbaes_1d(JI) = ZP_LWNETSNOW(1)
qns_ice_1d(JI) = qlw_ice_isbaes_1d(JI) - qla_ice_isbaes_1d(JI) - qsb_ice_isbaes_1d(JI) !ZP_LES3L(1) + ZP_LEL3L(1) + ZP_HSNOW(1) + ZP_LWNETSNOW(1)

qsr_ice_1d(JI) = ZP_SWNETSNOWS(1)
qemp_ice_1d(JI) = ZP_DELHEAT_SNWFL(1) * r1_Dt_ice

!  qcn_ice_top_1d(JI) = ZP_RESTOREN(JI)
!  qtr_ice_top_1d(JI) =  ZP_SW_RAD(JI) * ZP_SNOWALB     (JI)

!
END SUBROUTINE CALL_MODEL
!

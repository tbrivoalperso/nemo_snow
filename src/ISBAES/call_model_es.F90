!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
!================================================================
SUBROUTINE CALL_MODEL(KT,JI,KSIZE2,PTSTEP, ZP_A_S_FRA, ZP_SNOWBLOW, ZP_PA, ZP_RADXS, ZP_Q_REMA,  ZP_EVAP_REMA, ZP_BDG)
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
USE daymod

IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KT
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
REAL(wp), DIMENSION(1,KSIZE2) :: ZP_SCOND_ES

REAL(wp),DIMENSION(1)         :: ZP_SNOWALB
REAL(wp),DIMENSION(1)         :: ZP_SWNETSNOW
REAL(wp),DIMENSION(1)         :: ZP_SWNETSNOWS
REAL(wp),DIMENSION(1)         :: ZP_LWNETSNOW
REAL(wp),DIMENSION(1)         :: ZP_PS
REAL(wp),DIMENSION(1)         :: ZP_SRSNOW
REAL(wp),DIMENSION(1)         :: ZP_RRSNOW
REAL(wp),DIMENSION(1)         :: ZP_PSN3L
REAL(wp),DIMENSION(1)         :: ZP_TA
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
REAL(wp),DIMENSION(1)         :: ZP_PSN
REAL(wp),DIMENSION(1)         :: ZP_PSN_GFLXCOR
REAL(wp),DIMENSION(1)         :: ZP_WORK
REAL(wp),DIMENSION(1) :: ZP_DELHEAT_SNWFL, ZP_DELHEAT_SUB, ZP_DELHEAT_MLT, ZP_DELHEAT_DIF 
REAL(wp),DIMENSION(1)         :: ZP_TSUN, ZP_AZIMSOL

!

REAL(wp) :: wfs_res_bef, hfx_res_bef


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


! -------------------------------------------------------------------------------------------------------------------------------------------
! Model options 
! -------------------------------------------------------------------------------------------------------------------------------------------

HSNOWRES = 'RIL' !'DEF'
HIMPLICIT_WIND = 'OLD'
OMEB = .false.
OSI3 = .true.
OSNOWDRIFT = 'DFLT' !'DFLT' ! 'NONE' 
OSNOWDRIFT_SUBLIM = .false.

! -------------------------------------------------------------------------------------------------------------------------------------------
! Initialise constants (in modd_csts) 
! -------------------------------------------------------------------------------------------------------------------------------------------

! Constants in modd_csts
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
ZP_UREF         = rn_zu_isbaes !10. ! atm. level for wind => zu => No flux computation no need 
ZP_DIRCOSZW     = 1. ! Cosine of the angle between the normal to the surface and the vertical => = 1 (bertrand) 
ZP_ZREF         = rn_zqt_isbaes !10. ! atm. level for temp. and humidity 

ZP_LVTT    (1)     = XLVTT  ! Fourni par modd_csts 
ZP_LSTT    (1)     = XLSTT  ! Fourni par modd_csts 
ZP_EMISNOW (1)     = 0.99 ! Snow Emissivity 

ZP_VEGTYPE(1)  = 1. ! fraction of permanet snow/ice => ???? ! Mettre 1 => bertrand 
ZP_FOREST(1)   = 0. ! => ?????

!
! Initialize:
!
ZP_PSN_GFLXCOR  = 0.
ZP_WORK         = 0.
ZP_SOILD        = 0.

!---------------------------------------------------------------------------------------------
! Initial mass and heat - used for conservation checks 
!---------------------------------------------------------------------------------------------

zq_ini = SUM( e_s_1d(JI,1:nlay_s))
zm_ini = SUM(rho_s_1d(JI,1:nlay_s) * dv_s_1d(JI,1:nlay_s))


!---------------------------------------------------------------------------------------------
!                   Pack variables (= SI3 variables to ISBA-es variables)
! ---------------------------------------------------------------------------------------------
!   - Snow 3D variables 
! ---------------------------------------------------------------------------------------------

DO JWRK=1,KSIZE2
     IF((dh_s_1d(JI,JWRK) .eq. 0._wp) .OR. (a_i_1d(JI) .eq. 0._wp)) THEN
        ZP_SNOWSWE (1,JWRK) = 0. 
        ZP_SNOWRHO (1,JWRK) = 330. 
        ZP_SNOWTEMP(1,JWRK) = 273.15 ! t_su_1d(JI) !273.15 
        ZP_SNOWAGE (1,JWRK) = 10.
        ZP_SNOWLIQ (1,JWRK) = 0.
        ZP_SNOWDZ  (1,JWRK) = 0. 
        ZP_SNOWHEAT(1,JWRK) = 0.
     ELSE
        ZP_SNOWDZ  (1,JWRK) = dh_s_1d(JI,JWRK) 
        ZP_SNOWSWE (1,JWRK) = rho_s_1d(JI,JWRK) * dh_s_1d(JI,JWRK) ! Snow layer(s) liquid Water Equivalent (SWE:kg m-2) 
        ZP_SNOWRHO (1,JWRK) = rho_s_1d(JI,JWRK) ! Snow layer(s) averaged density (kg/m3) 
        ZSCAP     = SNOW3LSCAP(ZP_SNOWRHO(1,JWRK))
        ZP_SNOWTEMP(1,JWRK) = t_s_1d(JI,JWRK)   ! Snow temperature 
        ZP_SNOWAGE (1,JWRK) = 1. !o_s_1d (JI,JWRK)  ! Snow age 
        ZP_SNOWLIQ (1,JWRK) = lwc_s_1d(JI,JWRK) ! Snow liquid water content (diagnostic)
        
        ! Enthalpy: divide by the ice concentration because the enthalpy of SI3
        ! is in J/m2 per unit area + negative factor because of ISBA-es convention
        ZP_SNOWHEAT(1,JWRK) = - e_s_1d(JI,JWRK) / a_i_1d(JI) ! Snow enthalpy in J/m2  
        ! 
     ENDIF
ENDDO
!

!---------------------------------------------------------------------------------------------
!   - Surface snow variables
!---------------------------------------------------------------------------------------------

ZP_SNOWALB (1)     = albs_isbaes_1d (JI) ! Snow albedo
ZP_PSN3L   (1)     = 1. ! Snow fraction within ISBA-ES (not used in our case...  
                        ! if a snow fraction is used, it should be applied to the variables before and after ISBA-es) 
ZP_PSN           = ZP_PSN3L ! 

!---------------------------------------------------------------------------------------------
!   - Below-snow (sea-ice) variables
!---------------------------------------------------------------------------------------------

ZP_TG          (1) = t_i_1d(JI,1)       !* ZP_EXNS(1) ! Ground T° => 1st ice level T° (K)
ZP_ALB         (1) = 0. ! Ice albedo, set to zero for now to be consistent with SI3 (set to albi_isbaes_1d(JI)) 
ZP_SOILCOND    (1) = cnd_i_isbaes_1d(JI) ! Conductivity of 1st ice layer 
ZP_D_G          = h_i_1d(JI) * r1_nlay_i ! First sea-ice layer thickness (m) 

!---------------------------------------------------------------------------------------------
!   - Surface Heat & Mass fluxes & atmospheric variables 
!           * In forced mode: Downward heat fluxes and atmospheric variables should be set.
!           * In coupled mode, or in the case of flux forcing, we should give the net fluxes instead.      
!---------------------------------------------------------------------------------------------

!------------------------------
! CASE 1: Forced mode

! Atmospheric variables
ZP_Z0NAT   (1)     = 2.3e-04 ! Values from surfex
ZP_Z0HNAT  (1)     = 2.3e-04 ! Values from surfex
ZP_Z0EFF   (1)     = ZP_Z0NAT(1)  ! effective roughness length for momentum => equal to z0  
ZP_PS          (1) = slp_isbaes_1d(JI)      ! pressure at the surface (Pa) 
ZP_QA          (1) = qair_isbaes_1d(JI) ! air humidity at atm. level (kg/kg) 
ZP_VMOD        (1) = wndm_isbaes_1d(JI) ! module of the horizontal wind (m/s) 
ZP_RHOA        (1) = rho_air_isbaes_1d(JI) ! Air density (kg/m3) 
ZP_EXNS        (1) = 1. !(ZP_PS(1)/XP00)**(XRD/XCPD) ! Exner function at sea surface 
ZP_EXNA        (1) = 1. !(ZP_PA/XP00)**(XRD/XCPD) ! Exner function at atm level 
ZP_TA          (1) = tair_isbaes_1d(JI) !* ZP_EXNA(1) ! Absolute air temperature (K)   

! Downward heat fluxes
ZP_SW_RAD      (1) = qsr_ice_isbaes_1d(JI) ! Incoming shortwave radiation (W/m2) 
ZP_LW_RAD      (1) = qlwdwn_ice_isbaes_1d(JI) ! Incoming longwave radiation (W/m2)


!------------------------------
! CASE 2: Coupled or flux forced mode 

ZP_RNSNOW  (1)     = 0. ! Net radiative flux (W/m2): SW + LW 
ZP_HSNOW   (1)     = 0. ! Sensible heat flux (W/m2) 
ZP_HPSNOW  (1)     = 0. ! Heat release from rainfall (usually set to zero, because it implies removing heat from the atmosphere) 
ZP_LES3L       (1) = 0. ! Evaporation heat flux from snow (W/m2) 
ZP_LEL3L       (1) = 0. ! Sublimation heat flux from snow (W/m2) 
ZP_EVAP        (1) = 0. ! ZP_LES3L + ZP_LEL3L
ZP_SWNETSNOW   (1) = qsr_ice_isbaes_1d(JI) * (1. - albs_isbaes_1d(JI)) ! Net shortwave radiation entering top of snowpack (W/m2) !
! (qsr_ice_isbaes_1d(JI) * (1. - albs_isbaes_1d(JI)))
ZP_SWNETSNOWS  (1) = qtr_ice_top_1d(JI) ! Net shortwave radiation in uppermost layer of snowpack (qtr_ice_top_1d(JI))
ZP_LWNETSNOW    (1)= 0. ! Net longwave flux (W/m2) (qlw_ice_isbaes_1d(JI) )

!------------------------------
! Mass fluxes, should be set for both cases: 
ZP_SRSNOW      (1) = (snow_isbaes_1d(JI) * ZP_SNOWBLOW / at_i_1d(JI)) ! Snow rate (kg/m2/s)
ZP_RRSNOW      (1) = rain_isbaes_1d(JI) ! Rain rate over snow (kg/m2/s) 

!---------------------------------------------------------------------------------------------
! Geographic variables
!---------------------------------------------------------------------------------------------

ZP_LAT         (1) = gphit_1d(JI) ! Latitude 
ZP_LON         (1) = glamt_1d(JI) ! Londitude

! Zenital angle
CALL SUNPOS (nyear, nmonth, nday, nsec_day, ZP_LON(1), ZP_LAT(1), ZP_TSUN, ZP_ZENITH, ZP_AZIMSOL)

!---------------------------------------------------------------------------------------------
! Variables needed because of the implicit budget formulation
!---------------------------------------------------------------------------------------------

ZP_PEW_A_COEF  (1) = 0.
ZP_PEW_B_COEF  (1) = ZP_VMOD(1)
ZP_PET_A_COEF  (1) =  0.
ZP_PET_B_COEF  (1) =  ZP_TA(1) !/ (ZP_PA/XP00)**(XRD/XCPD)
ZP_PEQ_A_COEF  (1) =  0.
ZP_PEQ_B_COEF  (1) =  ZP_QA(1)

!---------------------------------------------------------------------------------------------
! Dummy variables
!---------------------------------------------------------------------------------------------

DO JWRK=1,KSIZE2
   ZP_SNOWGRAN1(1,JWRK) = XUNDEF ! Not used
   ZP_SNOWGRAN2(1,JWRK) = XUNDEF ! Not used
   ZP_SNOWHIST (1,JWRK) = XUNDEF ! Not used
ENDDO
!

!---------------------------------------------------------------------------------------------
! Initialisation of output vaiables 
!---------------------------------------------------------------------------------------------

ZP_GRNDFLUX    (1) = 0. ! qcn_snw_bot_1d(JI) ! snow-ground flux before correction (W m-2)

ZP_DELHEATG    (1) = 0. ! ground heat content change (diagnostic) (W/m2) 
ZP_DELHEATG_SFC(1) = 0. ! ground heat content change in sfc only (diagnostic) (W/m2) 
ZP_DELHEATN    (1) = 0. ! total snow heat content change in the surface layer (W m-2)
ZP_DELHEATN_SFC(1) = 0. ! total snow heat content change during the timestep (W m-2)
ZP_SNOWSFCH    (1) = 0. ! snow surface layer pseudo-heating term owing to changes in grid thickness (W m-2) ! 
ZP_MELTSTOT     (1)= 0.
ZP_SNREFREEZ    (1)= 0.
ZP_DELPHASEN    (1)= 0.
ZP_DELPHASEN_SFC(1)= 0.

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
           ZP_DELHEAT_SNWFL, ZP_DELHEAT_SUB, ZP_DELHEAT_MLT, ZP_DELHEAT_DIF, ZP_SCOND_ES)


! -------------------------------------------------------------------------------------------------------------------------------------------
! Unpack variables: from ISBA-ES to SI3 variables 
! -------------------------------------------------------------------------------------------------------------------------------------------
! Nb: The enthalpy is defined in J/m2 in ISBA-ES whereas it should be defined in J/m2 per unit area in SI3,
! thus we have to multiply the enthalpy per the ice area to be consistent with SI3 conventions

! ---------------------------------------------------------------------------------------------
!   - Snow 3D variables 
! ---------------------------------------------------------------------------------------------

DO JWRK=1,KSIZE2

     swe_s_1d(JI,JWRK) = ZP_SNOWSWE  (1,JWRK) ! snow water equivalent (kg/m2)
     rho_s_1d(JI,JWRK) = ZP_SNOWRHO  (1,JWRK) ! density (kg/m3)
     !rho_s_1d(JI,JWRK) = 330.

     e_s_1d(JI,JWRK)   =  - ZP_SNOWHEAT(1,JWRK) * a_i_1d(JI) ! Enthalpy (converted in J/m2 per unit area)
     o_s_1d(JI,JWRK)   = ZP_SNOWAGE  (1,JWRK) ! Snow age (s)
     t_s_1d(JI,JWRK)   = ZP_SNOWTEMP (1,JWRK) ! Temperature (K)
     lwc_s_1d(JI,JWRK) = ZP_SNOWLIQ  (1,JWRK) ! Liquid water content (m) 
     dh_s_1d(JI,JWRK)  = ZP_SNOWSWE(1,JWRK)/ZP_SNOWRHO(1,JWRK) ! Snow layer thicknesses (m) 
     dv_s_1d(JI,JWRK)  = dh_s_1d(JI,JWRK) * a_i_1d(JI) ! Snow layer volume (m per unit area)
     rhov_s_1d(JI,JWRK)= rho_s_1d(JI,JWRK) * dv_s_1d(JI,JWRK) ! Snow mass (kg / m2 per unit area)
     ov_s_1d(JI,JWRK)  = o_s_1d(JI,JWRK) * dv_s_1d(JI,JWRK) 
ENDDO


! ---------------------------------------------------------------------------------------------
!   - Snow 2D variables 
! ---------------------------------------------------------------------------------------------
albs_isbaes_1d(JI)   = ZP_SNOWALB(1) ! Snow albedo
v_s_1d(JI) = SUM(dv_s_1d(JI,:)) ! Snowpack volume (m per unit area)
h_s_1d(JI) = SUM(dh_s_1d(JI,:)) ! Snowpack thickness (m) 
t_su_1d(JI) = t_s_1d(JI,1)      ! Surface temperature (K)

! Average snowpack conductivity (diagnostic)
IF(h_s_1d(JI) >  0.000001) THEN
   cnd_s_isbaes_1d(JI) = MIN(100._wp, SUM(ZP_SCOND_ES(1,:) * dh_s_1d(JI,:)) / h_s_1d(JI))
ELSE
   cnd_s_isbaes_1d(JI) = 0.
ENDIF

! Diagnotics 
Cd_ice_isbaes_1d(JI) = ZP_CDSNOW(1)
Ch_ice_isbaes_1d(JI) = ZP_CHSNOW(1)

! ---------------------------------------------------------------------------------------------
! Snow to ice fluxes
! ---------------------------------------------------------------------------------------------
! Conductive heat flux at the snow/ice interface, used to force the sea-ice thermodynamics (W/m2)
qcn_snw_bot_1d(JI)  = ZP_GRNDFLUX(1) ! 

! Remaining heat and mass fluxes if the whole snowpack melts
ZP_Q_REMA   =  ZP_GFLXCOR(1) * rDt_ice * a_i_1d(JI) ! Remaining heat J/m2
ZP_EVAP_REMA      =  (ZP_EVAPCOR(1) + ZP_SOILCOR(1)) * rDt_ice * a_i_1d(JI) !  kg/m2
 
! ---------------------------------------------------------------------------------------------
!   - Surface heat fluxes 
! ---------------------------------------------------------------------------------------------
! Surface heat fluxes to be updated:
!    * qns_ice and qsr_ice should be updated because they are (for now) recomputed by isba-es 
!    * After vertical thermodynamic processes, qns_tot and qsr_tot are updated in icethd_zdf
qns_ice_1d(JI) = ZP_GFLUXSNOW(1) - ZP_SWNETSNOW(1)  ! Non-solar surface heat flux (W/m2)
qsr_ice_1d(JI) = ZP_SWNETSNOW(1)                    ! Net short-wave heat flux (W/m2)
qemp_ice_1d(JI) = ZP_SNOWHMASS(1) * r1_Dt_ice       ! E-P heat flux (W/m2) (precip only for now) 

! Other surface heat fluxes (for diagnotics)
qla_ice_isbaes_1d(JI) = ZP_LES3L(1) + ZP_LEL3L(1)    ! Latent heat flux (W/m2)
qsb_ice_isbaes_1d(JI) = ZP_HSNOW(1)                  ! Sensible heat flux (W/m2)
qlw_ice_isbaes_1d(JI) = ZP_LWNETSNOW(1)              ! Net longwave heat flux (W/m2)

! -------------------------------------------------------------------------------------------------------------------------------------------
! Heat fluxes for budget diagnostics
! -------------------------------------------------------------------------------------------------------------------------------------------

! Compute the variation of enthalpy and mass for the budget / conservation diagnostics

zdq = SUM( e_s_1d(JI,1:nlay_s)) - zq_ini  !  * dh_s_1d(JI,1:nlay_s) )
zdm = - zm_ini + SUM(rho_s_1d(JI,1:nlay_s) * dv_s_1d(JI,1:nlay_s))

! Compute the heat budget
ZP_BDG(1) = - zdq * r1_Dt_ice - (ZP_GFLUXSNOW(1) + ZP_SNOWHMASS(1)* r1_Dt_ice - ZP_GRNDFLUX(1) - ZP_GFLXCOR(1) - ZP_RADXS(1)) * a_i_1d(JI)

!Budget   =  d(e_s)/ dt - sum of heat fluxes (IN and OUT of the snowpack)
! with the sum of the heat fluxes being the sum of the contributions from:
! - The surface heat flux at the snow / atmosphere interface (ZP_GFLUXSNOW)
! - The heat flux from snowfall (i.e: the enthalpy added to the system through snowfall) (ZP_SNOWHMASS(1)* r1_Dt_ice)
! - The conduction flux at the snow-ice interface (ZP_GRNDFLUX(1))
! - The excess heat flux if all snow vanishes (ZP_GFLXCOR(1))
! - The remaining radiative flux flux at the snow-ice interface (ZP_RADXS(1))   
!
! Nb: the sum of the heat fluxes (from ISBA-ES) should be multiplied by the ice area a_i as SI3 enthalpy convention is J/m2 per unit
! area

! Sublimation heat flux (W/M2)
hfx_sub_1d(JI)  = hfx_sub_1d(JI)  - ZP_DELHEAT_SUB(1) * a_i_1d(JI) * r1_Dt_ice ! Minus factor + multiplied per the ice area following SI3 conventions

! Surface precipitation heat flux (W/M2)
hfx_spr_1d(JI)  = hfx_spr_1d(JI)  - ZP_SNOWHMASS  (1) * a_i_1d(JI) * r1_Dt_ice ! Minus factor + multiplied per the ice area following SI3 conventions 

! Heat flux from other snow processes (W/m2) 
hfx_snw_1d(JI)  = hfx_snw_1d(JI)  + ((ZP_GFLUXSNOW(1) - ZP_GRNDFLUX(1) - ZP_GFLXCOR(1) - ZP_RADXS(1)) & 
        + ZP_DELHEAT_SUB(1) * r1_Dt_ice ) * a_i_1d(JI)  ! Minus factor + multiplied per the ice area following SI3 conventions

! Residual heat flux (W/m2, should be very small)
hfx_res_1d(JI)  = hfx_res_1d(JI)  + ZP_BDG(1) ! ZP_BDG already in J/m2 per unit area

! -------------------------------------------------------------------------------------------------------------------------------------------
! Mass fluxes for budget diagnostics
! -------------------------------------------------------------------------------------------------------------------------------------------

! Precipitaion mass flux (kg/m2/s)
wfx_spr_1d    (JI)   = wfx_spr_1d    (JI) -(ZP_SRSNOW(1) + ZP_RRSNOW(1) ) * a_i_1d(JI)  ! METTRE EVAP A SUB 

! Snow sublimation mass flux (kg/m2/s)
wfx_snw_sub_1d   (JI)   =  wfx_snw_sub_1d   (JI) + ((ZP_PSN3L(1)*ZP_LES3L(1)/XLSTT))* a_i_1d(JI) - (ZP_EVAPCOR(1) + ZP_SOILCOR(1)) * a_i_1d(JI)

! Liquid water flux leaving the snowpack (kg/m2/s) 
wfx_snw_sum_1d(JI)   = wfx_snw_sum_1d(JI) + ZP_THRUFAL(1) * a_i_1d(JI)  ! rate that liquid water leaves snow pack (kg/(m2 s))

drhov_s_mlt_1d(JI)   = ZP_THRUFAL(1) * a_i_1d(JI) * rDt_ice 

! Residual (kg/m2/s) (Should be very small)
wfx_res_1d(JI)   = wfx_res_1d(JI) - zdm * r1_Dt_ice - (((ZP_PSN3L(1)*ZP_LES3L(1)/XLSTT) - (ZP_EVAPCOR(1) + ZP_SOILCOR(1))) &
       -(ZP_SRSNOW(1) + ZP_RRSNOW(1) ) + ZP_THRUFAL(1) ) * a_i_1d(JI)

!
END SUBROUTINE CALL_MODEL
!

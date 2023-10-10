!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
!================================================================
SUBROUTINE CALL_MODEL(KSIZE1,KSIZE2,KSIZE3,KMASK,PTSTEP, za_s_fra, ZP_RADXS, zq_rema, zevap_rema)
!
USE MODD_CSTS !,       ONLY : XLMTT, XLSTT
USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_MEB_PAR
USE ice1D ! SI3 1D variables
USE ice   ! SI3 variables 
USE sbc_phy        ! Catalog of functions for physical/meteorological parameters in the marine boundary layer
USE par_kind 
USE snwvar
USE MODI_SNOW3L
USE MODE_SNOW3L
USE sbcblk

IMPLICIT NONE
!
REAL(wp), INTENT(IN)                    :: PTSTEP
!                                      PTSTEP    = time step of the integration
INTEGER, INTENT(IN) :: KSIZE1
INTEGER, INTENT(IN) :: KSIZE2
INTEGER, INTENT(IN) :: KSIZE3
INTEGER, DIMENSION(KSIZE1), INTENT(IN) :: KMASK
REAL, DIMENSION(KSIZE1), INTENT(out) ::   za_s_fra    ! ice fraction covered by snow
REAL, DIMENSION(KSIZE1), INTENT(out) ::   ZP_RADXS    ! Radiation transmited through the snow
REAL, DIMENSION(KSIZE1), INTENT(out) ::   zq_rema     ! remaining heat flux from snow melting       (J.m-2)
REAL, DIMENSION(KSIZE1), INTENT(out) ::   zevap_rema  ! remaining mass flux from snow sublimation   (kg.m-2)

!
REAL(wp), DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSWE
REAL(wp), DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWDZ
REAL(wp), DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWRHO
REAL(wp), DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWHEAT
REAL(wp), DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWTEMP
REAL(wp), DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWLIQ
REAL(wp), DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWGRAN1
REAL(wp), DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWGRAN2
REAL(wp), DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWHIST
REAL(wp), DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWAGE
REAL(wp), DIMENSION(KSIZE1)        :: ZP_SNOWALB
REAL(wp), DIMENSION(KSIZE1)        :: ZP_SWNETSNOW
REAL(wp), DIMENSION(KSIZE1)        :: ZP_SWNETSNOWS
REAL(wp), DIMENSION(KSIZE1)        :: ZP_LWNETSNOW
REAL(wp), DIMENSION(KSIZE1)        :: ZP_PS
REAL(wp), DIMENSION(KSIZE1)        :: ZP_SRSNOW
REAL(wp), DIMENSION(KSIZE1)        :: ZP_RRSNOW
REAL(wp), DIMENSION(KSIZE1)        :: ZP_PSN3L
REAL(wp), DIMENSION(KSIZE1)        :: ZP_TA
REAL(wp), DIMENSION(KSIZE1)        :: ZP_CT
REAL(wp), DIMENSION(KSIZE1) :: ZP_TG
REAL(wp), DIMENSION(KSIZE1) :: ZP_D_G
!REAL(wp), DIMENSION(KSIZE1,KSIZE3) :: ZP_SOILHCAPZ
REAL(wp), DIMENSION(KSIZE1)        :: ZP_SOILD
REAL(wp), DIMENSION(KSIZE1)        :: ZP_DELHEATG
REAL(wp), DIMENSION(KSIZE1)        :: ZP_DELHEATG_SFC
REAL(wp), DIMENSION(KSIZE1)        :: ZP_SW_RAD
REAL(wp), DIMENSION(KSIZE1)        :: ZP_QA
REAL(wp), DIMENSION(KSIZE1)        :: ZP_LVTT
REAL(wp), DIMENSION(KSIZE1)        :: ZP_LSTT
REAL(wp), DIMENSION(KSIZE1)        :: ZP_VMOD
REAL(wp), DIMENSION(KSIZE1)        :: ZP_LW_RAD
REAL(wp), DIMENSION(KSIZE1)        :: ZP_RHOA
REAL(wp), DIMENSION(KSIZE1)        :: ZP_UREF
REAL(wp), DIMENSION(KSIZE1)        :: ZP_EXNS
REAL(wp), DIMENSION(KSIZE1)        :: ZP_EXNA
REAL(wp), DIMENSION(KSIZE1)        :: ZP_DIRCOSZW
REAL(wp), DIMENSION(KSIZE1)        :: ZP_ZREF
REAL(wp), DIMENSION(KSIZE1)        :: ZP_Z0NAT
REAL(wp), DIMENSION(KSIZE1)        :: ZP_Z0HNAT
REAL(wp), DIMENSION(KSIZE1)        :: ZP_Z0EFF
REAL(wp), DIMENSION(KSIZE1)        :: ZP_ALB
REAL(wp), DIMENSION(KSIZE1)        :: ZP_SOILCOND
REAL(wp), DIMENSION(KSIZE1)        :: ZP_THRUFAL
REAL(wp), DIMENSION(KSIZE1)        :: ZP_GRNDFLUX
! REAL(wp), DIMENSION(KSIZE1)        :: ZP_FLSN_COR   ! Not used 
REAL(wp), DIMENSION(KSIZE1)        :: ZP_RESTOREN
REAL(wp), DIMENSION(KSIZE1)        :: ZP_EVAPCOR
REAL(wp), DIMENSION(KSIZE1)        :: ZP_SOILCOR
REAL(wp), DIMENSION(KSIZE1)        :: ZP_GFLXCOR
REAL(wp), DIMENSION(KSIZE1)        :: ZP_RNSNOW
REAL(wp), DIMENSION(KSIZE1)        :: ZP_HSNOW
REAL(wp), DIMENSION(KSIZE1)        :: ZP_GFLUXSNOW
REAL(wp), DIMENSION(KSIZE1)        :: ZP_DELHEATN
REAL(wp), DIMENSION(KSIZE1)        :: ZP_DELHEATN_SFC
REAL(wp), DIMENSION(KSIZE1)        :: ZP_DELPHASEN
REAL(wp), DIMENSION(KSIZE1)        :: ZP_DELPHASEN_SFC
REAL(wp), DIMENSION(KSIZE1)        :: ZP_MELTSTOT
REAL(wp), DIMENSION(KSIZE1)        :: ZP_SNREFREEZ

REAL(wp), DIMENSION(KSIZE1)        :: ZP_SNOWSFCH
REAL(wp), DIMENSION(KSIZE1)        :: ZP_HPSNOW
REAL(wp), DIMENSION(KSIZE1)        :: ZP_LES3L
REAL(wp), DIMENSION(KSIZE1)        :: ZP_LEL3L
REAL(wp), DIMENSION(KSIZE1)        :: ZP_EVAP
REAL(wp), DIMENSION(KSIZE1)        :: ZP_SNDRIFT
REAL(wp), DIMENSION(KSIZE1)        :: ZP_RI
REAL(wp), DIMENSION(KSIZE1)        :: ZP_QS
REAL(wp), DIMENSION(KSIZE1)        :: ZP_EMISNOW
REAL(wp), DIMENSION(KSIZE1)        :: ZP_CDSNOW
REAL(wp), DIMENSION(KSIZE1)        :: ZP_USTARSNOW
REAL(wp), DIMENSION(KSIZE1)        :: ZP_CHSNOW
REAL(wp), DIMENSION(KSIZE1)        :: ZP_SNOWHMASS
REAL(wp), DIMENSION(KSIZE1)        :: ZP_VEGTYPE
REAL(wp), DIMENSION(KSIZE1)        :: ZP_FOREST
REAL(wp), DIMENSION(KSIZE1)        :: ZP_PEW_A_COEF
REAL(wp), DIMENSION(KSIZE1)        :: ZP_PEW_B_COEF
REAL(wp), DIMENSION(KSIZE1)        :: ZP_PET_A_COEF
REAL(wp), DIMENSION(KSIZE1)        :: ZP_PET_B_COEF
REAL(wp), DIMENSION(KSIZE1)        :: ZP_PEQ_A_COEF
REAL(wp), DIMENSION(KSIZE1)        :: ZP_PEQ_B_COEF
REAL(wp), DIMENSION(KSIZE1)        :: ZP_ZENITH
REAL(wp), DIMENSION(KSIZE1)        :: ZP_LAT,ZP_LON
REAL(wp), DIMENSION(KSIZE1)        :: ZP_PSN_INV
REAL(wp), DIMENSION(KSIZE1)        :: ZP_PSN
REAL(wp), DIMENSION(KSIZE1)        :: ZP_PSN_GFLXCOR
REAL(wp), DIMENSION(KSIZE1)        :: ZP_WORK
!
REAL(wp), DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWDEND
REAL(wp), DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSPHER
REAL(wp), DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSIZE
REAL(wp), DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSSA
REAL(wp), DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWTYPEMEPRA
REAL(wp), DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWRAM
REAL(wp), DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSHEAR
REAL(wp), DIMENSION(KSIZE1) :: ZP_SNDPT_1DY
REAL(wp), DIMENSION(KSIZE1) :: ZP_SNDPT_3DY
REAL(wp), DIMENSION(KSIZE1) :: ZP_SNDPT_5DY
REAL(wp), DIMENSION(KSIZE1) :: ZP_SNDPT_7DY
REAL(wp), DIMENSION(KSIZE1) :: ZP_SNSWE_1DY
REAL(wp), DIMENSION(KSIZE1) :: ZP_SNSWE_3DY
REAL(wp), DIMENSION(KSIZE1) :: ZP_SNSWE_5DY
REAL(wp), DIMENSION(KSIZE1) :: ZP_SNSWE_7DY
REAL(wp), DIMENSION(KSIZE1) :: ZP_SNRAM_SONDE
REAL(wp), DIMENSION(KSIZE1) :: ZP_SN_WETTHCKN
REAL(wp), DIMENSION(KSIZE1) :: ZP_SN_REFRZNTHCKN

REAL(wp), DIMENSION(KSIZE1) :: h_s_bef
REAL(wp), DIMENSION(KSIZE1,KSIZE2) :: dh_s_bef
REAL(wp), DIMENSION(KSIZE1,KSIZE2) :: oh_s_1d


! Additional intermediate variables
REAL(wp), DIMENSION(KSIZE1) :: ZP_PA ! Atmospheric pressure at forcing level
!
REAL(wp), PARAMETER :: ZDEPTHABS = 0.60 ! m
REAL(wp) ::  ZSCAP
REAL(wp), DIMENSION(KSIZE1)          ::   zq_ini      ! diag errors on heat
REAL(wp), DIMENSION(KSIZE1)          ::   zm_ini      ! diag errors on mass
REAL(wp), DIMENSION(KSIZE1) ::   zdq             ! diag errors on heat
REAL(wp), DIMENSION(KSIZE1) ::   zdm             ! diag errors on mass


!
INTEGER :: JWRK, JJ, JI


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
OSNOWDRIFT = 'DFLT' 
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
XTTS   = XTT*(1-XICEC) + XTTSI*XICEC
XICEC  = 0.5
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

PRINT*,'H before snow3l', h_s_1d(1)
PRINT*,'T° before snow3l',t_s_1d(1,:)
!
! Initialize:
!
ZP_PSN_GFLXCOR(:)  = 0.
ZP_WORK(:)         = 0.
ZP_SOILD(:)        = 0.
!
! pack the variables
!
! Compute the T° from SI3 enthalpy 
!DO JWRK=1,KSIZE2
!   DO JJ=1,KSIZE1
!      IF( h_s_1d(JJ) > 0._wp ) THEN
!         t_s_1d(JJ,JWRK) = rt0 + ( - e_s_1d(JJ,JWRK) * (1 / rho_s_1d(JJ,JWRK)) * r1_rcpi + rLfus * r1_rcpi ) 
!      ELSE
!         t_s_1d(JJ,JWRK) = rt0
!      ENDIF
!   END DO
!END DO

! --- diag error on heat diffusion - PART 1 --- !
DO JI = 1, npti
   zq_ini(JI) = SUM( e_s_1d(JI,1:nlay_s))!  * dh_s_1d(JI,1:nlay_s) )!* r1_nlay_s 
   zm_ini(JI) = SUM(rho_s_1d(JI,1:nlay_s) * dh_s_1d(JI,1:nlay_s)) * a_i_1d(JI) 
END DO


DO JWRK=1,KSIZE2
   DO JJ=1,KSIZE1
      JI = KMASK(JJ)
      ZP_SNOWSWE (JI,JWRK) = rho_s_1d(JI,JWRK) * dh_s_1d(JI,JWRK) !swe_s_1d(JI,JWRK) ! Snow layer(s) liquid Water Equivalent (SWE:kg m-2) 
      ZP_SNOWRHO (JI,JWRK) = rho_s_1d(JI,JWRK) ! Snow layer(s) averaged density (kg/m3) 
      ZSCAP     = SNOW3LSCAP(ZP_SNOWRHO(JI,JWRK))
      ZP_SNOWTEMP(JI,JWRK) = t_s_1d(JI,JWRK)   ! Snow temperature => °C ou K ?
      ZP_SNOWAGE (JI,JWRK) = o_s_1d (JI,JWRK)  ! Snow age (verifier si c'est x area ou pas)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ZP_SNOWLIQ (JI,JWRK) = lwc_s_1d(JI,JWRK) ! Diagnostique => Pas besoin d'advecter 
      ZP_SNOWDZ  (JI,JWRK) = dh_s_1d(JI,JWRK)       ! Snow layer(s) thickness (m) (per layer 

      ! Compute the snow heat from SI3 T°
!      ZP_SNOWHEAT(JI,JWRK) = ZP_SNOWDZ(JI,JWRK)*( ZSCAP*(ZP_SNOWTEMP(JI,JWRK)-XTT)        &
!                   - XLMTT*ZP_SNOWRHO(JI,JWRK) ) + XLMTT*XRHOLW*ZP_SNOWLIQ(JI,JWRK)
      ZP_SNOWHEAT(JI,JWRK) = e_s_1d(JI,JWRK) !* dh_s_1d(JI,JWRK) ! * a_i_1d(JI) !- e_s_1d(JI,JWRK) * dh_s_1d(JI,JWRK)  * a_i_1d(JI)  ! Snow layer(s) heat content (J/m2) (=> verifier correspondance unités) 
   ENDDO
     
ENDDO
!
DO JWRK=1,KSIZE2
   DO JJ=1,KSIZE1
      ZP_SNOWGRAN1(JJ,JWRK) = XUNDEF ! Not used
      ZP_SNOWGRAN2(JJ,JWRK) = XUNDEF ! Not used
      ZP_SNOWHIST (JJ,JWRK) = XUNDEF ! Not used
   ENDDO
ENDDO
!  
DO JJ=1,KSIZE1
   JI = KMASK(JJ)
   ZP_D_G      (JI) = h_i_1d(JI) * r1_nlay_i ! Assumed first soil layer thickness (m) 
ENDDO
!


DO JJ=1,KSIZE1
   JI = KMASK(JJ)

   h_s_bef    (JI) = SUM(ZP_SNOWDZ  (JI,:)) ! Save height for later
   ZP_LVTT    (JI) = XLVTT  ! Fourni par modd_csts 
   ZP_LSTT    (JI) = XLSTT  ! Fourni par modd_csts 
   ZP_EMISNOW (JI) = 0.99 ! Snow Emissivity 
   ZP_SNOWALB (JI) = albs_isbaes_1d (JI) ! Snow albedo
   ZP_PSN3L   (JI) = za_s_fra      (JI) ! Total Snow fraction : za_s_fra => à mettre en intent=IN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ZP_Z0NAT   (JI) = 1e-03 ! Values from surfex
   ZP_Z0HNAT  (JI) = 1e-04 ! Values from surfex
   ZP_Z0EFF   (JI) = ZP_Z0NAT   (JI) ! effective roughness length for momentum => equal to z0  
   ZP_RNSNOW  (JI) = (1. - albs_isbaes_1d(JI)) * qsr_ice_isbaes_1d(JI) + qlw_ice_isbaes_1d(JI) ! net radiative flux from snow => qsr_ice ?? => verifier avec VIRGINIE !!!!!!!!!!!!!!!!!!!
   ZP_HSNOW   (JI) = qsb_ice_isbaes_1d(JI) 
   ZP_HPSNOW  (JI) = qprec_ice_1d(JI) ! heat release from rainfall 

   ZP_PS      (JI) = slp_isbaes_1d(JI)      ! pressure at the surface => slp !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ZP_SRSNOW  (JI) = snow_isbaes_1d(JI) * a_i_1d(JI)      ! Snow rate => sprecip_1d (Kg/m2/s) => 
   ZP_CT      (JI) = 1. / (rcpi * SUM(rho_s_1d(JI,:)) ) !inverse of the product of snow heat capacity and layer thickness [(m2 K)/J] 
   ZP_DELHEATG(JI) = 0. ! ground heat content change (diagnostic) (W/m2) JUST NEED TO DECLARE POINTER  
   ZP_DELHEATG_SFC(JI) = 0. ! ground heat content change in sfc only (diagnostic) (W/m2) JUST NEED TO DECLARE POINTER
   ZP_SW_RAD  (JI) = qsr_ice_isbaes_1d(JI) !* a_i_1d(JI) ! qtr_ice_top_1d(JI) !/ (1 - albs_isbaes_1d (JI)) !qsr_ice_isbaes_1d(JI) ! Incoming solar radiation 
   ZP_QA      (JI) = qair_isbaes_1d(JI) ! air humidity at atm. level (kg/kg) => pqair (specific Q in SI3, whatabout ISBA?)
   ZP_VMOD    (JI) = wndm_isbaes_1d(JI) ! module of the horizontal wind => wndm_ice 
   ZP_LW_RAD  (JI) = qlwdwn_ice_isbaes_1d(JI) ! NOT USED IF FLUX ARE PRESCRIBED qlw_ice_isbaes_1d(JI) 
   ZP_RHOA    (JI) = rho_air_isbaes_1d(JI) ! => rhoa => No flux computation no need
   ZP_UREF    (JI) = 10. ! atm. level for wind => zu => No flux computation no need 

   ZP_DIRCOSZW(JI) = 1. ! Cosine of the angle between the normal to the surface and the vertical => = 1 (bertrand) 
   ZP_ZREF    (JI) = 2. ! atm. level for temp. and humidity => zqt => No flux computation no nee
   ZP_PA     (JI)  = pres_temp(ZP_QA(JI), ZP_PS(JI), ZP_ZREF(JI), ptpot=tair_isbaes_1d(JI), l_ice=.true. ) ! Compute pressure at atmospheric level

   ZP_EXNS    (JI) = (ZP_PS(JI)/XP00)**(XRD/XCPD) ! Exner function at sea surface => No flux computation no need
   ZP_EXNA    (JI) = (ZP_PA(JI)/XP00)**(XRD/XCPD) ! Exner function at atm level => No flux computation no need
   ZP_TA      (JI) = tair_isbaes_1d(JI) *ZP_EXNA(JI)      ! DOIT ETRE LA TEMPERATURE ABSOLUE = multiplier tpot par fonction exner air temperature at atm. level 
   ZP_TG       (JI) = t_i_1d(JI,1) * ZP_EXNS(JI) ! Ground T° => 1st ice level 
   ZP_ALB     (JI) = albi_isbaes_1d(JI) ! green areas albedo => snow free albedo   

   ZP_RRSNOW  (JI) = rain_isbaes_1d(JI) * a_i_1d(JI) !* za_s_fra(JI) ! rain rate over snow [kg/(m2 s)] !!!!!!! MULTIPLIER PAR LA FRACTION DE NEIGE 
   ZP_SOILCOND(JI) = cnd_i_isbaes_1d(JI) ! Temporary heat conductivity of litter + soil => conductivity of 1st ice layer 
 
   !
   ZP_PEW_A_COEF(JI) = 0. ! XUNDEF !0. !ZP_VMOD(JI) ! 0. ! Coeff flux => No flux computation no need ! B COEFF 0 et A COEFF 1
   ZP_PEW_B_COEF(JI) = ZP_VMOD(JI) ! => No flux computation no need

   ZP_PET_A_COEF(JI) =  0. ! ZP_TA(JI) / (ZP_PA(JI)/XP00)**(XRD/XCPD) 
   ZP_PET_B_COEF(JI) =  ZP_TA(JI) / (ZP_PA(JI)/XP00)**(XRD/XCPD) 
   ZP_PEQ_A_COEF(JI) =  0. !ZP_QA(JI)
   ZP_PEQ_B_COEF(JI) =  ZP_QA(JI)
   !
   ZP_LAT  (JI)      = gphit_1d(JI) 
   ZP_LON  (JI)      = glamt_1d(JI) ! => glamt_1d

   ZP_ZENITH(JI)     = 0. ! solar zenith angle
!
   ZP_GRNDFLUX    (JI) = 0. !qcn_snw_bot_1d(JI) ! snow-ground flux before correction (W m-2) 
   ZP_DELHEATN    (JI) = 0. ! total snow heat content change in the surface layer (W m-2)
   ZP_DELHEATN_SFC(JI) = 0. ! total snow heat content change during the timestep (W m-2)
   ZP_SNOWSFCH    (JI) = 0. ! snow surface layer pseudo-heating term owing to changes in grid thickness (W m-2) ! DIAG => PAS BESOIN
   !!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!! KEZAKO
   ZP_LES3L       (JI) = 0. ! qla_ice_isbaes_1d(JI) ! hfx_sub_1d(JI) ! ATTENTION INVERSE AVEC SUBLIMATION !evaporation heat flux from snow (W/m2) !!!!!!!!!!!!!!!! EAU LIQUIDE QUI SEVAPORE ????????????
   ZP_LEL3L       (JI) = 0. ! ATTENTION INVERSE !! hfx_sub_1d(JI) ! sublimation (W/m2) ! A
   ZP_EVAP        (JI) = 0. !ZP_LES3L(JI) + ZP_LEL3L(JI) ! total evaporative flux (kg/m2/s) LES + LEL ?????
   !
   ZP_SWNETSNOW   (JI) = qsr_ice_isbaes_1d(JI) * (1. - albs_isbaes_1d(JI)) ! net shortwave radiation entering top of snowpack (W m-2)
   ZP_SWNETSNOWS  (JI) = qsr_ice_isbaes_1d(JI) * (1. - albs_isbaes_1d(JI)) ! net shortwave radiation in uppermost layer of snowpack
   !/!\ /!\/!\/!\/!\/!\/!\/!\/!\ SWNETSNOWS must take into account solar penetration => for now it is qsr * (1 -albedo)
   ZP_LWNETSNOW   (JI) = qlw_ice_isbaes_1d(JI) ! net longwave radiation entering top of snowpack
   ZP_MELTSTOT    (JI) = 0.0
   ZP_SNREFREEZ   (JI) = 0.0
   ZP_DELPHASEN   (JI) = 0.0 
   ZP_DELPHASEN_SFC(JI)= 0.0

ENDDO
!ZP_PA     (JJ)  = pres_temp(ZP_QA(JJ), ZP_PS(JJ), ZP_ZREF(JJ), ptpot=ZP_TA(JJ), l_ice=.true. ) ! Compute pressure at atmospheric level
!
DO JJ=1,KSIZE1
   ZP_VEGTYPE (JJ) = 1. ! fraction of permanet snow/ice => ???? ! Mettre 1 => bertrand 
   ZP_FOREST  (JJ) = 0. ! => ?????
ENDDO
!
!
! ===============================================================
! conversion of snow heat from J/m3 into J/m2
!WHERE(ZP_SNOWSWE(:,:)>0.) &
!  ZP_SNOWHEAT(:,:) = ZP_SNOWHEAT(:,:) / ZP_SNOWRHO (:,:) * ZP_SNOWSWE (:,:) 
! ===============================================================
!
ZP_PSN_INV(:)       = 0.
ZP_PSN(:)           = ZP_PSN3L(:) 
!
!
PRINT*,'ZP_SNOWRHO bef',ZP_SNOWRHO(1,:)
  CALL SNOW3L(HSNOWRES, TPTIME, OMEB, OSI3, HIMPLICIT_WIND,                   &
             ZP_PEW_A_COEF, ZP_PEW_B_COEF,                                 &
             ZP_PET_A_COEF, ZP_PEQ_A_COEF,ZP_PET_B_COEF, ZP_PEQ_B_COEF,    &
             ZP_SNOWSWE, ZP_SNOWRHO, ZP_SNOWHEAT, ZP_SNOWALB,              &
             ZP_SNOWGRAN1, ZP_SNOWGRAN2, ZP_SNOWHIST, ZP_SNOWAGE, PTSTEP,  &
             ZP_PS, ZP_SRSNOW, ZP_RRSNOW, ZP_PSN3L, ZP_TA, ZP_TG(:),     &
             ZP_SW_RAD, ZP_QA, ZP_VMOD, ZP_LW_RAD, ZP_RHOA, ZP_UREF,       &
             ZP_EXNS, ZP_EXNA, ZP_DIRCOSZW, ZP_ZREF, ZP_Z0NAT, ZP_Z0EFF,   &
             ZP_Z0HNAT, ZP_ALB, ZP_SOILCOND, ZP_D_G(:),                  &
             ZP_LVTT, ZP_LSTT, ZP_SNOWLIQ,                                 &
             ZP_SNOWTEMP, ZP_SNOWDZ, ZP_THRUFAL, ZP_MELTSTOT, ZP_SNREFREEZ,&
             ZP_GRNDFLUX, ZP_EVAPCOR, ZP_SOILCOR, ZP_GFLXCOR, ZP_SNOWSFCH, &
             ZP_DELHEATN, ZP_DELHEATN_SFC, ZP_DELPHASEN, ZP_DELPHASEN_SFC, &
             ZP_SWNETSNOW, ZP_SWNETSNOWS, ZP_LWNETSNOW, ZP_RESTOREN,       &
             ZP_RNSNOW, ZP_HSNOW, ZP_GFLUXSNOW, ZP_HPSNOW, ZP_LES3L,       &
             ZP_LEL3L, ZP_EVAP, ZP_SNDRIFT, ZP_RI,                         &
             ZP_EMISNOW, ZP_CDSNOW, ZP_USTARSNOW,                          &
             ZP_CHSNOW, ZP_SNOWHMASS, ZP_QS, ZP_RADXS, ZP_VEGTYPE,  ZP_FOREST,       &
             ZP_ZENITH, ZP_LAT, ZP_LON, OSNOWDRIFT,OSNOWDRIFT_SUBLIM)
PRINT*,'ZP_SNOWRHO aft',ZP_SNOWRHO(1,:)

! unpack variables
!
!h_s_1d(JJ) = 0.
DO JJ=1,KSIZE1
h_s_1d(JJ) = 0.
  DO JWRK=1,KSIZE2
    swe_s_1d(JJ,JWRK) = ZP_SNOWSWE  (JJ,JWRK) ! Should we advect it ???
    rho_s_1d(JJ,JWRK) = ZP_SNOWRHO  (JJ,JWRK)
    e_s_1d(JJ,JWRK) =  ZP_SNOWHEAT (JJ,JWRK) !/ (ZP_SNOWDZ   (JJ,JWRK) )!* a_i_1d(JJ))
    o_s_1d(JJ,JWRK) = ZP_SNOWAGE  (JJ,JWRK)
    t_s_1d(JJ,JWRK)   = ZP_SNOWTEMP (JJ,JWRK)
    lwc_s_1d(JJ,JWRK)   = ZP_SNOWLIQ  (JJ,JWRK) ! No need because it is a diagnostic
    dh_s_1d(JJ,JWRK)   = ZP_SNOWSWE(JJ,JWRK)/ZP_SNOWRHO(JJ,JWRK) !ZP_SNOWDZ   (JJ,JWRK) 
    !h_s_1d(JJ) = h_s_1d(JJ) + ZP_SNOWDZ   (JJ,JWRK)
    oh_s_1d(JJ,JWRK) = o_s_1d(JJ,JWRK) * dh_s_1d(JJ,JWRK)

  ENDDO
    h_s_1d(JJ) = SUM(dh_s_1d(JJ,:))
ENDDO

!CALL snw_var_enthalpy ! Compute e_s_1d from T°

!DO JWRK=1,KSIZE3
!   DO JJ=1,KSIZE1
!      JJ              = KMASK          (JJ)
!      PTG    (JJ,JWRK)= ZP_TG        (JJ,JWRK)
!   ENDDO
!ENDDO
!
DO JI = 1, npti
   zdq(JI) = - zq_ini(JI) + SUM( e_s_1d(JI,1:nlay_s))!  * dh_s_1d(JI,1:nlay_s) )
   zdm(JI) = - zm_ini(JI) + SUM(rho_s_1d(JI,1:nlay_s) * dh_s_1d(JI,1:nlay_s)) * a_i_1d(JI) 
ENDDO

DO JJ=1,KSIZE1
!  JJ                  = KMASK          (JJ)
  albs_isbaes_1d(JJ)   = ZP_SNOWALB     (JJ)
!   PEK%TSNOW%EMIS(JJ)  = ZP_EMISNOW     (JJ) ! Constante
!  DMK%XCDSNOW   (JJ)  = ZP_CDSNOW      (JJ) ! OUT
!  DMK%XUSTARSNOW(JJ)  = ZP_USTARSNOW   (JJ) ! OUT
!  DMK%XCHSNOW   (JJ)  = ZP_CHSNOW      (JJ) ! OUT
!  DMK%XSNOWHMASS(JJ)  = ZP_SNOWHMASS   (JJ) 
!  DMK%XRNSNOW   (JJ)  = ZP_RNSNOW      (JJ)
!  DMK%XHSNOW    (JJ)  = ZP_HSNOW       (JJ)
!  DMK%XHPSNOW  (JJ)   = ZP_HPSNOW      (JJ)
!  DMK%XGFLUXSNOW(JJ)  = ZP_GFLUXSNOW   (JJ)
  !
!  PDELHEATG    (JJ)   = ZP_DELHEATG    (JJ)
!  PDELHEATG_SFC(JJ)   = ZP_DELHEATG_SFC(JJ)

!  PRI          (JJ)   = ZP_RI          (JJ)
!  PQS          (JJ)   = ZP_QS          (JJ)
  qcn_snw_bot_1d(JJ)  = ZP_GRNDFLUX    (JJ) + ZP_GFLXCOR     (JJ)! Somme des flux radiatifs et convectif ??
  
  !  PFLSN_COR     (JJ)  = ZP_FLSN_COR    (JJ) ! Not used
!  PDELHEATN    (JJ)   = ZP_DELHEATN    (JJ)
!  PDELHEATN_SFC(JJ)   = ZP_DELHEATN_SFC(JJ)
!  PSNOWSFCH    (JJ)   = ZP_SNOWSFCH    (JJ) ! Check avec bertrand
!  PGSFCSNOW    (JJ)   = ZP_RESTOREN    (JJ)  !heat flux between the surface and sub-surface 
!                                                  snow layers (W/m2)
!  PLES3L       (JJ)   = ZP_LES3L       (JJ)
!  PLEL3L       (JJ)   = ZP_LEL3L       (JJ)
!  PEVAP        (JJ)   = ZP_EVAP        (JJ)
!  ZSOILCOR     (JJ)   = ZP_SOILCOR     (JJ)  ! PSOILCOR => C'est quoi la diff avec EVAPCOR => Faire somme des deux pour zevap_ram ??
  !
!  qsr_ice_1d(JJ) = ZP_SWNETSNOW   (JJ)
   !qns_ice_1d(JJ) = ZP_LWNETSNOW   (JJ)
!  ZSWNET_NS  (JJ) = ZP_SWNETSNOWS  (JJ)
!  ZLWNET_N   (JJ) = ZP_LWNETSNOW   (JJ)
   
! Heat fluxes for budget diagnostics
  ! We put 0 to all heat fluxes except from hfx_snw that we now consider as the total heat content change in the snow
!  hfx_sub_1d(JJ) = 0. !ZP_LES3L(JJ) + ZP_LEL3L(JJ) ! Sublimation + liquid water Evaporation (for now) (W/m2)
!  hfx_spr_1d(JJ) = 0. ! ZP_SNOWHMASS(JJ)! heat release from rainfall (W/m2)
!  hfx_res_1d(JJ) = 0. ! 
  hfx_snw_1d(JJ) = hfx_snw_1d(JJ) - zdq(JJ) * r1_Dt_ice !* a_i_1d(JJ)  !ZP_DELHEATN(JJ) ! total heat content change in the snow 
  zq_rema       (JJ)   = (ZP_GFLXCOR     (JJ) + ZP_RADXS (JJ)) * rDt_ice ! En J / m2
   
  ! We put 0 to all heat fluxes except from wfx_snw_sum that we now consider as the total mass change in snow
  wfx_spr_1d    (JJ)   =  wfx_spr_1d    (JJ) -(ZP_SRSNOW  (JJ) + ZP_RRSNOW(JJ) ) * a_i_1d(JJ) ! METTRE EVAP A SUB 
  wfx_snw_sub_1d   (JJ)   = wfx_snw_sub_1d   (JJ) + ZP_EVAP(JJ) * a_i_1d(JJ)
!  IF((h_s_bef(JJ) .eq. 0.) .AND. (h_s_1d(JJ) .eq. 0.)) THEN 
!     ! IF all snow is melted, all precip are considered as melted 
!     PRINT*,'WE ARE HEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEERE'
!     wfx_snw_sum_1d(JJ)   = wfx_spr_1d    (JJ)  
!  ELSE
     wfx_snw_sum_1d(JJ)   = wfx_snw_sum_1d(JJ) + ZP_THRUFAL     (JJ) * a_i_1d(JJ)! rate that liquid water leaves snow pack (kg/(m2 s)): 
                                                ! partitioned into soil infiltration/runoff by ISBA
!  ENDIF
  zevap_rema    (JJ)   = (ZP_EVAPCOR(JJ) + ZP_SOILCOR(JJ)) * rDt_ice
  ! Surface total flux => Qtot = qns_ice_1d(ji) + qsr_ice_1d(ji) - qtr_ice_top_1d(ji) - qcn_ice_top_1d(ji)                                            
  qns_ice_1d(JJ) = ZP_LES3L(JJ) + ZP_LEL3L(JJ) + ZP_HSNOW(JJ) + ZP_LWNETSNOW(JJ)
!  qcn_ice_top_1d(JJ) = ZP_RESTOREN(JJ)
!  qtr_ice_top_1d(JJ) =  ZP_SW_RAD(JJ) * ZP_SNOWALB     (JJ)

ENDDO

PRINT*,'zevap rema',zevap_rema
PRINT*,'H after snow3l', h_s_1d(1)
PRINT*,'T° after snow3l',t_s_1d(1,:)
PRINT*,'Mass difference', zdm *a_i_1d * r1_Dt_ice
! Other inout variables that we do not need (à priori??)
! PRNSNOW, ! Not needed (à priori)
! PHSNOW, ! 
! PLES3L, 
! PLEL3L, 
! PHPSNOW, PEVAP,  PGRNDFLUX, PEMISNOW

 !  ZP_GRNDFLUX    (JJ) = qcn_snw_bot_1d(JJ) ! snow-ground flux before correction (W m-2) 

!
!
END SUBROUTINE CALL_MODEL
!

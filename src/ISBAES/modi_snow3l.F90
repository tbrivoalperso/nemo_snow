!depfile:snow3l.F90
MODULE MODI_SNOW3L
INTERFACE
      SUBROUTINE SNOW3L(HSNOWRES, TPTIME, OMEB,OSI3, HIMPLICIT_WIND,           &
                PPEW_A_COEF, PPEW_B_COEF,                                 &
                PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,       &
                PSNOWSWE,PSNOWRHO,PSNOWHEAT,PSNOWALB,                     &
                PSNOWGRAN1,PSNOWGRAN2,PSNOWHIST,PSNOWAGE,                 &
                PTSTEP,PPS,PSR,PRR,PPSN3L,                                &
                PTA,PTG,PSW_RAD,PQA,PVMOD,PLW_RAD, PRHOA,                 &
                PUREF,PEXNS,PEXNA,PDIRCOSZW,                              &
                PZREF,PZ0,PZ0EFF,PZ0H,PALB,                               &
                PSOILCOND,PD_G,PLVTT,PLSTT,                               &
                PSNOWLIQ,PSNOWTEMP,PSNOWDZ,                               &
                PTHRUFAL,PSNOWMELT,PSNREFREEZ,                            &
                PGRNDFLUX,PEVAPCOR,PSOILCOR,                              &
                PGFLXCOR,PSNOWSFCH,PDELHEATN,PDELHEATN_SFC,               &
                PDELPHASEN, PDELPHASEN_SFC,                               &
                PSWNETSNOW,PSWNETSNOWS,PLWNETSNOW,PRESTOREN,              &
                PRNSNOW,PHSNOW,PGFLUXSNOW,                                &
                PHPSNOW,PLES3L,PLEL3L,PEVAP,PSNDRIFT,PRI,                 &
                PEMISNOW,PCDSNOW,PUSTAR,PCHSNOW,PSNOWHMASS,PQS,ZRADXS,    &
                PPERMSNOWFRAC,PFORESTFRAC,PZENITH,PXLAT,PXLON,            &
                HSNOWDRIFT,OSNOWDRIFT_SUBLIM                              )
                
                
USE MODD_CSTS,     ONLY : XTT, XRHOLW, XLMTT, XCL, XDAY
USE MODD_SURF_PAR, ONLY : XUNDEF
!
USE MODD_SNOW_METAMO, ONLY : XSNOWDZMIN
USE MODD_SNOW_PAR,    ONLY : XSNOWDMIN, NSPEC_BAND_SNOW
!
USE MODE_SNOW3L
USE MODE_THERMOS,  ONLY : QSATI
!
USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!                                      PTSTEP    = time step of the integration
TYPE(DATE_TIME), INTENT(IN)         :: TPTIME      ! current date and time
!
 CHARACTER(LEN=*),     INTENT(IN)    :: HSNOWRES
!                                      HSNOWRES  = ISBA-SNOW3L turbulant exchange option
!                                      'DEF' = Default: Louis (ISBA: Noilhan and Mahfouf 1996)
!                                      'RIL' = Limit Richarson number under very stable
!                                              conditions (currently testing)
!
LOGICAL, INTENT(IN)                 :: OMEB       ! True = coupled to MEB. This means surface fluxes ae IMPOSED
!                                                 ! as an upper boundary condition to the explicit snow schemes. 
!                                                 ! If = False, then energy
!                                                 ! budget and fluxes are computed herein.
LOGICAL, INTENT(IN)                 :: OSI3       ! True = coupled with SI3 
!
 CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
!                                                     ! 'OLD' = direct
!                                                     ! 'NEW' = Taylor serie, order 1
!
REAL, DIMENSION(:), INTENT(IN)    :: PPS, PTA, PSW_RAD, PQA,                       &
                                         PVMOD, PLW_RAD, PSR, PRR  
!                                      PSW_RAD = incoming solar radiation (W/m2)
!                                      PLW_RAD = atmospheric infrared radiation (W/m2)
!                                      PRR     = rain rate [kg/(m2 s)]
!                                      PSR     = snow rate (SWE) [kg/(m2 s)]
!                                      PTA     = atmospheric temperature at level za (K)
!                                      PVMOD   = modulus of the wind parallel to the orography (m/s)
!                                      PPS     = surface pressure
!                                      PQA     = atmospheric specific humidity
!                                                at level za
!
REAL, DIMENSION(:), INTENT(IN)    :: PSOILCOND, PD_G, PPSN3L
!                                      PSOILCOND = soil thermal conductivity [W/(m K)]
!                                      PD_G      = Assumed first soil layer thickness (m)
!                                                  Used to calculate ground/snow heat flux
!                                      PPSN3L    = snow fraction
!
REAL, DIMENSION(:), INTENT(IN)    :: PZREF, PUREF, PEXNS, PEXNA, PDIRCOSZW, PRHOA, PZ0, PZ0EFF, &
                                       PALB, PZ0H, PPERMSNOWFRAC, PFORESTFRAC 
!                                      PZ0EFF    = roughness length for momentum
!                                      PZ0       = grid box average roughness length
!                                      PZ0H      = grid box average roughness length for heat
!                                      PZREF     = reference height of the first
!                                                  atmospheric level
!                                      PUREF     = reference height of the wind
!                                      PRHOA     = air density
!                                      PEXNS     = Exner function at surface
!                                      PEXNA     = Exner function at lowest atmos level
!                                      PDIRCOSZW = Cosinus of the angle between the
!                                                  normal to the surface and the vertical
!                                      PALB      = soil/vegetation albedo
!                                      PPERMSNOWFRAC  = fraction of permanet snow/ice
!                                      PFORESTFRAC = fraction of forest
!
REAL, DIMENSION(:), INTENT(IN)      :: PPEW_A_COEF, PPEW_B_COEF,                   &
                                         PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF,      &
                                         PPEQ_B_COEF  
!                                      PPEW_A_COEF = wind coefficient (m2s/kg)
!                                      PPEW_B_COEF = wind coefficient (m/s)
!                                      PPET_A_COEF = A-air temperature coefficient
!                                      PPET_B_COEF = B-air temperature coefficient
!                                      PPEQ_A_COEF = A-air specific humidity coefficient
!                                      PPEQ_B_COEF = B-air specific humidity coefficient
!
REAL, DIMENSION(:), INTENT(IN)    :: PTG
!                                      PTG       = Surface soil temperature (effective
!                                                  temperature the of layer lying below snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLVTT, PLSTT ! = latent heats for hydrology
REAL, DIMENSION(:), INTENT(INOUT) :: PSNOWALB
!                                      PSNOWALB = Prognostic surface snow albedo
!                                                 (does not include anything but
!                                                 the actual snow cover)
!
REAL, DIMENSION(:,:), INTENT(INOUT):: PSNOWHEAT, PSNOWRHO, PSNOWSWE
!                                      PSNOWHEAT = Snow layer(s) heat content (J/m2)
!                                      PSNOWRHO  = Snow layer(s) averaged density (kg/m3)
!                                      PSNOWSWE  = Snow layer(s) liquid Water Equivalent (SWE:kg m-2)
!
REAL, DIMENSION(:,:), INTENT(INOUT):: PSNOWGRAN1, PSNOWGRAN2, PSNOWHIST
!                                      PSNOWGRAN1 = Snow layers grain feature 1
!                                      PSNOWGRAN2 = Snow layer grain feature 2
!                                      PSNOWHIST  = Snow layer grain historical
!                                                   parameter (only for non
!                                                   dendritic snow)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWAGE  ! Snow grain age
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PRNSNOW, PHSNOW, PLES3L, PLEL3L, &
                                       PHPSNOW, PEVAP,  PGRNDFLUX, PEMISNOW
!                                      PLES3L      = evaporation heat flux from snow (W/m2)
!                                      PLEL3L      = sublimation (W/m2)
!                                      PHPSNOW     = heat release from rainfall (W/m2)
!                                      PRNSNOW     = net radiative flux from snow (W/m2)
!                                      PHSNOW      = sensible heat flux from snow (W/m2)
!                                      PHSNOW      = sensible heat flux from snow (W/m2)
!                                      PEVAP       = total evaporative flux (kg/m2/s)
!                                      PGRNDFLUX   = soil/snow interface heat flux (W/m2)
!                                      PEMISNOW    = snow surface emissivity
!
REAL, DIMENSION(:), INTENT(OUT)     :: PGFLUXSNOW
!                                      PGFLUXSNOW  = net heat flux from snow (W/m2)
!
REAL, DIMENSION(:), INTENT(INOUT) :: PSWNETSNOW, PLWNETSNOW, PSWNETSNOWS
!                                      PSWNETSNOW = net shortwave radiation entering top of snowpack 
!                                                  (W m-2) Imposed if MEB=T, diagnosed herein if MEB=F
!                                      PSWNETSNOWS= net shortwave radiation in uppermost layer of snowpack 
!                                                  (W m-2) Imposed if MEB=T, diagnosed herein if MEB=F
!                                                   Used for surface energy budget diagnostics
!                                      PLWNETSNOW = net longwave radiation entering top of snowpack 
!                                                  (W m-2) Imposed if MEB=T, diagnosed herein if MEB=F
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PUSTAR, PCDSNOW, PCHSNOW, PRI
!                                      PCDSNOW    = drag coefficient for momentum over snow (-)
!                                      PUSTAR     = friction velocity over snow (m/s)
!                                      PCHSNOW    = drag coefficient for heat over snow (-)
!                                      PRI        = Richardson number (-)
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWTEMP
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWLIQ, PSNOWDZ
!                                      PSNOWLIQ  = Snow layer(s) liquid water content (m)
!                                      PSNOWTEMP = Snow layer(s) temperature (m)
!                                      PSNOWDZ   = Snow layer(s) thickness (m)
!
REAL, DIMENSION(:), INTENT(OUT)     :: PTHRUFAL, PEVAPCOR, PSOILCOR, PGFLXCOR, &
                                       PRESTOREN, PSNOWSFCH, PDELHEATN, PDELHEATN_SFC, &
                                       PDELPHASEN, PDELPHASEN_SFC, PSNOWMELT, PSNREFREEZ
!                                      PTHRUFAL  = rate that liquid water leaves snowpack :
!                                                  paritioned into soil infiltration/runoff
!                                                  by ISBA [kg/(m2 s)]
!                                      PEVAPCOR  = evaporation/sublimation correction term:
!                                                  extract any evaporation exceeding the
!                                                  actual snow cover (as snow vanishes)
!                                                  and apply it as a surface soil water
!                                                  sink. [kg/(m2 s)]
!                                      PSOILCOR = for vanishingy thin snow cover,
!                                                 allow any excess evaporation
!                                                 to be extracted from the soil
!                                                 to maintain an accurate water
!                                                 balance [kg/(m2 s)]
!                                      PGFLXCOR  = flux correction to underlying soil for vanishing snowpack
!                                                  (to put any energy excess from snow to soil) (W/m2)
!                                      PRESTOREN = heat flux between the surface and sub-surface 
!                                                  snow layers (W/m2)
!                                      PSNOWSFCH = snow surface layer pseudo-heating term owing to
!                                                  changes in grid thickness            (W m-2)
!                                      PDELHEATN = total snow heat content change  (W m-2)
!                                      PDELHEATN_SFC = surface layer snow heat content change during the timestep (W m-2)
!                                      PDELPHASEN = latent heating due to snow melt/freeze  (W m-2)
!                                      PDELPHASEN_SFC = latent heating due to surface layer snow melt/freeze  (W m-2)
!                                      PSNOWMELT = snowmelt in the snowpack  [kg/(m2 s)]
!                                      PSNREFREEZ = refreezing of water in the snowpack  [kg/(m2 s)]
!
REAL, DIMENSION(:), INTENT(OUT)     :: PSNDRIFT
!                                      PSNDRIFT    = blowing snow sublimation (kg/m2/s)
!
REAL, DIMENSION(:), INTENT(OUT)   ::   PSNOWHMASS
!                                      PSNOWHMASS  = heat content change due to mass
!                                                    changes in snowpack (J/m2): for budget
!                                                    calculations only.
!
REAL, DIMENSION(:), INTENT(OUT)   :: PQS
!                                    PQS = surface humidity
REAL, DIMENSION(:), INTENT(OUT)   :: ZRADXS
!                                      ZRADXS   = shortwave radiation absorbed by soil surface
!                                                 (for thin snow sover) (W m-2)
!
REAL, DIMENSION(:), INTENT(IN)    :: PZENITH ! solar zenith angle
REAL, DIMENSION(:), INTENT(IN)    :: PXLAT,PXLON ! LAT/LON after packing
!
CHARACTER(4), INTENT(IN)          :: HSNOWDRIFT  ! Snowdrift scheme :
                                                 ! 'NONE': No snowdrift scheme
                                                 !  'DFLT':  Snowdrift scheme activated
                                                 !  Other options are available in Crocus

LOGICAL, INTENT(IN)               ::  OSNOWDRIFT_SUBLIM ! activate snowdrift, sublimation during drift


END SUBROUTINE SNOW3L

END INTERFACE
END MODULE MODI_SNOW3L

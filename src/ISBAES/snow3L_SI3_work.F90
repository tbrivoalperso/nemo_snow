!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
!================================================================
SUBROUTINE SNOW3L_SI3(KSIZE1,KSIZE2,KSIZE3,PTSTEP, za_s_fra, isnow, ZP_RADXS, zq_rema, zevap_rema) 
!
USE MODD_CSTS !,       ONLY : XLMTT, XLSTT
USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_SNOW_PAR,   ONLY : XRHOSMAX_ES, XSNOWDMIN, XRHOSMIN_ES, XEMISSN
USE ice1D ! SI3 1D variables
USE ice   ! SI3 variables 
IMPLICIT NONE
!
REAL, INTENT(IN)                    :: PTSTEP
!                                      PTSTEP    = time step of the integration
INTEGER, INTENT(IN) :: KSIZE1
INTEGER, INTENT(IN) :: KSIZE2
INTEGER, INTENT(IN) :: KSIZE3
REAL, DIMENSION(KSIZE1), INTENT(IN) ::   za_s_fra    ! ice fraction covered by snow
REAL, DIMENSION(KSIZE1), INTENT(out) ::   isnow       ! snow presence (1) or not (0)
REAL, DIMENSION(KSIZE1), INTENT(out) ::   ZP_RADXS    ! Radiation transmited through the snow
REAL, DIMENSION(KSIZE1), INTENT(out) ::   zq_rema     ! remaining heat flux from snow melting       (J.m-2)
REAL, DIMENSION(KSIZE1), INTENT(out) ::   zevap_rema  ! remaining mass flux from snow sublimation   (kg.m-2)

!*      0.2    declarations of local variables
!
REAL, PARAMETER                     :: ZCHECK_TEMP = 50.0
!                                      Limit to check suspicious low temperature (K)
!
INTEGER                             :: JWRK, JJ ! Loop control
!
INTEGER                             :: INLVLS   ! maximum number of snow layers
INTEGER                             :: INLVLG   ! number of ground layers
!
REAL, DIMENSION(KSIZE1)          :: ZRRSNOW, ZSOILCOND, ZSNOW, ZSNOWFALL,  &
                                       ZSNOWABLAT_DELTA, ZSNOWSWE_1D, ZSNOWD, &
                                       ZSNOWH, ZSNOWH1, ZGRNDFLUXN, ZPSN,     &
                                       ZSOILCOR, ZSNOWSWE_OUT, ZTHRUFAL,      &
                                       ZSNOW_MASS_BUDGET, ZWGHT, ZWORK, ZC2
!                                      ZSOILCOND    = soil thermal conductivity [W/(m K)]
!                                      ZRRSNOW      = rain rate over snow [kg/(m2 s)]
!                                      ZSNOW        = snow depth (m) 
!                                      ZSNOWFALL    = minimum equivalent snow depth
!                                                     for snow falling during the
!                                                     current time step (m)
!                                      ZSNOWABLAT_DELTA = FLAG =1 if snow ablates completely
!                                                     during current time step, else=0
!                                      ZSNOWSWE_1D  = TOTAL snowpack SWE (kg m-2)
!                                      ZSNOWD       = snow depth
!                                      ZSNOWH       = snow total heat content (J m-2)
!                                      ZSNOWH1      = snow surface layer heat content (J m-2)
!                                      ZGRNDFLUXN   = corrected snow-ground flux (if snow fully ablated during timestep)
!                                      ZPSN         = snow fraction working array
!                                      ZSOILCOR = for vanishingy thin snow cover,
!                                                 allow any excess evaporation
!                                                 to be extracted from the soil
!                                                 to maintain an accurate water
!                                                 balance [kg/(m2 s)]
!                                      ZSNOW_MASS_BUDGET = snow water equivalent budget (kg/m2/s)
!                                      ZWGHT        = MEB surface layer weight for distributing energy
!                                                     between litter and ground layers for the case
!                                                     of total ablation during a timestep (-).
!                                      ZWORK        = local working variable (*)
!                                      ZC2          = sub-surface heat capacity [(K m2)/J]

INTEGER                            :: ISIZE_SNOW ! number of points where computations are done
INTEGER, DIMENSION(KSIZE1)      :: NMASK      ! indices correspondance between arrays
!


! - Snow and rain falling onto the 3-L grid space:
!
!  DMK%XSRSFC(:) = 0.0
!
!snow_isbaes_1d(:) = 0.
  DO JJ=1,KSIZE1
    ZRRSNOW(JJ)        = rain_isbaes_1d(JJ) * za_s_fra(JJ) 
!    DMK%XRRSFC(JJ)    = PRR(JJ) - ZRRSNOW(JJ)
    ZSNOWFALL(JJ)      = snow_isbaes_1d(JJ)*PTSTEP/XRHOSMAX_ES    ! maximum possible snowfall depth (m)
  ENDDO
!
! Calculate preliminary snow depth (m)

  ZSNOW(:)      =0.
  ZSNOWH(:)     =0.
!  ZSNOWSWE_1D(:)=0.
!  ZSNOWH1(:)    = PEK%TSNOW%HEAT(:,1)*PEK%TSNOW%WSNOW(:,1)/PEK%TSNOW%RHO(:,1) ! sfc layer only
!
    DO JJ=1,KSIZE1
!      ZSNOWSWE_1D(JJ) = ZSNOWSWE_1D(JJ) + rho_s_1d(JJ,JWRK) * h_s_1d(JJ) * r1_nlay_s 
      ZSNOW(JJ)       = ZSNOW(JJ)       + h_s_1d(JJ) 
    END DO

  DO JWRK=1,KSIZE2
    DO JJ=1,KSIZE1
      ZSNOWH(JJ)      = ZSNOWH(JJ)      + e_s_1d(JJ,JWRK)*h_s_1d(JJ)
    END DO
  ENDDO
!
! ===============================================================
! === Packing: Only call snow model when there is snow on the surface
!              exceeding a minimum threshold OR if the equivalent
!              snow depth falling during the current time step exceeds 
!              this limit.
!
! counts the number of points where the computations will be made
!
!
  ISIZE_SNOW = 0
  NMASK(:) = 0
!
  DO JJ=1,KSIZE1
    IF (ZSNOW(JJ) >= XSNOWDMIN .OR. ZSNOWFALL(JJ) >= XSNOWDMIN) THEN
      isnow(JJ) = 1.
      ISIZE_SNOW = ISIZE_SNOW + 1
      NMASK(ISIZE_SNOW) = JJ
    ELSE
      isnow(JJ) = 0.      
    ENDIF
  ENDDO
  IF (ISIZE_SNOW>0) CALL CALL_MODEL(KSIZE1,KSIZE2,KSIZE3,NMASK,PTSTEP, za_s_fra, ZP_RADXS, zq_rema, zevap_rema) 
!
END SUBROUTINE SNOW3L_SI3
!

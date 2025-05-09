#include "SEAICE_OPTIONS.h"
#ifdef ALLOW_EXF
# include "EXF_OPTIONS.h"
#endif
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

CBOP
C !ROUTINE: SEAICE_MODEL

C !INTERFACE: ==========================================================
      SUBROUTINE SEAICE_MODEL( myTime, myIter, myThid )

C !DESCRIPTION: \bv
C     *===========================================================*
C     | SUBROUTINE SEAICE_MODEL                                   |
C     | o Time stepping of a dynamic/thermodynamic sea ice model. |
C     |  Dynamics solver: Zhang/Hibler, JGR, 102, 8691-8702, 1997 |
C     |  Thermodynamics:        Hibler, MWR, 108, 1943-1973, 1980 |
C     |  Rheology:              Hibler, JPO,   9,  815- 846, 1979 |
C     |  Snow:          Zhang et al.  , JPO,  28,  191- 217, 1998 |
C     |  Parallel forward ice model written by Jinlun Zhang PSC/UW|
C     |  & coupled into MITgcm by Dimitris Menemenlis (JPL) 2/2001|
C     |  zhang@apl.washington.edu / menemenlis@jpl.nasa.gov       |
C     | o The code has been rewritten substantially to use the    |
C     |   MITgcm C-grid, see Losch et al. OM, 33,  129- 144, 2010 |
C     *===========================================================*
      IMPLICIT NONE
C \ev

C !USES: ===============================================================
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DYNVARS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "SEAICE_SIZE.h"
#include "SEAICE_PARAMS.h"
#include "SEAICE.h"
#include "SEAICE_TRACER.h"
#ifdef ALLOW_EXF
# include "EXF_FIELDS.h"
#endif
#ifdef ALLOW_AUTODIFF_TAMC
# include "tamc.h"
#endif

C !INPUT PARAMETERS: ===================================================
C     myTime :: Current time in simulation
C     myIter :: Current iteration number in simulation
C     myThid :: my Thread Id number
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

C !LOCAL VARIABLES: ====================================================
C     i,j,bi,bj :: Loop counters
      INTEGER i, j
      INTEGER bi, bj
#ifdef ALLOW_EXF
      INTEGER grpDiag
#endif
#ifdef ALLOW_SITRACER
      INTEGER iTr
#endif
#ifdef SEAICE_BGRID_DYNAMICS
      _RL uLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#endif
CEOP
#ifdef ALLOW_DIAGNOSTICS
      LOGICAL  DIAGNOSTICS_IS_ON
      EXTERNAL DIAGNOSTICS_IS_ON
      LOGICAL  diag_SIenph_isOn
      _RL mTmp(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL recip_deltaTtherm
      _RL addSnow
#endif
#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_ENTER( 'SEAICE_MODEL', myThid )
#endif

      IF ( usingPCoords ) THEN
C     In z-coordinates, phiHydLow is just a diagnostics and hence its
C     overlaps are never filled properly
C     in p-coordinates, phiHydLow is the sea level elevation that is required
C     for the sea level tilt term in the momentum equation, so we need
C     to fill the overlaps properly here
C     (or elsewhere, not sure where it would be best)
       _EXCH_XY_RL( phiHydLow, myThid )
      ENDIF
#ifdef ALLOW_THSICE
      IF ( useThSice ) THEN
C--   Map thSice-variables to HEFF, HSNOW, and AREA, because they are
C     needed in S/R SEAICE_DYNSOLVER
       CALL SEAICE_MAP_THSICE( myTime, myIter, myThid )
      ENDIF
#endif /* ALLOW_THSICE */
#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
       diag_SIenph_isOn = DIAGNOSTICS_IS_ON('SIenph  ',myThid)
       IF ( diag_SIenph_isOn ) THEN
        IF ( SEAICEaddSnowMass ) THEN
         addSnow=1. _d 0
        ELSE
         addSnow=0. _d 0
        ENDIF
        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            mTmp(i,j,bi,bj)=SEAICE_rhoIce*HEFF(i,j,bi,bj)
     &                     +addSnow*SEAICE_rhoSnow*HSNOW(i,j,bi,bj)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDIF
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#ifdef ALLOW_EXF
      IF ( useEXF ) THEN
C--   Winds are from pkg/exf, which does not update edges.
       CALL EXCH_UV_AGRID_3D_RL( uwind, vwind, .TRUE., 1, myThid )
       IF ( useDiagnostics ) THEN
C-    Fill-in EXF wind-stess diags, weighted by open-ocean fraction
        grpDiag = -1
        IF ( SEAICEuseDYNAMICS ) grpDiag = 1
        CALL EXF_WEIGHT_SFX_DIAGS(
     I                  AREA, grpDiag, myTime, myIter, myThid )
       ENDIF
      ENDIF
#endif /* ALLOW_EXF */

#ifdef ALLOW_AUTODIFF
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          uIceNm1(i,j,bi,bj) = 0. _d 0
          vIceNm1(i,j,bi,bj) = 0. _d 0
# ifdef ALLOW_SITRACER
          DO iTr = 1, SItrMaxNum
           SItrBucket(i,j,bi,bj,iTr) = 0. _d 0
          ENDDO
# endif
         ENDDO
        ENDDO
       ENDDO
      ENDDO
#endif /* ALLOW_AUTODIFF */
#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE heff  = comlev1, key=ikey_dynamics, kind=isbyte
CADJ STORE area  = comlev1, key=ikey_dynamics, kind=isbyte
CADJ STORE hsnow = comlev1, key=ikey_dynamics, kind=isbyte
# ifdef SEAICE_CGRID
CADJ STORE fu    = comlev1, key=ikey_dynamics, kind=isbyte
CADJ STORE fv    = comlev1, key=ikey_dynamics, kind=isbyte
#  ifdef SEAICE_ALLOW_EVP
CADJ STORE seaice_sigma1  = comlev1, key=ikey_dynamics, kind=isbyte
CADJ STORE seaice_sigma2  = comlev1, key=ikey_dynamics, kind=isbyte
CADJ STORE seaice_sigma12 = comlev1, key=ikey_dynamics, kind=isbyte
#  endif /* SEAICE_ALLOW_EVP */
# endif /* SEAICE_CGRID */
# ifdef ALLOW_SITRACER
CADJ STORE siceload = comlev1, key=ikey_dynamics, kind=isbyte
CADJ STORE sitracer = comlev1, key=ikey_dynamics, kind=isbyte
# endif
#endif /* ALLOW_AUTODIFF_TAMC */

C--   Solve ice momentum equations and calculate ocean surface stress.
C     The surface stress always needs to be updated, even if neither B-
C     or C-grid dynamics are compiled, and SEAICEuseDYNAMICS = F.
#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_CALL( 'SEAICE_DYNSOLVER', myThid )
#endif
#ifdef SEAICE_BGRID_DYNAMICS
      CALL TIMER_START('DYNSOLVER          [SEAICE_MODEL]',myThid)
      CALL DYNSOLVER ( myTime, myIter, myThid )
      CALL TIMER_STOP ('DYNSOLVER          [SEAICE_MODEL]',myThid)
#else /* use default C-grid solver */
      CALL TIMER_START('SEAICE_DYNSOLVER   [SEAICE_MODEL]',myThid)
      CALL SEAICE_DYNSOLVER ( myTime, myIter, myThid )
      CALL TIMER_STOP ('SEAICE_DYNSOLVER   [SEAICE_MODEL]',myThid)
#endif

C--   Apply ice velocity open boundary conditions
#ifdef ALLOW_OBCS
# ifndef DISABLE_SEAICE_OBCS
      IF ( useOBCS ) CALL OBCS_ADJUST_UVICE( uice, vice, myThid )
# endif /* DISABLE_SEAICE_OBCS */
#endif /* ALLOW_OBCS */

#if ( defined ALLOW_AUTODIFF_TAMC &&  defined SEAICE_CGRID )
CADJ STORE uice  = comlev1, key=ikey_dynamics, kind=isbyte
CADJ STORE vice  = comlev1, key=ikey_dynamics, kind=isbyte
C     Note: Storing u/vice **after** seaice_dynsolver (and obcs_adjust_uvice)
C     has the effect that seaice_dynsolver is not called from seaice_model_ad
C     anymore. This is important because with the numerous tricks to avoid
C     complicated code in the backward integration (see pkg/autodiff),
C     the extra call would update the ice velocities with the wrong
C     set of flag values.
#endif

#ifdef ALLOW_THSICE
      IF ( useThSice ) THEN
#ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL( 'THSICE_DO_ADVECT', myThid )
#endif
       CALL THSICE_DO_ADVECT( 0, 0, myTime, myIter, myThid )
      ELSE
#endif /* ALLOW_THSICE */
C--   Only call advection of heff, area, snow, and salt and
C--   growth for the generic 0-layer thermodynamics of seaice
C--   if useThSice=.false., otherwise the 3-layer Winton thermodynamics
C--   (called from DO_OCEANIC_PHYSICS) take care of this
C NOW DO ADVECTION and DIFFUSION
       IF ( SEAICEadvHeff .OR. SEAICEadvArea .OR. SEAICEadvSnow
     &        .OR. SEAICEadvSalt ) THEN
#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_CALL( 'SEAICE_ADVDIFF', myThid )
#endif
C--   There always needs to be advection, even if neither B- or C-grid
C     dynamics are compiled.
#ifdef SEAICE_BGRID_DYNAMICS
        CALL SEAICE_ADVDIFF( uLoc, vLoc, myTime, myIter, myThid )
#else /* default C-grid advection */
        CALL SEAICE_ADVDIFF( uIce, vIce, myTime, myIter, myThid )
#endif
       ENDIF
#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE heff  = comlev1, key=ikey_dynamics, kind=isbyte
CADJ STORE area  = comlev1, key=ikey_dynamics, kind=isbyte
CADJ STORE hsnow = comlev1, key=ikey_dynamics, kind=isbyte
# ifdef SEAICE_VARIABLE_SALINITY
CADJ STORE hsalt = comlev1, key=ikey_dynamics, kind=isbyte
# endif
C     Note: This store has the effect that seaice_advdiff is not called
C     from seaice_model_ad anymore. Instead, the stored values are used.
#endif /* ALLOW_AUTODIFF_TAMC */

C     After advection, the sea ice variables may have unphysical values
C     e.g., < 0, that are regularized here. Concentration as a special case
C     may be > 1 in convergent motion and a ridging algorithm redistributes
C     the ice to limit the concentration to 1.
       CALL SEAICE_REG_RIDGE( myTime, myIter, myThid )

#ifdef ALLOW_EXF
       IF ( useEXF .AND. useDiagnostics ) THEN
C-    Fill-in EXF surface flux diags, weighted by open-ocean fraction
        grpDiag = -2
        IF ( usePW79thermodynamics ) grpDiag = 2
        CALL EXF_WEIGHT_SFX_DIAGS(
     I                  AREA, grpDiag, myTime, myIter, myThid )
       ENDIF
#endif /* ALLOW_EXF */

#ifdef DISABLE_SEAICE_GROWTH
       IF ( .TRUE. ) THEN
#else /* DISABLE_SEAICE_GROWTH */
C     thermodynamics growth
C     must call growth after calling advection
C     because of ugly time level business
       IF ( usePW79thermodynamics ) THEN
# ifdef SEAICE_USE_GROWTH_ADX
#  ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_CALL( 'SEAICE_GROWTH_ADX', myThid )
#  endif
        CALL SEAICE_GROWTH_ADX( myTime, myIter, myThid )
# else /* SEAICE_USE_GROWTH_ADX */
#  ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_CALL( 'SEAICE_GROWTH', myThid )
#  endif
        CALL SEAICE_GROWTH( myTime, myIter, myThid )
# endif /* SEAICE_USE_GROWTH_ADX */
       ELSE
#endif /* DISABLE_SEAICE_GROWTH */
        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
          DO j=1,sNy
           DO i=1,sNx
            sIceLoad(i,j,bi,bj) = HEFF(i,j,bi,bj)*SEAICE_rhoIce
     &                         + HSNOW(i,j,bi,bj)*SEAICE_rhoSnow
           ENDDO
          ENDDO
c#ifdef SEAICE_CAP_ICELOAD
c         sIceTooHeavy = rhoConst*drF(1) / 5. _d 0
c         DO j=1,sNy
c          DO i=1,sNx
c           sIceLoad(i,j,bi,bj) = MIN( sIceLoad(i,j,bi,bj),
c    &                                 sIceTooHeavy )
c          ENDDO
c         ENDDO
c#endif
         ENDDO
        ENDDO
       ENDIF

#ifdef ALLOW_SITRACER
# ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE sitracer  = comlev1, key=ikey_dynamics, kind=isbyte
# endif
       CALL SEAICE_TRACER_PHYS ( myTime, myIter, myThid )
#endif

C--   Apply ice tracer open boundary conditions
#ifdef ALLOW_OBCS
# ifndef DISABLE_SEAICE_OBCS
       IF ( useOBCS ) CALL OBCS_APPLY_SEAICE( myThid )
# endif /* DISABLE_SEAICE_OBCS */
#endif /* ALLOW_OBCS */

C--   Update overlap regions for a bunch of stuff
       _EXCH_XY_RL( HEFF,  myThid )
       _EXCH_XY_RL( AREA,  myThid )
       _EXCH_XY_RL( HSNOW, myThid )
#ifdef SEAICE_ITD
       CALL EXCH_3D_RL( HEFFITD,  nITD, myThid )
       CALL EXCH_3D_RL( AREAITD,  nITD, myThid )
       CALL EXCH_3D_RL( HSNOWITD, nITD, myThid )
#endif
#ifdef SEAICE_VARIABLE_SALINITY
       _EXCH_XY_RL( HSALT, myThid )
#endif
#ifdef ALLOW_SITRACER
       DO iTr = 1, SItrNumInUse
        _EXCH_XY_RL( SItracer(1-OLx,1-OLy,1,1,iTr),myThid )
       ENDDO
#endif
       _EXCH_XY_RS(EmPmR, myThid )
       _EXCH_XY_RS(saltFlux, myThid )
       _EXCH_XY_RS(Qnet , myThid )
#ifdef SHORTWAVE_HEATING
       _EXCH_XY_RS(Qsw  , myThid )
#endif /* SHORTWAVE_HEATING */
#ifdef ATMOSPHERIC_LOADING
       IF ( useRealFreshWaterFlux )
     &      _EXCH_XY_RS( sIceLoad, myThid )
#endif

#ifdef ALLOW_OBCS
C--   In case we use scheme with a large stencil that extends into overlap:
C     no longer needed with the right masking in advection & diffusion S/R.
c      IF ( useOBCS ) THEN
c       DO bj=myByLo(myThid),myByHi(myThid)
c        DO bi=myBxLo(myThid),myBxHi(myThid)
c          CALL OBCS_COPY_TRACER( HEFF(1-OLx,1-OLy,bi,bj),
c    I                            1, bi, bj, myThid )
c          CALL OBCS_COPY_TRACER( AREA(1-OLx,1-OLy,bi,bj),
c    I                            1, bi, bj, myThid )
c          CALL OBCS_COPY_TRACER( HSNOW(1-OLx,1-OLy,bi,bj),
c    I                            1, bi, bj, myThid )
#ifdef SEAICE_VARIABLE_SALINITY
c          CALL OBCS_COPY_TRACER( HSALT(1-OLx,1-OLy,bi,bj),
c    I                            1, bi, bj, myThid )
#endif
c        ENDDO
c       ENDDO
c      ENDIF
#endif /* ALLOW_OBCS */

#ifdef ALLOW_DIAGNOSTICS
       IF ( useDiagnostics ) THEN
C     diagnostics for "non-state variables" that are modified by
C     the seaice model ...
        CALL DIAGNOSTICS_FILL_RS(EmPmR,'SIempmr ',0,1 ,0,1,1,myThid)
        CALL DIAGNOSTICS_FILL_RS(Qnet ,'SIqnet  ',0,1 ,0,1,1,myThid)
        CALL DIAGNOSTICS_FILL_RS(Qsw  ,'SIqsw   ',0,1 ,0,1,1,myThid)
C     ... and energy diagnostic
        IF ( diag_SIenph_isOn ) THEN
         recip_deltaTtherm = 1. _d 0 / SEAICE_deltaTtherm
         DO bj=myByLo(myThid),myByHi(myThid)
          DO bi=myBxLo(myThid),myBxHi(myThid)
           DO j=1,sNy
            DO i=1,sNx
             mTmp(i,j,bi,bj) = 0.25 _d 0 * (
     &              uIce (i,j,bi,bj)**2 + uIce (i+1,j,  bi,bj)**2
     &            + vIce (i,j,bi,bj)**2 + vIce (i,  j+1,bi,bj)**2 )
     &            * (SEAICE_rhoIce*HEFF(i,j,bi,bj)
     &            + addSnow*SEAICE_rhoSnow*HSNOW(i,j,bi,bj)
     &            - mTmp(i,j,bi,bj))
     &             *recip_deltaTtherm
            ENDDO
           ENDDO
          ENDDO
         ENDDO
         CALL DIAGNOSTICS_FILL(mTmp,'SIenph  ',0,1,0,1,1,myThid)
        ENDIF
       ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#ifdef ALLOW_THSICE
C     endif .not.useThSice
      ENDIF
#endif /* ALLOW_THSICE */
CML   This has already been done in seaice_ocean_stress/ostres, so why repeat?
CML   CALL EXCH_UV_XY_RS(fu,fv,.TRUE.,myThid)

#ifdef ALLOW_EXF
# ifdef ALLOW_AUTODIFF
#  ifdef ALLOW_AUTODIFF_MONITOR
      CALL EXF_ADJOINT_SNAPSHOTS( 3, myTime, myIter, myThid )
#  endif
# endif
#endif

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_LEAVE( 'SEAICE_MODEL', myThid )
#endif

      RETURN
      END

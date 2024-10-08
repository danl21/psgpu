! -*- mode:F90 -*-
!  
      module GLOBAL_PARAMS
      IMPLICIT NONE
!
! ------------------------------------------------------------------------
!     PARAMETERS
! ------------------------------------------------------------------------
!
      INTEGER, PARAMETER  :: rk = 8       ! floating point precision      
!
      REAL(rk), PARAMETER :: alpha = 1.0_rk
!
      INTEGER, PARAMETER  :: NY = 256
      INTEGER, PARAMETER  :: NX = INT(NY/alpha)
      INTEGER, PARAMETER ::  NX2=NX/2+1
!
      INTEGER, PARAMETER  :: KTX = NX/3
      INTEGER, PARAMETER  :: KTY = NY/3
      INTEGER, PARAMETER  :: IKTX = KTX+1
      INTEGER, PARAMETER  :: IKTY = KTY+KTY+1
!
!$$$      INTEGER, PARAMETER  :: STORESIZE = 1000 !CHANGE 333 format() too
!$$$      INTEGER, PARAMETER  :: NTIMES = 100000 !CHANGE 333 format() too
!
      INTEGER, PARAMETER  :: GMRESMAX = 2000 
      INTEGER, PARAMETER  :: NEWTONMAX = 100
      INTEGER, PARAMETER  :: LWORKsvd = 7*GMRESMAX + 4*GMRESMAX*GMRESMAX
      INTEGER, PARAMETER  :: LRWORKsvd= 7*GMRESMAX + 5*GMRESMAX*GMRESMAX
!
      INTEGER, PARAMETER  :: FFTW_EXHAUSTIVE = 8       !see fftw3.f
      INTEGER, PARAMETER  :: FFTW_PATIENT = 32         !see fftw3.f
      INTEGER, PARAMETER  :: FFTW_MEASURE = 0          !see fftw3.f
      INTEGER, PARAMETER  :: FFTW_ESTIMATE = 64        !see fftw3.f
      INTEGER, PARAMETER  :: FFTW_FLAG = FFTW_ESTIMATE

!
!  -----------------------------------------------------------------------
!  CONSTANTS
!
      COMPLEX(rk) :: ZI 
      REAL(rk)    :: TWOPI 
      REAL(rk)    :: ROOT2 
!
! ------------------------------------------------------------------------
!  MASK
! 
      INTEGER, DIMENSION(IKTX,IKTY) :: L  ! stencil for wavenumber plane
      LOGICAL, DIMENSION(IKTX,IKTY) :: LL  ! stencil for wavenumber plane
!
! ---------------------------------------------------------------------------
!  FFTs
!
      INTEGER(KIND=8) :: PlanK2ZR,PlanK2UR,PlanK2VR,PlanNZR2K,PlanUR2K
!
! ----------------------------------------------------------------------------
!  PHYSICAL PARAMETERS
!
      REAL(rk), DIMENSION(IKTX,IKTY) :: ROSSBY,NU
      REAL(rk) ::  V1,V2,AMPFOR,TSTEP,DELT,BETA,AMPV,ROBERT
      INTEGER  ::  KF,NSTOP,NOUT,NOUTV
!
!--------------------------------------------------------------------------------
!  FLAGS
!
      INTEGER :: myISEED,ISEED,GRFLAG,TRFLAG,RCFLAG,STFLAG,ICFLAG,IREST
      INTEGER :: WNFLAG,SubSpaceFLAG,lamFLAG,UPOflag,RsymFLAG
!
!--------------------------------------------------------------------------------
!  GMRES PARAMS
!
      INTEGER :: IGUESSEND
      REAL(rk) :: maxPERIOD,minPERIOD,gmresResT,SNEWTX,SNEWTY, timeUNDER
!
! --------------------------------------------------------------------------
!  COUNTERS
!
      INTEGER ::  izkout
!
!----------------------------------------------------------------------------
! TTGIRD
!
      REAL(rk)    :: ResidualThreshold
      INTEGER     :: NSHIFT
!
!$$$      INTEGER :: POINTER,ITIME
!$$$!
!$$$      REAL(KIND=4), DIMENSION(STORESIZE,NTIMES) :: TtGRID_BIGT,TtGRID_DZ
!$$$      REAL(rk),     DIMENSION(NTIMES)        :: TtGRID_TIME
!$$$      REAL(rk),     DIMENSION(STORESIZE)        :: TIMESTORE
!$$$!
!$$$      COMPLEX(rk) :: ZSTORE(IKTX,IKTY,STORESIZE)
!
! ------------------------------------------------------------------------------
!  SPECTRA
!
      INTEGER      :: NS(KTY)
      REAL(KIND=8) :: SPZ(KTY),SPE(KTY)
!
! ------------------------------------------------------------------------------
!  TIMING
!
      REAL(KIND=4) :: time1(2),time2(2),TIME3,T(2)
      REAL(KIND=4) :: TIC,TOC,MYTIME
!
! ---------------------------------------------------------------------------------
!  ALC 
      REAL(rk) :: delta_v2,FixedPsize
      LOGICAL :: halveFLAG,doubleFLAG,changeDELT
      LOGICAL :: arclcFLAG,arclcFLAG2,CONVERGED
      LOGICAL :: FirstRun,SecondRun
      LOGICAL :: FirstOrder,FullSecondOrder,SecondOrderPredOnly
      LOGICAL :: drN_reset,FixedP
!
!-----------------------------------------------------------------------
! Arnoldi
      LOGICAL ::  ArnoldiFLAG  
!
! ---------------------------------------------------------------------------------
!  SYMMETRIC SUBSPACE
!
      INTEGER :: vertex,rotateYFLAG
      REAL(rk) :: initialSHIFTX,initialSHIFTY
!
! ---------------------------------------------------------------------------------
!     INFILE
!
      character(len=256) :: infile
      character(len=256) :: InfileType
!
!  ---------------------------------------------------------------------------------
      end module GLOBAL_PARAMS
!

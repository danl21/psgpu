!  
      module GLOBAL_PARAMS
      IMPLICIT NONE
!
! ------------------------------------------------------------------------
!     BASIC PARAMETERS
! ------------------------------------------------------------------------
!
      INTEGER, PARAMETER  :: rk = 8      ! floating point precision      
!
      REAL(rk), PARAMETER :: alpha1 = 0.5_rk ! streamwise aspect ratio of the box Ly/Lx or Lz/Lx
!
! Resolution
      INTEGER, PARAMETER  :: NZ = 64 
      INTEGER, PARAMETER  :: NY = NZ
      INTEGER, PARAMETER  :: NX = INT(NZ/alpha1)
      INTEGER, PARAMETER  :: NX2 = NX/2+1
!
! dealiased resolution
      INTEGER, PARAMETER  :: KTZ = NZ/3
      INTEGER, PARAMETER  :: KTX = NX/3
      INTEGER, PARAMETER  :: KTY = NY/3
      INTEGER, PARAMETER  :: IKTX = KTX+1
      INTEGER, PARAMETER  :: IKTY = KTY+KTY+1
      INTEGER, PARAMETER  :: IKTZ = KTZ+KTZ+1

      INTEGER, PARAMETER  :: NR  = NY*NX*NZ 
      INTEGER, PARAMETER  :: NKT = IKTY*IKTX*IKTZ 

!
!
!  -----------------------------------------------------------------------
!  CONSTANTS (i,2*pi,sqrt(2))
!
      COMPLEX(rk) :: ZI 
      REAL(rk)    :: TWOPI 
      REAL(rk)    :: ROOT2 
!
! ------------------------------------------------------------------------
!  MASK
! 
      INTEGER, DIMENSION(NKT) :: L  ! stencil for wavenumber plane
      LOGICAL, DIMENSION(NKT) :: LL  ! stencil for wavenumber plane
      LOGICAL, DIMENSION(4*NKT) :: L3  ! stencil for full state vector
! ----------------------------------------------------------------------------
!  PHYSICAL PARAMETERS
!
      REAL(rk) ::  V2,AMPFOR,TSTEP,DELT,AMPV,Ri,Theta,Sc,Dtarget,alpha
      INTEGER  ::  KF,NSTOP,NOUT,NOUTV
!
!--------------------------------------------------------------------------------
!  MPI
!
      INTEGER :: NPROCS,RANK,IERROR
!--------------------------------------------------------------------------------
!  FLAGS
!
      INTEGER :: ISEED,IREST
      INTEGER :: STATSFLAG,RCFLAG,UPOflag,ISTOP,ADAPTFLAG
!
!
! --------------------------------------------------------------------------
!  RECURRENCE
!
      REAL(rk) ::  ResidualThreshold
!
!
! ---------------------------------------------------------------------------------
!  TIMING                                                                                                                                                                                                                                   
!                                                                                                                                                                                                                                 
      REAL(KIND=4) :: time1(2),time2(2),TIME3,T(2)
      REAL(KIND=4) :: TIC,TOC,MYTIME
!                                                                                                                                                                                                  
! --------------------------------------------------------------------------------- 
!     input file and working directory strings
!
      character(len=256) :: infile
      character(len=256) :: InfileType
      character(len=256) :: simdir
!
!  ---------------------------------------------------------------------------------
      end module GLOBAL_PARAMS
!

!  
      module GLOBAL_PARAMS
      IMPLICIT NONE
!
! ------------------------------------------------------------------------
!     PARAMETERS
! ------------------------------------------------------------------------
!
      INTEGER, PARAMETER  :: rk = 8      ! floating point precision      
!
      REAL(rk),PARAMETER :: alpha1 = 0.5_rk
!
      INTEGER, PARAMETER  :: NZ = 64
      INTEGER, PARAMETER  :: NY = NZ
      INTEGER, PARAMETER  :: NX = INT(NZ/alpha1)
      INTEGER, PARAMETER  :: NX2 = NX/2+1
!
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
!  CONSTANTS
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
      LOGICAL, DIMENSION(4*NKT) :: L2  ! stencil for full state vector
! ----------------------------------------------------------------------------
!  PHYSICAL PARAMETERS
!
      REAL(rk) ::  v1,V2,AMPFOR,TSTEP,DELT,AMPV,Ri,Theta,Sc,alpha,scale
      INTEGER  ::  KF,NSTOP,NOUT,NOUTV,OUTPLANE
!
!-------------------------------------------------------------------------------- 
!  MPI
!
      INTEGER :: NPROCS,RANK,IERROR
!--------------------------------------------------------------------------------       
!                                                                         
      INTEGER, PARAMETER  :: GMRESMAX = 2000                               
      INTEGER, PARAMETER  :: NEWTONMAX = 150                               
      INTEGER, PARAMETER  :: LWORKsvd = 7*GMRESMAX + 4*GMRESMAX*GMRESMAX  
      INTEGER, PARAMETER  :: LRWORKsvd= 7*GMRESMAX + 5*GMRESMAX*GMRESMAX  
!                                                                         
! ---------------------------------------------------------------------------
!  FFTs                                                                      
!                                                                            
      INTEGER(KIND=8) :: PlanK2XIR,PlanK2ETR,PlanK2ZTR,PlanK2UR,PlanK2VR,PlanK2WR
      INTEGER(KIND=8) ::PlanUR2K,PlanVR2K,PlanWR2K       
!                                                                            
! ----------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!  FLAGS
!       
      INTEGER :: myISEED,ISEED,GRFLAG,TRFLAG,RCFLAG,STATSFLAG,ICFLAG,IREST
      INTEGER :: WNFLAG,SubSpaceFLAG,lamFLAG,UPOflag,ISTOP,ADAPTFLAG,ParFLAG             
!                                                                      
!--------------------------------------------------------------------------------
!  GMRES PARAMS                                                                  
!                                                                                
      INTEGER :: IGUESSEND                                                       
      REAL(rk) :: maxPERIOD,gmresResT,SNEWTX,SNEWTZ,SNEWTY
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
      REAL(rk) :: initialSHIFTX,initialSHIFTY,initialSHIFTZ
!
! ---------------------------------------------------------------------------------
!     INFILE
!
      character(len=256) :: infile
      character(len=256) :: InfileType
      character(len=256) :: simdir
!
!  ---------------------------------------------------------------------------------
      end module GLOBAL_PARAMS
!

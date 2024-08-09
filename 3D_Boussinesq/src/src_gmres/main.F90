! -------------------------------------------------------------------
     PROGRAM MAIN 
! -------------------------------------------------------------------
!    Main program for the Newton-GMRES-hookstep algorithm
!
!    Original code written by Gary Chandler in Bristol
!    See Chandler & Kerswell, J. Fluid Mech.  2013
!    "Invariant recurrent solutions embedded in a turbulent 
!    two-dimensional Kolmogorov flow"
!    For details of the method and for referencing.
!
!    Adapted for 3D and variable density flow by Dan Lucas
!
! -------------------------------------------------------------------
!
      use GLOBAL_PARAMS
      use PARAMETERS
      use io_arpack
      IMPLICIT NONE 
! -------------------------------------------------------------------
      COMPLEX(rk), DIMENSION(NKT) :: XII,ETA,ZET,RHO
      COMPLEX(rk), DIMENSION(4*NKT) :: ZN,ZO,UK,ZR
      COMPLEX(rk), DIMENSION(4*NKT) :: ZNEWT,ZNN,ZN2,ZNSHIFT
      COMPLEX(rk), DIMENSION(4*NKT) :: DIFFZ,JFDIFF
      COMPLEX(rk), DIMENSION(4*NKT) :: Zplus,Zplus2,FNplus,FNplus2
      COMPLEX(rk), DIMENSION(4*NKT) :: FN,FNpred,FNpredC
      COMPLEX(rk), DIMENSION(4*NKT) :: dZNdt,dzdt,v
      COMPLEX(rk), DIMENSION(4*NKT) :: dZNdSx,dzdSx
      COMPLEX(rk), DIMENSION(4*NKT) :: dZNdSz,dzdSz
      COMPLEX(rk), DIMENSION(4*NKT) :: dZrdr,dZrdrO2,d2Zrdr2
      COMPLEX(rk), DIMENSION(4*NKT) :: Zr0,Zr1,Zr2
      COMPLEX(rk),DIMENSION(:,:),ALLOCATABLE :: q
      COMPLEX(rk) :: TERMZ,TZ,C1,AVZ,SHIFT
! -------------------------------------------------------------------
      REAL(rk), DIMENSION(NKT)      :: KX,KY,KZ
      REAL(rk), DIMENSION(NR)       :: UR,VR,WR,XIR,ETR,ZTR
      REAL(rk), DIMENSION(GMRESMAX+1,GMRESMAX) :: H,Hcopy
      REAL(rk), DIMENSION(GMRESMAX)            :: Ssvd,Xh
      REAL(rk), DIMENSION(GMRESMAX+1)          :: P,Yh,qN1,qN2,qN3,qN4
      REAL(rk) :: Dtarget,TEMP,NORMZ,FK,RANNO,TSTART,TIME,kkx,kky,kkz,wk
      REAL(rk) :: Usvd(GMRESMAX+1,GMRESMAX+1),VHsvd(GMRESMAX,GMRESMAX)
      REAL(rk) :: WORKsvd(LWORKsvd)
      REAL(rk) :: NER,NERold,GMR,const,FN1pred,FN2pred,FN3pred,FN4pred
      REAL(rk) :: NORMFN,NORMFNplus,NORMFNplus2,NORMFNpred,NORMFNpredC
      REAL(rk) :: NETOL,GMTOL,epsilon1,epsilon2,epsilon,DELTA,GMRESdelta
      REAL(rk) :: vN1,vN2,vN3,vN4,normX,MU1,MU2,MUmid
      REAL(rk) :: TNEWT,DIFFT,Tplus,Tplus2,TRUEDELT,TrueTRUEDELT
      REAL(rk) :: DIFFSx,Sxplus,Sxplus2,DIFFSz,Szplus,Szplus2
      REAL(rk) :: TSTARTin,TNEWTin,SNEWTXin,SNEWTYin,SNEWTZin,NERin,v2in
      REAL(rk) :: delta_r1,delta_r2,Ndelta_r,Ndelta_rmax,drrat
      REAL(rk) :: NORMb,v2N,v2plus,v2plus2,DIFFv2,DIFFv3
      REAL(rk) :: Tr0,Tr1,Tr2,Sxr0,Sxr1,Sxr2,Szr0,Szr1,Szr2,v2r0,v2r1,v2r2
      REAL(rk) :: dTrdr,dSxrdr,dSzrdr,dv2rdr,dTrdrO2,dSxrdrO2,dSzrdrO2,dv2rdrO2
      REAL(rk) :: d2Trdr2,d2Sxrdr2,d2Szrdr2,d2v2rdr2
      REAL(rk) :: NORMbplus,NORMbpredC,R3,Narc,Narcplus
! -------------------------------------------------------------------!
      INTEGER ::  IKF(NX2*NZ*NY), IKN(NKT)
      INTEGER ::  IKX,IKY,IKZ,jw,KFX1,KFX2,KFY1,KFY2,KFZ2,KFZ1,KFT1,KFT2
      INTEGER ::  NT,IKK,IK,IT,NERcount,IKT,IYY,IZZ
      INTEGER ::  i,j,ix,jx,m,ic,jc,biggest,id,jd,nvort
      INTEGER ::  INE,IGM,IORTH,INFOsvd,IWORKsvd(8*GMRESMAX)
      INTEGER ::  IGUESS,GUESSARRAY(2),iuposout,IGUESSSTART,convergeFLAG
      INTEGER ::  IHOOK,IHOOKMAX,GMREScount,HOOKcount,AllocateStatus
! -------------------------------------------------------------------!
      EQUIVALENCE (XII,XIR), (ETA,ETR), (ZET,ZTR)
!
!     SET UP VARIOUS SPECIAL VARIABLES FOR ARNOLDI & ARPACK
      integer(kind=4), parameter ::                                 &
    &       nev = 30,                                               &
    &       ncv = 120,                                              &
    &       lworkl = 3*ncv*ncv + 6*ncv
!
      real(rk), parameter :: TOL = 0.000001_rk
!
      integer(kind=4) :: ido, info, cnt, switch, outfreqarn
      character(len=2) :: which
!
      integer(kind=4), dimension(11) :: iparam
      integer(kind=4), dimension(14) :: ipntr
      logical, dimension(ncv) :: sel
!
      real(rk),DIMENSION(:),ALLOCATABLE :: resid,VarnC,workd
      real(rk),DIMENSION(:,:),ALLOCATABLE :: Varn,myeigvecs
      real(rk), dimension(lworkl) :: workl
!
      real(rk), dimension(3*ncv) :: workev
      real(rk), dimension(nev+1) :: DR, DI
      real(rk), dimension(ncv) :: myritzr,myritzi,myritzerr
      real(rk), dimension(ncv) :: lambdar,lambdai
      real(rk), dimension(ncv) :: WWR,WWI
      real(rk), dimension(5*ncv) :: dgeevWORK
      real(rk), dimension(ncv,ncv) :: Hess,dgeevVL,dgeevVR
      real(rk) :: lambdarsum     
!
      integer :: ArnoldiEXIT, Iarn, Iarn2, dgeevINFO
!
!--------------------------------------------------------------------------
! ARPACK trace output

      integer(kind=4)  logfil, ndigit, mgetv0, &
           msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd, &
           mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd, &
           mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/ &
           logfil, ndigit, mgetv0, &
           msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd, &
           mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd, &
           mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      ndigit = -6
      mnaup2 = 2
      mnaitr = 3
!-----------------------------------------------------------------------

      call set_gpu(0)
! -------------------------------------------------------------------
!
!
!     Read in input parameters from file and set constants.
      CALL params_read_txt()   
      TNEWT = TSTEP   
      SNEWTX = 0._rk   
      SNEWTY = 0._rk
      SNEWTZ = 0._rk
      TrueTRUEDELT = DELT   
      STATSFLAG=0
      RCFLAG=0
      ADAPTFLAG=0
      RANK=0
      Dtarget=1._rk
      FK = REAL(KF,rk)
      write(simdir,"(A3)") ""
! -------------------------------------------------------------------
!     GMRES-HOOKSTEP PARAMETERS
!
      GUESSARRAY = (/ 4,14 /)
      NETOL  = 10._rk**(-10._rk)
      GMTOL  = 10._rk**(-7._rk)
      epsilon = 10._rk**(-7._rk)
      GMRESdelta = 0.05_rk
      const =  10._rk**(-4._rk)
      IHOOKMAX = 50
      changeDELT = .FALSE.
      convergeFLAG = 0
! -------------------------------------------------------------------
!
      CALL paramsEXTRA_read_txt()
      NORMb = 1._rk
      NORMbplus = 0._rk

      ALLOCATE ( q(4*NKT,GMRESMAX+1),STAT=AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory for q ***"

      IF(ArnoldiFLAG) THEN
         ALLOCATE ( resid(SzV),STAT=AllocateStatus)
         IF (AllocateStatus /= 0) STOP "*** Not enough memory for resid ***"
         ALLOCATE ( VarnC(SzV),STAT=AllocateStatus)
         IF (AllocateStatus /= 0) STOP "*** Not enough memory for VarnC ***"
         ALLOCATE ( Varn(SzV,ncv),STAT=AllocateStatus)
         IF (AllocateStatus /= 0) STOP "*** Not enough memory for Varn ***"
         ALLOCATE ( myeigvecs(SzV,ncv),STAT=AllocateStatus)
         IF (AllocateStatus /= 0) STOP "*** Not enough memory for myeigvecs ***"
         ALLOCATE ( workd(3*SzV),STAT=AllocateStatus)
         IF (AllocateStatus /= 0) STOP "*** Not enough memory for workd ***"
      ENDIF
!
! -------------------------------------------------------------------------------
!
!     Create files for outputting UPO information
      open (150,file='UPOinfo.dat',form='formatted')                    
      WRITE(150,1501)   'Guess No : Start time : guess Re :guess Period &
     &:guess Shiftx: guess Shifty: guess Shiftz: guess Res  :      Status     : UPO N&
     &o :   Re   : Period : Shift x : Shift y : Shift z : Residual     : #NEWTO&
     &N : #GMRES : #HOOK '        
      open (999,file='UPOs.out',form='unformatted')
      write(999) 1,IKTX,IKTY,IKTZ       
      iuposout=0
!
!
! ---------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------- 
! --------------- Main Loop that cycles through different UPO guesses -----------------------    
! -------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------
!
!     Set starting guess number
      IGUESSSTART = IREST
      IF(arclcFLAG) IGUESSEND = IREST + IGUESSEND  ! IGUESSEND in this case is maximum number of tracing steps                 
!      IF(ArnoldiFLAG) IGUESSEND = IREST  ! Only 1 UPO is put through ARPACK
!
!     Start loop 
      DO IGUESS=IGUESSSTART,IGUESSEND
!
         IREST = IGUESS
!        Reset IREST for tracing                                                                                                      
         IF(arclcFLAG) THEN
            IF(FirstRun) THEN
               IREST = IGUESSSTART
            ELSE
               IREST = 0
            END IF
         END IF
!
!     Reset some paramters
         DELTA = GMRESdelta
         TRUEDELT = TrueTRUEDELT
!
! Most of the following standard initialization stuff is the same as in the DNS code
! ------------------------------------------------------------------------------------
!
! INITIALISE STUFF, GET THE FFT READY, ETC.
!
      TIME  = RANNO(myISEED) !Initializes random number generator 
!
      IF (NOUT.EQ.0) NOUT = 1
      IF (NOUTV.EQ.0) NOUTV = 1
!
      open (99,file='Zk.out',form='unformatted')
! -------------------------------------------------------------------------------
!     Initialize variables
      NT = 0
      MYTIME =   0._rk
      ZN  = CMPLX(0._rk,0._rk,rk)
      ZNN  = CMPLX(0._rk,0._rk,rk)
      ZN2  = CMPLX(0._rk,0._rk,rk) 
      ZNSHIFT  = CMPLX(0._rk,0._rk,rk)
      ZNEWT = CMPLX(0._rk,0._rk,rk)
      ZO  = CMPLX(0._rk,0._rk,rk)
      UK  = CMPLX(0._rk,0._rk,rk)
!
      DIFFT = 0._rk
      FN     = CMPLX(0._rk,0._rk,rk)
      FNplus = CMPLX(0._rk,0._rk,rk)
      DIFFZ  = CMPLX(0._rk,0._rk,rk)
      JFDIFF = CMPLX(0._rk,0._rk,rk)
      dZNdt  = CMPLX(0._rk,0._rk,rk)
      dzdt   = CMPLX(0._rk,0._rk,rk)
      dZNdsx  = CMPLX(0._rk,0._rk,rk)
      dzdsx   = CMPLX(0._rk,0._rk,rk)
      dZNdsz  = CMPLX(0._rk,0._rk,rk)
      dzdsz   = CMPLX(0._rk,0._rk,rk)
      v      = CMPLX(0._rk,0._rk,rk)
      q      = CMPLX(0._rk,0._rk,rk)
      halveFLAG = .FALSE.
      doubleFLAG = .FALSE.
!
      TSTARTin = 0._rk 
      TNEWTin = FixedPsize
      SNEWTXin = 0._rk 
      SNEWTYin = 0._rk 
      SNEWTZin = 0._rk 
      NERin = 0._rk 
      NER = 0._rk 
      NERold = 0._rk
      NERcount = 0
!     
      if (irest.eq.0) TSTART=0._rk
! -----------------------------------------------------------------------------  
!   
!     Read in starting state and initialize arrays 
      CALL INIT (XII,ETA,ZET,RHO,KX,KY,KZ,TSTART,TNEWTin,SNEWTXin,SNEWTYin,SNEWTZin,NERin)
      
      if(parFLAG.eq.1) v2=v1
      if(parFLAG.eq.2) v2=Ri
      if(parFLAG.eq.3) v2=alpha
      normX = SQRT( SUM( REAL(XII*CONJG(XII)), LL ) + SUM( REAL(ETA*CONJG(ETA)), LL )  &
           & + SUM( REAL(ZET*CONJG(ZET)), LL )  + SUM( REAL(RHO*CONJG(RHO)), LL ))
      print*, 'normX = ',normX
      
      scale = 1._rk
!      scale = sqrt(Ri)
      if(scale .eq. 0._rk) scale=1._rk
      print*, 'scale = ',scale

      DO IKT = 1,NKT
               ZO(IKT)       = XII(IKT)
               ZO(NKT+IKT)   = ETA(IKT)
               ZO(2*NKT+IKT) = ZET(IKT)
               ZO(3*NKT+IKT) = scale*RHO(IKT)
      ENDDO
      NORMX = SQRT( SUM( REAL(ZO*CONJG(ZO)),L3) )
      print*, 'normZ3 = ',normX

      !     Store index array for padding before KR FFT
      DO IKX=1,IKTX
         DO IKY=1,KTY+1
            IYY = IKY-1+KTY
            DO IKZ=1,KTZ+1
               IZZ = IKZ-1+KTZ
               IK = ((IKZ-1)*NY +IKY-1)*NX2 + IKX

               IKF(IK) = IKX-1+IKTX*(IYY+IKTY*IZZ)
            ENDDO
            DO IKZ=KTZ+2,NZ-KTZ
               IK = ((IKZ-1)*NY +IKY-1)*NX2 + IKX
               IKF(IK) = -1
            ENDDO
            DO IKZ=NZ-KTZ+1,NZ
               IZZ = IKZ-1+KTZ-NZ
               IK = ((IKZ-1)*NY +IKY-1)*NX2 + IKX

               IKF(IK) = IKX-1+IKTX*(IYY+IKTY*IZZ)
            ENDDO
         ENDDO
         DO IKY=KTY+2,NY-KTY
            DO IKZ=1,NZ
               IK = ((IKZ-1)*NY +IKY-1)*NX2 + IKX
               IKF(IK) = -1
            ENDDO
         ENDDO
         DO IKY=NY-KTY+1,NY
            IYY = IKY-1+KTY-NY
            DO IKZ=1,KTZ+1
               IZZ = IKZ-1+KTZ
               IK = ((IKZ-1)*NY +IKY-1)*NX2 + IKX

               IKF(IK) = IKX-1+IKTX*(IYY+IKTY*IZZ)
            ENDDO
            DO IKZ=KTZ+2,NZ-KTZ
               IK = ((IKZ-1)*NY +IKY-1)*NX2 + IKX
               IKF(IK) = -1
            ENDDO
            DO IKZ=NZ-KTZ+1,NZ
               IZZ = IKZ-1+KTZ-NZ
               IK = ((IKZ-1)*NY +IKY-1)*NX2 + IKX

               IKF(IK) = IKX-1+IKTX*(IYY+IKTY*IZZ)
            ENDDO
         ENDDO
      ENDDO

!
      DO IKZ=1,NZ
         DO IKY=1,NY
            DO IKX=IKTX+1,NX/2+1
               IK = ((IKZ-1)*NY +IKY-1)*NX2 + IKX
               IKF(IK) = -1
            ENDDO
         ENDDO
      ENDDO

!     Store index array for rewriting and normalising after RK FFT
       DO IKX =1,IKTX
          DO IKZ =1,KTZ
             IZZ = NZ+IKZ-1-KTZ
             DO IKY =1,KTY
                IYY = NY+IKY-1-KTY

                IK  = ((IKZ-1)*IKTY+IKY-1)*IKTX +IKX
                IKK = (IZZ*NY+IYY)*NX2 + IKX-1

                IKN(IK) = IKK
             ENDDO
             DO IKY =KTY+1,IKTY
                IYY=IKY-1-KTY

                IK  = ((IKZ-1)*IKTY+IKY-1)*IKTX +IKX
                IKK = (IZZ*NY+IYY)*NX2 + IKX-1

                IKN(IK) = IKK
             ENDDO
          ENDDO
          DO IKZ =KTZ+1,IKTZ
             IZZ=IKZ-1-KTZ
             DO IKY =1,KTY
                IYY = NY+IKY-1-KTY

                IK  = ((IKZ-1)*IKTY+IKY-1)*IKTX +IKX
                IKK = (IZZ*NY+IYY)*NX2 + IKX-1

                IKN(IK) = IKK
             ENDDO
             DO IKY =KTY+1,IKTY
                IYY=IKY-1-KTY

                IK  = ((IKZ-1)*IKTY+IKY-1)*IKTX +IKX
                IKK = (IZZ*NY+IYY)*NX2 + IKX-1

                IKN(IK) = IKK
             ENDDO

          ENDDO
       ENDDO
! --------------------------------------------------------------------------------
!   *****************************************************************
!   *       GMRES Guess Filters                                      *
!   *****************************************************************
!
!     Decide whether or not to attempt guess 
      IF (UPOflag.eq.1) THEN
         Print*,'   '
         Print*,'PERIODin = ',TNEWTin
         Print*,'SHIFTXin = ',SNEWTXin
         Print*,'SHIFTYin = ',SNEWTYin
         Print*,'SHIFTZin = ',SNEWTZin
         Print*,'Residual = ',NERin
!
         IF( .NOT.(arclcFLAG .OR. ArnoldiFLAG)) THEN  ! only apply guess filters if alc and stability are off
!        *****************************************************************************
            IF (TNEWTin.GT.maxPERIOD .OR. NERin.GT.gmresResT) THEN
!        *****************************************************************************
               print*,'Not trying Guess number ',IREST
               WRITE(150,1502) IREST,'    ',                                 &
     &              TSTARTin,v2,TNEWTin,SNEWTXin,SNEWTYin,SNEWTZin,NERin,   &
     &   ' :   "NOT TRIED"   : ',0,v2,TNEWTin,SNEWTXin,SNEWTYin,   &
     &               NERin 
               CYCLE
            ENDIF   
         END IF
!        Initialize Guess variables
         TNEWT = TNEWTin
         SNEWTX = SNEWTXin
         SNEWTY = SNEWTYin
         SNEWTZ = SNEWTZin
         NER = NERin
         IF (FixedP) THEN
            SNEWTX = FixedPsize*SNEWTX/TNEWT
            SNEWTZ = FixedPsize*SNEWTZ/TNEWT
            TNEWT = FixedPsize
         END IF
      ENDIF
! -----------------------------------------------------------------------------------------
!
!     Setup Arnoldi variables
      IF(ArnoldiFLAG) THEN
         IF( SzV .NE. 2*COUNT(L3) ) THEN
            Print*,'Arnoldi error: SzV=',SzV,' is not equal to 2*COUNT(L3)=',2*COUNT(L3)
            STOP
         ENDIF
         ArnoldiEXIT = 0
         which = 'LM'       ! specifies which nev eigenvalues are to be computed.
!                           ! 'LM' = Eigenvalues of largest magnitude
!                           ! 'LR' = Eigenvalues of largest real part
         outfreqarn = 20
         cnt = 0
         ido = 0
         info = 0               ! ARPACK creates starting vector
!         info = 1               ! use supplied vector in "resid" as starting vector
         iparam(1) = 1          ! "exact" shifting
         iparam(3) = 1000       ! max number of implicit restarts
         iparam(7) = 1          ! mode
         switch = 0
!         resid = 1._rk
      ENDIF
!
! ----------------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------------
!                                                                                   
      IF(arclcFLAG) THEN
!                                                                                                                                   
!        ---------------------------------------------------------------------------                                                
!        --------------------------------------------------------------------------                                                 
         IF(FirstRun) THEN
!        -------------------------------------------------------------------------                                                  
!        --------------------------------------------------------------------------                                                 
!                                                                                                                                    
!           Do standard continuation step with fixed change in Re                                                                    
            IF(CONVERGED) THEN
               v2 = v2 - delta_v2
            ELSE
               v2 = v2 + delta_v2
               delta_v2 = delta_v2 / 2._rk
               v2 = v2 - delta_v2
            END IF
            WHERE (L3) Zr2 = ZO
            IF(.NOT.FixedP) Tr2 = TNEWT
            Sxr2 = SNEWTX
            Szr2 = SNEWTZ
            v2r2 = v2 + delta_v2
            CONVERGED = .FALSE.
!                                                                                                                                    
!        ----------------------------------------------------------------------------                                                
!        ----------------------------------------------------------------------------                                                
         ELSE
!        ---------------------------------------------------------------------------                                                
!        ----------------------------------------------------------------------------                                                
!                                                                                                                                    
!           Code should hopefully realize when UPOs are equilibrium points or T-Ws                                                   
            IF (Tr2 .LT. FixedPsize) THEN
               print*,'Period has decreased below',FixedPsize
               print*,'UPO is probably an equilibrium point or travellin&                                                            
     &g wave'
               print*,'Setting Fixed Period flag to .TRUE.'
               FixedP = .TRUE.
               Sxr2 = FixedPsize*Sxr2/Tr2
               Szr2 = FixedPsize*Szr2/Tr2
               Tr2 = FixedPsize
               Sxr1 = FixedPsize*Sxr1/Tr1
               Szr1 = FixedPsize*Szr1/Tr1
               Tr1 = FixedPsize
               Sxr0 = FixedPsize*Sxr0/Tr0
               Szr0 = FixedPsize*Szr0/Tr0
               Tr0 = FixedPsize
            END IF
!        ----------------------------------------------------------------------------                                               
!                                                                                                                                    
            Print*,'delta_r2 from previous step: ',delta_r2
!           calculate delta_r (ALC step size)                                                                                        
            delta_r2 = SQRT( SUM( REAL((Zr2-Zr1)*CONJG(Zr2-Zr1)), L3 )  &
     &            + (Sxr2-Sxr1)**2+ (Szr2-Szr1)**2 + (v2r2-v2r1)**2 )
            IF(.NOT.FixedP) delta_r2 = SQRT(delta_r2**2 + (Tr2-Tr1)**2)
            Print*,'actual delta_r2: ',delta_r2
!                                                                                                                                   
!           approximate derivates wrt r                                                                                             
            WHERE (L3) dZrdr = (Zr2-Zr1)/delta_r2
            IF(.NOT.FixedP) dTrdr = (Tr2-Tr1)/delta_r2
            dSxrdr = (Sxr2-Sxr1)/delta_r2
            dSzrdr = (Szr2-Szr1)/delta_r2
            dv2rdr = (v2r2-v2r1)/delta_r2
            Print*,'set 1st Order derivatives for correction step'
!        -------------------------------------------------------------------------                                                   
!                                                                                                                                    
            IF(SecondRun) THEN
!            set new delta_r for third run                                                                                           
               IF(CONVERGED) Ndelta_r = 2._rk*delta_r2
               IF(.NOT.CONVERGED) Ndelta_r = Ndelta_r/2._rk
               Ndelta_rmax = Ndelta_r
            ELSE
               IF(CONVERGED) THEN
                  IF(delta_r2 .LT. Ndelta_rmax) THEN
                     Ndelta_r = delta_r2
                  ELSE
                     Ndelta_r = Ndelta_rmax
                  END IF
                  IF(1.1_rk*Ndelta_r .LT. Ndelta_rmax) THEN
                     Ndelta_r = 1.1_rk*Ndelta_r
                  END IF
               END IF
               IF(.NOT.CONVERGED) Ndelta_r = Ndelta_r/2._rk
            END IF
            CONVERGED = .FALSE.
!
!           ---------------------------------------------------------------------------                                              
            IF(SecondRun .OR. FirstOrder) THEN
!           -----------------------------------------------------------------------                                                  
!              DO first Order Prediction of variables                                                                                
!                                                                                                                                    
!              Extrapolate starting state                                                                                            
               WHERE (L3) ZO = Zr2 + dZrdr*Ndelta_r
               IF(FixedP)  THEN
                  TNEWT = FixedPsize
               ELSE
                  TNEWT = Tr2  + dTrdr*Ndelta_r
               END IF
               SNEWTX = Sxr2  + dSxrdr*Ndelta_r
               SNEWTZ = Szr2  + dSzrdr*Ndelta_r
               v2    = v2r2 + dv2rdr*Ndelta_r
!                                                                                                                                    
               Print*,'done 1st Order prediction'
!                                                                                                                                    
!           ---------------------------------------------------------------------------                                             
            ELSE
!           ------------------------------------------------------------------------                                                
!                                                                                                                                    
!              Do second order prediction of variables                                                                              
!                                                                                                                                   
               drrat = delta_r2 / (delta_r2 + delta_r1)
!                                                                                                                                   
!            approximate first derivatives                                                                                          
             WHERE (L3)                                                   &
     &     dZrdrO2 = ( (1._rk-drrat*drrat)*Zr2 - Zr1 + drrat*drrat*Zr0 )  &
     &               / ( delta_r2 * (1._rk-drrat) )
           IF(.NOT.FixedP) dTrdrO2 = ( (1._rk-drrat*drrat)*Tr2 - Tr1      &
     &   + drrat*drrat*Tr0 ) / ( delta_r2 * (1._rk-drrat) )
           dSxrdrO2 = ( (1._rk-drrat*drrat)*Sxr2 - Sxr1 + drrat*drrat*Sxr0 )  &
     &               / ( delta_r2 * (1._rk-drrat) )
           dSzrdrO2 = ( (1._rk-drrat*drrat)*Szr2 - Szr1 + drrat*drrat*Szr0 )  &
     &               / ( delta_r2 * (1._rk-drrat) )
         dv2rdrO2 = ( (1._rk-drrat*drrat)*v2r2 -v2r1 +drrat*drrat*v2r0 )  &
     &               / ( delta_r2 * (1._rk-drrat) )
!
!           ---------------------------------------------------------------------------                                             
!            approximate second derivatives                                                                                         
               WHERE (L3)                                               &
     &         d2Zrdr2 = ( (drrat-1._rk)*Zr2 + Zr1 - drrat*Zr0 )        &
     &               / ( 0.5_rk*delta_r2*delta_r2*(1._rk-1._rk/drrat) )
      IF(.NOT.FixedP) d2Trdr2 = ( (drrat-1._rk)*Tr2 + Tr1 - drrat*Tr0 ) &
     &               / ( 0.5_rk*delta_r2*delta_r2*(1._rk-1._rk/drrat) )
               d2Sxrdr2 = ( (drrat-1._rk)*Sxr2 + Sxr1 - drrat*Sxr0 )        &
     &               / ( 0.5_rk*delta_r2*delta_r2*(1._rk-1._rk/drrat) )
               d2Szrdr2 = ( (drrat-1._rk)*Szr2 + Szr1 - drrat*Szr0 )        &
     &               / ( 0.5_rk*delta_r2*delta_r2*(1._rk-1._rk/drrat) )
               d2v2rdr2 = ( (drrat-1._rk)*v2r2 + v2r1 - drrat*v2r0 )    &
     &               / ( 0.5_rk*delta_r2*delta_r2*(1._rk-1._rk/drrat) )
!            -------------------------------------------------------------------------                                              
!                                                                                                                                   
!     Extrapolate starting state                                                                                                    
               WHERE (L3) ZO = Zr2 + dZrdrO2*Ndelta_r                     &
     &                             + 0.5_rk*Ndelta_r*Ndelta_r*d2Zrdr2
               IF(FixedP) THEN
                  TNEWT = FixedPsize
               ELSE
                  TNEWT = Tr2  + dTrdrO2*Ndelta_r                         &
     &                      + 0.5_rk*Ndelta_r*Ndelta_r*d2Trdr2
               END IF
               SNEWTX = Sxr2  + dSxrdrO2*Ndelta_r                            &
     &                      + 0.5_rk*Ndelta_r*Ndelta_r*d2Sxrdr2
               SNEWTZ = Szr2  + dSzrdrO2*Ndelta_r                            &
     &                      + 0.5_rk*Ndelta_r*Ndelta_r*d2Szrdr2
               v2    = v2r2 + dv2rdrO2*Ndelta_r                           &
     &                      + 0.5_rk*Ndelta_r*Ndelta_r*d2v2rdr2
!                                                                                                                                   
               Print*,'done 2nd Order prediction'
!                                                                                                                                   
               IF(FullSecondOrder) THEN
!                 Reset derivatives for correction step using 2nd Order calcs                                                       
                  WHERE (L3) dZrdr =  dZrdrO2
                  IF(.NOT.FixedP) dTrdr = dTrdrO2
                  dSxrdr =  dSxrdrO2
                  dSzrdr =  dSzrdrO2
                  dv2rdr = dv2rdrO2
             Print*,'reset derivatives for correction step to 2nd Order'
               END IF
!                                                                                                                                    
!           ---------------------------------------------------------------------------                                             
            END IF
!           ---------------------------------------------------------------------------
!                                                                                                                                   
!           If requested recalculate ALC step size based on new predicted values                                                  
            Print*,'Ndelta_r = ',Ndelta_r
            IF (drN_reset) THEN
              Ndelta_r = SQRT( SUM( REAL((ZO-Zr2)*CONJG(ZO-Zr2)), L3 )  &
     &            + (SNEWTX-Sxr2)**2+ (SNEWTZ-Szr2)**2 + (v2-v2r2)**2 )
           IF(.NOT.FixedP) Ndelta_r = SQRT(Ndelta_r**2 + (TNEWT-Tr2)**2)
              Print*,'Reset Ndelta_r to: ',Ndelta_r
            END IF
!        ------------------------------------------------------------------------                                                    
!           Reset GMRES variables                                                                                                    
            TNEWTin = TNEWT
            SNEWTXin = SNEWTX
            SNEWTZin = SNEWTZ
!        ------------------------------------------------------------------------                                                    
!                                                                                                                                    
!           Calculate initial Narc                                                                                                   
            Narc = SUM( REAL(dZrdr*CONJG(ZO-Zr2)), L3 )                 &
     &     + dSxrdr*(SNEWTX-Sxr2)+ dSzrdr*(SNEWTZ-Szr2) + dv2rdr*(v2-v2r2)                       &
     &    - Ndelta_r
            IF(.NOT.FixedP) Narc = Narc + dTrdr*(TNEWT-Tr2)
            print*,'initial Narc = ',Narc
!                                                                                                                                    
!        -----------------------------------------------------------------------                                                     
!        ----------------------------------------------------------------------                                                      
         END IF
!        -----------------------------------------------------------------------                                                     
!        -----------------------------------------------------------------------                                    
      END IF
!              
! ------------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------------
!
!     Output to Screen
      Print*,'   '
      print*,'        alpha=',alpha
      print*,'           Re=',1._rk/v2
      print*,'          ktx=',ktx
      print*,'          kty=',kty
      print*,'          ktz=',ktz
      print*,'          NSTOP=',NSTOP
      print*,'          NOUT=',NOUT
      print*,'          NOUTV=',NOUTV
      print*,'       epsilon =',epsilon
      print*,'         delta =',delta
!
      CALL FLUSH(6)
!
      print*,'Number of degrees of freedom = ',2*COUNT(L3)
      print*,'Max coefficient of initial condition = ',                 &
     &                                               MAXVAL(ABS(ZO),L3)
!  ---------------------------------------------------------------------------
!
!     OUTPUT INITIAL STATE 
!
      TIME = TNEWT
      izkout = 0
!      print*,'Wrote to Zk.out number ',izkout
      print*,'Period = ',TIME,'X shift =',SNEWTX,'Z shift =',SNEWTZ
!
!     -----------------------------------------------------------------------
      CALL etime(time2,tic)
!
      CALL FLUSH(6)
!
!     ------------------------------------------------------------------------
!     ------------------------------------------------------------------------
!
!     Calculate UPO and check error
!
!        If this is a guess from DNS then, due to numerical errors
!        the period stored may not be precisely the period used in the DNS.
!        The following 'IF' statement corrects for this, whilst hopefully
!        ignoring any partly converged guesses.
      DELT = TRUEDELT
      IF (ABS(TNEWT-NINT(TNEWT/DELT)*DELT).LT.10._rk**(-10._rk)) THEN
         print*,'Guess period is assumed to be from DNS'
         print*,'Adjusting period so it is exactly equivalent to DNS'
         TNEWT = NINT(TNEWT/DELT)*DELT
      ENDIF
!

!        Set ZNEWT = ZO (ZO gets updated during timestepping)
         WHERE (L3) ZNEWT = ZO  
         v2N = v2
!        From here onwards v2N mirrors ZNEWT and v2 mirrors Z0, i.e. v2N and ZNEWT
!        store the information for the current Newton step, whereas v2 and Z0 are 
!        temporary work variables that are used directly by the timestepping routine.
!
!      -------------------------------------------------------------------------
!        Timestep once round orbit
         v2in = v2
         TSTART = 0._rk
         DELT = TRUEDELT
         NSTOP  =   INT(TNEWT/DELT)  !Note this rounds down
         print*,'Calculating UPO'
         Print*,'Period = ',TNEWT,', dt = ',DELT,', NSTOP = ',NSTOP
         IF(parFLAG .eq.1)THEN
            Print*,'Re = ',1._rk / v2N
            v1=v2
            
         ELSE IF(parFLAG .eq.2) THEN
            Print*,'Ri = ', v2N
            Ri = v2
         ELSE IF(parFLAG .eq.3) THEN
            Print*,'alpha = ', v2N
            alpha = v2
            call makeK(KX,KY,KZ)
         ENDIF
         

         DO IKT =1,NKT
            IF (LL(IKT))THEN
               XII(IKT) = ZO(IKT)
               ETA(IKT) = ZO(NKT+IKT)
               ZET(IKT) = ZO(2*NKT+IKT)
               RHO(IKT) = ZO(3*NKT+IKT)/scale
            ENDIF
         ENDDO
!
      CALL TIMESTEP_CUDA(XII,ETA,ZET,RHO,KX,KY,KZ,TIME,TSTART,TNEWT,AMPFOR,DELT,ResidualThreshold,&
     &                   FK,v1,Ri,Sc,Theta,alpha,Dtarget,IKF,IKN,L,IKTX,IKTY,IKTZ,KTZ,NKT,NX,NY,NZ,&
     &                   NOUT,NOUTV,STATSFLAG,RCFLAG,ADAPTFLAG,RANK)

!     ----------------------------------------------------------------------------
!
!        Store final Z state 
         DO IKT=1,NKT
            IF(LL(IKT))THEN
               ZN(IKT)       = XII(IKT)
               ZN(NKT+IKT)   = ETA(IKT)
               ZN(2*NKT+IKT) = ZET(IKT)
               ZN(3*NKT+IKT) = RHO(IKT)*scale
               ZNN(IKT)       = XII(IKT)
               ZNN(NKT+IKT)   = ETA(IKT)
               ZNN(2*NKT+IKT) = ZET(IKT)
               ZNN(3*NKT+IKT) = RHO(IKT)*scale
            ENDIF
         ENDDO
         epsilon1=epsilon
         epsilon2=epsilon

         print*, 'epsilons = ',epsilon1,epsilon2
!
         !        Shift final state back 
         DO IKT = 1, NKT

                  IF(.NOT.LL(IKT))     CYCLE
                  SHIFT = EXP(-ZI*KY(IKT)*SNEWTY)*EXP(-ZI*KZ(IKT)*SNEWTZ)*EXP(-ZI*KX(IKT)*SNEWTX)
                  ZNSHIFT(IKT) = SHIFT*ZN(IKT)
                  ZNSHIFT(NKT+IKT) = SHIFT*ZN(NKT+IKT)
                  ZNSHIFT(2*NKT+IKT) = SHIFT*ZN(2*NKT+IKT)
                  ZNSHIFT(3*NKT+IKT) = SHIFT*ZN(3*NKT+IKT)
         ENDDO

!     -------------------------------------------------------------------------------     
         NORMZ = SQRT( SUM( REAL(ZN*CONJG(ZN)), L3 ) +SUM( REAL(ZNEWT*CONJG(ZNEWT)), L3 ))
!    -------------------------------------------------------------------------------
!      Calculate 1st residual0
!      Calculate F = Z(z,T) - z
         WHERE (L3) FN = ZNSHIFT - ZNEWT       
!
         NORMFN = SQRT( SUM( REAL(FN*CONJG(FN)), L3 ) )
!        Use NormFN to check for convergence 
         NER = NORMFN/NORMZ
!    ---------------------------------------------------------------------------------
!
         Print*,'-----------------------------------------------------'
         PRINT*,'Newton convergence check: ||F|| / ||Z|| = ',NER
         Print*,'-----------------------------------------------------'
         Print*,'-----------------------------------------------------'
!
!        Check if calculated residuals match the ones from timestepping
!        These should be an exact match, if everything is correct!
         IF (UPOflag.eq.1) THEN
            IF (NER.ne.NERin) THEN
               print*,'WARNING: calculated NER does not ',              &
     &                       'exactly match NER from Zk.in'
               Print*,'Setting starting NER to calculated NER'
               NERin = NER
            ENDIF
         ENDIF
!  -----------------------------------------------------------------------------
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ----------------------------------------------------------------------------
!     NEWTON ITERATIONS
! ----------------------------------------------------------------------------    
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!         
      GMREScount = 0
      HOOKcount = 0
      CALL CPU_TIME(tic)
      CALL FLUSH(6)

      DO INE=1,NEWTONMAX
!
!        qN1 and qN2 are zero on first GMRES iteration 
!        and so do not have to be added to NORMFN
         NORMb = NORMFN
!        But qN3 is not always zero!                                                                                                
         IF(arclcFLAG2) NORMb = SQRT(NORMFN*NORMFN + Narc*Narc)
!
!
         ! SET q(1) = b/|b| (in this case b = - FN )
         WHERE (L3) q(:,1) = - FN/NORMb  
         qN1(1) = 0._rk
         qN2(1) = 0._rk
         qN3(1) = 0._rk
         IF(arclcFLAG2) qN4(1) = - Narc/NORMb
!
         IF(.NOT.FixedP) THEN
            !        Calculate dZN/dt before GMRES iteration
            DO IKT =1,NKT
               IF (LL(IKT))THEN
                  XII(IKT) = ZNSHIFT(IKT)
                  ETA(IKT) = ZNSHIFT(NKT+IKT)
                  ZET(IKT) = ZNSHIFT(2*NKT+IKT)
                  RHO(IKT) = ZNSHIFT(3*NKT+IKT)/scale
               ENDIF
            ENDDO
            !                                                                                                                      
            CALL TIMEDERIV_CUDA(XII,ETA,ZET,RHO,KX,KY,KZ,TIME,TSTART,TEMP,AMPFOR,DELT,ResidualThreshold,&
                 &                   FK,v1,Ri,Sc,Theta,alpha,Dtarget,IKF,IKN,L,IKTX,IKTY,IKTZ,KTZ,NKT,NX,NY,NZ,&
                 &                   NOUT,NOUTV,STATSFLAG,RCFLAG,ADAPTFLAG,RANK)
!     ----------------------------------------------------------------------------                                     
            DO IKT=1,NKT
               IF(LL(IKT))THEN
                  dZNdt(IKT)       = XII(IKT)
                  dZNdt(NKT+IKT)   = ETA(IKT)
                  dZNdt(2*NKT+IKT) = ZET(IKT)
                  dZNdt(3*NKT+IKT) = RHO(IKT)*scale
               ENDIF
            ENDDO
            print*,'||dZNdt|| = ', SQRT( SUM( REAL(dZNdt*CONJG(dZNdt)), L3 ) )

!
!        Calculate dz/dt before GMRES iteration
            DO IKT =1,NKT
               IF (LL(IKT))THEN
                  XII(IKT) = ZNEWT(IKT)
                  ETA(IKT) = ZNEWT(NKT+IKT)
                  ZET(IKT) = ZNEWT(2*NKT+IKT)
                  RHO(IKT) = ZNEWT(3*NKT+IKT)/scale
               ENDIF
            ENDDO
            !                                                                                                                      
            CALL TIMEDERIV_CUDA(XII,ETA,ZET,RHO,KX,KY,KZ,TIME,TSTART,TEMP,AMPFOR,DELT,ResidualThreshold,&
                 &                   FK,v1,Ri,Sc,Theta,alpha,Dtarget,IKF,IKN,L,IKTX,IKTY,IKTZ,KTZ,NKT,NX,NY,NZ,&
                 &                   NOUT,NOUTV,STATSFLAG,RCFLAG,ADAPTFLAG,RANK)

!     ----------------------------------------------------------------------------                                     
!                                                                                                                      
            DO IKT=1,NKT
               IF(LL(IKT))THEN
                  dzdt(IKT)       = XII(IKT)
                  dzdt(NKT+IKT)   = ETA(IKT)
                  dzdt(2*NKT+IKT) = ZET(IKT)
                  dzdt(3*NKT+IKT) = RHO(IKT)*scale
               ENDIF
            ENDDO
            print*,'||dzdt|| = ', SQRT( SUM( REAL(dzdt*CONJG(dzdt)), L3 ) )
         END IF
!
!        Calculate dZN/dS and dz/dS before GMRES iteration
         DO IKT = 1, NKT
            IF(.NOT.LL(IKT))     CYCLE
            dZNdSx(IKT) = -ZI*KX(IKT)*ZNSHIFT(IKT)
            dzdSx(IKT) = -ZI*KX(IKT)*ZNEWT(IKT)

            dZNdSx(NKT+IKT) = -ZI*KX(IKT)*ZNSHIFT(NKT+IKT)
            dzdSx(NKT+IKT) = -ZI*KX(IKT)*ZNEWT(NKT+IKT)

            dZNdSx(2*NKT+IKT) = -ZI*KX(IKT)*ZNSHIFT(2*NKT+IKT)
            dzdSx(2*NKT+IKT) = -ZI*KX(IKT)*ZNEWT(2*NKT+IKT)

            dZNdSx(3*NKT+IKT) = -ZI*KX(IKT)*ZNSHIFT(3*NKT+IKT)
            dzdSx(3*NKT+IKT) = -ZI*KX(IKT)*ZNEWT(3*NKT+IKT)

            dZNdSz(IKT) = -ZI*KZ(IKT)*ZNSHIFT(IKT)
            dzdSz(IKT) = -ZI*KZ(IKT)*ZNEWT(IKT)

            dZNdSz(NKT+IKT) = -ZI*KZ(IKT)*ZNSHIFT(NKT+IKT)
            dzdSz(NKT+IKT) = -ZI*KZ(IKT)*ZNEWT(NKT+IKT)

            dZNdSz(2*NKT+IKT) = -ZI*KZ(IKT)*ZNSHIFT(2*NKT+IKT)
            dzdSz(2*NKT+IKT) = -ZI*KZ(IKT)*ZNEWT(2*NKT+IKT)

            dZNdSz(3*NKT+IKT) = -ZI*KZ(IKT)*ZNSHIFT(3*NKT+IKT)
            dzdSz(3*NKT+IKT) = -ZI*KZ(IKT)*ZNEWT(3*NKT+IKT)

         ENDDO
            print*,'||dzdSx|| = ', SQRT( SUM( REAL(dzdSx*CONJG(dzdSx)), L3 ) )
            print*,'||dZNdSx|| = ', SQRT( SUM( REAL(dZNdSx*CONJG(dZNdSx)), L3 ) )
            print*,'||dzdSz|| = ', SQRT( SUM( REAL(dzdSz*CONJG(dzdSz)), L3 ) )
            print*,'||dZNdSz|| = ', SQRT( SUM( REAL(dZNdSz*CONJG(dZNdSz)), L3 ) )
!
!     Initialize GMRES Hessenberg matrix to zeros
         H = 0._rk 
         Hcopy = 0._rk 
         Ssvd = 0._rk
         Usvd = 0._rk
         VHsvd = 0._rk
         WORKsvd = 0._rk
         INFOsvd = 0
!
         Print*,'Starting new Newton step and GMRES Hookstep iteration.'
         Print*,'This is Newton step number', INE
!
!    __________________________________________________________________________________
!    --------------------------------------------------------------------------------
!        GMRES
!    _________________________________________________________________________________
!    ----------------------------------------------------------------------------------
!       
!        GMRES LOOP to solve Az=b
         DO IGM=1,GMRESMAX
            CALL FLUSH(6)
!
!      ********************************************************************************
!       All  ARPACK related stuff happens here and exits the Newton loops on completion 
!      *********************************************************************************            
            IF(ArnoldiFLAG) THEN
!
!             ===================================================================
               IF(ArnoldiEXIT .EQ. 0) THEN
!             ===================================================================
!               Do regular arnoldi step
                  cnt = cnt + 1
                  write(*,*) "cnt = ", cnt

                  call dnaupd(ido, 'I', SzV, which, nev, TOL, resid,       &
     &    ncv, Varn, SzV, iparam, ipntr, workd, workl, lworkl, info)
!
!                 Check output from dnaupd to see what to do next
                  if (ido .eq. -1 .or. ido .eq. 1) then
!                    carry on arnoldi iteration
                     call workd2q(workd(ipntr(1)), q(:,IGM))
!
                  else
!                    ARPACK is finished - start output process
                     print*,'ido =',ido
!                     open (532,file='arpackfinal.out',form='unformatted')
!                     write(532) nev, ncv, Varn, ipntr, workl
!                     close (532)
!
                     myritzr = workl(ipntr(6):ipntr(6)+ncv-1)
                     myritzi = workl(ipntr(7):ipntr(7)+ncv-1)                    
                     myritzerr = workl(ipntr(8):ipntr(8)+ncv-1)
                     lambdar = LOG(SQRT(myritzr**2 + myritzi**2))/TNEWT
                     lambdai = ATAN2(myritzi,myritzr)/TNEWT
!
                     lambdarsum = 0._rk
                     Iarn2 = 0
                     open (222,file='eigenvalues.dat',form='formatted',access='APPEND')  
                     DO Iarn=1,ncv
                        IF (myritzerr(Iarn) .LT. 100._rk*TOL) THEN
                           IF (lambdar(Iarn) .GT. 0._rk) THEN 
                             Iarn2 = Iarn2 + 1
                             write(222,2221) Iarn2,lambdar(Iarn),lambdai(Iarn)
                             lambdarsum = lambdarsum + lambdar(Iarn)
                           END IF
                        ENDIF
                     END DO
                     write(222,2222) '     ' 
                     write(222,2221) IREST,lambdarsum,TNEWT
                     close (222)
 2221                FORMAT(I3,1x,F13.5,F13.5)
 2222                FORMAT(A)
!
!
         call arpack_final_output(info,nev,ncv,which,iparam,sel, DR, DI,     &
     &   Varn, SzV, workev, TOL,resid, ipntr, workd, workl,lworkl,TNEWT)
                     IF (iparam(5) .NE. nev) EXIT
!                 Check leading eigenvalue/eigenvector      
                     DO Iarn=1,nev
                      call workd2q(Varn(:,Iarn), q(:,IGM+Iarn-1))               
                     END DO
                     Iarn = 1
                     ArnoldiEXIT = IGM+nev
!
                  endif
!
!            ========================================================================
               ELSE
!            ========================================================================
!               DO stability output routines and exit Newton iterations
!
                  IF (DR(Iarn).EQ.DR(Iarn+1)   .AND.                         &
     &                DI(Iarn).EQ.-DI(Iarn+1)) THEN
                     VarnC = VarnC -                                      & 
     &                 (DR(Iarn)*Varn(:,Iarn) - DI(Iarn)*Varn(:,Iarn+1))
                     Print*,'Checking Real part of eigenvector ',Iarn
                     Print*,' A x_r - ( mu_r x_r - mu_i x_i ) ' 
                  ELSE IF (DR(Iarn).EQ.DR(Iarn-1)   .AND.                   &
     &                     DI(Iarn).EQ.-DI(Iarn-1)) THEN
                     VarnC = VarnC -                                    & 
     &             (DI(Iarn-1)*Varn(:,Iarn-1) + DR(Iarn-1)*Varn(:,Iarn))
                     Print*,'Checking Imag part of eigenvector ',Iarn-1
                     Print*,' A x_i - ( mu_i x_r + mu_r x_i ) '
                  ELSE
                     VarnC = VarnC - DR(Iarn)*Varn(:,Iarn)
                     Print*,'Checking eigenvector ',Iarn
                     Print*,' A x_r -  mu_r x_r ' 
                  END IF
                  Print*,'Norm of error vector: ',SQRT(SUM(VarnC*VarnC))
                  Print*,'error vector max component: ',                &
     &                 SQRT(MAXVAL(VarnC*VarnC))
                  Iarn = Iarn + 1
                  IF (IGM .EQ. ArnoldiEXIT) EXIT
!
!             ===================================================================
               END IF
!            ====================================================================
!
            ENDIF
!
!     **************************************************************************************
!         Main ARPACK Stuff ends here. Regular loops carry out Aq using existing code below 
!     **************************************************************************************
!
!     --------------------------------------------------------------------------------------
!         CALCULATE Aq
!     -------------------------------------------------------------------------------------
!
!           Calculate the main part of Aq
            DO IKT =1,NKT
               IF (LL(IKT))THEN
                  ZO(IKT      ) =  ZNEWT(IKT      ) + epsilon*q(IKT     ,IGM)
                  ZO(IKT+  NKT) =  ZNEWT(IKT+  NKT) + epsilon*q(IKT+  NKT,IGM)
                  ZO(IKT+2*NKT) =  ZNEWT(IKT+2*NKT) + epsilon*q(IKT+2*NKT,IGM)
                  ZO(IKT+3*NKT) =  ZNEWT(IKT+3*NKT) + epsilon*q(IKT+3*NKT,IGM)
               ENDIF
            ENDDO

            IF(arclcFLAG2) THEN
               v2 = v2N + epsilon*qN4(IGM)
            ENDIF
            TSTART = 0._rk
            DELT = TRUEDELT
            IF(parFLAG .eq.1)THEN
               Print*,'Re = ',1._rk / v2
               v1=v2
               
            ELSE IF(parFLAG .eq.2) THEN
               Print*,'Ri = ', v2
               Ri = v2
            ELSE IF(parFLAG .eq.3) THEN
               Print*,'alpha = ', v2
               alpha = v2
               call makeK(KX,KY,KZ)
            ENDIF

            DO IKT =1,NKT
               IF (LL(IKT))THEN
                  XII(IKT) = ZO(IKT)
                  ETA(IKT) = ZO(NKT+IKT)
                  ZET(IKT) = ZO(2*NKT+IKT)
                  RHO(IKT) = ZO(3*NKT+IKT)/scale
               ENDIF
            ENDDO

      CALL TIMESTEP_CUDA(XII,ETA,ZET,RHO,KX,KY,KZ,TIME,TSTART,TNEWT,AMPFOR,DELT,ResidualThreshold,&
     &                   FK,v1,Ri,Sc,Theta,alpha,Dtarget,IKF,IKN,L,IKTX,IKTY,IKTZ,KTZ,NKT,NX,NY,NZ,&
     &                   NOUT,NOUTV,STATSFLAG,RCFLAG,ADAPTFLAG,RANK)

!        Store final Z state 
         DO IKT=1,NKT
            IF(LL(IKT))THEN
               ZN(IKT)       = XII(IKT)
               ZN(NKT+IKT)   = ETA(IKT)
               ZN(2*NKT+IKT) = ZET(IKT)
               ZN(3*NKT+IKT) = RHO(IKT)*scale
            ENDIF
         ENDDO

!    -------------------------------------------------------------------------------------
!
!     ************************************************************************************
!     That's all ARPACK needs from the GMRES routines
!      so do some extra Stability stuff and skip rest of GMRES loop      
            IF(ArnoldiFLAG) THEN
               DO IKT = 1, NKT
                  IF(.NOT.LL(IKT))     CYCLE
                  SHIFT = EXP(-ZI*KX(IKT)*SNEWTX)*EXP(-ZI*KY(IKT)*SNEWTY)*EXP(-ZI*KZ(IKT)*SNEWTZ)
                  v(IKT) = SHIFT*( ZN(IKT)-ZNN(IKT) )/epsilon
                  v(NKT+IKT) =SHIFT*( ZN(NKT+IKT)-ZNN(NKT+IKT) )/epsilon
                  v(2*NKT+IKT) =SHIFT*( ZN(2*NKT+IKT)-ZNN(2*NKT+IKT) )/epsilon
                  v(3*NKT+IKT) =SHIFT*( ZN(3*NKT+IKT)-ZNN(3*NKT+IKT) )/epsilon
               ENDDO
               IF(ArnoldiEXIT .EQ. 0) THEN
                  call v2workd(v, workd(ipntr(2)))
               ELSE
                  call v2workd(v, VarnC)
               ENDIF
               CYCLE
            END IF
!     *************************************************************************************
!
!   ------------------------------------------------------------------------------------

!     Calculate rows 1 to N of vector v in GMRES
            DO IKT = 1, NKT
               IF(.NOT.LL(IKT))     CYCLE
               SHIFT = EXP(-ZI*KX(IKT)*SNEWTX)*EXP(-ZI*KY(IKT)*SNEWTY)*EXP(-ZI*KZ(IKT)*SNEWTZ) 
               DO IT = 1,4
                  IK = (IT-1)*NKT +IKT
                  v(IK) = SHIFT*( ZN(IK)-ZNN(IK) )/epsilon - q(IK,IGM) &
                       &                                         + dZNdSx(IK)*qN2(IGM) &
                       &                                         + dZNdSz(IK)*qN3(IGM)
                  IF(.NOT.FixedP) v(IK) = v(IK) + dZNdt(IK)*qN1(IGM)
               ENDDO
            ENDDO
!
!           Calculate v(N+1), v(N+2) and  v(N+3)      
            IF(.NOT.FixedP) vN1 = SUM( REAL(CONJG(dzdt)*q(:,IGM)),L3 )
            vN2 = SUM( REAL(CONJG(dzdSx)*q(:,IGM)),L3 )
            vN3 = SUM( REAL(CONJG(dzdSz)*q(:,IGM)),L3 )
!
!           If arc-length continuation calculate v(N+3)                                                                              
            IF(arclcFLAG2) THEN
               vN4 = SUM( REAL(CONJG(dZrdr)*q(:,IGM)),L3 )            &
     &           + dTrdr*qN1(IGM) + dSxrdr*qN2(IGM)+ dSzrdr*qN3(IGM) + dv2rdr*qN4(IGM)
            ENDIF
!   ------------------------------------------------------------------------------------
!
!           Carry out Orthogonalization of v
            DO IORTH=1,IGM
              H(IORTH,IGM) = SUM( REAL(CONJG(q(:,IORTH)) * v ),L3 ) 
             IF(.NOT.FixedP) H(IORTH,IGM) = H(IORTH,IGM) +qN1(IORTH)*vN1
                             H(IORTH,IGM) = H(IORTH,IGM) +qN2(IORTH)*vN2
                             H(IORTH,IGM) = H(IORTH,IGM) +qN3(IORTH)*vN3
             IF(arclcFLAG2) H(IORTH,IGM) = H(IORTH,IGM) +qN4(IORTH)*vN4
!
              WHERE (L3)        v = v   - H(IORTH,IGM)*q(:,IORTH)
              IF(.NOT.FixedP) vN1 = vN1 - H(IORTH,IGM)*qN1(IORTH)
                              vN2 = vN2 - H(IORTH,IGM)*qN2(IORTH)
                              vN3 = vN3 - H(IORTH,IGM)*qN3(IORTH)
              IF(arclcFLAG2)  vN4 = vN4 - H(IORTH,IGM)*qN4(IORTH)
             ENDDO
!   -----------------------------------------------------------------------------------
!
!        This part calculates: H(IGM+1,IGM) =  ||[ v(:,:) vN1 ] ||
            H(IGM+1,IGM) =  SUM( REAL(v*CONJG(v)),L3 )
            IF(.NOT.FixedP) H(IGM+1,IGM) = H(IGM+1,IGM) + vN1*vN1 
                            H(IGM+1,IGM) = H(IGM+1,IGM) + vN2*vN2
                            H(IGM+1,IGM) = H(IGM+1,IGM) + vN3*vN3
            IF(arclcFLAG2)  H(IGM+1,IGM) = H(IGM+1,IGM) + vN4*vN4
            H(IGM+1,IGM) =  SQRT(H(IGM+1,IGM))
!  -----------------------------------------------------------------------------------
!
!           Set vector for next GMRES step
            WHERE (L3) q(:,IGM+1) = v/H(IGM+1,IGM) 
            IF(.NOT.FixedP) qN1(IGM+1) = vN1/H(IGM+1,IGM) 
                            qN2(IGM+1) = vN2/H(IGM+1,IGM) 
                            qN3(IGM+1) = vN3/H(IGM+1,IGM) 
            IF(arclcFLAG2)  qN4(IGM+1) = vN4/H(IGM+1,IGM)
!
! -------------------------------------------------------------------------------------   
!     Do SVD of matrix H to determine if GMRES should stop
!
!           SVD destroys H, so make a copy
            Hcopy(1:IGM+1,1:IGM) = H(1:IGM+1,1:IGM) 
!
!           Call LAPACK SVD routine (this is a thin SVD)
            CALL DGESDD('A',IGM+1,IGM,Hcopy(1:IGM+1,1:IGM),IGM+1,       &
     &  Ssvd(1:IGM),Usvd(1:IGM+1,1:IGM+1),IGM+1,VHsvd(1:IGM,1:IGM),IGM, &
     &  WORKsvd,LWORKsvd,IWORKsvd,INFOsvd) 
!
! --------------------------------------------------------------------------------------
!       Invert H using SVD decompositon
!     
!           Set P = ||b|| U' e1
            P(1:IGM+1) = NORMb*Usvd(1,1:IGM+1)
!
!           Calculate X
            Xh(1:IGM) = P(1:IGM)/Ssvd(1:IGM)
!  
!  ---------------------------------------------------------------------------------------
!    Calculate GMRES residual
!
!    Set GMR = || H y - ||b|| e1 || / ||b|| = || Ssvd * Xh - P || / ||b||
            GMR = SQRT( SUM( (Ssvd(1:IGM)*Xh(1:IGM)-P(1:IGM))           &
     &               *(Ssvd(1:IGM)*Xh(1:IGM)-P(1:IGM)) )                &
     &          +  P(IGM+1)*P(IGM+1)  )/NORMb 
!
            print*,'GMRES step ',IGM
            print*,'| Ssvd * Xh - P | / |b| = ',GMR
!
            Yh(1:IGM) = matmul(                                         &
     &             TRANSPOSE( VHsvd(1:IGM,1:IGM) )  , Xh(1:IGM)   )
            Yh(1:IGM+1) = matmul( H(1:IGM+1,1:IGM),Yh(1:IGM) )
            Yh(1) = Yh(1) - NORMb
            GMR = SQRT( SUM( Yh(1:IGM+1)                                &
     &                            *Yh(1:IGM+1) ) )/NORMb 
!
            print*,' | H y - |b| e1 | / |b| = ',GMR
!
! ----------------------------------------------------------------------------------------
!
!    If residual is low enough exit GMRES loop
            IF (GMR.LE.GMTOL)  THEN
               GMREScount = GMREScount + IGM
               EXIT
            END IF
!
         ENDDO
!   ________________________________________________________________________________________
!   ----------------------------------------------------------------------------------------
!         GMRES Loop has now ended
!   ________________________________________________________________________________________
!   -----------------------------------------------------------------------------------------
!
!       *********************************
!        very last bit of arnoldi stuff
         IF (ArnoldiFLAG) EXIT
!      *********************************
! 
! -------------------------------------------------------------------------------------
!       Check if GMRES converged, if not exit
!           
         IF (IGM.EQ.GMRESMAX+1) THEN
           Print*,'GMRES did not converge, exiting Newton iteration'
           GMREScount = GMREScount + IGM
           WRITE(150,1502) IREST,'    ',                                 &
     &              TSTARTin,v2in,TNEWTin,SNEWTXin,SNEWTYin,SNEWTZin,NERin, &
     & ' : "GMRES FAILED"  : ',iuposout+1,v2N,TNEWT,SNEWTX,SNEWTY,SNEWTZ, &
     &                 NER,INE,GMREScount,HOOKcount  
            EXIT
         END IF
!
         Print*,'-----------------------------------------------------'
         Print*,'GMRES converged on GMRES step ',IGM
!
! -------------------------------------------------------------------------------------
! |||||||||    HOOK STEP          ||||||||||||||||||||||||||||||||||||||||
! -------------------------------------------------------------------------------------
! 
         IHOOK = 0
         DO   ! this is a DO WHILE type of DO loop. The exit statements are at the end
!
            IHOOK = IHOOK + 1
            IF(IHOOK.GT.IHOOKMAX) THEN
               Print*,'HOOKSTEP failed to improve residual, '            
               Print*,'exiting Newton iteration'
               WRITE(150,1502) IREST,'    ',                             &
     &              TSTARTin,v2in,TNEWTin,SNEWTXin,SNEWTYin,SNEWTZin,NERin, &
     & ' : "H-STEP FAILED" : ',iuposout+1,v2N,TNEWT,SNEWTX,SNEWTY,SNEWTZ, &
     &                 NER,INE,GMREScount,HOOKcount    
               EXIT
            END IF
            Print*,'                   '
            Print*,'Performing new Hookstep ...'
! -------------------------------------------------------------------------------
!           FIND MU
!
!           Calculate X
            Xh(1:IGM) = P(1:IGM)/Ssvd(1:IGM)
!
!           First check ||x|| with mu = 0
            normX = SQRT( SUM( Xh(1:IGM)*Xh(1:IGM) ) )
            Print*,'Starting ||x|| = ',normX,', Delta = ',Delta
!          
!           Set MU = 0. for print statement below
            MUmid = 0._rk
!
            IF (normX.GT.root2*Delta) THEN 
!              Set intital value for MU1 = 0 and MU2 = maxval(Ssvd)
               MU1 = 0._rk
               MU2 = MAXVAL(Ssvd(1:IGM))
!
!              Recalculate ||x||
               Xh(1:IGM) = P(1:IGM)*Ssvd(1:IGM)                         &
     &              /( MU2 + Ssvd(1:IGM)*Ssvd(1:IGM) )
               normX = SQRT( SUM( Xh(1:IGM)*Xh(1:IGM) ) )
!
!              If ||x|| is still greater than delta then double MU and try again...
               DO 
               IF(normX.LT.Delta) EXIT
                  MU1 = MU2
                  MU2 = 2._rk*MU2
!                 Recalculate ||x||
                  Xh(1:IGM) = P(1:IGM)*Ssvd(1:IGM)                      &
     &                 /( MU2 + Ssvd(1:IGM)*Ssvd(1:IGM) )
                  normX = SQRT( SUM( Xh(1:IGM)*Xh(1:IGM) ) )
               ENDDO
!
!             normX has just been calculated with MU2 so set MUmid = MU2
               MUmid = MU2
!
!      Now MU1 should give ||x|| greater than delta and MU2 should give ||x|| less than delta
! 
!      Bisection root finding algorithm to find MU such that delta/root2 < ||x|| < root2*delta              
               DO
               IF(normX.LT.root2*Delta .AND. normX.GT.Delta/root2) EXIT
                  MUmid = 0.5_rk*(MU2+MU1)
!                 Recalculate ||x||
                  Xh(1:IGM) = P(1:IGM)*Ssvd(1:IGM)                      &
     &                 /( MUmid + Ssvd(1:IGM)*Ssvd(1:IGM) )
                  normX = SQRT( SUM( Xh(1:IGM)*Xh(1:IGM) ) )
                  IF (normX.GT.Delta) THEN
                     MU1 = MUmid
                  ELSE
                     MU2 = MUmid
                  ENDIF
               ENDDO
!
            END IF
!          
            print*,'Final ||X|| = ',normX,' Mu = ',MUmid
!  ------------------------------------------------------------------------
!        NOW MU has been found check if Delta was appropriate.
!
!        Use VHsvd to calculate y = Vx
         Yh(1:IGM) = matmul(                                            &
     &          TRANSPOSE( VHsvd(1:IGM,1:IGM) ) , Xh(1:IGM)   )
!
         DIFFZ = CMPLX(0._rk,0._rk,rk)
         DIFFT = 0._rk
         DIFFSx = 0._rk
         DIFFSz = 0._rk
         IF(arclcFLAG2) DIFFv2 = 0._rk
         DO IORTH=1,IGM
            WHERE (L3) DIFFZ = DIFFZ +                                  &
     &                        q(:,IORTH) * CMPLX(Yh(IORTH),0._rk,rk)
            IF(.NOT.FixedP) DIFFT = DIFFT + qN1(IORTH) * Yh(IORTH) 
            DIFFSx = DIFFSx + qN2(IORTH) * Yh(IORTH) 
            DIFFSz = DIFFSz + qN3(IORTH) * Yh(IORTH) 
            IF(arclcFLAG2) DIFFv2 = DIFFv2 + qN4(IORTH) * Yh(IORTH)
         ENDDO
!
! --------------------------------------------------------------------------------
!      Calculate J*DIFF
!
         Print*,'Calculating J*DIFF ...'
!           Calculate the main part of JF*DIFF
            DO IKT =1,NKT
               IF (LL(IKT))THEN
                  ZO(IKT      ) =  ZNEWT(IKT      ) + epsilon*DIFFZ(IKT)
                  ZO(IKT+  NKT) =  ZNEWT(IKT+  NKT) + epsilon*DIFFZ(IKT+  NKT)
                  ZO(IKT+2*NKT) =  ZNEWT(IKT+2*NKT) + epsilon*DIFFZ(IKT+2*NKT)
                  ZO(IKT+3*NKT) =  ZNEWT(IKT+3*NKT) + epsilon*DIFFZ(IKT+3*NKT)
               ENDIF
            ENDDO

            IF(arclcFLAG2) THEN
               v2 = v2N + epsilon*DIFFv2
            IF(parFLAG .eq.1)THEN
               Print*,'Re = ',1._rk / v2
               v1=v2
               
            ELSE IF(parFLAG .eq.2) THEN
               Print*,'Ri = ', v2
               Ri = v2
            ENDIF
            ELSE IF(parFLAG .eq.3) THEN
               Print*,'alpha = ', v2
               alpha = v2
               call makeK(KX,KY,KZ)
            ENDIF
!
            TSTART = 0._rk
            DELT = TRUEDELT

            TSTART = 0._rk
            DELT = TRUEDELT

            DO IKT =1,NKT
               IF (LL(IKT))THEN
                  XII(IKT) = ZO(IKT)
                  ETA(IKT) = ZO(NKT+IKT)
                  ZET(IKT) = ZO(2*NKT+IKT)
                  RHO(IKT) = ZO(3*NKT+IKT)/scale
               ENDIF
            ENDDO

      CALL TIMESTEP_CUDA(XII,ETA,ZET,RHO,KX,KY,KZ,TIME,TSTART,TNEWT,AMPFOR,DELT,ResidualThreshold,&
     &                   FK,v1,Ri,Sc,Theta,alpha,Dtarget,IKF,IKN,L,IKTX,IKTY,IKTZ,KTZ,NKT,NX,NY,NZ,&
     &                   NOUT,NOUTV,STATSFLAG,RCFLAG,ADAPTFLAG,RANK)

!        Store final Z state 
         DO IKT=1,NKT
            IF(LL(IKT))THEN
               ZN(IKT)       = XII(IKT)
               ZN(NKT+IKT)   = ETA(IKT)
               ZN(2*NKT+IKT) = ZET(IKT)
               ZN(3*NKT+IKT) = RHO(IKT)*scale
            ENDIF
         ENDDO

!     Calculate rows 1 to N of vector v in GMRES
            DO IKT = 1, NKT
               IF(.NOT.LL(IKT))     CYCLE
               SHIFT= EXP(-ZI*KX(IKT)*SNEWTX)*EXP(-ZI*KY(IKT)*SNEWTY)*EXP(-ZI*KZ(IKT)*SNEWTZ)

               DO IT = 1,4
                  IK = (IT-1)*NKT +IKT
                  
                  JFDIFF(IK) = SHIFT*( ZN(IK)-ZNN(IK) )/epsilon - DIFFZ(IK) &
                       &      + dZNdSx(IK)*DIFFSx &
                       &      + dZNdSz(IK)*DIFFSz
                  IF(.NOT.FixedP) JFDIFF(IK) = JFDIFF(IK) + dZNdt(IK)*DIFFT
               ENDDO
            ENDDO
!
!     Approimate gradient of NORMFN with respect to z
!     If this is small it suggests a minimum in NORMFN has been found 
            Print*,'|| J*DIFF || / || DIFF || = ',                          &
     &          SUM( REAL(JFDIFF*CONJG(JFDIFF)), L3 ) /                      &
     &          SUM( REAL(DIFFZ*CONJG(DIFFZ)), L3 )
!
            WHERE (L3) FNpred = FN + JFDIFF
            NORMFNpred = SUM( REAL(FNpred*CONJG(FNpred)), L3 ) 
            NORMFNpred = SQRT(NORMFNpred)
            Print*,'Predicted || F{N+1} || = || FN + J*DIFF || =',      &
     &                                                        NORMFNpred
!
            WHERE (L3) FNpredC = FN + const*JFDIFF
            NORMFNpredC = SUM( REAL(FNpredC*CONJG(FNpredC)), L3 ) 
            NORMFNpredC = SQRT(NORMFNpredC)
            Print*,'|| FN + const*J*DIFF || = ',NORMFNpredC
!
            Print*,'|| FN || = ',NORMFN 
!
            NORMbpredC = NORMFNpredC
!
            IF(arclcFLAG2) NORMbpredC =                                  &
     &                       SQRT(NORMFNpredC*NORMFNpredC + Narc*Narc)
!
! ------------------------------------------------------------------------------
!     Calculate UPO error with updated Z and T
!
!        Calculate potential new Z & T
         WHERE (L3) Zplus = ZNEWT + DIFFZ
         Tplus = TNEWT + DIFFT
         Sxplus = SNEWTX + DIFFSx
         Szplus = SNEWTZ + DIFFSz
         IF(arclcFLAG2)THEN
            v2plus = v2N + DIFFv2
         ENDIF
            !
!        Check if delta is appropriate
!
!        Set ZO = ZNEWT (ZO gets updated during timestepping)
         WHERE (L3) ZO = Zplus  
         IF(arclcFLAG2) THEN
            v2 = v2plus
            !        update NU                                                                                                     
            IF(parFLAG .eq.1)THEN
               Print*,'Re = ',1._rk / v2
               v1=v2
               
            ELSE IF(parFLAG .eq.2) THEN
               Print*,'Ri = ', v2
               Ri = v2
            ELSE IF(parFLAG .eq.3) THEN
               Print*,'alpha = ', v2
               alpha = v2
               call makeK(KX,KY,KZ)
            ENDIF
         ENDIF
!
         TSTART = 0._rk
         DELT = TRUEDELT
         print*,'Calculating UPO'
         Print*,'Period = ',Tplus,', dt = ',DELT


            DO IKT =1,NKT
               IF (LL(IKT))THEN
                  XII(IKT) = ZO(IKT)
                  ETA(IKT) = ZO(NKT+IKT)
                  ZET(IKT) = ZO(2*NKT+IKT)
                  RHO(IKT) = ZO(3*NKT+IKT)/scale
               ENDIF
            ENDDO

      CALL TIMESTEP_CUDA(XII,ETA,ZET,RHO,KX,KY,KZ,TIME,TSTART,Tplus,AMPFOR,DELT,ResidualThreshold,&
     &                   FK,v1,Ri,Sc,Theta,alpha,Dtarget,IKF,IKN,L,IKTX,IKTY,IKTZ,KTZ,NKT,NX,NY,NZ,&
     &                   NOUT,NOUTV,STATSFLAG,RCFLAG,ADAPTFLAG,RANK)

!        Store final Z state
         DO IKT=1,NKT
            IF(LL(IKT))THEN
               ZN(IKT)       = XII(IKT)
               ZN(NKT+IKT)   = ETA(IKT)
               ZN(2*NKT+IKT) = ZET(IKT)
               ZN(3*NKT+IKT) = RHO(IKT)*scale
            ENDIF
         ENDDO
!
!        Shift final state back in x-direction
         DO IKT = 1, NKT
            IF(.NOT.LL(IKT))     CYCLE
            SHIFT = EXP(-ZI*KX(IKT)*Sxplus)*EXP(-ZI*KY(IKT)*SNEWTY)*EXP(-ZI*KZ(IKT)*Szplus)
            DO IT = 1,4
               IK = (IT-1)*NKT +IKT
               ZNSHIFT(IK) =SHIFT*ZN(IK)
            ENDDO
         ENDDO
!
!        Calculate F = Z(z,T) - z
         WHERE (L3) FNplus = ZNSHIFT - Zplus       
!
         NORMFNplus = SQRT( SUM( REAL(FNplus*CONJG(FNplus)), L3 ) )
         Print*,'|| FN{N+1} || = ',NORMFNplus
!
         NORMbplus = NORMFNplus
         IF (arclcFLAG2) THEN
            Narcplus = SUM( REAL(dZrdr*CONJG(Zplus-Zr2)), L3 )           &
                 &    + dSxrdr*(Sxplus-Sxr2)+ dSzrdr*(Szplus-Szr2) + dv2rdr*(v2plus-v2r2)                     &
                 &    - Ndelta_r
            IF(.NOT.FixedP) Narcplus = Narcplus + dTrdr*(Tplus-Tr2)
!                                                                                                                  
            NORMbplus = SQRT(NORMFNplus*NORMFNplus + Narcplus*Narcplus)
         END IF
         Print*,'|| b || = ',NORMb   
         Print*,'|| b{N+1} || = ',NORMbplus  
!
! ----------------------------------------------------------------------------------------------
!       Now check new solution to see if it is better than the old one
!
        IF (NORMFNplus.NE.NORMFNplus .OR. NORMFNpred.NE.NORMFNpred) THEN                 
           Print*,'Time-stepping routine produced NAN'                                      
           Print*,'Decrease trust region (DELTA = DELTA/2 )',                &              
     &                                             ' and try again ...'            
           DELTA = DELTA / 2._rk                                                            
           halveFLAG = .TRUE.                                                               
           CYCLE                                                                            
        END IF                  
!
        IF (NORMbplus .GE. NORMbpredC) THEN
            Print*,'||b{N+1}|| is .GE. || b ||'
!           Residual has increased
!
!           If a satisfactory solution has already been found then 
!           revert back to it and exit Hook step
            IF (doubleFLAG) THEN
               Print*,'Reverting back to previous saved solution'
               WHERE (L3) Zplus = Zplus2
               Tplus = Tplus2
               Sxplus = Sxplus2
               Szplus = Szplus2
               IF(arclcFLAG2) v2plus = v2plus2
               WHERE (L3) ZN = ZN2
               WHERE (L3) FNplus = FNplus2
               NORMFNplus = NORMFNplus2
               DELTA = DELTA / 2._rk
               EXIT
            END IF
!               
!           If not decrease trust region 
            Print*,'Decrease trust region (DELTA = DELTA/2 )',                &
     &                                              ' and try again ...'
            DELTA = DELTA / 2._rk
            halveFLAG = .TRUE.
         ELSE 
            Print*,'||b{N+1}|| is .LT. || b ||'
!           Residual has decreased
!
            IF (halveFLAG) THEN
               print*,'Accept current solution'
               EXIT
            END IF
!           
            IF (doubleFLAG .AND. NORMFNplus.GT.NORMFNplus2) THEN 
               Print*,'New ||F|| is greater than previous ||F||'
               Print*,'Reverting back to previous saved solution'
               WHERE (L3) Zplus = Zplus2
               Tplus = Tplus2
               Sxplus = Sxplus2
               Szplus = Szplus2
               WHERE (L3) ZN = ZN2
               WHERE (L3) FNplus = FNplus2
               NORMFNplus = NORMFNplus2
               DELTA = DELTA / 2._rk
               EXIT
            END IF
!            
            Print*,'Considering whether to increase trust region ...'
            IF (MUmid .EQ. 0._rk) THEN
               Print*,'Mu = 0 => ',                                      &
     &              'Trust region is already larger than step size'
               Print*,'Accept current solution'
               EXIT
            ELSE IF ( ABS(NORMFNpred-NORMFNplus) .LT.                   &
     &                         0.1_rk*ABS(NORMFNplus-NORMFN) ) THEN
               Print*,'Prediction is too good, increasing trust region'
               DELTA = DELTA * 2._rk
               doubleFLAG = .TRUE.
               Print*,'Saving current solution before retrying'
               WHERE (L3) Zplus2 = Zplus
               Tplus2 = Tplus
               Sxplus2 = Sxplus               
               Szplus2 = Szplus               
               IF(arclcFLAG2) v2plus2 = v2plus
               WHERE (L3) ZN2 = ZN
               WHERE (L3) FNplus2 = FNplus
               NORMFNplus2 = NORMFNplus
            ELSE IF (NORMFNPLUS .LT. NORMFNpred) THEN
               Print*,'Negative curvature, increasing trust region'
               DELTA = DELTA * 2._rk
               doubleFLAG = .TRUE.
               Print*,'Saving current solution before retrying'
               WHERE (L3) Zplus2 = Zplus
               Tplus2 = Tplus
               Sxplus2 = Sxplus
               Szplus2 = Szplus
               IF(arclcFLAG2) v2plus2 = v2plus
               WHERE (L3) ZN2 = ZN
               WHERE (L3) FNplus2 = FNplus
               NORMFNplus2 = NORMFNplus
            ELSE
               Print*,'Accept current solution'
               EXIT
            END IF
         ENDIF
      ENDDO
      HOOKcount = HOOKcount + IHOOK
!
!! ________________________________________________________________________
! ????????????????????????????????????????????????????????????????????????
! |||||||||    HOOK STEP FINISHED   ||||||||||||||||||||||||||||||||||||||
! ________________________________________________________________________
!
!
      halveFLAG = .FALSE.
      doubleFLAG = .FALSE.
!
!     Check if Timestepping has blown up
      IF (NORMFNpred.NE.NORMFNpred .OR. NORMFNplus.NE.NORMFNplus) THEN 
         Print*,'Time-stepping routine produced NAN'
         IF (TRUEDELT .LT. 0.1_rk*TrueTRUEDELT) THEN
            Print*,'Exiting Newton iteration'
            WRITE(150,1502) IREST,'    ',                                &
     &              TSTARTin,v2in,TNEWTin,SNEWTXin,SNEWTYin,SNEWTZin,NERin, &
     & ' : "T-STEP FAILED" : ',iuposout+1,v2N,TNEWT,SNEWTX,         &
     &                 SNEWTY,SNEWTZ,NER,INE,GMREScount,HOOKcount 
            EXIT
         ENDIF
         Print*,'Reducing dt and retrying newton step.'
         TRUEDELT = TRUEDELT/2._rk
         changeDELT = .TRUE.
!        Shift final state back in x-direction
         DO IKT = 1, NKT
            IF(.NOT.LL(IKT))     CYCLE
            SHIFT = EXP(-ZI*KX(IKT)*SNEWTX)*EXP(-ZI*KY(IKT)*SNEWTY)*EXP(-ZI*KZ(IKT)*SNEWTZ)
            DO IT = 1,4
               IK = (IT-1)*NKT +IKT
               ZNSHIFT(IK) = SHIFT*ZNN(IK)
            ENDDO
         ENDDO
         CYCLE
      ENDIF
!
      IF(IHOOK.GT.IHOOKMAX) EXIT   
!
! --------------------------------------------------------------------------------
!
!        Store final Z state 
         WHERE (L3) ZNN = ZN
!
!        Update Z and T for next newton step 
         WHERE (L3) ZNEWT = Zplus
         TNEWT = Tplus
         SNEWTX = Sxplus
         SNEWTZ = Szplus
         IF(arclcFLAG2) v2N   = v2plus
!
!        Shift final state back in x-direction
        DO IKT = 1, NKT
            IF(.NOT.LL(IKT))     CYCLE
            SHIFT = EXP(-ZI*KX(IKT)*SNEWTX)*EXP(-ZI*KY(IKT)*SNEWTY)*EXP(-ZI*KZ(IKT)*SNEWTZ)
            DO IT = 1,4
               IK = (IT-1)*NKT +IKT
               ZNSHIFT(IK) = SHIFT*ZNN(IK)
            ENDDO
         ENDDO
!     
         NORMZ = SQRT( SUM( REAL(ZNN*CONJG(ZNN)), L3 ) +SUM( REAL(ZNEWT*CONJG(ZNEWT)), L3 ))
!
!      Calculate 1st residual
         WHERE (L3) FN = FNplus
         NORMFN = NORMFNplus
         NORMb = NORMFN
         Print*,'||FN|| = ',NORMFN
!
         IF(arclcFLAG2) THEN
!        Calculate new Narc                                                                                                          
            Narc = SUM( REAL(dZrdr*CONJG(ZNEWT-Zr2)), L3 )              &
     &   + dSxrdr*(SNEWTX-Sxr2)+ dSzrdr*(SNEWTZ-Szr2) + dv2rdr*(v2N-v2r2)                        &
     &   - Ndelta_r
            IF(.NOT.FixedP) Narc = Narc + dTrdr*(TNEWT-Tr2)
            Print*,'Narc = ',Narc
            NORMb = SQRT(NORMFN*NORMFN + Narc*Narc)
         END IF
         Print*,'||b|| = ',NORMb
!
!
!        Use NormFN to check for convergence 
         NER = NORMFN/NORMZ
         Print*,'  '
         Print*,'-----------------------------------------------------'
         PRINT*,'Newton convergence check: ||F|| / ||Z|| = ',NER
         Print*,'-----------------------------------------------------'
         Print*,'  '
         Print*,'||Z(T)|| = ',SQRT( SUM( REAL(ZNN*CONJG(ZNN)), L3 ) )
         Print*,'||Z(0)|| = ',SQRT( SUM( REAL(Zplus*CONJG(Zplus)), L3 ))
         Print*,'  '
            DO IKT =1,NKT
               IF (LL(IKT))THEN
                  XII(IKT) = ZNEWT(IKT)
                  ETA(IKT) = ZNEWT(NKT+IKT)
                  ZET(IKT) = ZNEWT(2*NKT+IKT)
                  RHO(IKT) = ZNEWT(3*NKT+IKT)/scale
               ENDIF
            ENDDO

!                                                                                                                                                                                                                                           
!                                                                                                                                                                                                                                            
         rewind(99)
         write(99) 0._rk,v1,Ri,alpha,TNEWT,SNEWTX,SNEWTY,SNEWTZ,NER,                          &
              &                         (XII(IKT),IKT=1,NKT), &
              &                         (ETA(IKT),IKT=1,NKT), &
              &                         (ZET(IKT),IKT=1,NKT), &
              &                         (RHO(IKT),IKT=1,NKT)
         print*,'period =',TNEWT,'X shift =',SNEWTX,'Z shift =',SNEWTZ
         IF(arclcFLAG) Print*,'Re = ',1._rk / v2N
         Print*,'-----------------------------------------------------'
         
         NERold = abs((NER-NERold)/NER)
         if(NERold .lt. 0.001_rk)then
            NERcount = NERcount+1

            if(NERcount .gt. 5)then
               WRITE(150,1502) IREST,'    ',                                &
     &             TSTARTin,v2in,TNEWTin,SNEWTXin,SNEWTYin,SNEWTZin,NERin, &
     & ' : "NEWTON STALLED" : ',iuposout+1,v2N,TNEWT,SNEWTX,SNEWTY,SNEWTZ, &
     &                 NER,INE,GMREScount,HOOKcount
               EXIT
            END IF
         ELSE
            NERcount = 0
         END IF
         
         NERold = NER
!
!        If Newton method has converged then reset everything and exit
         IF (NER.LE.NETOL)  THEN
            Print*,'-------------------------------------------------'
            Print*,'-------------------------------------------------'
            Print*,'  '
            PRINT*,'Newton method converged in',INE,' steps.'         
            WRITE(150,1502) IREST,'    ',                                &
     &              TSTARTin,v2in,TNEWTin,SNEWTXin,SNEWTYin,SNEWTZin,NERin, &
     & ' :   "CONVERGED"   : ',iuposout+1,v2N,TNEWT,SNEWTX,         &
     &                 SNEWTY,SNEWTZ,NER,INE,GMREScount,HOOKcount
            convergeFLAG = 1
            IF(arclcFLAG) THEN
               CONVERGED = .TRUE.
               IF(SecondRun) SecondRun = .FALSE.
               IF(FirstRun) THEN
                  FirstRun = .FALSE.
                  SecondRun = .TRUE.
                  arclcFLAG2 = .TRUE.
               END IF
               WHERE (L3) Zr0 = Zr1
               Tr0 = Tr1
               Sxr0 = Sxr1
               Szr0 = Szr1
               v2r0 = v2r1
               WHERE (L3) Zr1 = Zr2
               Tr1 = Tr2
               Sxr1 = Sxr2
               Szr1 = Szr2
               v2r1 = v2r2
               delta_r1 = delta_r2
               WHERE (L3) Zr2 = ZNEWT
               Tr2 = TNEWT
               Sxr2 = SNEWTX
               Szr2 = SNEWTZ
               v2r2 = v2N
               delta_r2 = Ndelta_r
            END IF
            EXIT
         END IF
!
!        If max newton steps have been reached then admit failure and move on to next guess
         IF (INE.EQ.NEWTONMAX)  THEN
            WRITE(150,1502) IREST,'    ',                                &
     &             TSTARTin,v2in,TNEWTin,SNEWTXin,SNEWTYin,SNEWTZin,NERin, &
     & ' : "NEWTON FAILED" : ',iuposout+1,v2N,TNEWT,SNEWTX,SNEWTY,SNEWTZ, &
     &                 NER,INE,GMREScount,HOOKcount 
         END IF
         !
      ENDDO
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ----------------------------------------------------------------------------
!     NEWTON ITERATIONS END
! ----------------------------------------------------------------------------    
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!  
!     Write out final state of UPO to file
      if( convergeFLAG .eq. 1) then 
            DO IKT =1,NKT
               IF (LL(IKT))THEN
                  XII(IKT) = ZNEWT(IKT)
                  ETA(IKT) = ZNEWT(NKT+IKT)
                  ZET(IKT) = ZNEWT(2*NKT+IKT)
                  RHO(IKT) = ZNEWT(3*NKT+IKT)/scale
               ENDIF
            ENDDO
         write(999) 0._rk,v1,Ri,alpha,TNEWT,SNEWTX,SNEWTY,SNEWTZ,NER,                          &
              &                         (XII(IKT),IKT=1,NKT), &
              &                         (ETA(IKT),IKT=1,NKT), &
              &                         (ZET(IKT),IKT=1,NKT), &
              &                         (RHO(IKT),IKT=1,NKT)
         iuposout = iuposout + 1
         print*,'Wrote to UPOs.out number ',iuposout
         flush(999)
         print*,'period =',TNEWT,'X shift =',SNEWTX,'Z shift =',SNEWTZ
         IF(arclcFLAG) Print*,'Re = ',1._rk / v2N
         Print*,'-----------------------------------------------------'
         convergeFLAG = 0
      endif
!________
!
! WRAP UP
!________
!
      CALL CPU_TIME(toc)
      time3 = toc-tic
      print*,'    '
      write(6,5000) time3, time3/60._rk, time3/3600._rk
 5000  format(1x,'CPU time required for main loop = ',F10.3,' s = ',      &
      &               F7.3,' m = ',F7.3,' h.')
      print*,'    '
!
!
      CALL FLUSH(150)
!
      ENDDO
!
      close(99)
      close(150)
      close(999)
!
      IF (changeDELT) THEN
         print*,'WARNING: dt was changed for 1 or more NEWTON attempts.'
      END IF
!
      print*,'Finished'
!
 1501 FORMAT(A,A,A,A,A,A,A,A,A,A,A,A)  
 1502 FORMAT(I5,A,E16.6,E16.6,F12.5,F12.5,F12.5,F12.5,E14.5,                        & 
     &                      A,I4,E16.6,F10.5,F10.5,F10.5,F10.5,E14.5,I6,I6,I6)     
!
      END PROGRAM MAIN
!

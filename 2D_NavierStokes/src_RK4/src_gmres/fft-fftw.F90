!   _______________________________________________________
!
      SUBROUTINE KR_FFTW(ZK,ZR,FF2,PlanK2R)
!
!      CALLS SPECTRAL -> GRID POINT TRANSFORMS.
!
      USE GLOBAL_PARAMS
      IMPLICIT NONE
!
      INTEGER         :: IKX,IKY,INKY
      INTEGER(KIND=8) :: PlanK2R
!     
      REAL(rk)     :: ZR(NX,NY)
      COMPLEX(rk)  :: ZK(IKTX,IKTY),FF2(NX/2+1,NY)
!
!      CALL etime(T,tic);
!
!      On Kx=0, set -Ky values to conjugate of +Ky values
      DO IKY=1,KTY
         INKY = - IKY + 2*KTY + 2
         ZK(1,IKY) = CONJG(ZK(1,INKY))
      ENDDO
!
!      Map -Ky values to -Ky+NY and fill trucated wave numbers with zeros
      DO IKY=1,KTY+1
         DO IKX=1,IKTX
            FF2(IKX,IKY) = ZK(IKX,IKY+KTY)
         ENDDO
      ENDDO
      DO IKY=KTY+2,NY-KTY
         DO IKX=1,IKTX
            FF2(IKX,IKY) = CMPLX(0._rk,0._rk,rk)
         ENDDO
      ENDDO
      DO IKY=NY-KTY+1,NY
         DO IKX=1,IKTX
            FF2(IKX,IKY) = ZK(IKX,IKY+KTY-NY)
         ENDDO
      ENDDO
!
!     fill remaining trucated wave numbers with zeros
      DO IKY=1,NY
         DO IKX=IKTX+1,NX/2+1 
            FF2(IKX,IKY) = CMPLX(0._rk,0._rk,rk)
         ENDDO
      ENDDO!
!             
      call dfftw_execute(PlanK2R)
!
!       CALL etime(T,toc);
!       MYTIME = MYTIME+TOC-TIC
!
!
      END SUBROUTINE KR_FFTW
!
!
!
! _______________________________________________________________
!
      SUBROUTINE RK_FFTW(ZR,ZK,FF2,PlanR2K)
!
!      CALLS GRID POINT -> SPECTRAL TRANSFORMS.
!
      USE GLOBAL_PARAMS
      IMPLICIT NONE
!
      INTEGER         :: IKX,IKY
      INTEGER(KIND=8) :: PlanR2K
!     
      REAL(rk)     :: ZR(NX,NY),iNxNy
      COMPLEX(rk)  :: ZK(IKTX,IKTY),FF2(NX/2+1,NY)
!              
!      CALL etime(T,tic);
!
!
      call dfftw_execute(PlanR2K)
!
!
      iNxNy = 1._rk/REAL(NX*NY,rk)
!
      DO IKY=1,KTY+1
         DO IKX=1,IKTX
            ZK(IKX,IKY+KTY)  = FF2(IKX,IKY)*iNxNy 
         ENDDO
      ENDDO
      DO IKY=NY-KTY+1,NY
         DO IKX=1,IKTX
            ZK(IKX,IKY+KTY-NY) = FF2(IKX,IKY)*iNxNy
         ENDDO
      ENDDO
!
!
!       CALL etime(T,toc);
!       MYTIME = MYTIME+TOC-TIC
!
!
      END SUBROUTINE RK_FFTW

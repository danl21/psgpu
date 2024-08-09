!
      SUBROUTINE PROBE(ZK,FF2,UK,VK,UR,VR,ZR,TIME)
!
!     CALCULATES UR,VR and ZR at x=2*PI/alpha/SQRT(2), y=2*PI/SQRT(3)
!     and stores them in file probe.dat
!
      USE GLOBAL_PARAMS
      IMPLICIT NONE 
!
      INTEGER     :: IKX,IKY,I,J
      INTEGER     :: IKX1,IKY1,IKX2,IKY2,IKX3,IKY3
!
      REAL(rk)    :: KX,KY,WK
      REAL(rk)    :: UR(NX,NY),VR(NX,NY)
      REAL(rk)    :: ZR(NX,NY),WEIGHT
      REAL(rk)    :: TIME
!
      COMPLEX(rk) :: ZK(IKTX,IKTY)
      COMPLEX(rk) :: UK(IKTX,IKTY),VK(IKTX,IKTY)
      COMPLEX(rk) :: FF2(NX/2+1,NY)
!
!
      DO IKY = 1, IKTY
         KY =  REAL(IKY - KTY - 1,rk)
         DO IKX = 1, IKTX
            KX =  REAL(IKX-1,rk)*alpha 
            WK = MAX( KX*KX+KY*KY ,0.001_rk ) 
            UK(IKX,IKY)  = + L(IKX,IKY)*ZI*KY*ZK(IKX,IKY)/WK
            VK(IKX,IKY)  = - L(IKX,IKY)*ZI*KX*ZK(IKX,IKY)/WK
         ENDDO
      ENDDO
!
      CALL KR_FFTW(ZK,ZR,FF2,PlanK2ZR)
      CALL KR_FFTW(UK,UR,FF2,PlanK2UR)
      CALL KR_FFTW(VK,VR,FF2,PlanK2VR)
!
!      IKX = NINT(NX/SQRT(2.))
!      IKY = NINT(NY/SQRT(3.))
!      WRITE(51,5045) Time,UR(IKX,IKY),VR(IKX,IKY),ZR(IKX,IKY)
!
      IKX1 = NINT(NX/SQRT(2.0))
      IKY1 = NINT(NY/SQRT(3.0))
      IKX2 = NINT(NX/SQRT(6.0))
      IKY2 = NINT(NY/SQRT(17.0))
      IKX3 = NINT(NX/SQRT(15.0))
      IKY3 = NINT(NY/SQRT(1.5))
      WRITE(51,5045) DELT,ZR(IKX1,IKY1),ZR(IKX2,IKY2),ZR(IKX3,IKY3)
!
5045  FORMAT(1X,E27.20,4X,E27.20,4X,E27.20,4X,E27.20)       
!
      END SUBROUTINE PROBE 
!
!
!

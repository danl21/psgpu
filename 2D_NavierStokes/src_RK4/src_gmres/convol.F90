!
      SUBROUTINE CONVOL(ZK,NZK,FF2,UK,VK,UR,VR,ZR,NZR)
!
!     CALCULATES CONVOLUTION SUMS, CALLS FFT'S, ETC.
!
      USE GLOBAL_PARAMS
      IMPLICIT NONE 
!
      INTEGER                           :: IKX,IKY,I,J
!
      REAL(rk)                          :: KX,KY,WK,WEIGHT
      REAL(rk), DIMENSION(NX,NY)        :: UR,VR,NZR,ZR
!
      COMPLEX(rk), DIMENSION(IKTX,IKTY) :: NZK,ZK,UK,VK
      COMPLEX(rk)                       :: C1,FF2(NX/2+1,NY)
!
      DO IKY = 1, IKTY
         KY =  REAL(IKY - KTY - 1,rk)
         DO IKX = 1, IKTX
            KX =  REAL(IKX-1,rk)*alpha 
            WK = MAX( KX*KX+KY*KY ,0.001_rk ) 
            NZK(IKX,IKY) = CMPLX(0._rk,0._rk,rk) 
            UK(IKX,IKY)  = + L(IKX,IKY)*ZI*KY*ZK(IKX,IKY)/WK
            VK(IKX,IKY)  = - L(IKX,IKY)*ZI*KX*ZK(IKX,IKY)/WK
         ENDDO
      ENDDO
!
      CALL KR_FFTW(ZK,ZR,FF2,PlanK2ZR)
      CALL KR_FFTW(UK,UR,FF2,PlanK2UR)
      CALL KR_FFTW(VK,VR,FF2,PlanK2VR)
!
      DO J=1,NY
         DO I=1,NX
            NZR(I,J) = UR(I,J) * ZR(I,J)
         ENDDO
      ENDDO
!
      CALL RK_FFTW(NZR,NZK,FF2,PlanNZR2K)
!
      DO IKX = 1, IKTX
         KX =  REAL(IKX-1,rk)*alpha  
         DO IKY = 1, IKTY
            NZK(IKX,IKY) = + ZI*KX * NZK(IKX,IKY) 
         ENDDO
      ENDDO
!
      DO J=1,NY   
         DO I=1,NX  
            UR(I,J) = VR(I,J) * ZR(I,J)
         ENDDO
      ENDDO
!
      CALL RK_FFTW(UR,UK,FF2,PlanUR2K)
!
      DO IKY = 1, IKTY
         DO IKX = 1, IKTX
            KY = IKY - KTY - 1
            WEIGHT       =    L(IKX,IKY)
            C1           =  + ZI*KY * UK(IKX,IKY)
            NZK(IKX,IKY) = (- NZK(IKX,IKY) - C1)*WEIGHT
         ENDDO
      ENDDO
!
      END SUBROUTINE CONVOL

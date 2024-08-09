!
      SUBROUTINE TIMEDERIV(ZO,RHZ,                                      & 
     &                     KFX2,KFY2,NZK,FF2,UK,VK,UR,VR,ZR,NZR)
!
!     CALCULATES CONVOLUTION SUMS, CALLS FFT'S, ETC.
!
      USE GLOBAL_PARAMS
      IMPLICIT NONE 
!
      INTEGER  :: KFX2,KFY2,IKX,IKY
!
      REAL(rk), DIMENSION(NX,NY) :: UR,VR,NZR,ZR
!
      COMPLEX(rk), DIMENSION(IKTX,IKTY) :: ZO,RHZ,NZK,UK,VK
      COMPLEX(rk), DIMENSION(NX/2+1,NY) :: FF2
!
!
!---------------------------------------------------------
      CALL CONVOL(ZO,NZK,FF2,UK,VK,UR,VR,ZR,NZR)
!---------------------------------------------------------
!
!        NZK(KFX1,KFY1) = NZK(KFX1,KFY1) + AMPFOR*CMPLX(1._rk,0._rk,rk)
        NZK(KFX2,KFY2) = NZK(KFX2,KFY2) -2._rk*CMPLX(1._rk,0._rk,rk)
!        NZK(KFX2,KTY+2) = NZK(KFX2,KTY+2) -0.5_rk*AMPFOR*CMPLX(1._rk,0._rk,rk)
!
        WHERE (LL) 
           RHZ =  - NU*ZO  + NZK 
        ELSEWHERE
           RHZ = CMPLX(0._rk,0._rk,rk)
        END WHERE     
!
      END SUBROUTINE TIMEDERIV

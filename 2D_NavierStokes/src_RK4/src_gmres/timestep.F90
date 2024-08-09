!
      SUBROUTINE TIMESTEP(ZO,ZN,TSTART,RHZ,                             &
     &                    KFX2,KFY2,NZK,FF2,UK,VK,UR,VR,ZR,NZR)
!
!     CALCULATES CONVOLUTION SUMS, CALLS FFT'S, ETC.
!     This scheme is a 2nd-order Predictor Corrector 
!     as used in A.P. Willis' Pipe code 
!     This is the implicit Crank-Nicholson scheme for the viscous terms and 
!     a 2nd-order explicit (Heun's method) for the nonlinear and forcing terms.
!
      USE GLOBAL_PARAMS
      IMPLICIT NONE 
!
      INTEGER  :: KFX2,KFY2
      INTEGER  :: NT,IKX,IKY,j
!
      REAL(rk), DIMENSION(NX,NY)     :: UR,VR,NZR,ZR
      REAL(rk)                       :: TSTART,TIME
      REAL(rk), DIMENSION(IKTX,IKTY) :: NUZN,NU1,NU2
!
      COMPLEX(rk), DIMENSION(IKTX,IKTY)   :: ZO,ZN,RHZ,NZK,UK,VK,FN
      COMPLEX(rk), DIMENSION(NX/2+1,NY)   :: FF2
!
      CHARACTER(LEN=5)  :: xstring
      CHARACTER(LEN=50) :: filename
!
      WHERE (LL)
         FN = CMPLX(0._rk,0._rk,rk)
         ZN = ZO   
         NUZN = (1._rk - 0.5_rk*NU*DELT) / DELT
         NU1 = DELT / (1._rk + 0.5_rk*NU*DELT)
         NU2 = NU1/2._rk
      END WHERE
!
      DO NT=1,NSTOP
!
!---------------------------------------------------------
         CALL CONVOL(ZN,NZK,FF2,UK,VK,UR,VR,ZR,NZR)
!---------------------------------------------------------
!
!         NZK(KFX1,KFY1) = NZK(KFX1,KFY1) + AMPFOR*CMPLX(1._rk,0._rk,rk)
         NZK(KFX2,KFY2) = NZK(KFX2,KFY2) + AMPFOR*CMPLX(1._rk,0._rk,rk)
!
         WHERE (LL) 
            RHZ = NZK
            ZO = (NUZN*ZN + NZK)*NU1  
         END WHERE 
!
!---------------------------------------------------------
         CALL CONVOL(ZO,NZK,FF2,UK,VK,UR,VR,ZR,NZR)
!---------------------------------------------------------
!
!         NZK(KFX1,KFY1) = NZK(KFX1,KFY1) + AMPFOR*CMPLX(1._rk,0._rk,rk)
         NZK(KFX2,KFY2) = NZK(KFX2,KFY2) + AMPFOR*CMPLX(1._rk,0._rk,rk)
!
         WHERE (LL) 
            ZN = ZO + (NZK-RHZ)*NU2  
         ELSEWHERE
            ZN = CMPLX(0._rk,0._rk,rk)
         END WHERE  
!
         IF (SubSpaceFLAG.EQ.1) THEN
            ! Set Imaginary parts to zero
            WHERE (LL) ZN = CMPLX(REAL(ZN),0._rk,rk)  
         ENDIF  
!
         TIME = TSTART + NT*DELT                
!
      ENDDO
!
      END SUBROUTINE TIMESTEP

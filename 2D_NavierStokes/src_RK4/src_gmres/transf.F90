!
      SUBROUTINE TRANSF (Z,NZ,TIME)
!
!     CALCULATES SPECTRA.
!
      USE GLOBAL_PARAMS
      IMPLICIT NONE 
!
      INTEGER                            :: IKX,IKY,J
      REAL(rk)                           :: KX,KY,WK,VZ,AVZ,AVE,TIME
      COMPLEX(rk),  DIMENSION(IKTX,IKTY) :: Z,NZ
!

!
      PRINT *,'   '
      PRINT *,'ENERGY & ENSTROPHY TRANSFER SPECTRA AT T = ',TIME
      AVZ = 0._rk
      AVE = 0._rk
!
      DO J=1,KTX 
         SPE(J)   = 0._8
         SPZ(J)   = 0._8
         NS(J)     = 0
      ENDDO
!
      DO IKX=1,IKTX
         KX = REAL(IKX-1,rk)*alpha
         DO IKY=1,IKTY
            KY = REAL(IKY-KTY-1,rk)
            WK = SQRT(KX*KX+KY*KY)
            IF(.NOT.LL(IKX,IKY))     CYCLE
            J = INT(WK+0.5_rk)
            IF(J.LE.0 .OR. J.GT.KTX) PRINT *,'SCREW-UP.'
            VZ = REAL( Z(IKX,IKY)*CONJG(NZ(IKX,IKY)) )
            SPZ(J) = SPZ(J) + 2._rk*VZ
            SPE(J) = SPE(J) + 2._rk*VZ/WK**2
            NS(J)   = NS(J) + 2 
         ENDDO
      ENDDO
!
      WRITE(23,*) '              '
      WRITE(23,*) TIME,'  = TIME'
!
      DO J=1,KTX-1
         AVZ = AVZ + SPZ(J)
         AVE = AVE + SPE(J)
! aqui       WRITE( 6,5000) FLOAT(J),SPZ(J),SPE(J),NS(J)
! aqui       WRITE(23,5000) FLOAT(J),SPZ(J),SPE(J),NS(J)
      ENDDO
!
! aqui      WRITE(6,5010) AVZ,AVE
!
! aqui 5000  FORMAT(1X,F4.0,4X,E15.8,4X,E15.8,10X,I4)
! aqui 5010  FORMAT(1X, 'TOTALS   ', E14.8 , '     ', E14.8)
!
      END SUBROUTINE TRANSF


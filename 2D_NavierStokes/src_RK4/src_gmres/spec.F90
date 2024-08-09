!
      SUBROUTINE SPEC (Z,TIME)
!
!     CALCULATES SPECTRA.
!
      USE GLOBAL_PARAMS
      IMPLICIT NONE 
!
      INTEGER      :: IKX,IKY,J
      COMPLEX(rk)  :: Z(IKTX,IKTY)
      REAL(rk)     :: KX,KY,WK,VZ(IKTX,IKTY),AVZ,AVE,TIME
!
!      PRINT *,'ENSTROPHY AND ENERGY SPECTRA AT T = ',TIME
      AVZ = 0.0_rk
      AVE = 0.0_rk
!
      DO J=1,KTX 
        SPZ(J)   = 0._8
        SPE(J)   = 0._8
        NS(J)     = 0
      ENDDO
!
      DO IKX=1,IKTX
         KX = REAL(IKX-1,rk)*alpha
         DO IKY=1,IKTY
            KY = REAL(IKY-KTY-1,rk)
            WK = SQRT(KX*KX+KY*KY)
            VZ(IKX,IKY)=0._rk
            IF(.NOT.LL(IKX,IKY))     CYCLE
            J = INT(WK+0.5_rk)
            IF(J.LE.0 .OR. J.GT.KTY) PRINT *,'SCREW-UP.',J
            VZ(IKX,IKY) = REAL( Z(IKX,IKY)*CONJG(Z(IKX,IKY)) )
            SPZ(J) = SPZ(J) + 2._rk*VZ(IKX,IKY)
            VZ(IKX,IKY) = VZ(IKX,IKY)/(WK*WK)
            SPE(J) = SPE(J) + VZ(IKX,IKY)
            NS(J)   = NS(J) + 2 
         ENDDO
      ENDDO
!
      WRITE(13,*) '              '
      WRITE(13,*) TIME,'  = TIME'
!
      DO J=1,KTX-1
         AVZ = AVZ + SPZ(J)
         AVE = AVE + SPE(J)
! aqui       WRITE( 6,5000) FLOAT(J),SPZ(J),SPE(J),NS(J)
         WRITE(13,5000) REAL(J,rk),SPZ(J),SPE(J),NS(J)
      ENDDO
!
      write(31) TIME, ((VZ(ikx,iky),ikx=1,IKTX),iky=1,IKTY)
!
! aqui      WRITE(6,5010) AVZ,AVE
!
5000  FORMAT(1X,F5.0,4X,E15.8,4X,E15.8,10X,I4)
! aqui 5010  FORMAT(1X, 'TOTALS   ', E14.8 , '     ', E14.8)
!
      END SUBROUTINE SPEC 
!
!


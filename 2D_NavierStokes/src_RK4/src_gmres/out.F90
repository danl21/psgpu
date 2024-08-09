!
      SUBROUTINE OUT (Z,TIME)
!
!     CONTROLS OUTPUT
!
      use GLOBAL_PARAMS
      IMPLICIT NONE 
      INTEGER     :: IKX,IKY
!
      COMPLEX(rk) :: Z(IKTX,IKTY),UKKF
!
      REAL(rk)    :: VZ,I,D
      REAL(rk)    :: KX,KY,TIME,WK
!
      REAL(KIND=8) :: E,V,P
!
      E = 0._8
      V = 0._8
      P = 0._8
      DO IKX=1,IKTX
         KX = REAL(IKX-1,rk)*alpha
         DO IKY=1,IKTY
            KY = REAL(IKY-KTY-1,rk)
            WK = KX*KX + KY*KY
            IF(.NOT.LL(IKX,IKY)) CYCLE
            VZ = REAL( Z(IKX,IKY)*CONJG(Z(IKX,IKY)) )
            E = E + VZ/WK  ! Energy
            V = V + VZ
!            P = P + VZ*WK
         ENDDO
      ENDDO
      V = 2._rk*V  ! Enstrophy
!      P = 2._rk*P  ! Palinstrophy
      D = V*v2  ! Dissipation
      IKY = KF+KTY+1  ! index corresponding to wavenumber KF (KY=IKY-KTY-1)
      UKKF  = + L(1,IKY)*ZI*KF*Z(1,IKY)/(KF*KF) !UK(0,KF)
      I = - AIMAG(UKKF)
!
!
!      PRINT*,'    ' 
!      PRINT*,'    ' 
!      PRINT*,'_______________________________________________________' 
!      PRINT*,'    ' 
! aqui      WRITE( 6,5043) 
!      WRITE( 6,5044) Time,E,D,I,nt
!      PRINT *,' '
!
      WRITE(15,5045) Time,E,D,I
!
!
! aqui 5043  FORMAT(16X,'Time',8X,'E',9X,'Z',21X,'P',16X,'nt')
! 5044  FORMAT(1X,4(E12.5,1X),I8)
! 5045  FORMAT(1X,E12.5,4X,E12.5,4X,E12.5,4X,E12.5,4x,i10)
 5045  FORMAT(1X,E12.5,4X,E12.5,4X,E12.5,4X,E12.5)      
!
      END SUBROUTINE OUT 
!
!
!
!








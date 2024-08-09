      SUBROUTINE MAKEK (KX,KY,KZ)
!
!     Recomputes wavevector for when alpha is updated

      use GLOBAL_PARAMS
      IMPLICIT NONE 
!
      INTEGER             :: IKX,IKY,IKZ,ntime
      INTEGER             :: IK,jw,kx1,ky1,kz1,klag,jlag
      REAL(rk) :: KKX,KKZ,KKY
!
      REAL(rk), DIMENSION(NKT) :: KX,KY,KZ

!
      EXTERNAL RANNO
!
      DO IKX=1,IKTX
         KKX = REAL(IKX - 1,rk)*alpha
         DO IKY=1,IKTY
            KKY = REAL(IKY - KTY - 1,rk)
            DO IKZ=1,IKTZ
               KKZ = REAL(IKZ - KTZ - 1,rk)

               IK = ((IKZ-1)*IKTY +IKY-1)*IKTX + IKX
               
               KX(IK) = KKX
               KY(IK) = KKY
               KZ(IK) = KKZ
            ENDDO
         ENDDO
      ENDDO
    END SUBROUTINE MAKEK

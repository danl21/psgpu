      module SYMMETRIC
!
      use GLOBAL_PARAMS
      implicit none
!
      contains
!
!--------------------------------------------------------------
!
      subroutine attemptSym(ZO,SHIFT)
!
      COMPLEX(rk), DIMENSION(IKTX,IKTY)   :: ZO,ZNSHIFT,ZTEMP
      REAL(rk)  :: normIM,normIMmin,SHIFT,SHIFTtemp,SHIFTdir
      REAL(rk)  :: xs,xstemp 
      INTEGER :: i,imin
!
!      
      normIMmin = 9999999._rk
      DO i=1,4
         IF (i.eq.1) THEN
            Print*,'Attempt 1 - just x shift'
            WHERE(LL) ZNSHIFT = ZO
            SHIFTtemp = 1._rk
         ENDIF
         IF (i.eq.2) THEN
            Print*,'Attempt 2 - y shift = pi/2 followed by x shift'
            call shiftY(ZO,ZNSHIFT,twopi/4._rk)
            SHIFTtemp = 1._rk
         ENDIF
         IF (i.eq.3) THEN
            Print*,'Attempt 3 - y shift = pi/4',                   &
     &           ' and rotate followed by x shift'
            call shiftY(ZO,ZNSHIFT,twopi/8._rk)
            call rotateY(ZNSHIFT)
            SHIFTtemp = -1._rk
         ENDIF
         IF (i.eq.4) THEN
            Print*,'Attempt 4 - y shift = -pi/4',                  &
     &           ' and rotate followed by x shift'
            call shiftY(ZO,ZNSHIFT,-twopi/8._rk)
            call rotateY(ZNSHIFT)
            SHIFTtemp = -1._rk
         ENDIF
 !
         call findXshift(ZNSHIFT,xstemp)
         normIM = SQRT( SUM( AIMAG(ZNSHIFT)**2 , LL ) )
         IF (normIM.LT.normIMmin) THEN
            normIMmin = normIM
            imin = i
            SHIFTdir = SHIFTtemp
            ZTEMP = ZNSHIFT
            xs = xstemp
         ENDIF
         IF (normIM.LT.0.00001_rk) THEN
            Print*,'Symmetric solution found with attempt ',i
            Print*,'normIM = ',normIM,', x-shift = ',xs
            EXIT
         ENDIF
         Print*,'Failed attempt ',i
         IF (i.eq.4) THEN
            Print*,'NO SYMMETRIC SOLUTION FOUND'
            Print*,'min normIM found with attempt ',imin 
            Print*,'normIM = ',normIMmin,', x-shift = ',xs
            Print*,'Setting SubSpaceFLAG to 0' 
            SubSpaceFLAG = 0
            ZNSHIFT = ZTEMP
         ENDIF
      ENDDO
!
!     Shift so selected symmtric point is in centre      
      Print*,'vertex = ',vertex
      IF (vertex.eq.0) THEN
         Print*,'reverting back to original state'
      ENDIF
      IF (vertex.eq.1) THEN
         Print*,'new shifted state accepted'
         WHERE(LL) ZO = ZNSHIFT
         SHIFT = SHIFTdir*SHIFT
      ENDIF
      IF (vertex.eq.2) THEN
         Print*,'Shift in X by pi'
         xs = xs + twopi/2._rk
         call shiftX(ZNSHIFT,ZO,twopi/2._rk)
         SHIFT = SHIFTdir*SHIFT
         normIM = SQRT( SUM( AIMAG(ZNSHIFT)**2 , LL ) )
      ENDIF
      IF (vertex.eq.3) THEN
         Print*,'Shift in Y by pi'
         call shiftY(ZNSHIFT,ZO,twopi/2._rk)
         SHIFT = SHIFTdir*SHIFT
         normIM = SQRT( SUM( AIMAG(ZNSHIFT)**2 , LL ) )
      ENDIF
      IF (vertex.eq.4) THEN
         Print*,'Shift in X and Y by pi'
         xs = xs + twopi/2._rk
         call shiftY(ZNSHIFT,ZO,twopi/2._rk)
         WHERE(LL) ZNSHIFT=ZO
         call shiftX(ZNSHIFT,ZO,twopi/2._rk)
         SHIFT = SHIFTdir*SHIFT
         normIM = SQRT( SUM( AIMAG(ZNSHIFT)**2 , LL ) )
      ENDIF
!
      end subroutine attemptSym
!-----------------------------------------------
!
      subroutine findXshift(Z,xsfinal)
!
      COMPLEX(rk), DIMENSION(IKTX,IKTY)   :: Z,ZS
      REAL(rk)  :: gr,xs(1:4),fun(1:4),xsfinal  
      INTEGER :: is
!     
!     This is the golden section search algorithm
      gr = (1._rk + SQRT(5._rk)) / 2._rk
      xs(1) = 0._rk
      xs(3) = twopi / 2._rk
!
      xs(2) = (xs(3) + xs(1)*gr)/(1._rk + gr)
      xs(4) = (xs(3) + xs(2)*gr)/(1._rk + gr)
!
      DO is = 1,4
         call shiftX(Z,ZS,xs(is))
         fun(is) = SQRT( SUM( AIMAG(ZS)**2 , LL ) )
      ENDDO
!
      DO is = 1,75
         IF (fun(2).GE.fun(4)) THEN
            xs(1) = xs(2)
            fun(1) = fun(2) 
            xs(2) = xs(4) 
            fun(2) = fun(4) 
            xs(4) = (xs(3) + xs(2)*gr)/(1._rk + gr)
            call shiftX(Z,ZS,xs(4))
            fun(4) = SQRT( SUM( AIMAG(ZS)**2 , LL ) )
         ELSE
            xs(3) = xs(4)
            fun(3) = fun(4) 
            xs(4) = xs(2) 
            fun(4) = fun(2) 
            xs(2) = (xs(3) + xs(1)*gr)/(1._rk + gr)
            call shiftX(Z,ZS,xs(2))
            fun(2) = SQRT( SUM( AIMAG(ZS)**2 , LL ) )
         ENDIF
      ENDDO
!     
      IF (fun(2).GE.fun(4)) THEN
         xsfinal = xs(4)
         call shiftX(Z,ZS,xs(4))
         fun(4) = SQRT( SUM( AIMAG(ZS)**2 , LL ) )
      ELSE
         xsfinal = xs(2)
         call shiftX(Z,ZS,xs(2))
         fun(2) = SQRT( SUM( AIMAG(ZS)**2 , LL ) )
      ENDIF
      WHERE(LL) Z=ZS
!
      end subroutine findXshift
!
! -------------------------------------------------------------
!
      subroutine shiftX(Z,ZS,SHIFT)
!
      COMPLEX(rk), DIMENSION(IKTX,IKTY)   :: Z,ZS
      REAL(rk)  :: SHIFT,KX
      INTEGER :: IKX,IKY
!
      DO IKY = 1, IKTY
         DO IKX = 1, IKTX
            IF(.NOT.LL(IKX,IKY))     CYCLE
            KX =  REAL(IKX-1,rk)*alpha
            ZS(IKX,IKY) = EXP(-ZI*KX*SHIFT)*Z(IKX,IKY)
         ENDDO
      ENDDO
!
      end subroutine shiftX
!
!-----------------------------------------------------------------
!
      subroutine shiftY(Z,ZS,SHIFT)
!
      COMPLEX(rk), DIMENSION(IKTX,IKTY)   :: Z,ZS
      REAL(rk)  :: SHIFT,KY
      INTEGER :: IKX,IKY
!
      DO IKY = 1, IKTY
         KY = REAL(IKY-KTY-1,rk)
         DO IKX = 1, IKTX
            IF(.NOT.LL(IKX,IKY))     CYCLE
            ZS(IKX,IKY) = EXP(-ZI*KY*SHIFT)*Z(IKX,IKY)
         ENDDO
      ENDDO
!
      end subroutine shiftY
!
!-----------------------------------------------------------------
!
      subroutine rotateY(Z)
!
      COMPLEX(rk), DIMENSION(IKTX,IKTY)   :: Z,ZS
      INTEGER :: IKX,IKY
!
      DO IKY = 1, IKTY
         IF(LL(1,IKY)) ZS(1,IKY) = Z(1,IKY)
         DO IKX = 2, IKTX !Note this is from KX=1 NOT KX=0
            IF(.NOT.LL(IKX,IKY))     CYCLE
            ZS(IKX,IKY) = CONJG(Z(IKX,IKTY+1-IKY))
         ENDDO
      ENDDO
!
      WHERE(LL) Z = -ZS
!
      end subroutine rotateY
!
!-----------------------------------------------------------------
      end module SYMMETRIC

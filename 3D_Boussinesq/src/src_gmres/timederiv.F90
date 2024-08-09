!
      SUBROUTINE TIMEDERIV(ZO,RHZ,KX,KY,KZ,IKF,IKN)
!
!     Alternative way to construct dz/dt by finite difference
!
      USE GLOBAL_PARAMS
      IMPLICIT NONE 
!
      INTEGER  :: IKT
      INTEGER  :: IKF(NX2*NY*NZ),IKN(NKT)
!
      REAL(rk), DIMENSION(NKT) :: KX,KY,KZ
      REAL(rk) :: FK,TSMALL,DDT,TSTART,TIME,norm,Dtarget
!
      COMPLEX(rk), DIMENSION(4*NKT) :: ZO,RHZ
      COMPLEX(rk), DIMENSION(NKT) :: XII,ETA,ZET,RHO
!
!
!---------------------------------------------------------
      DO IKT =1,NKT
         IF (LL(IKT))THEN
            XII(IKT) = ZO(IKT)
            ETA(IKT) = ZO(NKT+IKT)
            ZET(IKT) = ZO(2*NKT+IKT)
            RHO(IKT) = ZO(3*NKT+IKT)/scale
         ENDIF
      ENDDO
      !                   
      Dtarget=0._rk
      TIME=0._rk
      TSTART=0._rk
      DDT=DELT*0.0001_rk
      TSMALL=DDT
      FK = REAL(KF,rk)
      CALL TIMESTEP_CUDA(XII,ETA,ZET,RHO,KX,KY,KZ,TIME,TSTART,TSMALL,AMPFOR,DDT,ResidualThreshold,&
     &                   FK,v2,Ri,Sc,Theta,alpha,Dtarget,IKF,IKN,L,IKTX,IKTY,IKTZ,KTZ,NKT,NX,NY,NZ,&
     &                   NOUT,NOUTV,STATSFLAG,RCFLAG,ADAPTFLAG,RANK)
     !  print*, 'check ||ZT|| ', SQRT( SUM( REAL(XII*CONJG(XII))) &
     ! &      + SUM( REAL(ETA*CONJG(ETA))) &
     ! &      + SUM( REAL(ZET*CONJG(ZET))) &
     ! &      + SUM( REAL(RHO*CONJG(RHO)))  )

      DO IKT = 1,NKT
         IF (LL(IKT))THEN 
            RHZ(IKT) =  (XII(IKT)-ZO(IKT) )/DDT
            RHZ(NKT+IKT) =  (ETA(IKT)-ZO(IKT+NKT))/DDT
            RHZ(2*NKT+IKT) =  (ZET(IKT)-ZO(IKT+2*NKT))/DDT
            RHZ(3*NKT+IKT) =  scale*(RHO(IKT)-ZO(IKT+3*NKT))/DDT
         ELSE
            RHZ(IKT) = CMPLX(0._rk,0._rk,rk)
            RHZ(NKT+IKT) = CMPLX(0._rk,0._rk,rk)
            RHZ(2*NKT+IKT) = CMPLX(0._rk,0._rk,rk)
            RHZ(3*NKT+IKT) = CMPLX(0._rk,0._rk,rk)
         ENDIF
      ENDDO
      print*, 'check ||RHZ|| ', SQRT( SUM( REAL(RHZ*CONJG(RHZ))) )
      FLUSH(6)
!
      END SUBROUTINE TIMEDERIV

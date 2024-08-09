
      program write_spec
C
C write 1D spectra from full 3D spectral restart data
C
      use GLOBAL_PARAMS
      IMPLICIT NONE
      INTEGER :: ikx,iky,J,ikz,nt,nw,ifile,iopt,kopt,ntime,IK
      REAL(rk)     :: KKX,KKY,KKZ,WK,WKH,WKV,maxVort
      INTEGER :: jw,kx1,ky1,kz1
C ----------------------------------------------------------------
      COMPLEX(rk), DIMENSION(NKT) :: XII,ETA,ZET,RHO
      REAL(rk):: SP(KTX),SH(KTX),SV(KTX)
      REAL(rk)    :: timein,PERIODin,SHIFTXin,SHIFTYin,SHIFTZin,NERin
      REAL(rk):: time,dummy
c
c
      open(12,file='Zk.in',form='UNFORMATTED')!,access='stream')
      open(15,file='spectrum.dat')
c
      DO J=1,KTX 
        SP(J)   = 0._rk
        SH(J)   = 0._rk
        SV(J)   = 0._rk
      ENDDO
c
      nw=1
c      read(12) jw,kx1,ky1,kz1
      do ntime = 1, nw
         read(12) timein,v2,PERIODin,SHIFTXin,SHIFTYin,SHIFTZin,NERin,
     .        (XII(ikx),ikx=1,nkt),
     .        (ETA(ikx),ikx=1,nkt),
     .        (ZET(ikx),ikx=1,nkt),
     .        (RHO(ikx),ikx=1,nkt)
      enddo
C
      PRINT*, timein,v2,PERIODin,SHIFTXin,SHIFTYin,SHIFTZin,NERin
c
      maxVort = MAXVAL(REAL(XII(:)*CONJG(XII(:))) + 
     .              REAL(ETA(:)*CONJG(ETA(:))) + 
     .              REAL(ZET(:)*CONJG(ZET(:))))

      DO IKX=1,IKTX
         KKX = REAL(IKX - 1,rk)*alpha
         DO IKY=1,IKTY
            KKY = REAL(IKY - KTY - 1,rk)
            DO IKZ=1,IKTZ
               KKZ = REAL(IKZ - KTZ - 1,rk)

               IK = ((IKZ-1)*IKTY +IKY-1)*IKTX + IKX
               WK=SQRT(KKX**2 + KKY**2 + KKZ**2)
               WKH=SQRT(KKX**2 + KKZ**2)
               WKV=SQRT(KKY**2)
c               IF(.NOT.LL(IK))     CYCLE
               J= INT(WK+0.5_rk)
               dummy = REAL(XII(IK)*CONJG(XII(IK))) + 
     .              REAL(ETA(IK)*CONJG(ETA(IK))) + 
     .              REAL(ZET(IK)*CONJG(ZET(IK)))

               if(dummy .gt. 0.01*maxVort)then
                  write(77,*) KKX,KKY,KKZ,dummy
               endif
               SP(J) = SP(J) + dummy/WK**2._rk
               J= INT(WKV+0.5_rk)
               SV(J) = SV(J) + dummy/WK**2._rk
               J= INT(WKH+0.5_rk)
               SH(J) = SH(J) + dummy/WK**2._rk
            ENDDO
         ENDDO
      ENDDO
      DO J=1,KTX-1
         WRITE(15,5000) REAL(J,rk),SP(J),SH(J),SV(J)
      ENDDO

 5000 FORMAT(1X,F5.0,3(1X,E15.8))
      END
C
c

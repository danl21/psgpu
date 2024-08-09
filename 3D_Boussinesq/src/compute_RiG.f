      program compute_RiG
C
C     Computation of 3D gradient Richardson number from snapshots
C
      use GLOBAL_PARAMS
      IMPLICIT NONE
      INTEGER :: ikx,iky,ikz,nt,nw,nlx,ip,im
C ----------------------------------------------------------------
      REAL(rk) :: UR(NX,NY,NZ),VR(NX,NY,NZ),TR(NX,NY,NZ),dy,Navg,Savg
      REAL(rk) :: RIG(NX,NY,NZ),B
      character(len=256) :: infile1,infile2,infile3,outfile
c
      TWOPI = 4._rk*ASIN(1._rk)
      print*,'Which file?'
      read(5,*) nw
      print*,'What B?'
      read(5,*) B

      write(infile1,"(A1 I3.3 A4)") "U",nw,".raw"
      write(infile2,"(A1 I3.3 A4)") "V",nw,".raw"
      write(infile3,"(A3 I3.3 A4)") "rho",nw,".raw"
      write(outfile,"(A3 I3.3 A4)") "RiG",nw,".raw"

      print*, trim(infile1), trim(infile2),trim(infile3),trim(outfile)

c
      open(10,file=trim(infile1),access="STREAM",form='UNFORMATTED')
      open(11,file=trim(infile2),access="STREAM",form='UNFORMATTED')
      open(12,file=trim(infile3),access="STREAM",form='UNFORMATTED')
      open(13,file=trim(outfile),access="STREAM",form='UNFORMATTED')

      read(10) (((UR(ikx,iky,ikz),ikx=1,NX),iky=1,NY),ikz=1,NZ)
      read(11) (((VR(ikx,iky,ikz),ikx=1,NX),iky=1,NY),ikz=1,NZ)
      read(12) (((TR(ikx,iky,ikz),ikx=1,NX),iky=1,NY),ikz=1,NZ)

      dy = 2._rk*TWOPI/REAL(NZ,rk)
      do ikz=1,NZ
         Navg=0._rk
         Savg=0._rk
         ip = ikz+1
         if(ikz .eq. NZ) ip=1
         im = ikz-1
         if(ikz .eq. 1) im=NZ
         do ikx=1,NX
            do iky=1,NY
               Navg = (TR(ikx,iky,ip)-TR(ikx,iky,im))/dy
               Savg = 0.5_rk*(((UR(ikx,iky,ip)-UR(ikx,iky,im))/dy)**2+
     .              ((VR(ikx,iky,ip)-VR(ikx,iky,im))/dy)**2)
               RIG(ikx,iky,ikz)=B*(Navg+1.0)/Savg
            enddo
         enddo
      enddo
      write(13) (((RIG(ikx,iky,ikz),ikx=1,NX),iky=1,NY),ikz=1,NZ)
      
C
      END
C

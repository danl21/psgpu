      program read
C
C DIRECT SIMULATION WITH CENTRED DIFFERENCES AND CRANK-NICHOLSON DAMPING.
C
      use GLOBAL_PARAMS
      IMPLICIT NONE
      INTEGER :: ikx,iky,ikz,nt,nw,nlx,iopt,i,ikp,ikm
C ----------------------------------------------------------------
      REAL(rk) :: ZR(NX,NY,NZ),dy,TIME,du,u,duz,dz,maxSy, maxSz
c
      TWOPI = 4._rk*ASIN(1._rk)
      print*,'Which average? (1) rho, (2) U, (3) V,'
      print*, '(4) W, (5) Vrms, (6) Wrms, (7) Urms?'
      read(5,*) iopt
c
      open(12,file='avgs.dat',access="STREAM",form='UNFORMATTED')
      open(13,file='avg.raw',access="STREAM",form='UNFORMATTED')
      open(14,file='avg_xy.dat',form='FORMATTED')
      open(15,file='avg_xz.dat',form='FORMATTED')
      open(16,file='avg_yz.dat',form='FORMATTED')
      do i=1,iopt
         read(12) TIME
         read(12) (((ZR(ikx,iky,ikz),ikx=1,NX),iky=1,NY),ikz=1,NZ)
      enddo

      write(13) (((ZR(ikx,iky,ikz),ikx=1,NX),iky=1,NY),ikz=1,NZ)
      do iky=1,NY
         write(14,333) (ZR(ikx,iky,NZ/2),ikx=1,NX)
         write(14,*) '           '
      enddo
      do ikz=1,NZ
         write(15,333) (ZR(ikx,NY/2,ikz),ikx=1,NX)
         write(15,*) '           '
      enddo
      do ikz=1,NZ
         write(16,333) (ZR(NX/2,iky,ikz),iky=1,NY)
         write(16,*) '           '
      enddo
      dy=TWOPI/NY
      dz=TWOPI/NZ
      maxSy=0._rk
      maxSz=0._rk
      do iky=2,NY-1
         du=0._rk
         duz=0._rk
         u=0._rk
         do ikz=1,NZ
            ikp = ikz+1
            if(ikz.eq.NZ) ikp=1
            ikm = ikz-1
            if(ikm.eq.0) ikm=NZ
            do ikx=1,NX
               u=u+ZR(ikx,iky,ikz)
               du = du+
     .              (ZR(ikx,iky+1,ikz)-ZR(ikx,iky-1,ikz))/(2._rk*dy)
               duz = duz+
     .              (ZR(ikx,iky,ikp)-ZR(ikx,iky,ikm))/(2._rk*dz)
               maxSy = 
     .       max(maxSy,(ZR(ikx,iky+1,ikz)-ZR(ikx,iky-1,ikz))/(2._rk*dy))
         maxSz = max(maxSz,(ZR(ikx,iky,ikp)-ZR(ikx,iky,ikm))/(2._rk*dz))
            enddo
         enddo
         write(33,*)iky*dy,u/REAL(NX*NZ),du/REAL(NX*NZ),duz/REAL(NX*NZ)
      enddo

      print*, 'max shear: ',maxSy,maxSz

      print*, 'mean: ',SUM(ABS(ZR))/REAL(NZ*NY*NX)
      close(12)
      close(13)
      close(14)
      close(15)
999   stop
 333  format(1x,E12.5,1x)
 3333  format(2x,E12.5,1x)
      END
C
c
c
      subroutine big(zr,N,M,bg)
c
c  Just for scaling the data to make graphs prettier...
c
      implicit none
      integer i,j,N,i2,M
      real    zr(N,M),zmax,zmin,bg
c
        bg=4.
c
        zmax=+1.01*sqrt(bg)
        zmin=-zmax
c
      do j=1,M
      do i=1,N
        if (zr(i,j).ge.0.) zr(i,j)=sqrt(zr(i,j))
        if (zr(i,j).lt.0.) zr(i,j)=-sqrt(-zr(i,j))
        if (zr(i,j).lt.zmin) zr(i,j)=zmin
        if (zr(i,j).gt.zmax) zr(i,j)=zmax
      enddo
      enddo
c
      zr(2,2)=zmax
      zr(1,1)=zmin
c
      return
      end
c

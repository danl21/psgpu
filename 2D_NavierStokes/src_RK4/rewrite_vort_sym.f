      program read
C
C DIRECT SIMULATION WITH CENTRED DIFFERENCES AND CRANK-NICHOLSON DAMPING.
C
      use GLOBAL_PARAMS
      IMPLICIT NONE
      INTEGER :: ikx,iky,nt,nw
C ----------------------------------------------------------------
      REAL(rk) :: ZR(NX,NY),Z2(NX,NY),bg,time
      REAL(rk) :: s1,s2,s3,s4,s5,s6,d1,d2,d3
c
      print*,'Which output?'
      read(5,*) nw
      print*,'You have selected number ',nw
c
      open(12,file='vort.dat',access="STREAM",form='UNFORMATTED')
      open(13,file='vort_a.dat',form='FORMATTED')

      read(12) time,((ZR(ikx,iky),ikx=1,NX),iky=1,NY)
      read(12) time,((Z2(ikx,iky),ikx=1,NX),iky=1,NY)
      do iky=1,NY
         write(13,333) (abs(ZR(ikx,iky)-Z2(ikx,iky)),ikx=1,NX)
         write(13,*) '           '
      enddo
999   stop
 333  format(1x,E12.5,1x)
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

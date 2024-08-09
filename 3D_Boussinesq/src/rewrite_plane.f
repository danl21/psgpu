      program read
C
C     Program to rewrite mid-planes of 3D binary data
C
      use GLOBAL_PARAMS
      IMPLICIT NONE
      INTEGER :: ikx,iky,ikz,nt,nw,nlx,iopt
C ----------------------------------------------------------------
      REAL(rk) :: ZR(NX,NY,NZ),dy
c
      TWOPI = 4._rk*ASIN(1._rk)
      print*,'Which file?'
      read(5,*) infile
c     If this is the buoyancy field, add the background gradient
      print*,'Is this rho? (1) yes, (2) no;?'
      read(5,*) iopt
c
      open(12,file=trim(infile),access="STREAM",form='UNFORMATTED')
      open(13,file='rho_xy.dat',form='FORMATTED')
      open(14,file='rho_yz.dat',form='FORMATTED')
      open(15,file='rho_xz.dat',form='FORMATTED')

      read(12) (((ZR(ikx,iky,ikz),ikx=1,NX),iky=1,NY),ikz=1,NZ)

      do iky=1,NY
         write(13,333) (ZR(ikx,iky,NZ/2),ikx=1,NX)
         write(13,*) '           '
         if(iopt.eq.1)then
            write(14,333) (ZR(NX/2,ikz,iky)+iky*TWOPI/NY,ikz=1,NZ)
            write(15,333) (ZR(ikx,NY/2,iky)+iky*TWOPI/NY,ikx=1,NX)
         else
            write(15,333) (ZR(ikx,NY/4,iky),ikx=1,NX)
            write(14,333) (ZR(NX/4,ikz,iky),ikz=1,NZ)
         endif
         write(14,*) '           '
         write(15,*) '           '
      enddo
      
C
 333  format(1x,E12.5,1x)
      END
C
c

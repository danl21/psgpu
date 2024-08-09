      program RiG_time
!                                                                                                                                                                                                                                           
!     Average in z the ZT gradient Richardson number
!                                                                                                                                                                                                                                           
      use GLOBAL_PARAMS
      IMPLICIT NONE
      INTEGER :: ikx,iky,ikz,nt,nw,nlx,ip,im,it
C ----------------------------------------------------------------                                                                                                                                                                           
      REAL(rk) :: dy,Navg,Savg,avgS,meanRi,time,told,trun
      REAL(rk) :: NZT(10000,NZ),SZT(10000,NZ),avgRi,meanN,meanS,B,avgN
      TWOPI = 4._rk*ASIN(1._rk)

      print*, 'B? '!  Read in the stratification parameter
      read(*,*) B
      print*, 'NT? '!  Number of ouput steps
      read(*,*) NT

      open(14,file='RiG.dat',form='FORMATTED')
      open(15,file='RiG_t.dat',form='FORMATTED')

      meanRi=0.d0
      meanN=0.d0
      meanS=0.d0
      im=0
      told=0.d0
      trun=0.d0
 
      do it=1,NT
         avgRi=0._rk
         avgN=0._rk
         avgS=0._rk
         read(14,*) 
         do ikz=2,NZ
            read(14,*) time,iky, Navg,Savg
            NZT(it,ikz)=B*Navg
            SZT(it,ikz)=Savg
            avgRi=avgRi+B*Navg/Savg
            avgN=avgN+B*Navg
            avgS=avgS+Savg
         enddo
         if(time.ge.told)then
            trun=trun+(time-told)
         endif
         write(15,*) time,avgRi/NZ,avgN/avgS
         read(14,*)
         meanRi=meanRi+(time-told)*avgRi/NZ
         meanN=meanN+(time-told)*avgN/NZ
         meanS=meanS+(time-told)*avgS/NZ
         im=im+1
         told=time
      enddo
      
C                          
      print*, trun, meanRi/dble(trun),meanN/dble(trun),meanS/dble(trun)

      meanRi=0.d0
      avgRi=0._rk
      do ikz=2,NZ
         meanN=0.d0
         meanS=0.d0
         do it=2,NT
            meanRi=meanRi+NZT(it,ikz)/SZT(it,ikz)
            meanN=meanN+NZT(it,ikz)
            meanS=meanS+SZT(it,ikz)
         enddo
         avgRi=avgRi+meanN/meanS
      enddo
      print*, meanRi/REAL(im*NZ),avgRi/NZ
c
      END

! ------------------------------------------------------------------                                                                                                                                           
      SUBROUTINE INIT (XII,ETA,ZET,RHO,KX,KY,KZ,TSTART,PERIODin,SHIFTXin,SHIFTYin,SHIFTZin,NERin)
! -------------------------------------------------------------------
        ! initialisation subroutine for CUDA pseudospectral DNS code.            

        ! See README for details on usage.
        ! Dan Lucas dl549@cam.ac.uk 2016
! -------------------------------------------------------------------                                                                                                                                                                       
!                                                                            
      use GLOBAL_PARAMS
      IMPLICIT NONE 
!
      INTEGER             :: IKX,IKY,IKZ,ntime
      INTEGER             :: IK,jw,kx1,ky1,kz1,klag,jlag
!
      COMPLEX(rk), DIMENSION(NKT):: XII,ETA,ZET,RHO
!
      REAL(rk), DIMENSION(NKT) :: KX,KY,KZ
      REAL(rk)        :: ranno,TSTART,normRE,normIM
      REAL(rk)        :: KKX,KKY,KKZ,WK,EK
      REAL(rk)        :: PHASE,ENSTRO,R1,R2,VK,R3,PE

      REAL(rk)    :: timein,PERIODin,SHIFTXin,SHIFTYin,SHIFTZin,NERin
      REAL(rk)    :: dummy,NER2
      character(len=256) :: temp

! -------------------------------------------------------------------                                                                                                                                                                       
!
      EXTERNAL RANNO
!
      ZET   = CMPLX(0._rk,0._rk,rk)
      ETA   = CMPLX(0._rk,0._rk,rk)
      XII   = CMPLX(0._rk,0._rk,rk)
      RHO   = CMPLX(0._rk,0._rk,rk)
! -------------------------------------------------------------------
      ! Read from file, the cases switch between the different file types
      if (irest.ne.0) then
         print*,'          IREST = ',IREST
         print*,'RESTART FROM OUTPUT DATA. '
         kx1 = iktx
         ky1 = ikty
         kz1 = iktz
!
         Print*,'infile type = ',trim(InfileType)

         SELECT CASE (trim(InfileType))

            CASE ("A") ! Read from recurrence guess
               open(22,file=trim(simdir)//'/'//trim(infile),form='unformatted',access='stream')

               do ntime = 1, irest
                  read(22) timein,v2,PERIODin,SHIFTXin,SHIFTYin,SHIFTZin,NERin, &
                 & (XII(ikx),ikx=1,nkt),&
                 & (ETA(ikx),ikx=1,nkt),&
                 & (ZET(ikx),ikx=1,nkt),&
                 & (RHO(ikx),ikx=1,nkt)
!
                  print*,'Just read time= ',timein
               enddo
               close(22)
            CASE ("B") ! Read from recurrence guess, ignoring parameter v2
               open(22,file=trim(simdir)//'/'//trim(infile),form='unformatted',access='stream')

               do ntime = 1, irest
                  read(22) timein,dummy,PERIODin,SHIFTXin,SHIFTYin,SHIFTZin,NERin, &
                 & (XII(ikx),ikx=1,nkt),&
                 & (ETA(ikx),ikx=1,nkt),&
                 & (ZET(ikx),ikx=1,nkt),&
                 & (RHO(ikx),ikx=1,nkt)
!
                  print*,'Just read time= ',timein
               enddo
               close(22)
            CASE ("C") ! Read from raw data (no parameters)
               open(22,file=trim(simdir)//'/'//trim(infile),form='unformatted')

               do ntime = 1, irest
                  read(22) (XII(ikx),ikx=1,nkt),&
                 & (ETA(ikx),ikx=1,nkt),&
                 & (ZET(ikx),ikx=1,nkt),&
                 & (RHO(ikx),ikx=1,nkt)
!
                  print*,'Just read time= ',timein
               enddo
               close(22)
            CASE ("D") ! Read from Zk.out DNS output file
               open(22,file=trim(simdir)//'/'//trim(infile),form='unformatted')

               read(22) jw,kx1,ky1,kz1
               print*,'Number of records = ',jw
               print*,'FILE    Iktx,Ikty = ',kx1,ky1,kz1
               print*,'FILE      ktx,kty = ',kx1-1,(ky1-1)/2,(kz1-1)/2
               print*,'CURRENT Iktx,Ikty = ',IKTX,ikty,iktz
               do ntime = 1, irest
                  read(22) timein,v2,PERIODin,SHIFTXin,SHIFTYin,SHIFTZin,NERin, &
                 & (XII(ikx),ikx=1,nkt),&
                 & (ETA(ikx),ikx=1,nkt),&
                 & (ZET(ikx),ikx=1,nkt),&
                 & (RHO(ikx),ikx=1,nkt)
!
                  print*,'Just read time= ',timein
               enddo
               print*,'debug!', timein,v2,PERIODin,SHIFTXin,SHIFTYin,SHIFTZin,NERin
               close(22)
            CASE ("E")! Read from Zk.out DNS output file ignoring parameter v2
               open(22,file=trim(simdir)//'/'//trim(infile),form='unformatted')

               read(22) jw,kx1,ky1,kz1
               print*,'Number of records = ',jw
               print*,'FILE    Iktx,Ikty = ',kx1,ky1,kz1
               print*,'FILE      ktx,kty = ',kx1-1,(ky1-1)/2,(kz1-1)/2
               print*,'CURRENT Iktx,Ikty = ',IKTX,ikty,iktz
               if (kx1.gt.iktx .or. ky1.gt.ikty.or. kz1.gt.iktz) then
                  print*,'Do not know how to decrease resolution.'
                  stop
               endif
               
               klag = INT(kty) - (ky1-1)/2
               jlag = INT(ktz) - (kz1-1)/2
               do ntime = 1, irest
                  read(22) timein,dummy,PERIODin,SHIFTXin,SHIFTYin,SHIFTZin,NERin, &
                 & (XII(ikx),ikx=1,nkt),&
                 & (ETA(ikx),ikx=1,nkt),&
                 & (ZET(ikx),ikx=1,nkt),&
                 & (RHO(ikx),ikx=1,nkt)
!
                  print*,'Just read time= ',timein
               enddo
            CASE ("F") ! Read from Zk.out GMRES output file with ALL parameters 
               open(22,file=trim(simdir)//'/'//trim(infile),form='unformatted')

               do ntime = 1, irest
                  read(22) timein,v2,Ri,alpha,PERIODin,SHIFTXin,SHIFTYin,SHIFTZin,NERin, &
                 & (XII(ikx),ikx=1,nkt),&
                 & (ETA(ikx),ikx=1,nkt),&
                 & (ZET(ikx),ikx=1,nkt),&
                 & (RHO(ikx),ikx=1,nkt)
!
                  print*,'Just read time= ',timein
               enddo
               close(22)
            CASE ("G")! Read from Zk.out GMRES output file with ALL parameters and resolution
               open(22,file=trim(simdir)//'/'//trim(infile),form='unformatted')
                      read(22) jw,kx1,ky1,kz1

               print*,'Number of records = ',jw
               print*,'FILE    Iktx,Ikty = ',kx1,ky1,kz1
               print*,'FILE      ktx,kty = ',kx1-1,(ky1-1)/2,(kz1-1)/2
               print*,'CURRENT Iktx,Ikty = ',IKTX,ikty,iktz
               do ntime = 1, irest
                  read(22) timein,v2,Ri,alpha,PERIODin,SHIFTXin,SHIFTYin,SHIFTZin,NERin, &
                 & (XII(ikx),ikx=1,nkt),&
                 & (ETA(ikx),ikx=1,nkt),&
                 & (ZET(ikx),ikx=1,nkt),&
                 & (RHO(ikx),ikx=1,nkt)
!                                                                                                                      
                  print*,'Just read time= ',timein
                         enddo
               close(22)
            CASE ("H")! Read from Zk.out GMRES output file with ALL parameters, resolution ignoring alpha
               open(22,file=trim(simdir)//'/'//trim(infile),form='unformatted')
                      read(22) jw,kx1,ky1,kz1

               print*,'Number of records = ',jw
               print*,'FILE    Iktx,Ikty = ',kx1,ky1,kz1
               print*,'FILE      ktx,kty = ',kx1-1,(ky1-1)/2,(kz1-1)/2
               print*,'CURRENT Iktx,Ikty = ',IKTX,ikty,iktz
               do ntime = 1, irest
                  read(22) timein,v2,Ri,dummy,PERIODin,SHIFTXin,SHIFTYin,SHIFTZin,NERin, &
                 & (XII(ikx),ikx=1,nkt),&
                 & (ETA(ikx),ikx=1,nkt),&
                 & (ZET(ikx),ikx=1,nkt),&
                 & (RHO(ikx),ikx=1,nkt)
!                                                                                                                      
                  print*,'Just read time= ',timein
                         enddo
               close(22)
            CASE DEFAULT
               Print*,'infile type unknown, stopping code'
               STOP
!
          END SELECT

      endif

! ------------------------------------------------------------------                                                                                                                                           
      ! Set spectral filter for dealiasing 
! ------------------------------------------------------------------                                                                                                                                           
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

               L(IK)  = 0
               LL(IK) = .TRUE.
               WK = KKX*KKX + KKY*KKY + KKZ*KKZ 
               WK = SQRT( WK )

               IF (KKX .LT. 0._rk) LL(IK)=.FALSE.
               IF (KKX.EQ.0._rk .AND. KKY.LT.0._rk ) LL(IK) = .FALSE.
               IF (KKX.EQ.0._rk .AND. KKY.EQ.0._rk .AND. KKZ.LE.0._rk) LL(IK) = .FALSE.
!               IF (KKX.EQ.0._rk .AND. KKY.EQ.0._rk .AND. KKZ.EQ.0._rk) LL(IK) = .FALSE.
               IF (WK.GT.REAL(KTY,rk)-0.5_rk) LL(IK) = .FALSE.
               IF ( LL(IK) )THEN
                  L(IK)  = 1
                  L3(IK) = .TRUE.
                  L3(NKT+IK) = .TRUE.
                  L3(2*NKT+IK) = .TRUE.
                  L3(3*NKT+IK) = .TRUE.
               ELSE
                  L(IK)  = 0
                  L3(IK) = .FALSE.
                  L3(NKT+IK) = .FALSE.
                  L3(2*NKT+IK) = .FALSE.
                  L3(3*NKT+IK) = .FALSE.
               ENDIF

            ENDDO
         ENDDO
      ENDDO


! ------------------------------------------------------------------                                                                                                                                           
      ! If not reading data, set randomised initial condition
! ------------------------------------------------------------------                                                                                                                                           
      PE = 0._rk
      ENSTRO = 0.0
      DO IKX=1,IKTX
         DO IKY=1,IKTY
            DO IKZ=1,IKTZ
               IK = ((IKZ-1)*IKTY +IKY-1)*IKTX + IKX
               WK = SQRT(KX(IK)*KX(IK)+KY(IK)*KY(IK)+KZ(IK)*KZ(IK))
               !Initial Condition is set by setting the desired energy Spectrum
               ! ---------------------------------------------------------------- 
               IF ( WK.GT.2.5_rk .AND. WK.LE.9.5_rk ) THEN
                  EK = 0.1_rk !WK**(-1.6667_rk)
               ELSE
                  EK = 0.0_rk
               ENDIF

                IF (irest .eq. 0 ) THEN
                   PHASE = RANNO(0)
                   XII(IK) = EK*EXP(ZI*PHASE)
                   PHASE = RANNO(0)
                   ETA(IK) = EK*EXP(ZI*PHASE)
                   PHASE = RANNO(0)
                   ! MAKE INCOMPRESSIBLE
                   IF(KZ(IK) .ne. 0._rk)then
                      ZET(IK) = -(KX(IK)*XII(IK)+KY(IK)*ETA(IK))/KZ(IK)
                   else
                      XII(IK) = CMPLX(0._rk,0._rk)
                      ETA(IK) = CMPLX(0._rk,0._rk)
                      ZET(IK) = CMPLX(0._rk,0._rk)
                   endif
                ENDIF
!
               R1     = REAL(XII(IK),rk)**2._rk+ &
     &                  REAL(ETA(IK),rk)**2._rk+ &
     &                  REAL(ZET(IK),rk)**2._rk

               R2     = AIMAG(XII(IK))**2._rk+ &
     &                  AIMAG(ETA(IK))**2._rk+ &
     &                  AIMAG(ZET(IK))**2._rk

               ENSTRO = ENSTRO + R1 + R2
               PE = PE+REAL(RHO(IK)*CONJG(RHO(IK)))
            ENDDO
         ENDDO
      ENDDO
!
      ENSTRO = ENSTRO*2._rk
!
! -------------------------------------------------------
!     
      IF (irest.EQ.0 .and. ENSTRO .NE. 0._rk) THEN
!     NORMALISE AMPLITUDES:
      DO IKX=1,IKTX
         DO IKY=1,IKTY
            DO IKZ=1,IKTZ
               IK = ((IKZ-1)*IKTY +IKY-1)*IKTX + IKX
               WK = SQRT(KX(IK)*KX(IK)+KY(IK)*KY(IK)+KZ(IK)*KZ(IK))

               IF (L(IK).EQ.1) THEN
                  XII(IK) = XII(IK) * SQRT(AMPV/ENSTRO)
                  ETA(IK) = ETA(IK) * SQRT(AMPV/ENSTRO)
                  ZET(IK) = ZET(IK) * SQRT(AMPV/ENSTRO)
               ELSE
                  XII(IK) = CMPLX(0._rk,0._rk,rk)
                  ETA(IK) = CMPLX(0._rk,0._rk,rk)
                  ZET(IK) = CMPLX(0._rk,0._rk,rk)
               ENDIF

            ENDDO
         ENDDO
      ENDDO

   ENDIF
!
      END SUBROUTINE INIT
!
!
!
!
      FUNCTION ranno (i)
!
!     Controls random number generator.
!-----------------------------------------
!   - If argument i.ne.0 it performs
!     initialization with i=seed no.
!   - If argument i.eq.0 it draws a
!     random no.
!-----------------------------------------
      use GLOBAL_PARAMS    
      implicit none
      integer  :: i,junk,ihold,iv(32),iy
      real(rk) :: ranno,ran1
      common /random/ junk,iv,iy
!
      if (i.ne.0) then
        if (i.gt.0) i = - i
        junk  = i
        ranno = (ran1(i)-0.5_rk)*twopi
      else
        junk  = junk - 1
        ihold = junk
        ranno = (ran1(ihold)-0.5_rk)*twopi
      endif
!
      END FUNCTION ranno
!
!
!
      FUNCTION ran1(idum)
      use GLOBAL_PARAMS    
      implicit none
      integer  :: idum,ia,im,iq,ir,ntab,ndiv
      real(rk) :: ran1,am,eps,rnmx
      parameter (ia=16807,im=2147483647,am=1._rk/im,iq=127773,ir=2836,  & 
     &        ntab=32,ndiv=1+(im-1)/ntab,eps=1.2E-7_rk,rnmx=1._rk-eps)
      integer  :: j,k,iv(32),iy,junk
      common /random/ junk,iv,iy
!
      if (idum.le.0.or.iy.eq.0) then
         idum=max(-idum,1)
         do j=ntab+8,1,-1
            k=idum/iq
            idum=ia*(idum-k*iq)-ir*k
            if (idum.lt.0) idum=idum+im
            if (j.le.ntab) iv(j)=idum
         enddo
         iy=iv(1)
      endif
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      j=1+iy/ndiv
      iy=iv(j)
      iv(j)=idum
      ran1=min(am*iy,rnmx)
!
      END FUNCTION ran1
!
!
      block data
      parameter (ntab=32)
      integer   :: iv(32),iy,junk
      common /random/ junk,iv,iy
      data iv /ntab*0/, iy /0/
      end

! -------------------------------------------------------------------!                                                                            !
      module PARAMETERS
! -------------------------------------------------------------------
!      Module for the constant parameters for the simulation. 
      use GLOBAL_PARAMS 
      implicit none
!
      REAL(rk) :: Re,outfreq1,outfreq2
!
      contains
!------------------------------------------------
! Compute common groupings
!
      subroutine params_update()
      implicit none
!
      v2     =   1._rk / Re ! "viscosity"
      alpha = alpha1 ! working aspect ratio
      AMPFOR =  -0.5_rk*KF ! amplitude of forcing (when unthrottled)
      NSTOP  =   INT(TSTEP/DELT) !Total timesteps Note this rounds down
      NOUT   =   INT(outfreq1/DELT) ! high frequency outputs
      NOUTV  =   INT(outfreq2/DELT) ! low frequency outputs
!
      ZI = CMPLX(0._rk,1._rk,rk) ! imaginary number i
      TWOPI = 4._rk*ASIN(1._rk)  ! 2*pi
      ROOT2 = SQRT(2._rk) 
      THETA = THETA*TWOPI/360._rk ! convert angle from degrees to radians
!
      end subroutine params_update
!
!------------------------------------------------
!
      subroutine read_label(u, n)
      implicit none
!
      integer, intent(in) :: u
      character(len=*), intent(in) :: n
!
      integer :: i
      character :: t
!
      do i = 1,LEN(n)
         read(u,"(A1)",ADVANCE='NO') t
         if (t /= n(i:i)) then
            write(*,*) "bad label ", n
            stop
         endif
      enddo
!
      end subroutine read_label
!
!------------------------------------------------
!
      subroutine read_int(u, n, v)
      implicit none
!
      integer, intent(in) :: u
      character(len=*), intent(in) :: n
      integer, intent(out) :: v
!
      call read_label(u, n)
!
      read(u,"(I10)") v
      write(*,*) n, v
      end subroutine read_int
!
!------------------------------------------------
!
      subroutine read_logical(u, n, v)
      implicit none
!
      integer, intent(in) :: u
      character(len=*), intent(in) :: n
      logical, intent(out) :: v
!     
      call read_label(u, n)
!
      read(u,"(L1)") v
      write(*,*) n, v
      end subroutine read_logical
!
!------------------------------------------------
!
      subroutine read_float(u, n, v)
      implicit none
!
      integer, intent(in) :: u
      character(len=*), intent(in) :: n
      real(rk), intent(out) :: v
!      
      call read_label(u, n)
!      
      read(u,"(F12.4)") v
      write(*,*) n, v
      end subroutine read_float
!      
!------------------------------------------------
!
      subroutine read_string(u, n, v)
      implicit none
!      
      integer, intent(in) :: u
      character(len=*), intent(in) :: n
      character(len=*), intent(out) :: v
!      
      call read_label(u, n)
!      
      read(u,'(A256)') v
      write(*,*) n, v
      end subroutine read_string
!      
!------------------------------------------------
! Read parameters from a simple text file
!      
      subroutine params_read_txt()
      implicit none
!      
      open(55, FILE=trim(simdir)//'/params.txt', ACTION='READ')
!
      call read_int(55, 'ISEED = ', ISEED)
      call read_int(55, 'IREST = ', IREST)
      call read_string(55, 'infile = ', infile)
      call read_string(55, 'infile type = ', InfileType)
      call read_int(55, 'KF = ', KF)
      call read_float(55, 'Re = ', Re)
      call read_float(55, 'Ri = ', Ri)
      call read_float(55, 'Theta = ', Theta)
      call read_float(55, 'Sc = ', Sc)
      call read_float(55, 'TSTEP = ', TSTEP)
      call read_float(55, 'ResidualThreshold = ', ResidualThreshold)
      call read_float(55, 'outfreq1 = ', outfreq1)
      call read_float(55, 'outfreq2 = ', outfreq2)
      call read_float(55, 'DELT = ', DELT)
      call read_float(55, 'AMPV = ', AMPV)              
      call read_int(55, 'STATSFLAG = ', STATSFLAG)
      call read_int(55, 'RCFLAG = ', RCFLAG)
      call read_int(55, 'UPOflag = ', UPOflag)
      call read_int(55, 'ADAPTflag = ', ADAPTFLAG)
      call read_int(55, 'ISTOP = ', ISTOP)
      call read_float(55, 'Dtarget = ', Dtarget)              

!
      close(55)
!
      call params_update()
!
      end subroutine params_read_txt
!
!
!------------------------------------------------
!
      end module parameters
!

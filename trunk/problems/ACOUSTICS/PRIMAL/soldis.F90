!--------------------------------------------------------------------
!
!     routine name      - soldis
!
!--------------------------------------------------------------------
!
!     latest revision:  - Aug 17
!
!     purpose:          - display user-defined quantity
!
!     arguments:
!
!     in:
!             Mdle   - element (middle node) number 
!             Xi     - master element coordinates
!             X      - physical coordinates
!             Rn     - outward normal unit vector
!             ZsolH  - H1    sol
!             ZgradH - H1    grad
!             ZsolE  - Hcurl sol
!             ZcurlE - Hcurl curl
!             ZsolV  - Hdiv  sol
!             ZdivV  - Hdiv  div
!             ZsolQ  - L2    sol
!     out:
!             val   - quantity to display
!
!---------------------------------------------------------------------
! 
   subroutine soldis(Mdle,Xi,X,Rn,ZsolH,ZgradH,ZsolE,ZcurlE,ZsolV,ZdivV,ZsolQ, Val)
!      
   use common_prob_data, only : IEXACT_DISP, ICHOOSE_DISP
   use data_structure3D
! 
   implicit none
! 
   integer,                      intent(in)  :: Mdle
   real*8,dimension(3),          intent(in)  :: Xi,X,Rn
   complex*16,dimension(  MAXEQNH  ),intent(in)  :: ZsolH
   complex*16,dimension(  MAXEQNH,3),intent(in)  :: ZgradH
   complex*16,dimension(3,MAXEQNE  ),intent(in)  :: ZsolE
   complex*16,dimension(3,MAXEQNE  ),intent(in)  :: ZcurlE
   complex*16,dimension(3,MAXEQNV  ),intent(in)  :: ZsolV
   complex*16,dimension(  MAXEQNV  ),intent(in)  :: ZdivV
   complex*16,dimension(  MAXEQNQ  ),intent(in)  :: ZsolQ
   complex*16                                    :: zval
   real*8,                           intent(out) :: Val
! 
!..exact solution workspace
   integer                             :: icase
   complex*16,dimension(  MAXEQNH    ) ::   zvalH
   complex*16,dimension(  MAXEQNH,3  ) ::  zdvalH
   complex*16,dimension(  MAXEQNH,3,3) :: zd2valH
   complex*16,dimension(3,MAXEQNE    ) ::   zvalE
   complex*16,dimension(3,MAXEQNE,3  ) ::  zdvalE
   complex*16,dimension(3,MAXEQNE,3,3) :: zd2valE
   complex*16,dimension(3,MAXEQNV    ) ::   zvalV
   complex*16,dimension(3,MAXEQNV,3  ) ::  zdvalV
   complex*16,dimension(3,MAXEQNV,3,3) :: zd2valV
   complex*16,dimension(  MAXEQNQ    ) ::   zvalQ
   complex*16,dimension(  MAXEQNQ,3  ) ::  zdvalQ
   complex*16,dimension(  MAXEQNQ,3,3) :: zd2valQ
!   
!---------------------------------------------------------------------------------------
!
   select case (IEXACT_DISP)
   case(1)
      icase=0
      call exact(X,icase, zvalH,zdvalH,zd2valH, &
                          zvalE,zdvalE,zd2valE, &
                          zvalV,zdvalV,zd2valV, &
                          zvalQ,zdvalQ,zd2valQ)
      select case (ICHOOSE_DISP)
!
!  ...exact solution
      case(1)
         zval = ZvalH(1)
!
!  ...exact (normal) flux
      case(2)
          zval = ZvalV(1,1)*Rn(1) + ZvalV(2,1)*Rn(2) + ZvalV(3,1)*Rn(3) 
      end select
!
   case(2)
      select case (ICHOOSE_DISP)
!
!  ...approximate solution
      case(1)
         Val = ZsolH(1)
!
!  ...approximate (normal) flux
      case(2)
         Val = ZsolV(1,1)*Rn(1) + ZsolV(2,1)*Rn(2) + ZsolV(3,1)*Rn(3) 
      end select

   case(3)
      icase=0
      call exact(X,icase, zvalH,zdvalH,zd2valH, &
                          zvalE,zdvalE,zd2valE, &
                          zvalV,zdvalV,zd2valV, &
                          zvalQ,zdvalQ,zd2valQ)
      
      select case (ICHOOSE_DISP)
!
!  ...approximate solution
      case(1)
         Val = ZvalH(1) - ZsolH(1)
!
!  ...approximate (normal) flux
      case(2)
         Val = ZvalV(1,1)*Rn(1) + ZvalV(2,1)*Rn(2) + ZvalV(3,1)*Rn(3) -  &
              (ZsolV(1,1)*Rn(1) + ZsolV(2,1)*Rn(2) + ZsolV(3,1)*Rn(3))

      end select
   end select

!..choose between real and imaginary parts
   if (ICHOOSE_DISP.gt.0) then
      Val = dreal(zval)
   else
      Val = dimag(zval)
   endif
! 
!
   end subroutine soldis
!   
!--------------------------------------------------------------------
!
!     routine name      - soldis_select
!
!--------------------------------------------------------------------
!
!     latest revision:  - Aug 17
!
!     purpose:          - select quantities to display 
!
!---------------------------------------------------------------------
! 
   subroutine soldis_select
!
   use control, ONLY: NEXACT      
   use parameters, only : NSTD_OUT
   use common_prob_data, only : ICHOOSE_DISP, IEXACT_DISP
!   
   implicit none
! 
   integer :: iprev
   integer, save :: ivis=0 ! visitation flag
!   
!---------------------------------------------------------------------------------------
!
       if (ivis.eq.1) then
   10   write(*,*) 'USE PREVIOUS SELECTION? 1-Y,0-N'
        read(*,*) iprev
        select case(iprev)
        case(0)
        case(1)
          call disp_soldis(NSTD_OUT) ; return
        case default
          go to 10
        end select
      else
        ivis=1
      endif
!
      if ((NEXACT.eq.1).or.(NEXACT.eq.2)) then
   20   write(*,*) 'SET VARIABLE     1) EXACT 2) APPROX 3) ERROR'
        read(*,*) IEXACT_DISP
        select case(IEXACT_DISP)
        case(1,2,3)
        case default
          go to 20
        end select
      else
        IEXACT_DISP=2
      endif  
!
   30 write(*,*) 'SET VARIABLE: 1) p, 2) \hat p '
      write(*,*) 'USE NEGATIVE FLAGS FOR IMAGINARY PARTS'
      read(*,*) ICHOOSE_DISP
      select case(ICHOOSE_DISP)
      case(1,2)
      case(-1,-2)
      case default
        go to 30
      end select
!
      call disp_soldis(NSTD_OUT)

! 
   end subroutine soldis_select


   subroutine disp_soldis(Nstream)
!      
   use common_prob_data, only : IEXACT_DISP, ICHOOSE_DISP
! 
   implicit none
   integer, intent(in) :: Nstream

      write(Nstream,1005)
 1005 format('DISPLAY SETUP')
      write(Nstream,1000)
 1000 format('-----------------------')
      select case (IEXACT_DISP)
      case(1); write(Nstream,1010)
 1010 format('EXACT  SOLUTION SELECTED')
      case(2); write(Nstream,1020)
 1020 format('APPROX SOLUTION SELECTED')
      case(3); write(Nstream,1030)
 1030 format('ERROR SELECTED')
      end select
      if (ICHOOSE_DISP .eq. 1)  write(Nstream,1040)
 1040 format('REAL PART OF PRESSURE SELECTED')
      if (ICHOOSE_DISP .eq.-1)  write(Nstream,1041)
 1041 format('IMAG PART OF PRESSURE SELECTED')
      if (ICHOOSE_DISP .eq. 2)  write(Nstream,1050)
 1050 format('REAL PART OF TRACE SELECTED')
      if (ICHOOSE_DISP .eq.-2)  write(Nstream,1051)
 1051 format('IMAG PART OF TRACE SELECTED')
      write(*,1000)
!
!
!
   end subroutine disp_soldis
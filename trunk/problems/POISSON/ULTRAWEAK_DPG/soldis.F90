!---------------------------------------------------------------------------------------
!> Purpose : display user-defined quantity
!
!! @param[in] Mdle   - element (middle node) number
!! @param[in] Xi     - master element coordinates
!! @param[in] X      - physical coordinates
!! @param[in] Rn     - outward normal unit vector
!! @param[in] RsolH  - H1    sol
!! @param[in] RgradH - H1    grad
!! @param[in] RsolE  - Hcurl sol
!! @param[in] RcurlE - Hcurl curl
!! @param[in] RsolV  - Hdiv  sol
!! @param[in] RdivV  - Hdiv  div
!! @param[in] RsolQ  - L2    sol
!
!! @param[out] Val   - quantity to display
!---------------------------------------------------------------------------------------
!
#include "typedefs.h"
!
subroutine soldis(Mdle,Xi,X,Rn,RsolH,RgradH,RsolE,RcurlE,RsolV,RdivV,RsolQ, Val)
!
   use data_structure3D
   use common_prob_data
!
   implicit none
!
   integer,                       intent(in)  :: Mdle
   real(8),dimension(3),          intent(in)  :: Xi,X,Rn
   VTYPE  ,dimension(  MAXEQNH  ),intent(in)  :: RsolH
   VTYPE  ,dimension(  MAXEQNH,3),intent(in)  :: RgradH
   VTYPE  ,dimension(3,MAXEQNE  ),intent(in)  :: RsolE
   VTYPE  ,dimension(3,MAXEQNE  ),intent(in)  :: RcurlE
   VTYPE  ,dimension(3,MAXEQNV  ),intent(in)  :: RsolV
   VTYPE  ,dimension(  MAXEQNV  ),intent(in)  :: RdivV
   VTYPE  ,dimension(  MAXEQNQ  ),intent(in)  :: RsolQ
   real(8),                       intent(out) :: Val
!
#if DEBUG_MODE
   integer :: iprint
   iprint = 0
#endif
!
!---------------------------------------------------------------------------------------
!
   select case (IEXACT_DISP)
      case(0)
         select case (ICHOOSE_COMP)
            case(1)
               Val = (RsolQ(1))
            case(2)
               Val = (RsolQ(2))
            case(3)
               Val = (RsolQ(3))
            case(4)
               Val = (RsolQ(4))
         end select
   end select
!
#if DEBUG_MODE
   if (iprint.eq.1) then
      write(*,7003) IEXACT_DISP, ICHOOSE_COMP, Val
 7003 format('IEXACT_DISP, ICHOOSE_COMP, Val = ',2i2,e12.5)
      call pause
   endif
#endif
!
end subroutine soldis
!
subroutine soldis_select
   use control,    only: NEXACT    ! exact solution flag
   use parameters, only: NSTD_OUT  ! display file
   use common_prob_data
!
   implicit none
!
   integer :: iprev
!
!
   select case(ICHOOSE_COMP)
      case(1,2)
         write(*,*) 'USE PREVIOUS SELECTION? 1-Y,0-N'
         read(*,*) iprev
         if (iprev.eq.1) then
            call disp_soldis(NSTD_OUT) ; return
         endif
   end select
!
!
   IEXACT_DISP=0
   if ((NEXACT.eq.1).or.(NEXACT.eq.2)) then
   10 write(*,*) 'DISPLAY APPROXIMATE OR EXACT SOLUTION (0/1) ?'
      read(*,*) IEXACT_DISP
      if ((IEXACT_DISP.ne.0).and.(IEXACT_DISP.ne.1)) goto 10
   endif
!
!
   20 write(*,*) 'SET VARIABLE: u(1), sigma(2-4)'
   read(*,*) ICHOOSE_COMP
   if (ICHOOSE_COMP.lt.1) goto 20
!
   call disp_soldis(NSTD_OUT)
!
end subroutine soldis_select
!
!
!
subroutine disp_soldis(Nstream)
!
   use common_prob_data, only: IEXACT_DISP,ICHOOSE_COMP
!
   implicit none
   integer, intent(in) :: Nstream
!
   write(Nstream,1000)
   write(Nstream,3100)
   select case (IEXACT_DISP)
      case(1); write(Nstream,1010)
      case(2); write(Nstream,1020)
   end select
!
   select case (ICHOOSE_COMP)
      case(1); write(Nstream,3010)
   end select
!
   write(Nstream,3100)
!
   1000 format('DISPLAY SETUP')
   1010 format('DISPLAYING EXACT SOLUTION')
   1020 format('DISPLAYING APPROXIMATE SOLUTION')
   3010 format('DISPLAYING u')
   3100 format('-----------------------')
!
end subroutine disp_soldis

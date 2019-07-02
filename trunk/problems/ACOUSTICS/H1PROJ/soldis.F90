!---------------------------------------------------------------------------------------
!> Purpose : display user-defined quantity
!! @param[in] Mdle   - element (middle node) number
!! @param[in] Xi     - master element coordinates
!! @param[in] X      - physical coordinates
!! @param[in] Rn     - outward normal unit vector
!! @param[in] ZsolH  - H1    sol
!! @param[in] ZgradH - H1    grad
!! @param[in] ZsolE  - Hcurl sol
!! @param[in] ZcurlE - Hcurl curl
!! @param[in] ZsolV  - Hdiv  sol
!! @param[in] ZdivV  - Hdiv  div
!! @param[in] ZsolQ  - L2    sol
!!
!! @param[out] val   - quantity to display
!---------------------------------------------------------------------------------------
!
subroutine soldis(Mdle,Xi,X,Rn,ZsolH,ZgradH,ZsolE,ZcurlE,ZsolV,ZdivV,ZsolQ, Val)
!   use common_prob_data, only : IEXACT_DISP, ICHOOSE_DISP
!   use data_structure3D
! !---------------------------------------------------------------------------------------
!   implicit none
!   integer,                      intent(in)  :: Mdle
!   real*8,dimension(3),          intent(in)  :: Xi,X,Rn
!   complex*16,dimension(  MAXEQNH  ),intent(in)  :: ZsolH
!   complex*16,dimension(  MAXEQNH,3),intent(in)  :: ZgradH
!   complex*16,dimension(3,MAXEQNE  ),intent(in)  :: ZsolE
!   complex*16,dimension(3,MAXEQNE  ),intent(in)  :: ZcurlE
!   complex*16,dimension(3,MAXEQNV  ),intent(in)  :: ZsolV
!   complex*16,dimension(  MAXEQNV  ),intent(in)  :: ZdivV
!   complex*16,dimension(  MAXEQNQ  ),intent(in)  :: ZsolQ
!   complex*16,                       intent(out) :: Val
! !---------------------------------------------------------------------------------------
! ! exact solution workspace
!   integer                         :: icase
!   complex*16,dimension(  MAXEQNH    ) ::   zvalH
!   complex*16,dimension(  MAXEQNH,3  ) ::  zdvalH
!   complex*16,dimension(  MAXEQNH,3,3) :: zd2valH
!   complex*16,dimension(3,MAXEQNE    ) ::   zvalE
!   complex*16,dimension(3,MAXEQNE,3  ) ::  zdvalE
!   complex*16,dimension(3,MAXEQNE,3,3) :: zd2valE
!   complex*16,dimension(3,MAXEQNV    ) ::   zvalV
!   complex*16,dimension(3,MAXEQNV,3  ) ::  zdvalV
!   complex*16,dimension(3,MAXEQNV,3,3) :: zd2valV
!   complex*16,dimension(  MAXEQNQ    ) ::   zvalQ
!   complex*16,dimension(  MAXEQNQ,3  ) ::  zdvalQ
!   complex*16,dimension(  MAXEQNQ,3,3) :: zd2valQ
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!
  ! select case (IEXACT_DISP)
  ! case(1)
  !   icase=0
  !   call exact(X,icase, zvalH,zdvalH,zd2valH, &
  !                       zvalE,zdvalE,zd2valE, &
  !                       zvalV,zdvalV,zd2valV, &
  !                       zvalQ,zdvalQ,zd2valQ)
  !   ! ** Exact  solution - displacement
  !   tmpVal = zvalQ(1:3)
  ! case(2)
  !   ! ** Approximate solution - displacement
  !   tmpVal = ZsolQ(1:3)
  ! case(3)
  !   icase=0
  !   call exact(X,icase, zvalH,zdvalH,zd2valH, &
  !                       zvalE,zdvalE,zd2valE, &
  !                       zvalV,zdvalV,zd2valV, &
  !                       zvalQ,zdvalQ,zd2valQ)

  !   ! ** Solution error - displacement
  !   tmpVal = zvalQ(1:3)-ZsolQ(1:3)
  ! end select

  ! select case (ICHOOSE_DISP)
  ! case(1,2,3)
  !   Val = tmpVal(ICHOOSE_DISP)
  ! case(4)
  !   Val = dsqrt(tmpVal(1)**2 + tmpVal(2)**2 + tmpVal(3)**2)
  ! end select

!
end subroutine soldis
!
!
!
!---------------------------------------------------------------------------------------
!> Purpose : show the quantities to display
!---------------------------------------------------------------------------------------
subroutine soldis_select
!   use parameters, only : NSTD_OUT
!   use common_prob_data, only : ICHOOSE_DISP, IEXACT_DISP
! !---------------------------------------------------------------------------------------
!   implicit none
!   integer :: iprev
! !---------------------------------------------------------------------------------------
! !
!   write(*,*) 'USE PREVIOUS SELECTION? 1-Y,0-N'
!   read(*,*) iprev
!   if (iprev.eq.1) then
!     call disp_soldis(NSTD_OUT) ; return
!   endif
!   !
!   write(*,*) 'EXACT SOLUTION (1) or APPROXIMATE SOLUTION (2) or SIGNED ERROR (3)?'
!   read(*,*) IEXACT_DISP
!     write(*,200)
! 200 format('SET VARIABLE     1) u1    2) u2    3) u3    4) ||u||') 
!   read(*,*) ICHOOSE_DISP
!   select case(ICHOOSE_DISP)
!   case(1,2,3,4)
!   case default
!     write(*,200)
!   end select
!   !
!   call disp_soldis(NSTD_OUT)
! !
!
!---------------------------------------------------------------------------------------
!> Purpose : Soldis display interface
!---------------------------------------------------------------------------------------
end subroutine soldis_select

subroutine disp_soldis(Nstream)
!   use common_prob_data, only : IEXACT_DISP, ICHOOSE_DISP
! !---------------------------------------------------------------------------------------
!   implicit none
!   integer, intent(in) :: Nstream
! !---------------------------------------------------------------------------------------

!   write(Nstream,100)
!   write(Nstream,310)
!   select case (IEXACT_DISP)
!   case(1); write(Nstream,101)
!   case(2); write(Nstream,102)
!   case(3); write(Nstream,103)
!   end select
!   select case (ICHOOSE_DISP)
!   case(1); write(Nstream,301)
!   case(2); write(Nstream,302)
!   case(3); write(Nstream,303)
!   case(4); write(Nstream,304)
!   end select
!   write(Nstream,310)
! !
! 100 format('DISPLAY SETUP')
! 101 format('EXACT SOLUTION IS CHOSEN')
! 102 format('APPROX SOLUTION IS CHOSEN')
! 103 format('(EXACT SOLUTION - APPROX SOLUTION) IS CHOSEN')
! 301 format('u1 IS CHOSEN')
! 302 format('u2 IS CHOSEN')
! 303 format('u3 IS CHOSEN')
! 304 format('||u|| IS CHOSEN')
! 310 format('-----------------------')
!
end subroutine disp_soldis
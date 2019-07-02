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
  use lapl
  use data_structure3D
!---------------------------------------------------------------------------------------
  implicit none
  integer,                      intent(in)  :: Mdle
  real*8,dimension(3),          intent(in)  :: Xi,X,Rn
  real*8,dimension(  MAXEQNH  ),intent(in)  :: ZsolH
  real*8,dimension(  MAXEQNH,3),intent(in)  :: ZgradH
  real*8,dimension(3,MAXEQNE  ),intent(in)  :: ZsolE
  real*8,dimension(3,MAXEQNE  ),intent(in)  :: ZcurlE
  real*8,dimension(3,MAXEQNV  ),intent(in)  :: ZsolV
  real*8,dimension(  MAXEQNV  ),intent(in)  :: ZdivV
  real*8,dimension(  MAXEQNQ  ),intent(in)  :: ZsolQ
  real*8,                       intent(out) :: Val
!---------------------------------------------------------------------------------------
!
  select case (IEXACT_DISP)
  case(2)
     ! ** Approximation solution
     select case (ICHOOSE_DISP)
     case(1)
        Val = ZsolH(ICHOOSE_DISP)
     case (2,3,4)
        Val = ZgradH(1,ICHOOSE_DISP-1) 
     end select
  end select
!  
end subroutine soldis
!
!
!
!---------------------------------------------------------------------------------------
!> Purpose : show the quantities to display 
!---------------------------------------------------------------------------------------
subroutine soldis_select
  use lapl
!---------------------------------------------------------------------------------------
  implicit none
  integer :: iprev
!---------------------------------------------------------------------------------------
!
  write(*,*) 'USE PREVIOUS SELECTION? 1-Y,0-N'
  read(*,*) iprev
  if (iprev.eq.1) then
    call disp_soldis(NSTD_OUT) ; return
  endif
  !
  write(*,*) 'SET VARIABLE     1) u     2) u_x   3) u_y   4) u_z'
  read(*,*) iprev
  ICHOOSE_DISP = iprev
  !
  call disp_soldis(NSTD_OUT)
!  
end subroutine soldis_select

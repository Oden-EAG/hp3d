!----------------------------------------------------------------------------------
!> Purpose : calculate dirichlet boundary condition
!! @param[in]  X       - physical coordinates of a point
!! @param[in]  Icase   - node case
!!
!! @param[out] ZvalH   - value of the H1 solution
!! @param[out] ZdvalH  - H1 corresponding first derivatives
!! @param[out] ZvalE   - value of the H(curl) solution
!! @param[out] ZdvalE  - H(curl) corresponding first derivatives
!! @param[out] ZvalV   - value of the H(div) solution
!! @param[out] ZdvalV  - H(div) corresponding first derivatives
!----------------------------------------------------------------------------------
#include "implicit_none.h"
!
subroutine dirichlet(X,Icase, ZvalH,ZdvalH,ZvalE,ZdvalE,ZvalV,ZdvalV)
!
  use control    , only : NEXACT
  use parameters , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ, ZERO
!----------------------------------------------------------------------------------
  implicit none
  real*8, dimension(3),          intent(in)  :: X
  integer,                       intent(in)  :: Icase
! exact solution  
  VTYPE,dimension(  MAXEQNH    ),intent(out) ::   zvalH
  VTYPE,dimension(  MAXEQNH,3  ),intent(out) ::  zdvalH
  VTYPE,dimension(  MAXEQNH,3,3)             :: zd2valH
  VTYPE,dimension(3,MAXEQNE    ),intent(out) ::   zvalE
  VTYPE,dimension(3,MAXEQNE,3  ),intent(out) ::  zdvalE
  VTYPE,dimension(3,MAXEQNE,3,3)             :: zd2valE
  VTYPE,dimension(3,MAXEQNV    ),intent(out) ::   zvalV
  VTYPE,dimension(3,MAXEQNV,3  ),intent(out) ::  zdvalV
  VTYPE,dimension(3,MAXEQNV,3,3)             :: zd2valV
  VTYPE,dimension(  MAXEQNQ    )             ::   zvalQ
  VTYPE,dimension(  MAXEQNQ,3  )             ::  zdvalQ
  VTYPE,dimension(  MAXEQNQ,3,3)             :: zd2valQ
!----------------------------------------------------------------------------------
!
  select case(NEXACT)
  case(0);   
    ZvalH = ZERO; ZdvalH = ZERO
    ZvalE = ZERO; ZdvalE = ZERO
    ZvalV = ZERO; ZdvalV = ZERO
!
!**********************************************************************************
!   Y O U R    D I R I C H L E T    D A T A    H E R E
!**********************************************************************************
!
  case(1)
    call exact(X,Icase, zvalH,zdvalH,zd2valH, zvalE,zdvalE,zd2valE,   &
                        zvalV,zdvalV,zd2valV, zvalQ,zdvalQ,zd2valQ)
  case default
     write(*,*)'dirichlet: UNKNOWN EXACT SOLUTION FLAG', NEXACT
     stop 1
  end select
!
!
end subroutine dirichlet

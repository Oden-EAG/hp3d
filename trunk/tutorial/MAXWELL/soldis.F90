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
!
#include"implicit_none.h"
  use em
  use data_structure3D
!---------------------------------------------------------------------------------------
  implicit none
  integer,                     intent(in)  :: Mdle
  VTYPE,dimension(3),          intent(in)  :: Xi,X,Rn
  VTYPE,dimension(  MAXEQNH  ),intent(in)  :: ZsolH
  VTYPE,dimension(  MAXEQNH,3),intent(in)  :: ZgradH
  VTYPE,dimension(3,MAXEQNE  ),intent(in)  :: ZsolE
  VTYPE,dimension(3,MAXEQNE  ),intent(in)  :: ZcurlE
  VTYPE,dimension(3,MAXEQNV  ),intent(in)  :: ZsolV
  VTYPE,dimension(  MAXEQNV  ),intent(in)  :: ZdivV
  VTYPE,dimension(  MAXEQNQ  ),intent(in)  :: ZsolQ
  real*8,                      intent(out) :: Val
!
  VTYPE,dimension(3)                       :: zval
!---------------------------------------------------------------------------------------
!
  select case (IEXACT_DISP)
! Approximation solution
  case(2)
!
     select case (ITANGENT_DISP)
     case(1)
       zval(1:3) = zsolE(1:3,1)
     case(2)
       zval(1) =   Rn(2)*zsolE(3,1) - Rn(3)*zsolE(2,1)
       zval(2) = - Rn(1)*zsolE(3,1) + Rn(3)*zsolE(1,1)
       zval(3) =   Rn(1)*zsolE(2,1) - Rn(2)*zsolE(1,1)
     endselect
!
     select case (ICHOOSE_DISP)
     case (1,2,3)
       if (ICOMPLEX_DISP.eq.1) then
         Val = real(zval(ICHOOSE_DISP))
       else
         Val = aimag(zval(ICHOOSE_DISP))
       endif
     case (4)
       Val = sqrt(abs(zval(1)**2 + zval(2)**2 + zval(3)**2)) 
     end select
! 
  end select
!  
!
end subroutine soldis
!
!
!
!---------------------------------------------------------------------------------------
!> Purpose : show the quantities to display 
!---------------------------------------------------------------------------------------
subroutine soldis_select
  use em
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
  write(*,*) 'SET COMPONENT TYPE         1) FIELD  2) TANGENT COMP.'
  read(*,*) iprev
  if (iprev.eq.0) then
     call disp_soldis(NSTD_OUT); return
  endif
  ITANGENT_DISP = iprev
!
  write(*,*) 'SET VARIABLE               1) Ex     2) Ey    3) Ez   4) ABS(E)'
  read(*,*) iprev
  if (iprev.eq.0) then
     call disp_soldis(NSTD_OUT); return
  endif
  ICHOOSE_DISP = iprev
!
  write(*,*) 'SET REAL OR COMPLEX        1) Real   2) Imag'
  read(*,*) iprev
  if (iprev.eq.0) then
     call disp_soldis(NSTD_OUT); return
  endif
  ICOMPLEX_DISP = iprev
!
  call disp_soldis(NSTD_OUT)
!  
end subroutine soldis_select

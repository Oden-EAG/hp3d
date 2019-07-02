!----------------------------------------------------------------------
!> Purpose : display user-defined quantity
!! @param[in] Mdle   - element (middle node) number
!! @param[in] Xi     - master element coordinates
!! @param[in] X     - physical coordinates
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
#include "implicit_none.h"
subroutine soldis(Mdle,Xi,X,Rn, &
     ZsolH,ZgradH, &
     ZsolE,ZcurlE, &
     ZsolV,ZdivV, &
     ZsolQ, &
     Val)
  !
  use thermo_elast
  use data_structure3D
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

  VTYPE :: &
       zvalH(    MAXEQNH    ), &
       zdvalH(   MAXEQNH,3  ), &
       zd2valH(  MAXEQNH,3,3), &
       zvalE(  3,MAXEQNE    ), &
       zdvalE( 3,MAXEQNE,3  ), &
       zd2valE(3,MAXEQNE,3,3), &
       zvalV(  3,MAXEQNV    ), &
       zdvalV( 3,MAXEQNV,3  ), &
       zd2valV(3,MAXEQNV,3,3), &
       zvalQ(    MAXEQNQ    ), &
       zdvalQ(   MAXEQNQ,3  ), &
       zd2valQ(  MAXEQNQ,3,3)

  integer :: iprint
  !----------------------------------------------------------------------
  !
  iprint=0
  !
  Val = ZsolH(ICHOOSE_DISP)
  !
end subroutine soldis

subroutine soldis_select
  use thermo_elast
  implicit none
  integer :: iprev
  !----------------------------------------------------------------------
  write(*,*) 'USE PREVIOUS SELECTION (0/1) '
  read(*,*) iprev
  
  if (iprev.eq.1) then
     call disp_soldis(NSTD_OUT)
     return
  endif

  write(*,*) 'SET DISP COMPONENT         1) u_x   2) u_y   3) u_z   4) temp'
  read(*,*) iprev
  if (iprev.eq.0) then
     call disp_soldis(NSTD_OUT); return
  endif
  ICHOOSE_DISP = iprev
  
end subroutine soldis_select


!----------------------------------------------------------------------------------
!> Purpose : provide dirichlet boundary condition data
!!
!! @param[in]  Mdle  - an element (middle node) number
!! @param[in]  X     - physical coordinates of a point
!! @param[in]  Icase - node case
!! @param[in]  Bcond - node BC flags
!!
!! @param[out] ValH  - value of the H1 solution
!! @param[out] DvalH - H1 corresponding first derivatives
!! @param[out] ValE  - value of the H(curl) solution
!! @param[out] DvalE - H(curl) corresponding first derivatives
!! @param[out] ValV  - value of the H(div) solution
!! @param[out] DvalV - H(div) corresponding first derivatives
!----------------------------------------------------------------------------------
!
#include "typedefs.h"
!
      subroutine dirichlet(Mdle,X,Icase,Bcond, ValH,DvalH,ValE,DvalE,ValV,DvalV)
!
      use control    , only : NEXACT,GEOM_TOL
      use parameters , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ, ZERO
!
      implicit none
!
      integer,                       intent(in)  :: Mdle
      real(8),                       intent(in)  :: X(3)
      integer,                       intent(in)  :: Icase,Bcond
!
      VTYPE,dimension(  MAXEQNH    ),intent(out) ::   ValH
      VTYPE,dimension(  MAXEQNH,3  ),intent(out) ::  DvalH
      VTYPE,dimension(  MAXEQNH,3,3)             :: d2valH
      VTYPE,dimension(3,MAXEQNE    ),intent(out) ::   ValE
      VTYPE,dimension(3,MAXEQNE,3  ),intent(out) ::  DvalE
      VTYPE,dimension(3,MAXEQNE,3,3)             :: d2valE
      VTYPE,dimension(3,MAXEQNV    ),intent(out) ::   ValV
      VTYPE,dimension(3,MAXEQNV,3  ),intent(out) ::  DvalV
      VTYPE,dimension(3,MAXEQNV,3,3)             :: d2valV
      VTYPE,dimension(  MAXEQNQ    )             ::   valQ
      VTYPE,dimension(  MAXEQNQ,3  )             ::  dvalQ
      VTYPE,dimension(  MAXEQNQ,3,3)             :: d2valQ
!
#if DEBUG_MODE
!  ...printing flag : 0 - silent ; 1 - verbose
      integer :: iprint = 0
#endif
!
!----------------------------------------------------------------------------------
!
!  ...initialize
      ValH = ZERO; DvalH = ZERO
      ValE = ZERO; DvalE = ZERO
      ValV = ZERO; DvalV = ZERO
!
      select case(NEXACT)
!
!  ...exact solution UNKNOWN: 
      case(0)
        write(*,*) 'dirichlet: NEXACT = ',NEXACT,' not supported'
        stop 1
!
!  ...exact solution KNOWN
      case(1,2)
!
!  .....use the exact solution to determine Dirichlet data
        call exact(X,Mdle, ValH,DvalH,d2valH,ValE,DvalE,d2valE,  &
                           ValV,DvalV,d2valV,valQ,dvalQ,d2valQ)
!
      case default
        write(*,7010) NEXACT
 7010   format(' dirichlet: NEXACT = ',i2,'. stop.')
        stop 1
!
!  ...end select NEXACT
      end select
!
#if DEBUG_MODE
      if (iprint.eq.1) then
        write(*,7020) X(1:3), ValH(1:MAXEQNH)
 7020   format(' dirichlet: X,ValH = ',3(e12.5,2x),2x,10(e12.5,2x))
      endif
#endif
!
!
      end subroutine dirichlet

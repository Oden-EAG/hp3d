!----------------------------------------------------------------------------------
!> Purpose : calculate dirichlet boundary condition
!!
!! @param[in]  Mdle  - an element (middle node) number
!! @param[in]  X     - physical coordinates of a point
!! @param[in]  Icase - node case
!!
!! @param[out] ValH  - value of the H1 solution
!! @param[out] DvalH - H1 corresponding first derivatives
!! @param[out] ValE  - value of the H(curl) solution
!! @param[out] DvalE - H(curl) corresponding first derivatives
!! @param[out] ValV  - value of the H(div) solution
!! @param[out] DvalV - H(div) corresponding first derivatives
!----------------------------------------------------------------------------------
!
#include "implicit_none.h"
!
subroutine dirichlet(Mdle,X,Icase, ValH,DvalH, ValE,DvalE,ValV,DvalV)
!
   use control    , only : NEXACT, GEOM_TOL
   use CommonParam
   use LaserParam
   use parameters , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ, ZERO
!
   implicit none
!
   integer,                       intent(in)  :: Mdle
   real*8, dimension(3),          intent(in)  :: X
   integer,                       intent(in)  :: Icase
   VTYPE,dimension(  MAXEQNH    ),intent(out) ::   ValH
   VTYPE,dimension(  MAXEQNH,3  ),intent(out) ::  DvalH
   VTYPE,dimension(  MAXEQNH,3,3)             :: d2valH
   VTYPE,dimension(3,MAXEQNE    ),intent(out) ::   ValE
   VTYPE,dimension(3,MAXEQNE,3  ),intent(out) ::  DvalE
   VTYPE,dimension(3,MAXEQNE,3,3)             :: d2valE
   VTYPE,dimension(3,MAXEQNV    ),intent(out) ::   ValV
   VTYPE,dimension(3,MAXEQNV,3  ),intent(out) ::  DvalV
!
!..exact solution (UNUSED)
   VTYPE,dimension(3,MAXEQNV,3,3)             :: d2valV
   VTYPE,dimension(  MAXEQNQ    )             ::   valQ
   VTYPE,dimension(  MAXEQNQ,3  )             ::  dvalQ
   VTYPE,dimension(  MAXEQNQ,3,3)             :: d2valQ
!
!..printing flag
   integer :: iprint
!
!----------------------------------------------------------------------------------
!
!..printing flag : 0 - silent ; 1 - verbose
   iprint=0
!
!..initialize
   ValH=ZERO ; DvalH=ZERO ; ValE=ZERO ; DvalE=ZERO ; ValV=ZERO ; DvalV=ZERO
!
   select case(NEXACT)
!
!  ...exact solution UNKNOWN: solving homogeneous equation
      case(0)
!     ...Select case for problem number: between Maxwell (signal/pump) and heat
         select case(NO_PROBLEM)
!        ...first check if we are solving the Maxwell (signal/pump case)
            case(3,4)
!           ...Solving Maxwell problem with co-pump configuration of fiber
!           ...no need to check for X(3) = 0 for both signal and pump, since
!           ...both are launched at X(3) = 0
            if(COPUMP.eq.1) then
               if((X(3).lt.GEOM_TOL)) then
!              ...Efield value
                  call exact(X,Mdle, ValH,DvalH,d2valH, ValE,DvalE,d2valE,  &
                             ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
               endif
!        ...Solving Maxwell problem with counter-pump configuration of fiber
!        ...signal launched at X(3) = 0, pump launched at X(3) = ZL
            elseif(COPUMP.eq.0) then
               select case(NO_PROBLEM)
!              ...First check for signal and launch at X(3) = 0
                  case(3)
                     if((X(3).lt.GEOM_TOL)) then
!                    ...Efield value
                        call exact(X,Mdle, ValH,DvalH,d2valH, ValE,DvalE,d2valE,  &
                                   ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
                     endif
!              ...Next check for pump and launch at X(3) = ZL
                  case(4)
                     if((X(3).ge.(ZL-GEOM_TOL))) then
!                    ...Efield value
                        call exact(X,Mdle, ValH,DvalH,d2valH, ValE,DvalE,d2valE,  &
                                   ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
                     endif
!           ...end select for NO_PROBLEM launch in counter-pump case
               endselect
!        ...endif for COPUMP check
            endif
!        ...next check if we are solving the Heat system
            case(2)
!           ...do nothing since for heat loop, we have dirichlet BC = 0.d0 on fiber boundary
!           ...should not be running NO_PROBLEM =1 with NEXACT = 0
            case(1)
               write(*,*) 'from dirichlet: cannot run NEXACT = 0 with NO_PROBLEM = 1'
               stop
!           ...end select for NO_PROBLEM for Heat/Maxwell case
         endselect
!     ...exact solution KNOWN
      case(1,2)
!     ...use the exact solution to determine Dirichlet data
         call exact(X,Mdle, ValH,DvalH,d2valH,ValE,DvalE,d2valE,  &
                    ValV,DvalV,d2valV,valQ,dvalQ,d2valQ)
!
      case default
         write(*,1000) NEXACT
 1000    format(' dirichlet: NEXACT = ',i2)
         stop
   endselect
!
      if (iprint == 1) then
         write(*,1001)X(1:3),ValH(1:MAXEQNH)
 1001    format(' dirichlet: X,ValH = ',3(e12.5,2x),2x,10(e12.5,2x))
      endif
!
!
end subroutine dirichlet

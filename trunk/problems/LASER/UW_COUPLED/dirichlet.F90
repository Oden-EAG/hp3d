!-------------------------------------------------------------------------------
! REMARK 1: THIS ROUTINE MUST BE OMP THREAD-SAFE
!           DUE TO OMP PARALLELIZATION OF UPDATE_DDOF->DIRICHLET
!
! REMARK 2: In LASER problem, this dirichlet routine computes valid Dirichlet
!           values only for one particular field, depending on NO_PROBLEM:
!           NO_PROBLEM = 2 --> HEAT DOFS
!           NO_PROBLEM = 3 --> SIGNAL DOFS
!           NO_PROBLEM = 4 --> PUMP DOFS
!           The remaining fields are assumed not to be written by update_Ddof,
!           since the corresponding attributes are deactivated via PHYSAm.
!-------------------------------------------------------------------------------
!> @brief calculate Dirichlet boundary condition
!> last mod: June 2021
!!
!> @param[in]  Mdle  - an element (middle node) number
!> @param[in]  X     - physical coordinates of a point
!> @param[in]  Icase - node case
!!
!> @param[out] ValH  - value of the H1 solution
!> @param[out] DvalH - H1 corresponding first derivatives
!> @param[out] ValE  - value of the H(curl) solution
!> @param[out] DvalE - H(curl) corresponding first derivatives
!> @param[out] ValV  - value of the H(div) solution
!> @param[out] DvalV - H(div) corresponding first derivatives
!-------------------------------------------------------------------------------
!
#include "typedefs.h"
!
subroutine dirichlet(Mdle,X,Icase, ValH,DvalH,ValE,DvalE,ValV,DvalV)
!
   use control    , only : NEXACT,GEOM_TOL
   use commonParam
   use laserParam
   use parameters , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ, ZERO
!
   implicit none
!
   integer,                       intent(in)  :: Mdle
   real(8),                       intent(in)  :: X(3)
   integer,                       intent(in)  :: Icase
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
#if HP3D_DEBUG
   integer :: iprint
   iprint=0
#endif
!
!-------------------------------------------------------------------------------
!
!..initialize
   ValH = ZERO; DvalH = ZERO
   ValE = ZERO; DvalE = ZERO
   ValV = ZERO; DvalV = ZERO
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
!           ...both signal and pump are launched at X(3) = 0
            if(COPUMP.eq.1) then
               if((X(3).lt.GEOM_TOL)) then
!              ...E-trc value
                  call exact(X,Mdle, ValH,DvalH,d2valH,ValE,DvalE,d2valE,  &
                                     ValV,DvalV,d2valV,valQ,dvalQ,d2valQ)
               endif
!
!           ...Solving Maxwell problem with counter-pump configuration of fiber
!           ...signal launched at X(3) = 0, pump launched at X(3) = ZL
            elseif(COPUMP.eq.0) then
               select case(NO_PROBLEM)
!              ...First check for signal and launch at X(3) = 0
                  case(3)
                     if((X(3).lt.GEOM_TOL)) then
!                    ...E-trc value
                        call exact(X,Mdle, ValH,DvalH,d2valH,ValE,DvalE,d2valE,  &
                                           ValV,DvalV,d2valV,valQ,dvalQ,d2valQ)
                     endif
!              ...Next check for pump and launch at X(3) = ZL
                  case(4)
                     if((X(3).gt.(ZL-GEOM_TOL))) then
!                    ...E-trc value
                        call exact(X,Mdle, ValH,DvalH,d2valH,ValE,DvalE,d2valE,  &
                                           ValV,DvalV,d2valV,valQ,dvalQ,d2valQ)
                     endif
!           ...end select NO_PROBLEM
               end select
!        ...endif COPUMP
            endif
!        ...check if we are solving the Heat problem
            case(2)
!           ...do nothing since for heat loop, we have dirichlet BC = 0.d0 on fiber boundary
!              (i.e., homogeneous BC either for the temperature or for the heat flux)
!           ...should not be running linear heat problem NO_PROBLEM = 1 with NEXACT = 0
!           ...we could but we should then set non-zero BC, IC, or source term
            case(1)
               write(*,*) 'dirichlet: cannot run NEXACT = 0 with linear ', &
                                     'heat problem (NO_PROBLEM = 1). stop.'
               stop
!     ...end select NO_PROBLEM
         end select
!  ...exact solution KNOWN
      case(1,2)
!     ...use the exact solution to determine Dirichlet data
         call exact(X,Mdle, ValH,DvalH,d2valH,ValE,DvalE,d2valE,  &
                            ValV,DvalV,d2valV,valQ,dvalQ,d2valQ)
!
      case default
         write(*,1000) NEXACT
 1000    format(' dirichlet: NEXACT = ',i2,'. stop.')
         stop
!
!..end select NEXACT
   end select
!
#if HP3D_DEBUG
   if (iprint == 1) then
      write(*,1001)X(1:3),ValH(1:MAXEQNH)
 1001 format(' dirichlet: X,ValH = ',3(e12.5,2x),2x,10(e12.5,2x))
   endif
#endif
!
!
end subroutine dirichlet

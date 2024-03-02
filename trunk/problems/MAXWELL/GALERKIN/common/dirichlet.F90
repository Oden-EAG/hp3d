!----------------------------------------------------------------------
!
!     routine name      - dirichlet
!
!----------------------------------------------------------------------
!
!     latest revision:  - July 2019
!
!> @brief         - return dirichlet data at a point
!
!     arguments:
!
!     in:
!             Mdle      - middle node number
!             X         - a point in physical space
!             Icase     - node case (specifies what variables are supported)
!     out:
!             ValH      - value of the H1 solution
!             DvalH     - corresponding first derivatives
!             ValE      - value of the H(curl) solution
!             DvalE     - corresponding first derivatives
!             ValV      - value of the H(div) solution
!             DvalV     - corresponding first derivatives
!
!----------------------------------------------------------------------
#include "typedefs.h"
subroutine dirichlet(Mdle,X,Icase, ValH,DvalH,ValE,DvalE,ValV,DvalV)
!
   use control    , only : NEXACT
   use parameters , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ, ZERO
   use common_prob_data
!
   implicit none
!
   real(8), intent(in)  :: X(3)
   integer, intent(in)  :: Icase,Mdle
!..exact solution
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
!--------------------------------------------------------------------
!
!..initialize
   ValH = ZERO; DvalH = ZERO
   ValE = ZERO; DvalE = ZERO
   ValV = ZERO; DvalV = ZERO
!
!  NEXACT:
!  0: unknown exact solution
!  1:   known exact solution
   select case(NEXACT)
      case(0)
         write(*,*) 'dirichlet: missing data; unknown exact solution. stop.'
         stop
      case(1,2)
!     ...use the exact solution to determine Dirichlet data
         call exact(X,Icase, ValH,DvalH,d2valH,ValE,DvalE,d2valE, &
                             ValV,DvalV,d2valV,valQ,dvalQ,d2valQ)
      case default
         write(*,*) 'dirichlet: unknown exact solution flag. stop.'
         stop
   end select
!
end subroutine dirichlet

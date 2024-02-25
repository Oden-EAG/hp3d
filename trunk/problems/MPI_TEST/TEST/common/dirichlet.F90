!
#include "typedefs.h"
!
!----------------------------------------------------------------------
!
!     routine name      - dirichlet
!
!----------------------------------------------------------------------
!
!     latest revision:  - July 17
!
!     purpose:          - return dirichlet data at a point
!
!     arguments:
!
!     in:
!             X         - a point in physical space
!             Icase     - node case (specifies what variables are supported)
!     out:
!             ValH      - value of the H1 solution
!             DvalH     - corresponding first derivatives
!             DvalE     - value of the H(curl) solution
!             DdvalE    - corresponding first derivatives
!             DvalV     - value of the H(div) solution
!             DdvalV    - corresponding first derivatives
!
!----------------------------------------------------------------------
!
   subroutine dirichlet(Mdle,X,Icase, ValH,DvalH,ValE,DvalE,ValV,DvalV)
!
   use parameters, only : MAXEQNH,MAXEQNE,MAXEQNV,ZERO
!
   implicit none
!
   real(8), intent(in)  :: X(3)
   integer, intent(in)  :: Icase,Mdle
!
   VTYPE, dimension(  MAXEQNH    ) ::   ValH
   VTYPE, dimension(  MAXEQNH,3  ) ::  DvalH
   VTYPE, dimension(3,MAXEQNE    ) ::   ValE
   VTYPE, dimension(3,MAXEQNE,3  ) ::  DvalE
   VTYPE, dimension(3,MAXEQNV    ) ::   ValV
   VTYPE, dimension(3,MAXEQNV,3  ) ::  DvalV
!
!--------------------------------------------------------------------
!
!..initialize
   ValH = ZERO; DvalH = ZERO
   ValE = ZERO; DvalE = ZERO
   ValV = ZERO; DvalV = ZERO
!
   end subroutine dirichlet

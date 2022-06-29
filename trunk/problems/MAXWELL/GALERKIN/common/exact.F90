!------------------------------------------------------------------------------
!
!     routine name      - exact
!
!------------------------------------------------------------------------------
!
!     latest revision:  - Oct 2021
!
!     purpose:          - return exact solution value at a point X
!
!     arguments:
!        in:
!             X         - a point in physical space
!             Icase     - node case (specifies what variables are supported)
!        out:
!             ValH      - value of the H1 solution
!             DvalH     - corresponding first derivatives
!             D2valH    - corresponding second derivatives
!             DvalE     - value of the H(curl) solution
!             DdvalE    - corresponding first derivatives
!             Dd2valE   - corresponding second derivatives
!             DvalV     - value of the H(div) solution
!             DdvalV    - corresponding first derivatives
!             Dd2valV   - corresponding second derivatives
!             DvalQ     - value of the L2 solution
!             DdvalQ    - corresponding first derivatives
!             Dd2valQ   - corresponding second derivatives
!
!------------------------------------------------------------------------------
#include "typedefs.h"
subroutine exact(X,Icase, ValH,DvalH,D2valH, &
                          ValE,DvalE,D2valE, &
                          ValV,DvalV,D2valV, &
                          ValQ,DvalQ,D2valQ)
!
   use data_structure3D
!
   implicit none
!
!------------------------------------------------------------------------------
!
   real(8), intent(in)  :: X(3)
   integer, intent(in)  :: Icase
!
   VTYPE,dimension(  MAXEQNH    ), intent(out) ::   ValH
   VTYPE,dimension(  MAXEQNH,3  ), intent(out) ::  DvalH
   VTYPE,dimension(  MAXEQNH,3,3), intent(out) :: D2valH
   VTYPE,dimension(3,MAXEQNE    ), intent(out) ::   ValE
   VTYPE,dimension(3,MAXEQNE,3  ), intent(out) ::  DvalE
   VTYPE,dimension(3,MAXEQNE,3,3), intent(out) :: D2valE
   VTYPE,dimension(3,MAXEQNV    ), intent(out) ::   ValV
   VTYPE,dimension(3,MAXEQNV,3  ), intent(out) ::  DvalV
   VTYPE,dimension(3,MAXEQNV,3,3), intent(out) :: D2valV
   VTYPE,dimension(  MAXEQNQ    ), intent(out) ::   ValQ
   VTYPE,dimension(  MAXEQNQ,3  ), intent(out) ::  DvalQ
   VTYPE,dimension(  MAXEQNQ,3,3), intent(out) :: D2valQ
!
!------------------------------------------------------------------------------
!
!..initialize exact solution
   ValH = ZERO ; DvalH = ZERO ; D2valH = ZERO
   ValE = ZERO ; DvalE = ZERO ; D2valE = ZERO
   ValV = ZERO ; DvalV = ZERO ; D2valV = ZERO
   ValQ = ZERO ; DvalQ = ZERO ; D2valQ = ZERO
!
!..Field values (polynomial solution with homogeneous tangential BC field)
   ValE(1,1) =  X(2)*(1.d0-X(2))*X(3)*(1.d0-X(3))
   ValE(2,1) =  X(2)*X(1)*(1.d0-X(1))*X(3)*(1.d0-X(3))
   ValE(3,1) =  X(1)*(1.d0-X(1))*X(2)*(1.d0-X(2))
!
!..1st order derivatives
!  grad(E_x)
   DvalE(1,1,1) = ZERO !...derivative wrt x
   DvalE(1,1,2) = (1.d0-2.d0*X(2))*X(3)*(1.d0-X(3)) !...derivative wrt y
   DvalE(1,1,3) = X(2)*(1.d0-X(2))*(1.d0-2.d0*X(3)) !...derivative wrt z
!  grad(E_y)
   DvalE(2,1,1) = X(2)*(1.d0-2.d0*X(1))*X(3)*(1.d0-X(3))
   DvalE(2,1,2) = X(1)*(1.d0-X(1))*X(3)*(1.d0-X(3))
   DvalE(2,1,3) = X(2)*X(1)*(1.d0-X(1))*(1.d0-2.d0*X(3))
!  grad(E_z)
   DvalE(3,1,1) = (1.d0-2.d0*X(1))*X(2)*(1.d0-X(2))
   DvalE(3,1,2) = X(1)*(1.d0-X(1))*(1.d0-2.d0*X(2))
   DvalE(3,1,3) = ZERO
!
!..2nd order derivatives
!  grad grad E_x
   D2valE(1,1,1,1) = ZERO !...derivative wrt xx
   D2valE(1,1,1,2) = ZERO !...derivative wrt xy
   D2valE(1,1,1,3) = ZERO !...derivative wrt xz
   D2valE(1,1,2,1) = D2valE(1,1,1,2) !...derivative wrt yx
   D2valE(1,1,2,2) = -2.d0*X(3)*(1.d0-X(3))  !...derivative wrt yy
   D2valE(1,1,2,3) = (1.d0-2.d0*X(2))*(1.d0-2.d0*X(3)) !...derivative wrt yz
   D2valE(1,1,3,1) = D2valE(1,1,1,3) !...derivative wrt zx
   D2valE(1,1,3,2) = D2valE(1,1,2,3) !...derivative wrt zy
   D2valE(1,1,3,3) = X(2)*(1.d0-X(2))*(-2.d0) !...derivative wrt zz
!  grad grad E_y
   D2valE(2,1,1,1) = X(2)*(-2.d0)*X(3)*(1.d0-X(3))
   D2valE(2,1,1,2) = (1.d0-2.d0*X(1))*X(3)*(1.d0-X(3))
   D2valE(2,1,1,3) = X(2)*(1.d0-2.d0*X(1))*(1.d0-2.d0*X(3))
   D2valE(2,1,2,1) = D2valE(2,1,1,2)
   D2valE(2,1,2,2) = ZERO
   D2valE(2,1,2,3) = X(1)*(1.d0-X(1))*(1.d0-2.d0*X(3))
   D2valE(2,1,3,1) = D2valE(2,1,1,3)
   D2valE(2,1,3,2) = D2valE(2,1,2,3)
   D2valE(2,1,3,3) = X(2)*X(1)*(1.d0-X(1))*(-2.d0)
!  grad grad E_z
   D2valE(3,1,1,1) = (-2.d0)*X(2)*(1.d0-X(2))
   D2valE(3,1,1,2) = (1.d0-2.d0*X(1))*(1.d0-2.d0*X(2))
   D2valE(3,1,1,3) = ZERO
   D2valE(3,1,2,1) = D2valE(3,1,1,2)
   D2valE(3,1,2,2) = X(1)*(1.d0-X(1))*(-2.d0)
   D2valE(3,1,2,3) = ZERO
   D2valE(3,1,3,1) = D2valE(3,1,1,3)
   D2valE(3,1,3,2) = D2valE(3,1,2,3)
   D2valE(3,1,3,3) = ZERO
!
!
end subroutine exact

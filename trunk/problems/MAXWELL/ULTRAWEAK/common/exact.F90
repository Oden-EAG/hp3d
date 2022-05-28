!----------------------------------------------------------------------
!
!     routine name      - exact
!
!----------------------------------------------------------------------
!
!     latest revision:  - May 2020
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
!----------------------------------------------------------------------
#include "typedefs.h"
subroutine exact(X,Icase, ValH,DvalH,D2valH, &
                          ValE,DvalE,D2valE, &
                          ValV,DvalV,D2valV, &
                          ValQ,DvalQ,D2valQ)
!
   use data_structure3D
   use parameters       , only : ZIMG
   use common_prob_data , only : OMEGA, MU
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
   VTYPE                   :: E
   VTYPE, dimension(3)     :: dE
   VTYPE, dimension(3,3)   :: d2E
!
   integer :: icomp

!------------------------------------------------------------------------------
!
!..initialize exact solution
   ValH = ZERO ; DvalH = ZERO ; D2valH = ZERO
   ValE = ZERO ; DvalE = ZERO ; D2valE = ZERO
   ValV = ZERO ; DvalV = ZERO ; D2valV = ZERO
   ValQ = ZERO ; DvalQ = ZERO ; D2valQ = ZERO
!
!..get exact solution
   call solution(X, E, dE, d2E)
!
!..non zero component
   icomp = 1

!..1st H(curl) (Electric field)

!..value
   ValE(icomp,1) =  E ! first component non zero
!   
!..1st order derivatives 
   DvalE(icomp,1,1) = dE(1)  !...derivative wrt x
   DvalE(icomp,1,2) = dE(2)  !...derivative wrt y
   DvalE(icomp,1,3) = dE(3)  !...derivative wrt z
!
!..2nd order derivatives
   D2valE(icomp,1,1,1) = d2E(1,1)  !...derivative wrt xx
   D2valE(icomp,1,1,2) = d2E(1,2)  !...derivative wrt xy
   D2valE(icomp,1,1,3) = d2E(1,3)  !...derivative wrt xz
   D2valE(icomp,1,2,1) = d2E(2,1)  !...derivative wrt yx
   D2valE(icomp,1,2,2) = d2E(2,2)  !...derivative wrt yy
   D2valE(icomp,1,2,3) = d2E(2,3)  !...derivative wrt yz
   D2valE(icomp,1,3,1) = d2E(3,1)  !...derivative wrt zx
   D2valE(icomp,1,3,2) = d2E(3,2)  !...derivative wrt zy
   D2valE(icomp,1,3,3) = d2E(3,3)  !...derivative wrt zz
!

!..2nd H(curl) (Magnetic field)  H = - curl(E) / (i \omega \mu)
!
!..Value
   ValE(1,2)   =  DvalE(3,1,2) - DvalE(2,1,3)
   ValE(2,2)   =  DvalE(1,1,3) - DvalE(3,1,1)
   ValE(3,2)   =  DvalE(2,1,1) - DvalE(1,1,2)
   ValE(1:3,2) = -ValE(1:3,2)/(ZIMG*OMEGA*MU)   ! Magnetic field

!..1st order derivatives
   DvalE(1,2,1) = D2valE(3,1,2,1) - D2valE(2,1,3,1)
   DvalE(1,2,2) = D2valE(3,1,2,2) - D2valE(2,1,3,2)
   DvalE(1,2,3) = D2valE(3,1,2,3) - D2valE(2,1,3,3)
!
   DvalE(2,2,1) = D2valE(1,1,3,1) - D2valE(3,1,1,1)
   DvalE(2,2,2) = D2valE(1,1,3,2) - D2valE(3,1,1,2)
   DvalE(2,2,3) = D2valE(1,1,3,3) - D2valE(3,1,1,3)
!
   DvalE(3,2,1) = D2valE(2,1,1,1) - D2valE(1,1,2,1)
   DvalE(3,2,2) = D2valE(2,1,1,2) - D2valE(1,1,2,2)
   DvalE(3,2,3) = D2valE(2,1,1,3) - D2valE(1,1,2,3)
!
   DvalE(1:3,2,1:3) = -DvalE(1:3,2,1:3)/(ZIMG*OMEGA*MU)

!..2nd order derivatives (not needed)

!..L2 components, derivatives not needed
   ValQ(1:3) = ValE(1:3,1)
   ValQ(4:6) = ValE(1:3,2)

!
end subroutine exact

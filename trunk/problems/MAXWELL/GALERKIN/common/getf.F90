!-----------------------------------------------------------------------------------
!
!     routine name      - getf
!
!-----------------------------------------------------------------------------------
!
!     latest revision:  - Oct 2021
!
!> @brief         - return source term value at a point X
!
!     arguments:
!
!     in:
!             Mdle      - element (middle node) number
!             X         - a point in physical space
!     out:
!             zJ        - Value of source term (impressed current) at the point X
!
!-----------------------------------------------------------------------------------
#include "typedefs.h"
subroutine getf(Mdle,X, zJ)
!
   use data_structure3D
   use control         , only : NEXACT
   use parameters      , only : ZERO
   use common_prob_data
!
   implicit none
!-----------------------------------------------------------------------------------
   integer   , intent(in)  :: Mdle
   real(8)   , intent(in)  :: X(3)
   complex(8), intent(out) :: zJ(3)
!-----------------------------------------------------------------------------------
   VTYPE,dimension(  MAXEQNH    ) ::   valH
   VTYPE,dimension(  MAXEQNH,3  ) ::  dvalH
   VTYPE,dimension(  MAXEQNH,3,3) :: d2valH
   VTYPE,dimension(3,MAXEQNE    ) ::   valE
   VTYPE,dimension(3,MAXEQNE,3  ) ::  dvalE
   VTYPE,dimension(3,MAXEQNE,3,3) :: d2valE
   VTYPE,dimension(3,MAXEQNV    ) ::   valV
   VTYPE,dimension(3,MAXEQNV,3  ) ::  dvalV
   VTYPE,dimension(3,MAXEQNV,3,3) :: d2valV
   VTYPE,dimension(  MAXEQNQ    ) ::   valQ
   VTYPE,dimension(  MAXEQNQ,3  ) ::  dvalQ
   VTYPE,dimension(  MAXEQNQ,3,3) :: d2valQ
!
   integer    :: icase
   complex(8) :: za(3),zb(3)
   complex(8) :: CE_x(3), CE_y(3), CE_z(3), CC_E(3)
!
!-----------------------------------------------------------------------------------
!
   zJ = ZERO
!
!..make sure exact solution is available
   if (NEXACT==0) then
      write(*,*) 'getf: source term cannot be computed;', &
                     '  exact solution is unknown. stop.'
      stop
   endif
!
!..compute exact solution
   icase = 1
   call exact(X,icase, valH,dvalH,d2valH,valE,dvalE,d2valE, &
                       valV,dvalV,d2valV,valQ,dvalQ,d2valQ)
!
!
!..1st order derivatives of curl(E)
!  grad(curl(E)_x)
   CE_x(1) = D2valE(3,1,2,1) - D2valE(2,1,3,1) ! (curl E)_1,x
   CE_x(2) = D2valE(3,1,2,2) - D2valE(2,1,3,2) ! (curl E)_1,y
   CE_x(3) = D2valE(3,1,2,3) - D2valE(2,1,3,3) ! (curl E)_1,z
!  grad(curl(E)_y)
   CE_y(1) = D2valE(1,1,3,1) - D2valE(3,1,1,1) ! (curl E)_2,x
   CE_y(2) = D2valE(1,1,3,2) - D2valE(3,1,1,2) ! (curl E)_2,y
   CE_y(3) = D2valE(1,1,3,3) - D2valE(3,1,1,3) ! (curl E)_2,z
!  grad(curl(E)_z)
   CE_z(1) = D2valE(2,1,1,1) - D2valE(1,1,2,1) ! (curl E)_3,x
   CE_z(2) = D2valE(2,1,1,2) - D2valE(1,1,2,2) ! (curl E)_3,y
   CE_z(3) = D2valE(2,1,1,3) - D2valE(1,1,2,3) ! (curl E)_3,z
!
!..curl curl E
   CC_E(1) = CE_z(2) - CE_y(3)
   CC_E(2) = CE_x(3) - CE_z(1)
   CC_E(3) = CE_y(1) - CE_x(2)
!
!..compute impressed current from exact solution
!  -iωJ = curl((1/μ)curl E) - (ω^2ε-iωσ)E
   za = CC_E / MU
   zb = (OMEGA*OMEGA*EPSILON - ZI*OMEGA*SIGMA) * ValE(1:3,1)
   zJ = (za-zb) / (-ZI * OMEGA)
!
end subroutine getf

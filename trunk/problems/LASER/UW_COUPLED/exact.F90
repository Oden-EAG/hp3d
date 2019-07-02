!------------------------------------------------------------------------------
!> Purpose : exact (manufactured, no miracle!) solution
!!
!> @param[in]  X      - a point in physical space
!> @param[in]  Mdle   - element (middle node) number
!> @param[out] ValH   - value of the H1 solution
!> @param[out] DvalH  - corresponding first derivatives
!> @param[out] D2valH - corresponding second derivatives
!> @param[out] ValE   - value of the H(curl) solution
!> @param[out] DvalE  - corresponding first derivatives
!> @param[out] D2valE - corresponding second derivatives
!> @param[out] ValV   - value of the H(div) solution
!> @param[out] DvalV  - corresponding first derivatives
!> @param[out] D2valV - corresponding second derivatives
!> @param[out] ValQ   - value of the H(div) solution
!> @param[out] DvalQ  - corresponding first derivatives
!> @param[out] D2valQ - corresponding second derivatives
!------------------------------------------------------------------------------
!
#include "implicit_none.h"
!
subroutine exact(Xp,Mdle, ValH,DvalH,D2valH, ValE,DvalE,D2valE, &
                          ValV,DvalV,D2valV, ValQ,DvalQ,D2valQ)
!
  use data_structure3D
  use CommonParam
  use LaserParam
!
  implicit none
  real*8,dimension(3),            intent(in)  :: Xp
  integer                       , intent(in)  :: Mdle
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
!     Space for temporary solutions
!
  VTYPE                    :: u,E
  VTYPE,dimension(3)       :: gradu,dE
  VTYPE,dimension(3,3)     :: grad2u,d2E
  integer                  :: icomp

!
! initialize exact solution
  ValE=ZERO ; DvalE=ZERO ; D2valE=ZERO
  ValH=ZERO ; DvalH=ZERO ; D2valH=ZERO
  ValV=ZERO ; DvalV=ZERO ; D2valV=ZERO
  ValQ=ZERO ; DvalQ=ZERO ; D2valQ=ZERO
! set icomp
  icomp = ICOMP_EXACT

! initialize variables
  u = ZERO; E = ZERO
  gradu = ZERO; dE = ZERO
  grad2u = ZERO; d2E = ZERO
  call mfd_solutions(Xp, u,gradu,grad2u)
  call mfd_solutions(Xp, E,dE,d2E)

!
  ValH(1) = u                               ! H1 variable
  DvalH(1,1:3) = gradu(1:3)                 ! 1st der
  D2valH(1,1:3,1:3) = grad2u(1:3,1:3)       ! 2nd der
!
  ValV(1:3,1) = gradu(1:3)                  ! Hdiv
  DvalV(1:3,1,1:3) = grad2u(1:3,1:3)        ! 1st der

! Efield.....value
  ValE(icomp,1    ) =   E
  ValE(icomp,3    ) =   E
!
! 1st order derivatives
  DvalE(icomp,1,1  ) =  dE(1)
  DvalE(icomp,1,2  ) =  dE(2)
  DvalE(icomp,1,3  ) =  dE(3)

  DvalE(icomp,3,1  ) =  dE(1)
  DvalE(icomp,3,2  ) =  dE(2)
  DvalE(icomp,3,3  ) =  dE(3)
!
! 2nd order derivatives
  D2valE(icomp,1,1,1) = d2E(1,1)
  D2valE(icomp,1,1,2) = d2E(1,2)
  D2valE(icomp,1,1,3) = d2E(1,3)
  D2valE(icomp,1,2,1) = d2E(2,1)
  D2valE(icomp,1,2,2) = d2E(2,2)
  D2valE(icomp,1,2,3) = d2E(2,3)
  D2valE(icomp,1,3,1) = d2E(3,1)
  D2valE(icomp,1,3,2) = d2E(3,2)
  D2valE(icomp,1,3,3) = d2E(3,3)

  D2valE(icomp,3,1,1) = d2E(1,1)
  D2valE(icomp,3,1,2) = d2E(1,2)
  D2valE(icomp,3,1,3) = d2E(1,3)
  D2valE(icomp,3,2,1) = d2E(2,1)
  D2valE(icomp,3,2,2) = d2E(2,2)
  D2valE(icomp,3,2,3) = d2E(2,3)
  D2valE(icomp,3,3,1) = d2E(3,1)
  D2valE(icomp,3,3,2) = d2E(3,2)
  D2valE(icomp,3,3,3) = d2E(3,3)
!
! 2nd H(curl) ATTRIBUTE = curl of the first attribute/-i omega \mu
! value
  ValE(1,2) = DvalE(3,1,2) - DvalE(2,1,3)
  ValE(2,2) = DvalE(1,1,3) - DvalE(3,1,1)
  ValE(3,2) = DvalE(2,1,1) - DvalE(1,1,2)
  ValE(1:3,2) = ValE(1:3,2)/(-ZI*OMEGA*OMEGA_RATIO_SIGNAL*MU)

  ValE(1,4) = DvalE(3,3,2) - DvalE(2,3,3)
  ValE(2,4) = DvalE(1,3,3) - DvalE(3,3,1)
  ValE(3,4) = DvalE(2,3,1) - DvalE(1,3,2)
  ValE(1:3,4) = ValE(1:3,4)/(-ZI*OMEGA*OMEGA_RATIO_PUMP*MU)
!
! 1st order derivatives
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
  DvalE(1:3,2,1:3) = DvalE(1:3,2,1:3)/(-ZI*OMEGA*OMEGA_RATIO_SIGNAL*MU)

! For second set of H-curl traces
  DvalE(1,4,1) = D2valE(3,3,2,1) - D2valE(2,3,3,1)
  DvalE(1,4,2) = D2valE(3,3,2,2) - D2valE(2,3,3,2)
  DvalE(1,4,3) = D2valE(3,3,2,3) - D2valE(2,3,3,3)
!
  DvalE(2,4,1) = D2valE(1,3,3,1) - D2valE(3,3,1,1)
  DvalE(2,4,2) = D2valE(1,3,3,2) - D2valE(3,3,1,2)
  DvalE(2,4,3) = D2valE(1,3,3,3) - D2valE(3,3,1,3)
!
  DvalE(3,4,1) = D2valE(2,3,1,1) - D2valE(1,3,2,1)
  DvalE(3,4,2) = D2valE(2,3,1,2) - D2valE(1,3,2,2)
  DvalE(3,4,3) = D2valE(2,3,1,3) - D2valE(1,3,2,3)
!
  DvalE(1:3,4,1:3) = DvalE(1:3,4,1:3)/(-ZI*OMEGA*OMEGA_RATIO_PUMP*MU)
!
! fake 2nd order derivatives (not needed)
  D2valE(1:3,2,1:3,1:3) = ZERO
  D2valE(1:3,4,1:3,1:3) = ZERO
!
! L2 components, derivatives not needed
  ValQ(1:3) = ValE(1:3,1)
  ValQ(4:6) = ValE(1:3,2)
  ValQ(7:9) = ValE(1:3,3)
  ValQ(10:12) = ValE(1:3,4)
!
end subroutine exact

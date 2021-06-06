!
! REMARK: THIS ROUTINE MUST BE OMP THREAD-SAFE
!         DUE TO OMP PARALLELIZATION OF UPDATE_DDOF->DIRICHLET->EXACT
!------------------------------------------------------------------------------
!> Purpose : exact (manufactured) solution
!> last mod: May 21
!
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
#include "typedefs.h"
!
      subroutine exact(Xp,Mdle, ValH,DvalH,D2valH, ValE,DvalE,D2valE, &
                                ValV,DvalV,D2valV, ValQ,DvalQ,D2valQ)
!
      use data_structure3D
!
      implicit none
      real(8)                       , intent(in)  :: Xp(3)
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
!  ...Locals
      integer :: nsol
      real(8) :: x,y,z
!
!------------------------------------------------------------------------------
!
!  ...initialize exact solution
      ValH=ZERO ; DvalH=ZERO ; D2valH=ZERO
      ValE=ZERO ; DvalE=ZERO ; D2valE=ZERO
      ValV=ZERO ; DvalV=ZERO ; D2valV=ZERO
      ValQ=ZERO ; DvalQ=ZERO ; D2valQ=ZERO
!
      x = Xp(1); y = Xp(2); z = Xp(3)
      nsol=1
      select case(nsol)
!
!  ...a quadratic function
      case(1)
        ValH(1:3) = x**2 + 1.d0
        DvalH(1:3,1) = 2.d0*x
        D2valH(1:3,1,1) = 2.d0
!
!  ...a quadratic function with zero Dirichlet BC's 
      case(2)
        ValH(1:3) = x*(1.d0-x)
        DvalH(1:3,1) = -2.d0*x + 1.d0
        D2valH(1:3,1,1) = -2.d0

      end select
!
      end subroutine exact

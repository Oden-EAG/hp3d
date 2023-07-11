!
#include "typedefs.h"
!
!----------------------------------------------------------------------
!> @brief       Compute exact solution at a point
!!
!> @param[in]   X       - point in physical space
!> @param[in]   Icase   - node case (specifies what variables are supported)
!> @param[out]  ValH    - value of the H1 solution
!> @param[out]  DvalH   - corresponding first derivatives
!> @param[out]  D2valH  - corresponding second derivatives
!> @param[out]  DvalE   - value of the H(curl) solution
!> @param[out]  DdvalE  - corresponding first derivatives
!> @param[out]  Dd2valE - corresponding second derivatives
!> @param[out]  DvalV   - value of the H(div) solution
!> @param[out]  DdvalV  - corresponding first derivatives
!> @param[out]  Dd2valV - corresponding second derivatives
!> @param[out]  DvalQ   - value of the H(div) solution
!> @param[out]  DdvalQ  - corresponding first derivatives
!> @param[out]  Dd2valQ - corresponding second derivatives
!!
!> @date        July 2023
!----------------------------------------------------------------------
   subroutine exact(X,Icase, ValH,DvalH,D2valH, &
                             ValE,DvalE,D2valE, &
                             ValV,DvalV,D2valV, &
                             ValQ,DvalQ,D2valQ)
!
      use data_structure3D
      use common_prob_data_UW, ONLY: OMEGA
!
      implicit none
!
      real(8), intent(in) :: X(3)
      integer, intent(in) :: Icase
!
      VTYPE, intent(out) :: ValH  (  MAXEQNH    )
      VTYPE, intent(out) :: DvalH (  MAXEQNH,3  )
      VTYPE, intent(out) :: D2valH(  MAXEQNH,3,3)
      VTYPE, intent(out) :: ValE  (3,MAXEQNE    )
      VTYPE, intent(out) :: DvalE (3,MAXEQNE,3  )
      VTYPE, intent(out) :: D2valE(3,MAXEQNE,3,3)
      VTYPE, intent(out) :: ValV  (3,MAXEQNV    )
      VTYPE, intent(out) :: DvalV (3,MAXEQNV,3  )
      VTYPE, intent(out) :: D2valV(3,MAXEQNV,3,3)
      VTYPE, intent(out) :: ValQ  (  MAXEQNQ    )
      VTYPE, intent(out) :: DvalQ (  MAXEQNQ,3  )
      VTYPE, intent(out) :: D2valQ(  MAXEQNQ,3,3)
!  
      VTYPE                   :: p
      VTYPE, dimension(3)     :: gradp
      VTYPE, dimension(3,3)   :: grad2p
!
!------------------------------------------------------------------------------
!
!  ...initialize exact solution
      ValH = ZERO ; DvalH = ZERO ; D2valH = ZERO
      ValE = ZERO ; DvalE = ZERO ; D2valE = ZERO
      ValV = ZERO ; DvalV = ZERO ; D2valV = ZERO
      ValQ = ZERO ; DvalQ = ZERO ; D2valQ = ZERO

!  ...get solution for the pressure
      call acoustics_solution(X, p, gradp, grad2p)
!
      ValQ(1)           = p                             ! pressure
      DvalQ(1,  1:3)    = gradp                         ! gradient of the pressure
      D2valQ(1,1:3,1:3) = grad2p                        ! Hessian of the pressure
      ValQ(2:4)         = - gradp  / (ZIMG * OMEGA)     ! velocity
      DvalQ(2:4,1:3)    = - grad2p / (ZIMG * OMEGA)     ! Hessian of velocity
!
      ValH (1)          = p        !  p \hat      (H1)
      DvalH(1,1:3)      = gradp    !  p \hat      1st der
      D2valH(1,1:3,1:3) = grad2p   !  p \hat      2nd der
!
      ValV(1:3,1)       = -gradp  / (ZIMG * OMEGA)    ! u \hat     (Hdiv)
      DvalV(1:3,1,1:3)  = -grad2p / (ZIMG * OMEGA)    ! u \hat    1st der
! 
   end subroutine exact

!
#include "typedefs.h"
!
!-------------------------------------------------------------------------------
! REMARK 1: THIS ROUTINE MUST BE OMP THREAD-SAFE
!           DUE TO OMP PARALLELIZATION OF UPDATE_DDOF->DIRICHLET->EXACT
!-------------------------------------------------------------------------------
!> Purpose : exact (manufactured) solution
!> last mod: June 2021
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
!-------------------------------------------------------------------------------
   subroutine exact(Xp,Mdle, ValH,DvalH,D2valH, ValE,DvalE,D2valE, &
                             ValV,DvalV,D2valV, ValQ,DvalQ,D2valQ)
!
      use data_structure3D
      use commonParam
      use mpi_param
!
      implicit none
      real(8), intent(in)  :: Xp(3)
      integer, intent(in)  :: Mdle
      VTYPE,   intent(out) :: ValH  (  MAXEQNH    )
      VTYPE,   intent(out) :: DvalH (  MAXEQNH,3  )
      VTYPE,   intent(out) :: D2valH(  MAXEQNH,3,3)
      VTYPE,   intent(out) :: ValE  (3,MAXEQNE    )
      VTYPE,   intent(out) :: DvalE (3,MAXEQNE,3  )
      VTYPE,   intent(out) :: D2valE(3,MAXEQNE,3,3)
      VTYPE,   intent(out) :: ValV  (3,MAXEQNV    )
      VTYPE,   intent(out) :: DvalV (3,MAXEQNV,3  )
      VTYPE,   intent(out) :: D2valV(3,MAXEQNV,3,3)
      VTYPE,   intent(out) :: ValQ  (  MAXEQNQ    )
      VTYPE,   intent(out) :: DvalQ (  MAXEQNQ,3  )
      VTYPE,   intent(out) :: D2valQ(  MAXEQNQ,3,3)
!
!  ...space for temporary solutions
      VTYPE   :: E
      VTYPE   :: dE(3)
      VTYPE   :: d2E(3,3)
      integer :: icomp, idx
!
!------------------------------------------------------------------------------
!
!  ...initialize exact solution
      ValH=ZERO ; DvalH=ZERO ; D2valH=ZERO
      ValE=ZERO ; DvalE=ZERO ; D2valE=ZERO
      ValV=ZERO ; DvalV=ZERO ; D2valV=ZERO
      ValQ=ZERO ; DvalQ=ZERO ; D2valQ=ZERO
!
!  ...auxiliary value
      idx = 1
!
!  ...set icomp (the only non-zero component of the Etrc solution)
      icomp = ICOMP_EXACT
!
!  ...initialize variables
      E = ZERO; dE = ZERO; d2E = ZERO
!
!  ...get Electric field
      call mfd_solutions(Xp, E,dE,d2E)
!
!  ...E-field value
      ValE(icomp,idx) = E ! E-field trace; (icomp,idx+1) is H-field trace
!
!  ...E-field 1st order derivatives
      DvalE(icomp,idx,1:3) = dE(1:3)
!
!  ...E-field 2nd order derivatives
      D2valE(icomp,idx,1:3,1:3) = d2E(1:3,1:3)
!
!  ...2nd H(curl) ATTRIBUTE = curl of the first attribute/-i omega \mu
!     H-field value (H-field trace)
      ValE(1,idx+1) = DvalE(3,idx,2) - DvalE(2,idx,3)
      ValE(2,idx+1) = DvalE(1,idx,3) - DvalE(3,idx,1)
      ValE(3,idx+1) = DvalE(2,idx,1) - DvalE(1,idx,2)
!
      ValE(1:3,idx+1) = ValE(1:3,idx+1)/(-ZI*OMEGA*MU)
!
!  ...H-field 1st order derivatives
      DvalE(1,idx+1,1) = D2valE(3,idx,2,1) - D2valE(2,idx,3,1)
      DvalE(1,idx+1,2) = D2valE(3,idx,2,2) - D2valE(2,idx,3,2)
      DvalE(1,idx+1,3) = D2valE(3,idx,2,3) - D2valE(2,idx,3,3)
!
      DvalE(2,idx+1,1) = D2valE(1,idx,3,1) - D2valE(3,idx,1,1)
      DvalE(2,idx+1,2) = D2valE(1,idx,3,2) - D2valE(3,idx,1,2)
      DvalE(2,idx+1,3) = D2valE(1,idx,3,3) - D2valE(3,idx,1,3)
!
      DvalE(3,idx+1,1) = D2valE(2,idx,1,1) - D2valE(1,idx,2,1)
      DvalE(3,idx+1,2) = D2valE(2,idx,1,2) - D2valE(1,idx,2,2)
      DvalE(3,idx+1,3) = D2valE(2,idx,1,3) - D2valE(1,idx,2,3)
!
      DvalE(1:3,idx+1,1:3) = DvalE(1:3,idx+1,1:3)/(-ZI*OMEGA*MU)
!
!  ...2nd order derivatives (not needed)
!
!  ...L2 components, derivatives not needed
      ValQ(1:3) = ValE(1:3,1) ! E-field
      ValQ(4:6) = ValE(1:3,2) ! H-field
!
   end subroutine exact

!
#include "typedefs.h"
!
!-------------------------------------------------------------------------------
!> @remark     THIS ROUTINE MUST BE OMP THREAD-SAFE
!              DUE TO OMP PARALLELIZATION OF UPDATE_DDOF->DIRICHLET
!!
!> @brief      calculate Dirichlet boundary condition
!> @date       Sep 2023
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
   subroutine dirichlet(Mdle,X,Icase, ValH,DvalH,ValE,DvalE,ValV,DvalV)
!
      use control,         only : NEXACT,GEOM_TOL
      use commonParam
      use parameters,      only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ,ZERO
!
      implicit none
!
      integer, intent(in)  :: Mdle
      real(8), intent(in)  :: X(3)
      integer, intent(in)  :: Icase
!
      VTYPE, intent(out) :: ValH  (  MAXEQNH    )
      VTYPE, intent(out) :: DvalH (  MAXEQNH,3  )
      VTYPE              :: d2valH(  MAXEQNH,3,3)
      VTYPE, intent(out) :: ValE  (3,MAXEQNE    )
      VTYPE, intent(out) :: DvalE (3,MAXEQNE,3  )
      VTYPE              :: d2valE(3,MAXEQNE,3,3)
      VTYPE, intent(out) :: ValV  (3,MAXEQNV    )
      VTYPE, intent(out) :: DvalV (3,MAXEQNV,3  )
      VTYPE              :: d2valV(3,MAXEQNV,3,3)
      VTYPE              :: valQ  (  MAXEQNQ    )
      VTYPE              :: dvalQ (  MAXEQNQ,3  )
      VTYPE              :: d2valQ(  MAXEQNQ,3,3)
!
#if DEBUG_MODE
!  ...printing flag : 0 - silent ; 1 - verbose
      integer :: iprint = 0
#endif
!
!-------------------------------------------------------------------------------
!
!  ...initialize
      ValH = ZERO; DvalH = ZERO
      ValE = ZERO; DvalE = ZERO
      ValV = ZERO; DvalV = ZERO
!
      select case(NEXACT)
!
!     ...exact solution UNKNOWN: solving homogeneous equation
         case(0)
!
            continue
!
!     ...exact solution KNOWN
         case(1,2)
!        ...use the exact solution to determine Dirichlet data
            call exact(X,Mdle, ValH,DvalH,d2valH,ValE,DvalE,d2valE,  &
                            ValV,DvalV,d2valV,valQ,dvalQ,d2valQ)
!
         case default
            write(*,1000) NEXACT
 1000       format(' dirichlet: NEXACT = ',i2,'. stop.')
            stop
!
!  ...end select NEXACT
      end select
!
#if DEBUG_MODE
      if (iprint == 1) then
         write(*,1001)X(1:3),ValH(1:MAXEQNH)
 1001    format(' dirichlet: X,ValH = ',3(e12.5,2x),2x,10(e12.5,2x))
      endif
#endif
!
end subroutine dirichlet


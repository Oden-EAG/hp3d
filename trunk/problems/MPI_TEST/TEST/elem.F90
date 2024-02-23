!----------------------------------------------------------------------
!
!     routine name      - elem
!
!----------------------------------------------------------------------
!
!     latest revision:  - July 2019
!
!     purpose:          - driver for the element routine
!
!     arguments:
!        in:
!             Mdle      - an element middle node number, identified
!                         with the element
!        out:
!             Itest     - index for assembly
!             Itrial    - index for assembly
!
!----------------------------------------------------------------------
!
#include "typedefs.h"
!
subroutine elem(Mdle, Itest,Itrial)
!
   use physics, only: NR_PHYSA
!
   implicit none
!
   integer, intent(in) :: Mdle
   integer, dimension(NR_PHYSA), intent(out) :: Itest,Itrial
!
   Itest = 1; Itrial = 1
!
end subroutine elem

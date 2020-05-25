!
#include "implicit_none.h"
!
!----------------------------------------------------------------------
!
!   routine name       - dot_product
!
!----------------------------------------------------------------------
!
!   latest revision    - May 2020
!
!   purpose            - routine evaluates dot product  of two
!                        vectors in R^3
!
!   arguments :
!     in:
!           A,B        - vectors in R^3
!     out:
!           Prod       - value of the product
!
!----------------------------------------------------------------------
subroutine dot_product(A,B, Prod)
!
   implicit none
!
   real(8), intent(in)  :: A(3),B(3)
   real(8), intent(out) :: Prod
!
   Prod = 0.d0
   do i=1,3
      Prod = Prod + A(i)*B(i)
   enddo
!
end subroutine dot_product
!
!
!--------------------------------------------------------------------
!
!     routine name      - zdot_product
!
!-------------------------------------------------------------------- 
!
!     latest revision:  - May 2020
!
!     purpose:          - compute dot product of real and
!                         (possibly) complex valued vector
!
!---------------------------------------------------------------------
subroutine zdot_product(Rn,Za, Zprod)
!
   implicit none
!
   real(8), intent(in)  :: Rn(3)
   VTYPE  , intent(in)  :: Za(3)
   VTYPE  , intent(out) :: Zprod
!
!---------------------------------------------------------------------
!
   Zprod = Rn(1)*Za(1)+Rn(2)*Za(2)+Rn(3)*Za(3)
!
!
end subroutine zdot_product

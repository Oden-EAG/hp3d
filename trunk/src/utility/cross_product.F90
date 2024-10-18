!
#include "typedefs.h"
!
!--------------------------------------------------------------------
!
!     routine name      - cross_product
!
!     latest revision:  - Mar 2023
!
!> @brief         - compute 3D cross product of real-valued
!                         vectors
!
!--------------------------------------------------------------------
!
   subroutine cross_product(a,b, c)
!
      implicit none
!
      real(8), intent(in ) :: a(3),b(3)
      real(8), intent(out) :: c(3)
!
      c(1) =   a(2)*b(3) - a(3)*b(2)
      c(2) = - a(1)*b(3) + a(3)*b(1)
      c(3) =   a(1)*b(2) - a(2)*b(1)
!
   end subroutine cross_product
!
!
!--------------------------------------------------------------------
!
!     routine name      - zcross_product
!
!--------------------------------------------------------------------
!
!     latest revision:  - Mar 2023
!
!> @brief         - compute 3D cross product of real and
!                         (possibly) complex-valued vectors
!
!---------------------------------------------------------------------
!
   subroutine zcross_product(Rn,Za, Zcross)
!
      implicit none
!
      real(8), intent(in)  :: Rn(3)
      VTYPE  , intent(in)  :: Za(3)
      VTYPE  , intent(out) :: Zcross(3)
!
      Zcross(1) =   Rn(2)*Za(3) - Rn(3)*Za(2)
      Zcross(2) = - Rn(1)*Za(3) + Rn(3)*Za(1)
      Zcross(3) =   Rn(1)*Za(2) - Rn(2)*Za(1)
!
   end subroutine zcross_product
!
!
!--------------------------------------------------------------------
!
!     routine name      - zz_cross_product
!
!--------------------------------------------------------------------
!
!     latest revision:  - Mar 2023
!
!> @brief         - compute cross product of (possibly)
!                         complex-valued vectors
!
!               in:     - Za, Zb
!              out:     - Zcross = Za x Zb
!
!---------------------------------------------------------------------
!
   subroutine zz_cross_product(Za,Zb, Zcross)
!
      implicit none
!
      VTYPE , intent(in)  :: Za(3), Zb(3)
      VTYPE , intent(out) :: Zcross(3)
!
      Zcross(1) =   Za(2)*Zb(3) - Za(3)*Zb(2)
      Zcross(2) = - Za(1)*Zb(3) + Za(3)*Zb(1)
      Zcross(3) =   Za(1)*Zb(2) - Za(2)*Zb(1)
!
   end subroutine zz_cross_product
!
!--------------------------------------------------------------------
!
!     routine name      - cross_product2D
!
!--------------------------------------------------------------------
!
!     latest revision:  - Mar 2023
!
!> @brief         - compute 2D cross product of real-valued
!                         vectors
!
!---------------------------------------------------------------------
!
   subroutine cross_product2D(Vec1,Vec2, Vec_result)
!
      implicit none
!
      real(8), intent(in ) :: Vec1(2),Vec2(2)
      real(8), intent(out) :: Vec_result
!
      Vec_result=Vec1(1)*Vec2(2)-Vec1(2)*Vec2(1)
!
   end subroutine cross_product2D
!
!---------------------------------------------------------------------
!
!   routine name       - cross
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine evaluates cross product of two
!                        real-valued vectors in R^2 or R^3
!
!   arguments :
!     in:
!           N          - space dimension N=2,3
!           A,B        - vectors in R^N
!     out:
!           C          - scalar (R^2) or vector (R^3)
!
!----------------------------------------------------------------------
!
   subroutine cross(N,A,B, C)
!
      implicit none
!
      integer, intent(in)  :: N
      real(8), intent(in)  :: A(N),B(N)
      real(8), intent(out) :: C(2*N-3)
!
      if (N.eq.3) then
        C(1) =   A(2)*B(3) - A(3)*B(2)
        C(2) = - A(1)*B(3) + A(3)*B(1)
      endif
      C(2*N-3) = A(1)*B(2) - A(2)*B(1)
!
   end subroutine cross

!
#include "implicit_none.h"
!
!--------------------------------------------------------------------
!                                                                     
!     routine name      - cross_product
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
!-------------------------------------------------------------------- 
!
!..for backward compatibility
   subroutine cross_product3D(Vec1,Vec2, Vec)
      implicit none
      real(8), intent(in)  :: Vec1(3),Vec2(3)
      real(8), intent(out) :: Vec(3)
      call cross_product(Vec1,Vec2, Vec)
   end subroutine cross_product3D
!
!-------------------------------------------------------------------- 
!
!..for more backward compatibility
   subroutine cross_product3D_int(Vec1,Vec2, Vec)
      implicit none
      real(8), intent(in)  :: Vec1(3),Vec2(3)
      real(8), intent(out) :: Vec(3)
      call cross_product(Vec1,Vec2, Vec)
   end subroutine cross_product3D_int
!
!
!--------------------------------------------------------------------
!                                                                     
!     routine name      - zcross_product
!                                                                     
!-------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - Feb 2018
!                                                                     
!     purpose:          - compute cross product of real and
!                         (possibly) complex valued vectors
!                                                                    
!---------------------------------------------------------------------
!
   subroutine zcross_product(Rn,Za, Zcross)
!
      implicit none
!
      VTYPE   :: Za(3), Zcross(3)
      real(8) :: Rn(3)
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
!     latest revision:  - Aug 2018
!                                                                     
!     purpose:          - compute cross product of two
!                         (possibly) complex valued vectors
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

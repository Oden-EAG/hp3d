
module cross_product_module

implicit none

interface cross_product2D
   module procedure cross_product2D_int
   module procedure cross_product2D_real
   module procedure cross_product2D_complex
end interface

interface cross_product
   module procedure cross_product3D_int
   module procedure cross_product3D_real
   module procedure cross_product3D_complex
end interface

interface cross_product3D
   module procedure cross_product3D_int
   module procedure cross_product3D_real
   module procedure cross_product3D_complex
end interface


contains

  subroutine cross_product2D_int(Vec1,Vec2, Vec_result)
    implicit none
    integer, dimension(2), intent(in)  :: Vec1,Vec2
    integer,               intent(out) :: Vec_result
    Vec_result = Vec1(1)*Vec2(2) - Vec1(2)*Vec2(1)
  end subroutine cross_product2D_int

  subroutine cross_product2D_real(Vec1,Vec2, Vec_result)
    implicit none
    real(kind=kind(1.d0)), dimension(2), intent(in)  :: Vec1,Vec2
    real(kind=kind(1.d0)),               intent(out) :: Vec_result
    Vec_result = Vec1(1)*Vec2(2) - Vec1(2)*Vec2(1)
  end subroutine cross_product2D_real

  subroutine cross_product2D_complex(Vec1,Vec2, Vec_result)
    implicit none
    complex(kind=kind(1.d0)), dimension(2), intent(in)  :: Vec1,Vec2
    complex(kind=kind(1.d0)),               intent(out) :: Vec_result
    Vec_result = Vec1(1)*Vec2(2) - Vec1(2)*Vec2(1)
  end subroutine cross_product2D_complex

  subroutine cross_product3D_int(Vec1,Vec2,Vec_result)
    !useful to make "symbolic" cross product of standars basis e_i
    implicit none
    integer, dimension(3), intent(in)  :: Vec1, Vec2
    integer, dimension(3), intent(out) :: Vec_result
    Vec_result(1) =   Vec1(2)*Vec2(3) - Vec1(3)*Vec2(2)
    Vec_result(2) = - Vec1(1)*Vec2(3) + Vec1(3)*Vec2(1)
    Vec_result(3) =   Vec1(1)*Vec2(2) - Vec1(2)*Vec2(1)
  end subroutine cross_product3D_int

  subroutine cross_product3D_real(Vec1,Vec2,Vec_result)
    implicit none
    real(kind=kind(1.d0)), dimension(3), intent(in)  :: Vec1, Vec2
    real(kind=kind(1.d0)), dimension(3), intent(out) :: Vec_result
    Vec_result(1) =   Vec1(2)*Vec2(3) - Vec1(3)*Vec2(2)
    Vec_result(2) = - Vec1(1)*Vec2(3) + Vec1(3)*Vec2(1)
    Vec_result(3) =   Vec1(1)*Vec2(2) - Vec1(2)*Vec2(1)
  end subroutine cross_product3D_real

  subroutine cross_product3D_complex(Vec1,Vec2,Vec_result)
    implicit none
    complex(kind=kind(1.d0)), dimension(3), intent(in)  :: Vec1, Vec2
    complex(kind=kind(1.d0)), dimension(3), intent(out) :: Vec_result
    Vec_result(1) =   Vec1(2)*Vec2(3) - Vec1(3)*Vec2(2)
    Vec_result(2) = - Vec1(1)*Vec2(3) + Vec1(3)*Vec2(1)
    Vec_result(3) =   Vec1(1)*Vec2(2) - Vec1(2)*Vec2(1)
  end subroutine cross_product3D_complex

end module cross_product_module

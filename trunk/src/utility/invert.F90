!
#include "typedefs.h"
!
!--------------------------------------------------------------------
!> @brief Computes determinant of real-valued 3x3 matrix
!> @date Mar 2023
subroutine get_ddet(A, Adet)
   implicit none
   real(8), intent(in)  :: A(3,3)
   real(8), intent(out) :: Adet
   Adet = A(1,1)*A(2,2)*A(3,3) &
        + A(2,1)*A(3,2)*A(1,3) &
        + A(3,1)*A(1,2)*A(2,3) &
        - A(3,1)*A(2,2)*A(1,3) &
        - A(1,1)*A(3,2)*A(2,3) &
        - A(2,1)*A(1,2)*A(3,3)
end subroutine get_ddet
!
!--------------------------------------------------------------------
!> @brief Computes determinant of complex-valued 3x3 matrix
!> @date Mar 2023
subroutine get_zdet(A, Adet)
   implicit none
   VTYPE, intent(in)  :: A(3,3)
   VTYPE, intent(out) :: Adet
   Adet = A(1,1)*A(2,2)*A(3,3) &
        + A(2,1)*A(3,2)*A(1,3) &
        + A(3,1)*A(1,2)*A(2,3) &
        - A(3,1)*A(2,2)*A(1,3) &
        - A(1,1)*A(3,2)*A(2,3) &
        - A(2,1)*A(1,2)*A(3,3)
end subroutine get_zdet
!
!--------------------------------------------------------------------
!> @brief Computes inverse of real-valued 3x3 matrix
!> @date Mar 2023
subroutine dinvert(A, Ainv,Adet)
   implicit none
   real(8), dimension(3,3), intent(in)  :: A
   real(8), dimension(3,3), intent(out) :: Ainv
   real(8),                 intent(out) :: Adet
   real(8) :: det(3)
   logical, external :: dnear

   call get_ddet(A, Adet)

   if (dnear(Adet,0.d0)) then
      write(*,*) 'dinvert: Adet = 0'
      stop
   endif

   det(1) =   A(2,2)*A(3,3) - A(3,2)*A(2,3)
   det(2) = - A(2,1)*A(3,3) + A(3,1)*A(2,3)
   det(3) =   A(2,1)*A(3,2) - A(3,1)*A(2,2)

   Ainv(1,1) = det(1)/Adet
   Ainv(2,1) = det(2)/Adet
   Ainv(3,1) = det(3)/Adet

   det(1) =   A(3,2)*A(1,3) - A(1,2)*A(3,3)
   det(2) =   A(1,1)*A(3,3) - A(3,1)*A(1,3)
   det(3) = - A(1,1)*A(3,2) + A(3,1)*A(1,2)

   Ainv(1,2) = det(1)/Adet
   Ainv(2,2) = det(2)/Adet
   Ainv(3,2) = det(3)/Adet

   det(1) =   A(1,2)*A(2,3) - A(2,2)*A(1,3)
   det(2) = - A(1,1)*A(2,3) + A(2,1)*A(1,3)
   det(3) =   A(1,1)*A(2,2) - A(2,1)*A(1,2)

   Ainv(1,3) = det(1)/Adet
   Ainv(2,3) = det(2)/Adet
   Ainv(3,3) = det(3)/Adet

end subroutine dinvert
!
!--------------------------------------------------------------------
!> @brief Computes inverse of complex-valued 3x3 matrix
!> @date Mar 2023
subroutine zinvert(A, Ainv,Adet)
   use parameters, only: ZERO
   implicit none
   VTYPE, intent(in)  :: A(3,3)
   VTYPE, intent(out) :: Ainv(3,3)
   VTYPE, intent(out) :: Adet
   VTYPE :: det(3)
   logical, external :: znear

   call get_zdet(A, Adet)

   if (znear(Adet,ZERO)) then
      write(*,*) 'zinvert: Adet = 0'
      stop
   endif

   det(1) =   A(2,2)*A(3,3) - A(3,2)*A(2,3)
   det(2) = - A(2,1)*A(3,3) + A(3,1)*A(2,3)
   det(3) =   A(2,1)*A(3,2) - A(3,1)*A(2,2)

   Ainv(1,1) = det(1)/Adet
   Ainv(2,1) = det(2)/Adet
   Ainv(3,1) = det(3)/Adet

   det(1) =   A(3,2)*A(1,3) - A(1,2)*A(3,3)
   det(2) =   A(1,1)*A(3,3) - A(3,1)*A(1,3)
   det(3) = - A(1,1)*A(3,2) + A(3,1)*A(1,2)

   Ainv(1,2) = det(1)/Adet
   Ainv(2,2) = det(2)/Adet
   Ainv(3,2) = det(3)/Adet

   det(1) =   A(1,2)*A(2,3) - A(2,2)*A(1,3)
   det(2) = - A(1,1)*A(2,3) + A(2,1)*A(1,3)
   det(3) =   A(1,1)*A(2,2) - A(2,1)*A(1,2)

   Ainv(1,3) = det(1)/Adet
   Ainv(2,3) = det(2)/Adet
   Ainv(3,3) = det(3)/Adet

end subroutine zinvert

!
#include "typedefs.h"
!
!----------------------------------------------------------------------
!
!   function name      - dnear
!
!----------------------------------------------------------------------
!> @brief Returns true if the arguments are almost identical
!> @param[in] a - Real double-precision value
!> @param[in] b - Real double-precision value
!> @date  Mar 2023
!----------------------------------------------------------------------
logical function dnear(a,b)
!
   implicit none
!
   real(8) :: a,b
!
   real(8), parameter :: eps = 1.0d-15
!
   if (abs(a-b) .lt. eps) then
      dnear = .true.
   else
      dnear = .false.
   endif
!
end function dnear
!
!
!----------------------------------------------------------------------
!
!   function name      - znear
!
!----------------------------------------------------------------------
!> @brief Returns true if the arguments are almost identical
!> @param[in] a - Double-precision value of type VTYPE
!> @param[in] b - Double-precision value of type VTYPE
!> @date  Mar 2023
!----------------------------------------------------------------------
logical function znear(a,b)
!
   implicit none
!
   VTYPE   :: a,b
   real(8) :: diff_real,diff_imag
!
   real(8), parameter :: eps = 1.0d-15
!
   real(8), external :: dreal_part,dimag_part
!
   diff_real = dreal_part(a-b)
   diff_imag = dimag_part(a-b)
!
   if (abs(diff_real) .lt. eps .and. &
       abs(diff_imag) .lt. eps ) then
      znear = .true.
   else
      znear = .false.
   endif
!
end function znear

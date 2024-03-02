#include "typedefs.h"

!> @brief Returns real part of complex number
!> @date Feb 2023
real(8) function dreal_part(Z)
!
   implicit none
   VTYPE, intent(in) :: Z
!
#if HP3D_COMPLEX
   dreal_part = dreal(Z)
#else
   dreal_part = Z
#endif
!
end function dreal_part
!
!
!> @brief Returns imag part of complex number
!> @date Feb 2023
real(8) function dimag_part(Z)
!
   implicit none
   VTYPE, intent(in) :: Z
!
#if HP3D_COMPLEX
   dimag_part = dimag(Z)
#else
   dimag_part = 0.d0
#endif
!
end function dimag_part

!> @date Mar 2023
program real_imag
!
   implicit none
!
   integer :: NPASS
!
#if HP3D_COMPLEX
   complex(8) :: a
#else
   real(8)    :: a
#endif
!
   real(8), parameter :: eps = 1.0d-15
   real(8), external  :: dreal_part,dimag_part
!
   NPASS = 1
!
#if HP3D_COMPLEX
   a = (1.d0, 2.d0)
   if (dabs(dreal_part(a) - 1.d0) > eps) NPASS = 0
   if (dabs(dimag_part(a) - 2.d0) > eps) NPASS = 0
#else
   a = 1.d0
   if (dabs(dreal_part(a) - 1.d0) > eps) NPASS = 0
   if (dabs(dimag_part(a) - 0.d0) > eps) NPASS = 0
#endif
!
!
   if (NPASS.ne.1) stop 1
!
end program real_imag

!> @date Mar 2023
program real_imag
!
   implicit none
!
   integer :: NPASS
!
#if C_MODE
   complex(8) :: a
#else
   real(8)    :: a
#endif
!
   real(8), external :: dreal_part,dimag_part
!
   NPASS = 1
!
#if C_MODE
   a = (1.d0, 2.d0)
   if (dreal_part(a) .ne. 1.d0) NPASS = 0
   if (dimag_part(a) .ne. 2.d0) NPASS = 0
#else
   a = 1.d0
   if (dreal_part(a) .ne. 1.d0) NPASS = 0
   if (dimag_part(a) .ne. 0.d0) NPASS = 0
#endif
!
!
   if (NPASS.ne.1) stop 1
!
end program real_imag

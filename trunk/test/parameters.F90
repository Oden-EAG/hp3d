!> @date Mar 2023
program test_parameters
   use parameters
!
   implicit none
!
   integer :: NPASS
!
   integer :: nsum
!
   NPASS = 1
!
   call set_parameters(2,1)
!
   if (NRCOMS .ne. 2) NPASS = 0
   if (NRRHS  .ne. 1) NPASS = 0
!
   if (NPASS.ne.1) stop 1
!
end program test_parameters

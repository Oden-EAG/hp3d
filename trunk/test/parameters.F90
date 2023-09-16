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
   call set_parameters(2,1,3,4,5,6)
!
   nsum = NRCOMS + NRRHS + &
          MAXEQNH + MAXEQNE + MAXEQNV + MAXEQNQ
!
   if (nsum .ne. 21) NPASS = 0
!
!
   if (NPASS.ne.1) stop 1
!
end program test_parameters

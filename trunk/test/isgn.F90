!> @date Mar 2023
program test_isgn
!
   implicit none
!
   integer :: NPASS
!
   integer :: a
!
   integer, external :: isgn
!
   NPASS = 1
!
   a = 2
   if (isgn(a) .ne.  1) NPASS = 0
!
   a = -5
   if (isgn(a) .ne. -1) NPASS = 0
!
   a = 0
   if (isgn(a) .ne.  0) NPASS = 0
!
!
   if (NPASS.ne.1) stop 1
!
end program test_isgn

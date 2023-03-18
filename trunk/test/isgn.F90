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
   if (NPASS .ne. 1) then
      write(*,*) 'test_isgn FAILED.'
   else
      write(*,*) 'test_isgn PASSED.'
   endif
!
end program test_isgn

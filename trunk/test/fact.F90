!> @date Mar 2023
program test_fact
!
   implicit none
!
   integer :: NPASS
!
   integer :: a
!
   integer, external :: fact
!
   NPASS = 1
!
   a = 5
   if (fact(a) .ne. 120) NPASS = 0
!
   a = 0
   if (fact(a) .ne.   1) NPASS = 0
!
   a = -2
   if (fact(a) .ne.   0) NPASS = 0
!
   if (NPASS .ne. 1) then
      write(*,*) 'test_fact FAILED.'
   else
      write(*,*) 'test_fact PASSED.'
   endif
!
end program test_fact

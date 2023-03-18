!> @date Mar 2023
program test_locate
!
   implicit none
!
   integer :: NPASS
!
   integer :: a(10), loc
!
   NPASS = 1
!
   a = 0; a(7) = 5
!
   call locate(5,a,10, loc)
   if (loc .ne. 7) NPASS = 0
!
   call locate(4,a,10, loc)
   if (loc .ne. 0) NPASS = 0
!
   if (NPASS .ne. 1) then
      write(*,*) 'test_locate FAILED.'
   else
      write(*,*) 'test_locate PASSED.'
   endif
!
end program test_locate

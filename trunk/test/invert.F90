!> @date Mar 2023
program test_invert
!
   implicit none
!
   integer :: NPASS
!
   NPASS = 1
!
!  TODO: implement unit test
!
   if (NPASS .ne. 1) then
      write(*,*) 'test_invert FAILED.'
   else
      write(*,*) 'test_invert PASSED.'
   endif
!
end program test_invert

!> @date Mar 2023
program test_locate
   implicit none
!
   integer :: a(10), loc
!
   a = 0; a(7) = 5
!
   call locate(5, a, 10, loc)
!
   if (loc .ne. 7) then
      write(*,*) 'test_locate FAILED.'
      stop
   endif
   write(*,*) 'test_locate PASSED.'
!
end program test_locate

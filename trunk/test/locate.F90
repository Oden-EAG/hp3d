!> @date Mar 2023
program test_locate
   implicit none
!
   integer :: a(10), loc
!
   a(7) = 5
!
   call locate(5, a, 10, loc)
   write(*,*) 'loc ',loc
!
end program test_locate

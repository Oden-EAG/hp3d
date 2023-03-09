!> @date Mar 2023
program test_parameters
   use parameters
   implicit none
!
   integer :: sum
!
   call set_parameters(1,2,3,4,5,6)
!
   sum = NRCOMS + MAXNRHS + &
         MAXEQNH + MAXEQNE + MAXEQNV + MAXEQNQ
!
   write(*,*) 'sum (21) = ', sum
!
end program test_parameters

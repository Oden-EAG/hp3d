!> @date Mar 2023
program test_parameters
   use parameters
   implicit none
!
   integer :: nsum
!
   call set_parameters(1,2,3,4,5,6)
!
   nsum = NRCOMS + MAXNRHS + &
          MAXEQNH + MAXEQNE + MAXEQNV + MAXEQNQ
!
   if (nsum .ne. 21) then
      write(*,*) 'test_parameters FAILED.'
      stop
   endif
   write(*,*) 'test_parameters PASSED.'
!
end program test_parameters

      program test_parameters
      use parameters
      implicit none
      integer :: sum 
      call set_parameters(1,2,3,4,5,6,7)

      sum = NDIMEN() + NRCOMS() + MAXNRHS() + 
     .  MAXEQNS() + MAXEQNE() + MAXEQNV() + MAXEQNQ()

      write(*,*) 'sum (28) = ', sum
      
      
      end program test_parameters
      

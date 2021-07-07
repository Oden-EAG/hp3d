subroutine exact_error_l2(Mdle, Derr,Dnorm)
!
      use control
      use element_data
      use data_structure3D
!
      implicit none      
      integer,                  intent(in ) :: Mdle
      real*8,dimension(MAXEQNQ),intent(out) :: Derr
      real*8,dimension(MAXEQNQ),intent(out) :: Dnorm
!-------------------------------------------------------------      
!      
      Derr=0.d0 ; Dnorm=0.d0
!
!
endsubroutine exact_error_l2
 

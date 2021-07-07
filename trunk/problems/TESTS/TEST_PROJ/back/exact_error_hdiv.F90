subroutine exact_error_hdiv(Mdle, Derr,Dnorm)
!
      use control
      use element_data
      use data_structure3D
!
      implicit none      
      integer,                  intent(in ) :: Mdle
      real*8,dimension(MAXEQNV),intent(out) :: Derr
      real*8,dimension(MAXEQNV),intent(out) :: Dnorm
!-------------------------------------------------------------      
!      
      Derr=0.d0 ; Dnorm=0.d0
!
!
endsubroutine exact_error_hdiv
 

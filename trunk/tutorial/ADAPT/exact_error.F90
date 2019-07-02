!------------------------------------------------------------------------------------
!> Purpose : computes norm of the error and the exact solution      
!!
!! @param[out]  Err   - norm of the error      
!! @param[out]  Rnorm - norm of the exact solution
!------------------------------------------------------------------------------------
!
subroutine exact_error(Err,Rnorm)
  !
  use data_structure3D , only : MAXEQNH,NRELES
  !------------------------------------------------------------------------------------
  implicit none
  real*8,dimension(MAXEQNH),intent(out) :: Err
  real*8,dimension(MAXEQNH),intent(out) :: Rnorm
  !------------------------------------------------------------------------------------
  !  ...increments      
  real*8,dimension(MAXEQNH) :: derr
  real*8,dimension(MAXEQNH) :: dnorm
  !  ...miscellanea      
  integer :: mdle,iel,n
  integer :: iprint
  !------------------------------------------------------------------------------------
  !
  iprint=0
  !
  !  ...initialize
  Err = 0.d0; Rnorm = 0.d0
  !
  !  ...loop over ACTIVE elements      
  mdle=0
  do iel=1,NRELES
     !      
     !  .....active element iterator
     call nelcon(mdle, mdle)
     !  .....compute element error
     call exact_error_element(mdle, derr(1:MAXEQNH),dnorm(1:MAXEQNH))
     !  .....printing
     if (iprint.eq.1) then
        do n=1,MAXEQNH
           write(*,7001) n,mdle,derr(n),dnorm(n),derr(n)/dnorm(n)
7001       format(' compute_error: n,mdle,derr,dnorm,derr/dnorm = ',i4,i6,10(e12.5))
        enddo
     endif
     !  .....accumulate error^2
     Err(  1:MAXEQNH) = Err(  1:MAXEQNH) + derr( 1:MAXEQNH)
     Rnorm(1:MAXEQNH) = Rnorm(1:MAXEQNH) + dnorm(1:MAXEQNH)

     !  ...loop over ACTIVE elements        
  enddo
  !      
  !  ...compute square root
  do n=1,MAXEQNH
     Err(n)   = sqrt(Err(  n))
     Rnorm(n) = sqrt(Rnorm(n))
  enddo
  !
  !
end subroutine exact_error

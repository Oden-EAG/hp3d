!> Purpose : calculate exact error
!! @param[out] Err  
!! @param[out] Rnorm
subroutine exact_error(Err,Rnorm)
  use data_structure3D
  use proj

  implicit none
  ! total quantity
  real*8, dimension(MAXEQN_PROB), intent(inout) ::  Err, Rnorm
  real*8, dimension(MAXEQN_PROB)                :: derr, dnorm

  integer :: iprint, neq, mdle, iel, itype, n
  !-----------------------------------------------------------------------
  iprint = 0
  !
  Err = 0.d0; Rnorm = 0.d0
  neq = 0

  ! loop over ACTIVE elements      
  mdle=0
  do iel=1,NRELES
     call nelcon(mdle, mdle)
     neq = MAXEQNH
     call exact_error_element(mdle, derr,dnorm)

     if (iprint.eq.1) then
        write(*,7000) itype, mdle
        write(*,7001) 'error', derr
        write(*,7001) 'norm ', dnorm
        
7000    format('exact_error: problem(1 EM, 2 HEAT), mdle', i3,2x,i6)
7001    format('            ',a5, e12.5)
     endif

     ! accumulate        
     Err(1:neq)   = Err(1:neq)   + derr(1:neq)
     Rnorm(1:neq) = Rnorm(1:neq) + dnorm(1:neq)
  enddo
  !
  do n=1,neq
     Err(n) = sqrt(Err(n))
     Rnorm(n) = sqrt(Rnorm(n))
  enddo
  !
end subroutine exact_error

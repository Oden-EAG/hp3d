!> Purpose : calculate geometry error
!! @param[out] Err  
!! @param[out] Rnorm
subroutine geometry_error(Err,Rnorm)
  use data_structure3D
  implicit none

  ! total quantity
  real*8, intent(inout) ::  Err, Rnorm
  real*8                :: derr, dnorm

  integer :: iprint, mdle, iel
  !-----------------------------------------------------------------------
  iprint = 0

  Err = 0.d0; Rnorm = 0.d0

  mdle=0
  do iel=1,NRELES
     call nelcon(mdle, mdle)
     call geometry_error_element(mdle, derr,dnorm)

     if (iprint.eq.1) then
        write(*,7000) mdle
        write(*,7001) 'error', derr
        write(*,7001) 'norm ', dnorm
        
7000    format('geometry_error: mdle', i6)
7001    format('            ',a5, e12.5)
     endif

     ! accumulate        
     Err   = Err   + derr
     Rnorm = Rnorm + dnorm
  enddo

  Err   = sqrt(Err)
  Rnorm = sqrt(Rnorm)

end subroutine geometry_error

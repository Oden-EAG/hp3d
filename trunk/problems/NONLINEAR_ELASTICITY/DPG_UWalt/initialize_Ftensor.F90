! subroutine that sets the "previous" value of F as the identity matrix
subroutine initialize_Ftensor
use data_structure3D
implicit none
integer :: mdle,F_address2
real*8  :: vec_id(9)

! vectorized identity matrix
vec_id = 0.d0
vec_id(1) = 1.d0
vec_id(5) = 1.d0
vec_id(9) = 1.d0
! set the base index for the 2nd copy of components of F
F_address2 = NRQVAR + 12
! we assume we are at the initial mesh so we can write directly at NODES
do mdle=1,NRELIS
! we set the first degree of freedom (column of zdofQ),  
! which for L2 variables is the constant Legendre polynomial
  ! write(*,*) 'memory address', loc(NODES(mdle)%dof%zdofQ)
  NODES(mdle)%dof%zdofQ(F_address2+1:F_address2+9,1) = vec_id(1:9)
 !  write(*,1000) NODES(mdle)%dof%zdofQ(1:42,1)
 ! 1000 format(42(es8.1))
  ! write(*,*) lbound(NODES(mdle)%dof%zdofQ,2),ubound(NODES(mdle)%dof%zdofQ,2)
enddo


end subroutine
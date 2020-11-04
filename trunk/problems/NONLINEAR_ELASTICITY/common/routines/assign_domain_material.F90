! We use GMP domains to determine materials
! July 2020
subroutine assign_domain_material
use hyperelasticity, only : GMP_MAT
use GMP, only : NRDOMAIN
implicit none
integer :: n

allocate(GMP_MAT(NRDOMAIN))
! this is the default: 
! domain number equals material number
do n=1,NRDOMAIN
   GMP_MAT(n) = n
enddo

! For particular problems uncomment and fill the following
!
if (NRDOMAIN.gt.1) then
GMP_MAT = 1
GMP_MAT(NRDOMAIN) = 2
endif

end subroutine
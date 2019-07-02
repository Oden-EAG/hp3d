!-----------------------------------------------------------------------------
!> Purpose : routine determines the domain of an element possibly originated
!!           from refinements. Remark: array NODES must have been previously
!!           created. Routine domain_number can be used for initial mesh 
!!           elements, when nodes have not been created yet.
!!
!! @param[in]  Mdle    - an element middle node
!! @param[out] Ndomain - domain number 
!!
!! @revision Oct 12
!-----------------------------------------------------------------------------
!
subroutine find_domain(Mdle, Ndom)
  use GMP
  use data_structure3D
  implicit none
  integer, intent(in)  :: Mdle
  integer, intent(out) :: Ndom
!
  integer :: nfath, nel
!-----------------------------------------------------------------------------
!
! determine initial mesh ancestor
  nfath = NODES(Mdle)%father
  do while(nfath.gt.0)
     nfath = NODES(nfath)%father
  enddo
!
! determine domain number
  nel = abs(nfath)
  Ndom = Ndomain_gmp(ELEMS(nel)%GMPblock)
!
!
end subroutine find_domain




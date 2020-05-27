!-----------------------------------------------------------------------------------------
!> Purpose : routine determines domain number of an initial mesh element using
!!           the information contained in the ELEMS array (cf. subroutine
!!           find_domain to determine domain of a refined element).
!!
!! @param[in]  N_elem - element number
!! @param[out] N_dom  - domain number (= 0, if element N_elem does not exits, e.g., N_elem
!!                      is the non-existing neighbor of an element on the boundary)
!!
!! @revision Oct 12
!-----------------------------------------------------------------------------------------
!
subroutine domain_number(N_elem, N_dom)
!
  use data_structure3D
  use GMP
!
  IMPLICIT NONE
!
! DUMMY ARGUMENTS
  integer, intent(in)  :: N_elem
  integer, intent(out) :: N_dom
!
! LOCAL VARIABLES
  integer :: nbln,labn
!-----------------------------------------------------------------------------------------
!
    if (N_elem.eq.0) then
      N_dom=0 ; return
    endif
!
    call decode(ELEMS(N_elem)%GMPblock, nbln,labn)
    select case(labn)
    case(1) ; N_dom=PRISMS(  nbln)%Domain
    case(2) ; N_dom=HEXAS(   nbln)%Domain
    case(3) ; N_dom=TETRAS(  nbln)%Domain
    case(4) ; N_dom=PYRAMIDS(nbln)%Domain
    case default
      write(*,*)'domain_number: invalid block type!'
      write(*,7000)N_elem,labn
 7000 format(' N_elem,labn = ',i8,2x,i1)
      call result
      stop
    endselect
!
!
end subroutine domain_number

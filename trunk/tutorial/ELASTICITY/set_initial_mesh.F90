!------------------------------------------------------------------------------------
!> Purpose : define problem dependent data ( multiphysics, BC, approxiamtion )
!! @param[out] Nelem_order - order for initial mesh elements
!------------------------------------------------------------------------------------
!
subroutine set_initial_mesh(Nelem_order)
  use GMP
  use data_structure3D
!------------------------------------------------------------------------------------
  implicit none
  integer,dimension(NRELIS),intent(out) :: Nelem_order
!------------------------------------------------------------------------------------
! BC flags
  integer, dimension(6) :: ibc
! miscellanea  
  integer :: nel,ndom,nbl,lab
!------------------------------------------------------------------------------------
!
! loop over initial mesh elements
  do nel=1,NRELIS
!
!    set physics
     ELEMS(nel)%nrphysics = 1
     allocate(ELEMS(nel)%physics(1)) 
     ELEMS(nel)%physics(1) = 'elast'
!
!    set order of approximation
     select case(ELEMS(nel)%Type)
     case('pris'); Nelem_order(nel) = 33
     case('bric'); Nelem_order(nel) = 333
     case('tetr'); Nelem_order(nel) = 3
     case('pyra'); Nelem_order(nel) = 3
     end select
!
!    set BC flags: 0 - no BC ; 1 - Dirchlet ; 2 - Neumann
     ibc(1:6) = 0 ; ibc(1) = 1
     allocate(ELEMS(nel)%bcond(NR_PHYSA))
     call encodg(ibc,10,6, ELEMS(nel)%bcond(1))
!     
! end of loop over initial mesh elements
  enddo
!
end subroutine set_initial_mesh

!> Purpose : define problem dependent data ( multiphysics, BC, approxiamtion )
!! @param[out] Nelem_order - order for initial mesh elements
subroutine set_initial_mesh(Nelem_order)
  use GMP
  use data_structure3D

  integer, dimension(NRELIS), intent(out) :: Nelem_order
  integer, dimension(6) :: ibc
  integer :: nel, ndom, nbl, lab

  ! TODO :: This routine is ambiguous...
  !         User should not allocate or deallocate for the datastructure
  ! temporary setting single physics for electromagnetics

  do nel=1,NRELIS

     ELEMS(nel)%nrphysics = 1
     allocate(ELEMS(nel)%physics(1)) 
     ELEMS(nel)%physics(1) = 'poiss'
     !
     ! ** Assign order of approximation
     select case(ELEMS(nel)%Type)
     case('pris'); Nelem_order(nel) = 22
     case('bric'); Nelem_order(nel) = 222
     case('tetr'); Nelem_order(nel) = 1
     end select
     !
     ! ** Assign boundary condition
     ibc(1:6) = 0
     do i=1, nface(ELEMS(nel)%Type)
        if (ELEMS(nel)%neig(i).eq.0) then
           ibc(i) = 1
        endif
     enddo
     

!!!     ibc(1:6) = 0
     allocate(ELEMS(nel)%bcond(NR_PHYSA))
     call encodg(ibc,10,6, ELEMS(nel)%bcond(1))

  enddo
end subroutine set_initial_mesh

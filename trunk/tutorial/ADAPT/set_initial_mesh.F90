!> Purpose : define problem dependent data ( multiphysics, BC, approxiamtion )
!! @param[out] Nelem_order - order for initial mesh elements
subroutine set_initial_mesh(Nelem_order)
  use GMP
  use data_structure3D
  use lapl
  implicit none
  integer, dimension(NRELIS), intent(out) :: Nelem_order
  integer, dimension(6) :: ibc
  integer :: i, nel, neig, nblk_elem

  do nel=1,NRELIS

     ELEMS(nel)%nrphysics = 1
     allocate(ELEMS(nel)%physics(1)) 
     ELEMS(nel)%physics(1) = PHYSA(1)
     !
     ! ** Assign order of approximation
     nblk_elem = Ndomain_gmp(ELEMS(nel)%GMPblock)

     select case(ELEMS(nel)%Type)
     case('pris'); Nelem_order(nel) = IORDER_PRIS_PROB
     case('tetr'); Nelem_order(nel) = IORDER_TETR_PROB
     end select
     !
     ! ** Assign boundary condition
     ibc(1:6) = 0
     do i=1, nface(ELEMS(nel)%Type)
        neig = ELEMS(nel)%neig(i)
        select case (neig)
        case(0); ibc(i) = IBC_PROB 
        end select
     enddo
     !
     allocate(ELEMS(nel)%bcond(NR_PHYSA))
     call encodg(ibc,10,6, ELEMS(nel)%bcond(1))
     !
  enddo
end subroutine set_initial_mesh

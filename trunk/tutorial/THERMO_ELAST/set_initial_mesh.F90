!> Purpose            - routine sets problem dependent attributes
!!                      of the initial mesh for the projection
!!                      problem
!! @param[out] Nelem_order - order for initial mesh elements
subroutine set_initial_mesh(Nelem_order)
  !
  use GMP
  use data_structure3D
  !

  integer, dimension(NRELIS), intent(out) :: Nelem_order
  integer, dimension(6) :: ibcE, ibcT
  integer :: nel, ndom, nbl, lab
  !
  !  ...loop through initial mesh elements
  do nel=1,NRELIS
     ELEMS(nel)%nrphysics = 2
     allocate(ELEMS(nel)%physics(2)) 
     ibcE(1:6)=0; ibcT(1:6)=0
     !
     !  .....get the corresponding GMP domain
     call decode(ELEMS(nel)%GMPblock, nbl,lab)
     select case(lab)
     case(1); ndom = PRISMS(nbl)%Domain
     case(2); ndom = HEXAS(nbl)%Domain
     case(3); ndom = TETRAS(nbl)%Domain
     case(4); ndom = PYRAMIDS(nbl)%Domain
     case default
        write(*,*) 'set_initial_mesh: lab = ',lab
        stop1
     end select
     !
     !  .....set up the physics
     select case(ndom)
        !
        !  .....thermoealsticity
     case(1)
        ELEMS(nel)%physics(1) = 'displ'
        ELEMS(nel)%physics(2) = 'tempr'
     case default
        stop1
     end select
     !
     !  .....set up the element order
     select case(ELEMS(nel)%Type)
     case('pris'); Nelem_order(nel) = 44   
     case('bric'); Nelem_order(nel) = 323
     case('tetr'); Nelem_order(nel) = 4
     case('pyra'); Nelem_order(nel) = 4
     end select
     !
     !  .....set BC flags
     select case(nel)
     case(1)
        ibcE(1:6) = (/0,0,2,3,2,1/)
        ibcT(1:6) = (/0,0,2,3,2,1/)
     case(2)
        ibcE(1:6) = (/0,0,2,3,2,1/)
        ibcT(1:6) = (/0,0,2,3,2,1/)
     case(3)
        ibcE(1:6) = (/0,0,2,3,2,1/)
        ibcT(1:6) = (/0,0,2,3,2,1/)
     case(4)
        ibcE(1:6) = (/0,0,2,3,2,1/)
        ibcT(1:6) = (/0,0,2,3,2,1/)
     end select
     !
     !  .....encode the BC's
     allocate(ELEMS(nel)%bcond(2)) 

     call encodg(ibcE,10,6, ELEMS(nel)%bcond(1))
     call encodg(ibcT,10,6, ELEMS(nel)%bcond(2))
     !
     !  ...end of loop through initial mesh elements
  enddo
  !
  !
end subroutine set_initial_mesh

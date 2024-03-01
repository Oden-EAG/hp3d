!--------------------------------------------------------------------
!
!     routine name      - set_initial_mesh
!
!--------------------------------------------------------------------
!
!     latest revision:  - Oct 2021
!
!> @brief         - define problem dependent data
!                         (multiphysics, BC, approximation)
!
!     arguments:
!
!     in:
!     out:
!           Nelem_order - order for initial mesh elements
!
!---------------------------------------------------------------------
!
   subroutine set_initial_mesh(Nelem_order)
!
   use GMP
   use common_prob_data
   use data_structure3D
!
   implicit none
!
!----------------------------------------------------------------------
!
   integer, intent(out) :: Nelem_order(NRELIS)
!..BC flags
   integer :: ibc(6,NRINDEX)
!..miscellaneous
   integer :: ifc,iel,neig
!
!------------------------------------------------------------------------------------
!
!..check if exceeding the maximum order
   if (IP.gt.MAXP) then
      write(*,*) 'set_initial_mesh: IP, MAXP = ', IP,MAXP
      stop
   endif
!
!..set BC
!..loop over initial mesh elements
   do iel=1,NRELIS
!
!  ...set physics
      ELEMS(iel)%nrphysics = 1
      allocate(ELEMS(iel)%physics(1))
      ELEMS(iel)%physics(1) ='field'
!
!  ...set order of approximation
      if (IP.gt.0) then
!     ...uniform order of approximation
         select case(ELEMS(iel)%etype)
            case(BRIC); Nelem_order(iel) = 111*IP
            case(PRIS); Nelem_order(iel) = 11*IP
            case(TETR); Nelem_order(iel) = 1*IP
            case(PYRA); Nelem_order(iel) = 1*IP
         end select
      else
!     ...custom order of approximation (NOT IMPLEMENTED)
         write(*,1003) IP
1003     format('ERROR in set_initial_mesh: IP = ',i3)
         stop
      endif
!
!  ...set BC flags: 0 - no BC ; 1 - Dirichlet
!     (Maxwell Galerkin: NRINDEX = 1)
      ibc(1:6,1:NRINDEX) = 0
!
      do ifc=1,nface(ELEMS(iel)%etype)
         neig = ELEMS(iel)%neig(ifc)
         select case(neig)
            case(0); ibc(ifc,1) = 1
         end select
      enddo
!
!  ...allocate BC flags (one per attribute)
      allocate(ELEMS(iel)%bcond(1))
!  ...encode face BC into a single BC flag
      call encodg(ibc(1:6,1),10,6, ELEMS(iel)%bcond(1))
!
!..end of loop over initial mesh elements
   enddo
!
!
end subroutine set_initial_mesh

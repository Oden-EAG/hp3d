!------------------------------------------------------------------------------------
!
!     routine name      - set_initial_mesh
!
!------------------------------------------------------------------------------------
!
!     latest revision:  - Oct 2021
!
!     purpose:          - define problem dependent data
!                         (multiphysics, BC, approximation)
!
!     arguments:
!
!     in:
!     out:
!           Nelem_order - order for initial mesh elements
!
!------------------------------------------------------------------------------------
!
   subroutine set_initial_mesh(Nelem_order)
!
   use GMP
   use common_prob_data
   use data_structure3D
!
   implicit none
!
!------------------------------------------------------------------------------------
!
   integer,dimension(NRELIS),intent(out) :: Nelem_order
!..BC flags
   integer, dimension(6,NRINDEX) :: ibc
!..miscellaneous
   integer :: i,ifc,iel,neig
!
!------------------------------------------------------------------------------------
!
!..check if have not exceeded the maximum order
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
      ELEMS(iel)%nrphysics = 2
      allocate(ELEMS(iel)%physics(2))
      ELEMS(iel)%physics(1) ='field'
      ELEMS(iel)%physics(2) ='trace'
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
      ibc(1:6,1:NRINDEX) = 0
!
      select case(IBC_PROB)
!
!     ...uniform BC
         case(BC_DIRICHLET)
!        ...if exterior face, set boundary condition to IBC_PROB
            do ifc=1,NFACE(ELEMS(iel)%etype)
               neig = ELEMS(iel)%neig(ifc)
               select case(neig)
                  case(0); ibc(ifc,1) = 1 ! Dirichlet BC (H1 field variable)
               end select
            enddo
!
         case default
            write(*,*) 'set_initial_mesh: problem currently supports Dirichlet BC only. stop.'
            stop
      end select
!
!  ...allocate BC flags (one per attribute component)
      allocate(ELEMS(iel)%bcond(NRINDEX))
!
!  ...for each component, encode face BC into a single BC flag
      do i=1,NRINDEX
         call encodg(ibc(1:6,i),10,6, ELEMS(iel)%bcond(i))
      enddo
!
!..end of loop over initial mesh elements
   enddo
!
!
end subroutine set_initial_mesh

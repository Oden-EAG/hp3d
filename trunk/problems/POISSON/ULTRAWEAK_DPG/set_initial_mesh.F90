!------------------------------------------------------------------------------------
!
!     routine name      - set_initial_mesh
!
!------------------------------------------------------------------------------------
!
!     latest revision:  - May 2023
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
!..0 for cube with Dirichlet boundary condition,
!..1 for Fichera corner (make sure to change default file-geometry in common/set_enviroment.F90)
!..2 for the cube but arc tan test case
   integer :: geomtype 
!
!------------------------------------------------------------------------------------
!
!  setting the parameter for choosing BC associated with a given geometry
   geomtype = 1
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
!..set physics
      ELEMS(iel)%nrphysics = 4
      allocate(ELEMS(iel)%physics(4))
      ELEMS(iel)%physics(1) ='trc_a'
      ELEMS(iel)%physics(2) ='trc_b'
      ELEMS(iel)%physics(3) ='field'
      ELEMS(iel)%physics(4) ='gradu'
      
!..set order of approximation
      if (IP.gt.0) then
!..uniform order of approximation
         select case(ELEMS(iel)%etype)
            case(BRIC); Nelem_order(iel) = 111*IP
            case(PRIS); Nelem_order(iel) = 11*IP
            case(TETR); Nelem_order(iel) = 1*IP
            case(PYRA); Nelem_order(iel) = 1*IP
         end select
      else
!..custom order of approximation (NOT IMPLEMENTED)
         write(*,1003) IP
1003     format('ERROR in set_initial_mesh: IP = ',i3)
         stop
      endif
!
!..set BC flags: 0 - no BC ; 1 - Dirichlet
      ibc(1:6,1:NRINDEX) = 0
!
      select case(IBC_PROB)
!
!..uniform BC
         case(BC_DIRICHLET)
!..if exterior face, set boundary condition to IBC_PROB
            if(geomtype .eq. 0) then
            ! all external faces have Dirichlet BC on H1 trace
               do ifc=1,nface(ELEMS(iel)%etype)
               neig = ELEMS(iel)%neig(ifc)
                  select case(neig)
                     case(0); ibc(ifc,1) = 1 ! Dirichlet BC (H1 trace variable)
                  end select
               enddo
!
            else if (geomtype .eq. 2) then
            ! three faces on co-ordinate planes have H1 trace Dirichlet BC
            ! and the other three have H-div trace Dirichlet BC
               ibc(1,1) = 1
               ibc(2,2) = 1
               ibc(3,1) = 1
               ibc(4,2) = 1
               ibc(5,2) = 1
               ibc(6,1) = 1
!
            else if (geomtype .eq. 1) then ! Fichera Corner BCs (Dirichlet BC for normal flux)
               if(iel .eq. 4) then
                  do ifc=1,nface(ELEMS(iel)%etype)
                     neig = ELEMS(iel)%neig(ifc)
                     select case(neig)
                        case(0)
                           if(ifc .eq. 2) then
                               ibc(ifc,1) = 1 ! Dirichlet BC (H1 trace variable)
                           else
                               ibc(ifc,2) = 1 ! Dirichlet BC (Hdiv trace variable)
                           endif
                     end select
                  enddo
               elseif (iel .eq. 6) then
                  do ifc=1,nface(ELEMS(iel)%etype)
                     neig = ELEMS(iel)%neig(ifc)
                     select case(neig)
                        case(0)
                           if(ifc .eq. 5) then
                               ibc(ifc,1) = 1 ! Dirichlet BC (H1 trace variable)
                           else
                               ibc(ifc,2) = 1 ! Dirichlet BC (Hdiv trace variable)
                           endif
                     end select
                  enddo
               elseif (iel .eq. 7) then
                  do ifc=1,nface(ELEMS(iel)%etype)
                     neig = ELEMS(iel)%neig(ifc)
                     select case(neig)
                        case(0)
                           if(ifc .eq. 4) then
                               ibc(ifc,1) = 1 ! Dirichlet BC (H1 trace variable)
                           else
                               ibc(ifc,2) = 1 ! Dirichlet BC (Hdiv trace variable)
                           endif
                     end select
                  enddo
               else 
                  do ifc=1,nface(ELEMS(iel)%etype)
                     neig = ELEMS(iel)%neig(ifc)
                        select case(neig)
                           case(0); ibc(ifc,2) = 1 ! Dirichlet BC (Hdiv trace variable)
                        end select
                  enddo
               endif
            endif    
         case default
            write(*,*) 'set_initial_mesh: problem currently supports Dirichlet BC only. stop.'
            stop
      end select
!
!..allocate BC flags (one per attribute component)
      allocate(ELEMS(iel)%bcond(NRINDEX))
!
!..for each component, encode face BC into a single BC flag
      do i=1,NRINDEX
         call encodg(ibc(1:6,i),10,6, ELEMS(iel)%bcond(i))
      enddo
!
!..end of loop over initial mesh elements
   enddo
!
!
end subroutine set_initial_mesh

!-------------------------------------------------------------------------------
!
!     routine name      - set_initial_mesh
!
!-------------------------------------------------------------------------------
!
!     latest revision:  - May 2021
!
!     purpose:          - define problem dependent data
!                         (multiphysics, BC, approximation)
!
!     arguments:
!        out:
!           Nelem_order - order for initial mesh elements
!
!-------------------------------------------------------------------------------
!
subroutine set_initial_mesh(Nelem_order)
!
   use GMP
   use common_prob_data
   use data_structure3D
!
   implicit none
!
!..output argument
   integer, intent(out) :: Nelem_order(NRELIS)
!
!..BC flags
   integer :: ibc(6,NRINDEX)
!
!..miscellaneous
   integer :: ifc,iel,neig,ivar
!
#if HP3D_DEBUG
   integer :: iprint
   iprint=0
#endif
!
!-------------------------------------------------------------------------------
!
   Nelem_order = 0
!
!..check if exceeding maximum order
   if (IP.gt.MAXP) then
      write(*,*) 'set_initial_mesh: IP, MAXP = ', IP,MAXP
      stop 1
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
!
!     ...uniform order of approximation
         select case(ELEMS(iel)%etype)
            case(BRIC); Nelem_order(iel) = 111*IP
            case(PRIS); Nelem_order(iel) = 11*IP
            case(TETR); Nelem_order(iel) = 1*IP
            case(PYRA); Nelem_order(iel) = 1*IP
         end select
      else
!
!     ...custom order of approximation (NOT IMPLEMENTED)
         write(*,1003) IP
 1003    format('ERROR in set_initial_mesh: IP = ',i3)
         stop 1
      endif
!
!  ...set BC flags for each of 3 components:
!     0 - no BC ; 1 - Dirichlet ; 2 - Neumann
      ibc(1:6,1:NRINDEX) = 0
!
!  ...loop through the element faces
      do ifc=1,nface(ELEMS(iel)%etype)
         neig = ELEMS(iel)%neig(ifc)
!
!     ...no neighbor, set the BC flags
         if (neig .eq. 0) then
            select case(ifc)
               case(1) ! bottom
                  ibc(ifc,1) = 1; ibc(ifc,2) = 1; ibc(ifc,3) = 2
               case(2) ! top
                  ibc(ifc,1) = 2; ibc(ifc,2) = 2; ibc(ifc,3) = 1
               case(3,4,5,6)
                  ibc(ifc,1:3) = 2
            end select
         endif
      enddo
!
!  ...allocate BC flags (one per physical attribute component)
      allocate(ELEMS(iel)%bcond(NRINDEX))
!
!  ...encode face BCs into a single BC flag, one component at a time
      do ivar=1,NRINDEX
         call encodg(ibc(1:6,ivar),10,6, ELEMS(iel)%bcond(ivar))
#if HP3D_DEBUG
         if (iprint.eq.1) then
            write(*,*) 'ivar, ELEMS(iel)%bcond(ivar) = ', &
                        ivar,', ',ELEMS(iel)%bcond(ivar)
         endif
#endif
      enddo
!
#if HP3D_DEBUG
!  ...print order of approximation
      if (iprint.eq.1 .and. IP.gt.0) then
         write(*,*) '-- uniform order of approximation --'
         write(*,999) NRELIS
         select case(ELEMS(iel)%etype)
            case(BRIC); write(*,1002) IP,IP,IP
            case(PRIS); write(*,1001) IP,IP
            case(TETR); write(*,1000) IP
            case(PYRA); write(*,1000) IP
         end select
         write(*,*)
!
 999     format('Element# : ',i4)
1000     format(' p = ',       i3)
1001     format(' p, q = ',   2i3)
1002     format(' p, q, r = ',3i3)
      endif
#endif
!
!..end of loop over initial mesh elements
   enddo
!
!
end subroutine set_initial_mesh

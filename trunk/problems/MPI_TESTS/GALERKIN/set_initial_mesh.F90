!-------------------------------------------------------------------------------
!
!     routine name      - set_initial_mesh
!
!-------------------------------------------------------------------------------
!
!     latest revision:  - Jul 21
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
#if DEBUG_MODE
   integer :: iprint = 0
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
      ELEMS(iel)%nrphysics = 4
      allocate(ELEMS(iel)%physics(4))
      ELEMS(iel)%physics(1) ='ffldH'
      ELEMS(iel)%physics(2) ='ffldE'
      ELEMS(iel)%physics(3) ='ffldV'
      ELEMS(iel)%physics(4) ='ffldQ'
!
!  ...set order of approximation
      if (IP.gt.0) then
!
!     ...uniform order of approximation
         select case(ELEMS(iel)%Type)
            case('tetr'); Nelem_order(iel) = 1*IP
            case('pyra'); Nelem_order(iel) = 1*IP
            case('pris'); Nelem_order(iel) = 11*IP
            case('bric'); Nelem_order(iel) = 111*IP
         end select
      else
!
!     ...custom order of approximation (NOT IMPLEMENTED)
         write(*,1003) IP
 1003    format('ERROR in set_initial_mesh: IP = ',i3)
         stop 1
      endif
!
!  ...set BC flags 
!     0 - no BC 
      ibc(1:6,1:NRINDEX) = 0
!
!  ...allocate BC flags (one per physical attribute component)
      allocate(ELEMS(iel)%bcond(4))
!
!  ...encode face BCs into a single BC flag, one component at a time
      do ivar=1,4
        call encodg(ibc(1:6,ivar),10,6, ELEMS(iel)%bcond(ivar))
      enddo
!
#if DEBUG_MODE
!
!  ...print order of approximation
      if (iprint.eq.1 .and. IP.gt.0) then
         write(*,*) '-- uniform order of approximation --'
         write(*,999) NRELIS
         select case(ELEMS(iel)%Type)
            case('tetr'); write(*,1000) IP
            case('pyra'); write(*,1000) IP
            case('pris'); write(*,1001) IP,IP
            case('bric'); write(*,1002) IP,IP,IP
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

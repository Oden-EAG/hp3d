!------------------------------------------------------------------------------------
!> Purpose : define problem dependent data ( multiphysics, BC, approxiamtion )
!! @param[out] Nelem_order - order for initial mesh elements
!------------------------------------------------------------------------------------
!
subroutine set_initial_mesh(Nelem_order)
!
   use GMP
   use data_structure3D
   use common_prob_data, only : IBC_PROB,IP
!
   implicit none
!
   integer,dimension(NRELIS),intent(out) :: Nelem_order
!
!..BC flag
   integer, dimension(6,NRINDEX) :: ibc
   integer :: ifc,iel,neig,i
!
   integer :: iprint = 0
!
!------------------------------------------------------------------------------------
!
   iprint=0
!
!..check if have not exceeded the maximum order
   if (IP.gt.MAXP) then
      write(*,*) 'set_initial_mesh: IP, MAXP = ', IP,MAXP
      stop
   endif
!
!..loop over initial mesh elements
   do iel=1,NRELIS
!
!     STEP 1 : set up order of approximation
!
!  ...set order of approximation
      if (IP.gt.0) then
!     ...uniform order of approximation
         select case(ELEMS(iel)%etype)
            case(TETR); Nelem_order(iel) = 1*IP
            case(PYRA); Nelem_order(iel) = 1*IP
            case(PRIS); Nelem_order(iel) = 11*IP
            case(BRIC); Nelem_order(iel) = 111*IP
         end select
      else
!     ...custom order of approximation (NOT IMPLEMENTED)
         write(*,1003) IP
1003     format('ERROR in set_initial_mesh: IP = ',i3)
         stop
      endif
!
!  STEP 2 : set up physics
!
!   ...set physics
      ELEMS(iel)%nrphysics = NR_PHYSA
      allocate(ELEMS(iel)%physics(NR_PHYSA))
      do i=1,NR_PHYSA
         ELEMS(iel)%physics(i) = PHYSA(i)
      enddo
!
!  STEP 3 : set BC flags
!              0 - no BC
!              1 - Dirchlet
!              2 - Neumann
!              3 - Robin
!
      ibc(1:6,1:NRINDEX) = 0
!
!  ...physics attributes components (3+3+3+6 = 15 components)
!     phys|  comp   | description
!     (1) |  1 - 3  | H1 diplacement trace
!     (2) |  4 - 6  | Hdiv traction trace
!     (3) |  7 - 9  | L2 displacement field
!     (4) | 10 - 15 | L2 traction field
!
      select case(IBC_PROB)
      case(1)
!     ...if exterior face, set boundary condition to IBC_PROB
         do ifc=1,nface(ELEMS(iel)%etype)
            neig = ELEMS(iel)%neig(ifc)
            select case(neig)
            case(0)
               ibc(ifc,1) = 1  ! TrDis (H1) - x
               ibc(ifc,2) = 1  ! TrDis (H1) - y
               ibc(ifc,3) = 1  ! TrDis (H1) - z
               ibc(ifc,4) = 0  ! TrStr (Hdiv) - x
               ibc(ifc,5) = 0  ! TrStr (Hdiv) - y
               ibc(ifc,6) = 0  ! TrStr (Hdiv) - z
            end select
         enddo
      case(2)
!   ...if exterior face, set boundary condition (traction on top + bottom, displacement on sides)
         do ifc=1,nface(ELEMS(iel)%etype)
            neig = ELEMS(iel)%neig(ifc)
            select case(neig)
            case(0)
               if (ifc.eq.1) then
                  ibc(ifc,1) = 0  ! TrDis (H1) - x
                  ibc(ifc,2) = 0  ! TrDis (H1) - y
                  ibc(ifc,3) = 0  ! TrDis (H1) - z
                  ibc(ifc,4) = 1  ! TrStr (Hdiv) - x
                  ibc(ifc,5) = 1  ! TrStr (Hdiv) - y
                  ibc(ifc,6) = 1  ! TrStr (Hdiv) - z
               else
                  ibc(ifc,1) = 1  ! TrDis (H1) - x
                  ibc(ifc,2) = 1  ! TrDis (H1) - y
                  ibc(ifc,3) = 1  ! TrDis (H1) - z
                  ibc(ifc,4) = 0  ! TrStr (Hdiv) - x
                  ibc(ifc,5) = 0  ! TrStr (Hdiv) - y
                  ibc(ifc,6) = 0  ! TrStr (Hdiv) - z
               endif
            end select
         enddo
      case(3)
         do ifc=1,nface(ELEMS(iel)%etype)
            neig = ELEMS(iel)%neig(ifc)
            select case(neig)
! TODO: fix this
            case(0)
               if (ifc.eq.1) then
                  ibc(ifc,1) = 0  ! TrDis (H1) - x
                  ibc(ifc,2) = 0  ! TrDis (H1) - y
                  ibc(ifc,3) = 0  ! TrDis (H1) - z
                  ibc(ifc,4) = 1  ! TrStr (Hdiv) - x
                  ibc(ifc,5) = 1  ! TrStr (Hdiv) - y
                  ibc(ifc,6) = 1  ! TrStr (Hdiv) - z
               else
                  ibc(ifc,1) = 1  ! TrDis (H1) - x
                  ibc(ifc,2) = 1  ! TrDis (H1) - y
                  ibc(ifc,3) = 1  ! TrDis (H1) - z
                  ibc(ifc,4) = 0  ! TrStr (Hdiv) - x
                  ibc(ifc,5) = 0  ! TrStr (Hdiv) - y
                  ibc(ifc,6) = 0  ! TrStr (Hdiv) - z
               endif
            end select
         enddo
      end select
!
      allocate(ELEMS(iel)%bcond(NRINDEX))
      do i=1,NRINDEX
         call encodg(ibc(1:6,i),10,6, ELEMS(iel)%bcond(i))
      enddo
!
#if DEBUG_MODE
!  ...print order of approximation
      if (iprint.eq.1 .and. IP.gt.0) then
!     ...uniform order of approximation
         write(*,*) '-- uniform order of approximation --'
         write(*,999) NRELIS
         select case(ELEMS(iel)%etype)
         case(TETR); write(*,1000) IP
         case(PYRA); write(*,1000) IP
         case(PRIS); write(*,1001) IP,IP
         case(BRIC); write(*,1002) IP,IP,IP
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
end subroutine set_initial_mesh

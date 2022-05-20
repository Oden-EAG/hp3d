!-------------------------------------------------------------------------------
!
!     routine name      - set_initial_mesh
!
!-------------------------------------------------------------------------------
!
!     latest revision   - May 2022
!
!     purpose:          - define problem dependent data
!                         (multiphysics, BC, approximation)
!
!     arguments:
!            out:       - Nelem_order
!                         (polynomial order of each initial mesh element)
!
!-------------------------------------------------------------------------------
!
subroutine set_initial_mesh(Nelem_order)
   use GMP
   use data_structure3D
   use common_prob_data, only : IBC_PROB,IP
   use mpi_param, only: RANK,ROOT
!------------------------------------------------------------------------------------
   implicit none
!
!..polynomial order for initial mesh elements
   integer, intent(out) :: Nelem_order(NRELIS)
!
!..element type
   character(len=4) :: etype
!
!..BC flags; dimension = num_faces * components
   integer :: ibc(6,NRINDEX)
!
!..miscellaneous
   integer :: ifc,iel,neig,i,max_order
!
!------------------------------------------------------------------------------------
!     I N I T I A L I Z E
!------------------------------------------------------------------------------------
!
   if (NRELIS.ne.8) then
      write(*,*) "set_initial_mesh : WRONG GEOMETRY FILE"
      stop 1
   endif
!
!..check if exceeding maximum order
   if (IP.gt.MAXP) then
      write(*,*) 'set_initial_mesh: IP, MAXP = ', IP,MAXP
      stop
   endif
!
!------------------------------------------------------------------------------------
!     S E T    B C s
!------------------------------------------------------------------------------------
!
! loop over initial mesh elements
   do iel=1,NRELIS

      etype = ELEMS(iel)%type
!
!   SET ORDER OF APPROXIMATION
!
!  ...uniform order of approximation
      if (IP.gt.0) then
         select case(etype)
         case('tetr'); Nelem_order(iel) = 1*IP
         case('pyra'); Nelem_order(iel) = 1*IP
         case('pris'); Nelem_order(iel) = 11*IP
         case('bric'); Nelem_order(iel) = 111*IP
         end select
      else
         write(*,1003) IP
1003     format('ERROR in set_initial_mesh: IP = ',i3)
         stop
      endif
!
!   SET PHYSICS
!
    select case(iel)
!   Primal
    case(1,3,6,8)    
      ELEMS(iel)%nrphysics = 2
      allocate(ELEMS(iel)%physics(2))
      ELEMS(iel)%physics(1) = 'TrDis'
      ELEMS(iel)%physics(2) = 'TrStr'
!   Ultra-weak
    case(2,4,5,7)
      ELEMS(iel)%nrphysics = 5
      allocate(ELEMS(iel)%physics(5))
      ELEMS(iel)%physics(1) = 'TrDis'
      ELEMS(iel)%physics(2) = 'TrStr'
      ELEMS(iel)%physics(3) = 'Displ'
      ELEMS(iel)%physics(4) = 'Stres'
      ELEMS(iel)%physics(5) = 'Omega'
    end select
!
!   SET BC FLAGS: 0 - no BC ; 1 - Dirchlet ; 2 - Neumann ; 3 - Robin ; >3 - Mixed
!
   ibc(1:6,1:NR_PHYSA) = 0
!
!   IBC_PROB : 0 - uniform traction ; 1 - clamped ends ; 2 - free ends ; 3 - periodic ends
!       ___ ___ _____________________________
!      /8__|__6\                             \
!     / /7_|_5\ \                             \
!    /_/_/   \_\_\  _  _  _  _  _  _  _  _  _  \
!    \ \4\_ _/2/ /                             /
!     \3\__|__/1/                             /
!      \___|___/_____________________________/
!
!        elements 1,3,6,8 : E = 1000, nu = 0.2
!        elements 2,4,5,7 : E = 1, nu = 0.5
!
    select case(IBC_PROB)
!
!   Prescribed displacement over entire boundary
    case(0)
      do ifc=1,nface(etype)
         neig = ELEMS(iel)%neig(ifc)
         select case(neig)
         case(0)
            ibc(ifc,1) = 1  ! 1st H1 component
            ibc(ifc,2) = 1  ! 2nd H1 component
            ibc(ifc,3) = 1  ! 3rd H1 component
            ! ! ibcflag      -> physics variable
            ! ibc(ifc,1) = 1  ! TrDis (H1)
            ! ibc(ifc,2) = 0  ! TrStr (H(div))
            ! ibc(ifc,3) = 0  ! Displ (L2)
            ! ibc(ifc,4) = 0  ! Stres (L2)
            ! ibc(ifc,5) = 0  ! Omega (L2)
            ! ! Note that L2 variables should not take values at the boundary,
            ! ! so their ibc should always be 0 (no BC).
         end select
      enddo
!
!   Prescribed displacement at ends of pipe
!   Pressure BCs at exterior faces of pipe
    case(1)

      ! exterior faces
      select case(iel)
      case(1,3,5,7); ifc=1
      case(2,4,6,8); ifc=2
      end select

      ibc(ifc,1) = 0  ! TrDis (H1)
      ibc(ifc,2) = 1  ! TrStr (H(div))
      ibc(ifc,3) = 0  ! Displ (L2)
      ibc(ifc,4) = 0  ! Stres (L2)
      ibc(ifc,5) = 0  ! Omega (L2)

      ! fixed ends
      do ifc=3,6
        neig = ELEMS(iel)%neig(ifc)
        select case(neig)
        case(0)
          ibc(ifc,1) = 1  ! TrDis (H1)
          ibc(ifc,2) = 0  ! TrStr (H(div))
          ibc(ifc,3) = 0  ! Displ (L2)
          ibc(ifc,4) = 0  ! Stres (L2)
          ibc(ifc,5) = 0  ! Omega (L2)
        end select
      enddo
!
!   Pressure BCs at exterior faces of pipe
    case(2)
      do ifc=1,nface(etype)
          neig = ELEMS(iel)%neig(ifc)
          select case(neig)
          case(0)
            ! ibcflag      -> physics variable
            ibc(ifc,1) = 0  ! TrDis (H1)
            ibc(ifc,2) = 1  ! TrStr (H(div))
            ibc(ifc,3) = 0  ! Displ (L2)
            ibc(ifc,4) = 0  ! Stres (L2)
            ibc(ifc,5) = 0  ! Omega (L2)
            ! Note that L2 variables should not take values at the boundary,
            ! so their ibc should always be 0 (no BC).
          end select
        enddo

!
!   Pressure BCs at exterior faces of pipe
!   In this case, it is enough to incorporate zero displacement conditions in the x-direction
!   so long as the pressure is symmetric in x
   !  case(3)

   !    ! exterior faces
   !    select case(iel)
   !    case(1,3,5,7); ifc=1
   !    case(2,4,6,8); ifc=2
   !    end select

   !    ibc(ifc,1) = 0  ! TrDis (H1)
   !    ibc(ifc,2) = 1  ! TrStr (H(div))
   !    ibc(ifc,3) = 0  ! Displ (L2)
   !    ibc(ifc,4) = 0  ! Stres (L2)
   !    ibc(ifc,5) = 0  ! Omega (L2)

   !    ! periodic ends
   !    do ifc=3,6
   !      neig = ELEMS(iel)%neig(ifc)
   !      select case(neig)
   !      case(0)
   !        ibc(ifc,1) = 6  ! TrDis (H1)
   !        ibc(ifc,2) = 6  ! TrStr (H(div))
   !        ibc(ifc,3) = 0  ! Displ (L2)
   !        ibc(ifc,4) = 0  ! Stres (L2)
   !        ibc(ifc,5) = 0  ! Omega (L2)
   !      end select
   !    enddo

      end select
!
!  ...allocate BC flags (one per attribute component)
      allocate(ELEMS(iel)%bcond(NRINDEX))
!  ...for each component, encode face BC into a single BC flag
      do i=1,NRINDEX
         call encodg(ibc(1:6,i),10,6, ELEMS(iel)%bcond(i))
      enddo
!
!..end of loop over initial mesh elements
   enddo
!
end subroutine set_initial_mesh

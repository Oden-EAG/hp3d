!------------------------------------------------------------------------------------
!> Purpose : define problem dependent data ( multiphysics, BC, approxiamtion )
!! @param[out] Nelem_order - order for initial mesh elements
!------------------------------------------------------------------------------------
!
subroutine set_initial_mesh(Nelem_order)
  use GMP
  use data_structure3D
  use common_prob_data, only : IBC_PROB,IP
!------------------------------------------------------------------------------------
  implicit none
  integer,dimension(NRELIS),intent(out) :: Nelem_order
!------------------------------------------------------------------------------------
! BC flags
  integer, dimension(6,NR_PHYSA) :: ibc
! miscellanea
  integer :: iprint,ifc,iel,neig,iat
!
!------------------------------------------------------------------------------------
!     I N I T I A L I Z E
!------------------------------------------------------------------------------------
!
  iprint=0
!
  if (NRELIS.ne.8) then
    write(*,*) "set_initial_mesh : WRONG GEOMETRY FILE"
    stop 1
  endif
!
! check if have not exceeded the maximum order
  if (IP.gt.MAXP) then
    write(*,*) 'set_initial_mesh: IP, MAXP = ', IP,MAXP
    stop
  endif
!
! add flag `6' to the list of dirichlet flags
  call add_dirichlet_to_list(6)
!
!------------------------------------------------------------------------------------
!     S E T    B C s
!------------------------------------------------------------------------------------
!
! loop over initial mesh elements
  do iel=1,NRELIS
!
!   SET ORDER OF APPROXIMATION
!
    if (IP.gt.0) then
!     uniform order of approximation
      select case(ELEMS(iel)%Type)
      case('tetr'); Nelem_order(iel) = 1*IP
      case('pyra'); Nelem_order(iel) = 1*IP
      case('pris'); Nelem_order(iel) = 11*IP
      case('bric'); Nelem_order(iel) = 111*IP
      end select
    else
!     custom order of approximation (NOT IMPLEMENTED)
      write(*,1003) IP
1003  format('ERROR in set_initial_mesh: IP = ',i3)
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
      do ifc=1,nface(ELEMS(iel)%Type)
        neig = ELEMS(iel)%neig(ifc)
        select case(neig)
        case(0)
          ! ibcflag      -> physics variable
          ibc(ifc,1) = 1  ! TrDis (H1)
          ibc(ifc,2) = 0  ! TrStr (H(div))
          ibc(ifc,3) = 0  ! Displ (L2)
          ibc(ifc,4) = 0  ! Stres (L2)
          ibc(ifc,5) = 0  ! Omega (L2)
          ! Note that L2 variables should not take values at the boundary,
          ! so their ibc should always be 0 (no BC).
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
      do ifc=1,nface(ELEMS(iel)%Type)
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
    case(3)

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

      ! periodic ends
      do ifc=3,6
        neig = ELEMS(iel)%neig(ifc)
        select case(neig)
        case(0)
          ibc(ifc,1) = 6  ! TrDis (H1)
          ibc(ifc,2) = 6  ! TrStr (H(div))
          ibc(ifc,3) = 0  ! Displ (L2)
          ibc(ifc,4) = 0  ! Stres (L2)
          ibc(ifc,5) = 0  ! Omega (L2)
        end select
      enddo

    end select
!
!   allocate BC flags (one per attribute)
    allocate(ELEMS(iel)%bcond(NR_PHYSA))
!   for each attribute, encode face BC into a single BC flag
    do iat=1,NR_PHYSA
      call encodg(ibc(1:6,iat),10,6, ELEMS(iel)%bcond(iat))
    enddo
!
!-----------------------------------------------------------------------------------
!     P R I N T    S T A T E M E N T S
!-----------------------------------------------------------------------------------
!
!   print order of approximation
    if (iprint.eq.1 .and. IP.gt.0) then
      ! uniform order of approximation
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
 999  format('Element# : ',i4)
1000  format(' p = ',       i3)
1001  format(' p, q = ',   2i3)
1002  format(' p, q, r = ',3i3)
    endif
!
! end of loop over initial mesh elements
  enddo
!
end subroutine set_initial_mesh

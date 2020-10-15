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
!  ...check if have not exceeded the maximum order
      if (IP.gt.MAXP) then
        write(*,*) 'set_initial_mesh: IP, MAXP = ', IP,MAXP
        stop
      endif
!
!  ...loop through initial mesh elements
      if (NRELIS.ne.1) then
        write(*,*) 'set_initial_mesh: WRONG NRELIS = ',NRELIS
        stop
      endif

!
! add flag `8' to the list of dirichlet flags
  call add_dirichlet_to_list(8)
!
!------------------------------------------------------------------------------------
!     S E T    B C s
!------------------------------------------------------------------------------------
!
! loop over initial mesh elements
  do iel=1,NRELIS
!
!   set physics
    ELEMS(iel)%nrphysics = NR_PHYSA
    allocate(ELEMS(iel)%physics(NR_PHYSA))
    do iat=1,NR_PHYSA
      ELEMS(iel)%physics(iat) = PHYSA(iat)
    enddo
!
!   set order of approximation
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
!   set BC flags: 0 - no BC ; 1 - Dirchlet ; 2 - Neumann ; 3 - Robin ; >3 - Mixed
!
    ibc(1:6,1:NR_PHYSA) = 0
    select case(IBC_PROB)
!
!   U N I F O R M     B C s
!
    case(1)
      ! if exterior face, set boundary condition to IBC_PROB
      do ifc=1,nface(ELEMS(iel)%Type)
        neig = ELEMS(iel)%neig(ifc)
        select case(neig)
        case(0)
          ibc(ifc,1) = 1
          ibc(ifc,2) = 0
          ibc(ifc,3) = 0
          ibc(ifc,4) = 0
          ibc(ifc,5) = 0
        end select
      enddo
!
!   M I X E D   B C s   (and not uniform)
!
!  NEEDS TO BE IMPROVED!!!
    case(8)
      ! if exterior face, set boundary condition (Mixed on top + bottom, Dirichlet on sides)
      do ifc=1,nface(ELEMS(iel)%Type)
        neig = ELEMS(iel)%neig(ifc)
        select case(neig)
        case(0)
          if (ifc.lt.3) then
            ibc(ifc,1) = IBC_PROB
            ibc(ifc,2) = 0
          else
            ibc(ifc,1) = 1
            ibc(ifc,1) = 0
          endif
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

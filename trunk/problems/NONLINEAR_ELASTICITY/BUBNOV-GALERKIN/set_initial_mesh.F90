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
  integer, dimension(6) :: ibc
! miscellanea
  integer :: iprint,i,iel,neig,iDirichlet
!------------------------------------------------------------------------------------
!
  iprint=0
!




! add flag `8' to the list of dirichlet flags
  call add_dirichlet_to_list(8)
!
! loop over initial mesh elements
  do iel=1,NRELIS
  !
  ! set physics
    ELEMS(iel)%nrphysics = 1
    allocate(ELEMS(iel)%physics(1))
    ELEMS(iel)%physics(1) = 'Displ'
  !
  ! set order of approximation
    if (IP.gt.0) then
      ! uniform order of approximation
      select case(ELEMS(iel)%Type)
      case('tetr'); Nelem_order(iel) = 1*IP
      case('pyra'); Nelem_order(iel) = 1*IP
      case('pris'); Nelem_order(iel) = 11*IP
      case('bric'); Nelem_order(iel) = 111*IP
      end select
    else
      ! custom order of approximation (NOT IMPLEMENTED)
      write(*,1003) IP
1003  format('ERROR in set_initial_mesh: IP = ',i3)
      stop
    endif
  !
  !   set BC flags: 0 - no BC ; 1 - Dirchlet ; 2 - Neumann ; 3 - Robin ; >3 - Mixed
  !
    ibc(1:6) = 0
    select case(IBC_PROB)
    !----------------
    ! UNIFORM BCs
    !----------------
    case(1,2,3)
      ! if exterior face, set boundary condition to IBC_PROB
      do i=1,nface(ELEMS(iel)%Type)
        neig = ELEMS(iel)%neig(i)
        select case(neig)
        case(0); ibc(i) = IBC_PROB
        end select
      enddo
    !----------------
    ! MIXED BCs (and not uniform)
    !----------------
    case(8)
      iDirichlet = 1
      if (iDirichlet.eq.1 .or. ELEMS(iel)%Type.ne.'bric' .or. NRELIS.ne.3) then
        ! If exterior face, set boundary condition (Mixed on top + bottom, Dirichlet all 6 sides)
        !
        !   ___D___
        !  | 1 : 2 |    {1,2,3} : element number (iel)
        ! D|___:...|    D = Dirichlet (ibc = 1)
        !    D | 3 |D   N = Neumann   (ibc = 2)
        !     D|___|
        !        D
        !
        ! OR just set Dirichlet boundary conditions if the wrong geometry file is loaded
        do i=1, nface(ELEMS(iel)%Type)
          neig = ELEMS(iel)%neig(i)
          select case(neig)
          case(0)
            if (i.lt.3) then
              ibc(i) = IBC_PROB
            else
              ibc(i) = 1
            endif
          end select
        enddo
      else
        ! If exterior face, set boundary condition (Mixed on top + bottom, Dirichlet on opening sides, Neumann on other 4 sides)
        !
        !   ___N___
        !  | 1 : 2 |    {1,2,3} : element number (iel)
        ! N|___:...|    D = Dirichlet (ibc = 1)
        !    D | 3 |N   N = Neumann   (ibc = 2)
        !     D|___|
        !        N
        !
        ibc(1) = 8
        ibc(2) = 8
        select case(iel)
        case(1)
          ibc(3) = 2
          ibc(4) = 1
          ibc(5) = 0
          ibc(6) = 2
        case(2)
          ibc(3) = 0
          ibc(4) = 0
          ibc(5) = 2
          ibc(6) = 2
        case(3)
          ibc(3) = 1
          ibc(4) = 2
          ibc(5) = 2
          ibc(6) = 0
        end select
      endif
    end select
    !
    !  store the boundary conditions
    allocate(ELEMS(iel)%bcond(NR_PHYSA))
    call encodg(ibc,10,6, ELEMS(iel)%bcond(1))
    !
    !  print order of approximation
    if (iprint.eq.1 .and. IP.gt.0) then
      !  uniform order of approximation
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

  !  end of loop over initial mesh elements
  enddo

end subroutine set_initial_mesh

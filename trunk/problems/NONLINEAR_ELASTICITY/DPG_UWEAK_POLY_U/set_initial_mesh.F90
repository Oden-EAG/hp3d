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
  integer :: iprint,ifc,iel,neig,iat,iDisplacement
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
          ! ibcflag      -> physics variable
          ibc(ifc,1) = 1  ! TrDis (H1)
          ibc(ifc,2) = 0  ! TrStr (H(div))
          ibc(ifc,3) = 0  ! Displ (L2)
          ibc(ifc,4) = 0  ! Stres (L2)
          ! Note that L2 variables should not take values at the boundary,
          ! so their ibc should always be 0 (no BC).
        end select
      enddo
    case(2)
      ! if exterior face, set boundary condition (traction on top + bottom, displacement on sides)
      do ifc=1,nface(ELEMS(iel)%Type)
        neig = ELEMS(iel)%neig(ifc)
        select case(neig)
        case(0)
          if (ifc.lt.3) then
            ! ibcflag      -> physics variable
            ibc(ifc,1) = 0  ! TrDis (H1)
            ibc(ifc,2) = 1  ! TrStr (H(div))
            ibc(ifc,3) = 0  ! Displ (L2)
            ibc(ifc,4) = 0  ! Stres (L2)
            ! Note that L2 variables should not take values at the boundary,
            ! so their ibc should always be 0 (no BC).
          else
            ! ibcflag      -> physics variable
            ibc(ifc,1) = 1  ! TrDis (H1)
            ibc(ifc,2) = 0  ! TrStr (H(div))
            ibc(ifc,3) = 0  ! Displ (L2)
            ibc(ifc,4) = 0  ! Stres (L2)
            ! Note that L2 variables should not take values at the boundary,
            ! so their ibc should always be 0 (no BC).
          endif
        end select
      enddo
    case(3)
      if (ELEMS(iel)%Type.ne.'bric' .or. NRELIS.ne.2) then
        write(*,*) ' set_initial_mesh : UNSUPPORTED GEOMETRY'
        stop
      endif
      ! set all boundary conditions to traction
      do ifc=1,nface(ELEMS(iel)%Type)
        neig = ELEMS(iel)%neig(ifc)
        select case(neig)
        case(0)
          ! ibcflag      -> physics variable
          ibc(ifc,1) = 0  ! TrDis (H1)
          ibc(ifc,2) = 1  ! TrStr (H(div))
          ibc(ifc,3) = 0  ! Displ (L2)
          ibc(ifc,4) = 0  ! Stres (L2)
        end select
      enddo
      !  reset BCs on bottom face
      if (iel.eq.1) ibc(1,1:4) = (/1,0,0,0/)
!
!   M I X E D     B C s    (and not uniform)
!
    case(8)
      iDisplacement = 0
      if (iDisplacement.eq.1 .or. ELEMS(iel)%Type.ne.'bric' .or. NRELIS.ne.3) then
        ! If exterior face, set boundary condition (Mixed on top + bottom, Dirichlet all 6 sides)
        !
        !   ___D___
        !  | 1 : 2 |    {1,2,3} : element number (iel)
        ! D|___:...|    D = Displacement (ibc = (1,0))
        !    D | 3 |D   T = Traction     (ibc = (0,1))
        !     D|___|
        !        D
        !
        ! OR just set Displacement boundary conditions if the wrong geometry file is loaded
        do ifc=1, nface(ELEMS(iel)%Type)
          neig = ELEMS(iel)%neig(ifc)
          select case(neig)
          case(0)
            if (ifc.lt.3) then
              ! ibcflag      -> physics variable
              ibc(ifc,1) = 8  ! TrDis (H1)
              ibc(ifc,2) = 8  ! TrStr (H(div))
              ibc(ifc,3) = 0  ! Displ (L2)
              ibc(ifc,4) = 0  ! Stres (L2)
              ! Note that L2 variables should not take values at the boundary,
              ! so their ibc should always be 0 (no BC).
            else
              ! ibcflag      -> physics variable
              ibc(ifc,1) = 1  ! TrDis (H1)
              ibc(ifc,2) = 0  ! TrStr (H(div))
              ibc(ifc,3) = 0  ! Displ (L2)
              ibc(ifc,4) = 0  ! Stres (L2)
              ! Note that L2 variables should not take values at the boundary,
              ! so their ibc should always be 0 (no BC).
            endif
          end select
        enddo
      else
        ! If exterior face, set boundary condition (Mixed on top + bottom, Displacement on opening sides, Traction on other 4 sides)
        !
        !   ___T___
        !  | 1 : 2 |    {1,2,3} : element number (iel)
        ! T|___:...|    D = Displacement (ibc = (1,0))
        !    D | 3 |T   T = Traction     (ibc = (0,1))
        !     D|___|
        !        T
        !
        ! Regardless of the element, top and bottom are mixed
        ! Element=1,2,3 - Face 1 - Mixed (this is assuming a very specific geometry file)
        ! ibcflag      -> physics variable
        ifc=1
        ibc(ifc,1) = 8  ! TrDis (H1)
        ibc(ifc,2) = 8  ! TrStr (H(div))
        ibc(ifc,3) = 0  ! Displ (L2)
        ibc(ifc,4) = 0  ! Stres (L2)
        ! Note that L2 variables should not take values at the boundary,
        ! so their ibc should always be 0 (no BC).
        ! Element=1,2,3 - Face 2 - Mixed (this is assuming a very specific geometry file)
        ! ibcflag      -> physics variable
        ifc=2
        ibc(ifc,1) = 8  ! TrDis (H1)
        ibc(ifc,2) = 8  ! TrStr (H(div))
        ibc(ifc,3) = 0  ! Displ (L2)
        ibc(ifc,4) = 0  ! Stres (L2)
        ! Note that L2 variables should not take values at the boundary,
        ! so their ibc should always be 0 (no BC).
        !
        ! Then depending on the element (provided the very specific geometry file)
        ! each face has different BC as pictured above
        select case(iel)
        case(1)
          ! Element=1 - Face 3 - Traction (this is assuming a very specific geometry file)
          ! ibcflag      -> physics variable
          ifc=3
          ibc(ifc,1) = 0  ! TrDis (H1)
          ibc(ifc,2) = 1  ! TrStr (H(div))
          ibc(ifc,3) = 0  ! Displ (L2)
          ibc(ifc,4) = 0  ! Stres (L2)
          ! Note that L2 variables should not take values at the boundary,
          ! so their ibc should always be 0 (no BC).
          ! Element=1 - Face 4 - Displacement (this is assuming a very specific geometry file)
          ! ibcflag    -> physics variable
          ifc=4
          ibc(ifc,1) = 1  ! TrDis (H1)
          ibc(ifc,2) = 0  ! TrStr (H(div))
          ibc(ifc,3) = 0  ! Displ (L2)
          ibc(ifc,4) = 0  ! Stres (L2)
          ! Note that L2 variables should not take values at the boundary,
          ! so their ibc should always be 0 (no BC).
          ! Element=1 - Face 5 - Inner face (this is assuming a very specific geometry file)
          ! ibcflag    -> physics variable
          ifc=5
          ibc(ifc,1) = 0  ! TrDis (H1)
          ibc(ifc,2) = 0  ! TrStr (H(div))
          ibc(ifc,3) = 0  ! Displ (L2)
          ibc(ifc,4) = 0  ! Stres (L2)
          ! Note that L2 variables should not take values at the boundary,
          ! so their ibc should always be 0 (no BC).
          ! Element=1 - Face 6 - Traction (this is assuming a very specific geometry file)
          ! ibcflag    -> physics variable
          ifc=6
          ibc(ifc,1) = 0  ! TrDis (H1)
          ibc(ifc,2) = 1  ! TrStr (H(div))
          ibc(ifc,3) = 0  ! Displ (L2)
          ibc(ifc,4) = 0  ! Stres (L2)
          ! Note that L2 variables should not take values at the boundary,
          ! so their ibc should always be 0 (no BC).
        case(2)
          ! Element=2 - Face 3 - Inner face (this is assuming a very specific geometry file)
          ! ibcflag      -> physics variable
          ifc=3
          ibc(ifc,1) = 0  ! TrDis (H1)
          ibc(ifc,2) = 0  ! TrStr (H(div))
          ibc(ifc,3) = 0  ! Displ (L2)
          ibc(ifc,4) = 0  ! Stres (L2)
          ! Note that L2 variables should not take values at the boundary,
          ! so their ibc should always be 0 (no BC).
          ! Element=2 - Face 4 - Inner face (this is assuming a very specific geometry file)
          ! ibcflag    -> physics variable
          ifc=4
          ibc(ifc,1) = 0  ! TrDis (H1)
          ibc(ifc,2) = 0  ! TrStr (H(div))
          ibc(ifc,3) = 0  ! Displ (L2)
          ibc(ifc,4) = 0  ! Stres (L2)
          ! Note that L2 variables should not take values at the boundary,
          ! so their ibc should always be 0 (no BC).
          ! Element=2 - Face 5 - Traction (this is assuming a very specific geometry file)
          ! ibcflag    -> physics variable
          ifc=5
          ibc(ifc,1) = 0  ! TrDis (H1)
          ibc(ifc,2) = 1  ! TrStr (H(div))
          ibc(ifc,3) = 0  ! Displ (L2)
          ibc(ifc,4) = 0  ! Stres (L2)
          ! Note that L2 variables should not take values at the boundary,
          ! so their ibc should always be 0 (no BC).
          ! Element=2 - Face 6 - Traction (this is assuming a very specific geometry file)
          ! ibcflag    -> physics variable
          ifc=6
          ibc(ifc,1) = 0  ! TrDis (H1)
          ibc(ifc,2) = 1  ! TrStr (H(div))
          ibc(ifc,3) = 0  ! Displ (L2)
          ibc(ifc,4) = 0  ! Stres (L2)
          ! Note that L2 variables should not take values at the boundary,
          ! so their ibc should always be 0 (no BC).
        case(3)
          ! Element=3 - Face 3 - Displacement (this is assuming a very specific geometry file)
          ! ibcflag      -> physics variable
          ifc=3
          ibc(ifc,1) = 1  ! TrDis (H1)
          ibc(ifc,2) = 0  ! TrStr (H(div))
          ibc(ifc,3) = 0  ! Displ (L2)
          ibc(ifc,4) = 0  ! Stres (L2)
          ! Note that L2 variables should not take values at the boundary,
          ! so their ibc should always be 0 (no BC).
          ! Element=3 - Face 4 - Traction (this is assuming a very specific geometry file)
          ! ibcflag    -> physics variable
          ifc=4
          ibc(ifc,1) = 0  ! TrDis (H1)
          ibc(ifc,2) = 1  ! TrStr (H(div))
          ibc(ifc,3) = 0  ! Displ (L2)
          ibc(ifc,4) = 0  ! Stres (L2)
          ! Note that L2 variables should not take values at the boundary,
          ! so their ibc should always be 0 (no BC).
          ! Element=3 - Face 5 - Traction (this is assuming a very specific geometry file)
          ! ibcflag    -> physics variable
          ifc=5
          ibc(ifc,1) = 0  ! TrDis (H1)
          ibc(ifc,2) = 1  ! TrStr (H(div))
          ibc(ifc,3) = 0  ! Displ (L2)
          ibc(ifc,4) = 0  ! Stres (L2)
          ! Note that L2 variables should not take values at the boundary,
          ! so their ibc should always be 0 (no BC).
          ! Element=3 - Face 6 - Inner face (this is assuming a very specific geometry file)
          ! ibcflag    -> physics variable
          ifc=6
          ibc(ifc,1) = 0  ! TrDis (H1)
          ibc(ifc,2) = 0  ! TrStr (H(div))
          ibc(ifc,3) = 0  ! Displ (L2)
          ibc(ifc,4) = 0  ! Stres (L2)
          ! Note that L2 variables should not take values at the boundary,
          ! so their ibc should always be 0 (no BC).
        end select
      endif
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

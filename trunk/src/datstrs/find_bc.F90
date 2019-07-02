!--------------------------------------------------------------------------------
!> Purpose            - routine determines BC flags for faces of an element
!!
!! @param[in]  Mdle   - middle node
!! @param[out] Ibc    - BC
!--------------------------------------------------------------------------------
!
subroutine find_bc(Mdle, Ibc)
!
  use element_data
  use data_structure3D
!--------------------------------------------------------------------------------
  implicit none
  ! ** Arguments
  integer, intent(in)    :: Mdle
  integer, intent(out)   :: Ibc(6,NR_PHYSA)
  ! ** Locals
  integer, dimension(27) :: nodesl, norientl
  integer :: i, j, nve, iprint
!--------------------------------------------------------------------------------
!
  iprint=0
  Ibc(6,NR_PHYSA) = 0

  call elem_nodes(Mdle, nodesl,norientl)
!
  nve = nvert(NODES(Mdle)%type) + nedge(NODES(Mdle)%type)
!
! printing
  if (iprint.eq.1) then
    write(*,9999)Mdle
9999  format(' find_bc: Mdle = ',i7)
  endif
!
! loop over faces of Mdle and copy BC from face nodes
  do i=1,nface(NODES(Mdle)%type)
    call decod(NODES(nodesl(nve+i))%bcond,10,NR_PHYSA, Ibc(i,1:NR_PHYSA))
!!!     Ibc(i) = NODES(nodesl(nve+i))%bcond
!
!   printing
    if (iprint.eq.1) then
    write(*,9998)i,Ibc(i,1:NR_PHYSA)
9998  format('   i,Ibc = ',i1,4x,5(i1,2x))
    endif
!
  enddo
!
!
end subroutine find_bc

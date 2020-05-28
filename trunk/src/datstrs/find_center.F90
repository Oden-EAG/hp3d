!> Purpose : find out center coordinate of middle node
!! @param[in]  Mdle - middle node
!! @param[out] X    - center coordinate
subroutine find_center(Mdle, X)
  use data_structure3D
  implicit none
  integer,               intent(in)  :: Mdle
  real(8), dimension(3), intent(out) :: X

  integer, dimension(27) :: nodesl,norientl
  integer                :: i, nv

  call elem_nodes(Mdle, nodesl, norientl)
  nv = nvert(NODES(Mdle)%Type)

  X(1:3) = 0.d0
  do i=1,nv
     X(1:3) = X(1:3) + NODES(nodesl(i))%dof%coord(1:3,1)
  enddo
  X(1:3) = X(1:3)/nv

end subroutine find_center

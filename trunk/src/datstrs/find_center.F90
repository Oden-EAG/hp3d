!> @brief      find out center coordinate of middle node
!> @param[in]  Mdle - middle node
!> @param[out] X    - center coordinate
!> @date       Feb 2023
subroutine find_center(Mdle, X)
   use data_structure3D
   implicit none
   integer, intent(in)  :: Mdle
   real(8), intent(out) :: X(3)
!
   integer :: nodesl(27),norientl(27)
   integer :: i, nv
!
   call elem_nodes(Mdle, nodesl,norientl)
   nv = nvert(NODES(Mdle)%ntype)
!
   X(1:3) = 0.d0
   do i=1,nv
      X(1:3) = X(1:3) + NODES(nodesl(i))%dof%coord(1:3,1)
   enddo
   X(1:3) = X(1:3)/nv
!
end subroutine find_center

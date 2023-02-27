!> @brief      find order of approximation for Mdle;
!!             remark: routines calls elem_nodes
!> @param[in]  Mdle         - middle node number
!> @param[out] Norder       - order of approximation
!> @param[out] Nedge_orient - edge orientation
!> @param[out] Nface_orient - face orientation
!> @date       Feb 2023
subroutine find_elem_nodes(Mdle, Norder,Nedge_orient,Nface_orient)
   use data_structure3D
   implicit none
   integer, intent(in)  :: Mdle
   integer, intent(out) :: Norder(19), Nedge_orient(12), Nface_orient(6)
!
   integer :: nodesl(27), norientl(27)
   integer :: ntype
!
!..initialize
   Nedge_orient = 0; Nface_orient = 0
!
   call elem_nodes(Mdle, nodesl,norientl)
!
   ntype = NODES(Mdle)%ntype
   call find_orient_from_list(ntype,norientl, Nedge_orient,Nface_orient)
   call find_order_from_list(ntype,nodesl, Norder)
!
end subroutine find_elem_nodes

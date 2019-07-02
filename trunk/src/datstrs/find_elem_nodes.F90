!> Purpose : find order of approximation for Mdle. It calls elem_nodes
!! @param[in]  Mdle   - middle node number
!! @param[out] Norder - order of approximation 
subroutine find_elem_nodes(Mdle, Norder, Nedge_orient,Nface_orient)
  use data_structure3D
  implicit none
  integer, intent(in)    :: Mdle
  integer, intent(out)   :: Norder(19), Nedge_orient(12), Nface_orient(6)

  integer, dimension(27) :: nodesl, norientl
  character(len=4)       :: type

  call elem_nodes(Mdle, nodesl, norientl)

  type = NODES(Mdle)%Type
  call find_orient_from_list(type, norientl, nedge_orient,nface_orient)
  call find_order_from_list(type, nodesl, norder)

end subroutine find_elem_nodes

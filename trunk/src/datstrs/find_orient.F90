!> Purpose : compute the orientation with a given Mdle. It calls elem_nodes
!! @param[in]  Mdle         - middle node number
!! @param[out] Nedge_orient - edge orientation
!! @param[out] Nface_orient - face orientation
subroutine find_orient(Mdle, Nedge_orient,Nface_orient)
  use data_structure3D
  implicit none
  integer, intent(in)    :: Mdle
  integer, intent(out)   :: Nedge_orient(12), Nface_orient(6)
  integer, dimension(27) :: nodesl, norientl

  call elem_nodes(Mdle, nodesl,norientl)
  call find_orient_from_list(NODES(Mdle)%Type, norientl, Nedge_orient, Nface_orient)
end subroutine find_orient

!> Purpose : compute the orientation from list
!! @param[in]  Type                       - middle node type
!! @param[in]  Nodesl, Norientl           - nodal connectivity
!! @param[out] Nedge_orient, Nface_orient - edge and face orientation
subroutine find_orient_from_list(Type, Norientl, Nedge_orient, Nface_orient)
  use data_structure3D
  implicit none
  character(len=4), intent(in)  :: Type
  integer,          intent(in)  :: Norientl(27) 
  integer,          intent(out) :: Nedge_orient(12), Nface_orient(6)
  integer :: nrv, nre, nrf

  nrv = NVERT(Type)
  nre = NEDGE(Type)
  nrf = NFACE(Type)

  Nedge_orient(1:nre) = Norientl(nrv+1:nrv+nre)
  Nface_orient(1:nrf) = Norientl(nrv+nre+1:nrv+nre+nrf)
end subroutine find_orient_from_list


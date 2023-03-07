!----------------------------------------------------------------------
!> @brief      compute the orientation with a given Mdle;
!!             remark: this routine calls elem_nodes
!! @param[in]  Mdle         - middle node number
!! @param[out] Nedge_orient - edge orientation
!! @param[out] Nface_orient - face orientation
!> @date       Feb 2023
subroutine find_orient(Mdle, Nedge_orient,Nface_orient)
   use data_structure3D
   implicit none
   integer, intent(in)    :: Mdle
   integer, intent(out)   :: Nedge_orient(12), Nface_orient(6)
   integer                :: nodesl(27), norientl(27)
!
   call elem_nodes(Mdle, nodesl,norientl)
   call find_orient_from_list(NODES(Mdle)%ntype,norientl, Nedge_orient,Nface_orient)
!
end subroutine find_orient
!
!----------------------------------------------------------------------
!> @brief      compute the orientation from list
!! @param[in]  Ntype         - middle node type
!! @param[in]  Nodesl        - nodal connectivity (node list)
!! @param[in]  Norientl      - nodal connectivity (orientations)
!! @param[out] Nedge_orient  - edge orientation
!! @param[out] Nface_orient  - face orientation
!> @date       Feb 2023
subroutine find_orient_from_list(Ntype,Norientl, Nedge_orient,Nface_orient)
   use data_structure3D
   implicit none
   integer, intent(in)  :: Ntype, Norientl(27)
   integer, intent(out) :: Nedge_orient(12), Nface_orient(6)
   integer :: nrv, nre, nrf
!
   Nedge_orient(1:12) = 0
   Nface_orient(1: 6) = 0
!
   nrv = NVERT(Ntype)
   nre = NEDGE(Ntype)
   nrf = NFACE(Ntype)
!
   Nedge_orient(1:nre) = Norientl(nrv+1:nrv+nre)
   Nface_orient(1:nrf) = Norientl(nrv+nre+1:nrv+nre+nrf)
!
end subroutine find_orient_from_list


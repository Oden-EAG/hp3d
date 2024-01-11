!> @date Feb 2023
subroutine elem_face(Mdle,Iface, Mdlf)
!
   use data_structure3D
   implicit none
   integer, intent(in)  :: Mdle,Iface
   integer, intent(out) :: Mdlf
!
   integer :: ntype,inode
   integer :: nodesl(27),nvoid(27)
!
!..check middle node type
   ntype=NODES(Mdle)%ntype
   select case(ntype)
      case(MDLB,MDLP,MDLN,MDLD)
      case default
        write(*,9999) Mdle,S_Type(ntype)
 9999   format(' elem_face: inconsistent node type! Mdle,type = ',i7,2x,a4)
        stop
   endselect
!
!..determine face node
   call elem_nodes(Mdle, nodesl,nvoid)
   inode = nvert(ntype) + nedge(ntype) + Iface
   mdlf  = nodesl(inode)
!
end subroutine elem_face

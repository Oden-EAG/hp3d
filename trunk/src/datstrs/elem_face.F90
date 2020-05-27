subroutine elem_face(Mdle,Iface, Mdlf)
!
      use data_structure3D
      implicit none
      integer, intent(in)  :: Mdle
      integer, intent(in)  :: Iface
      integer, intent(out) :: Mdlf
!
      character(len=4)      :: etype
      integer,dimension(27) :: nodesl,void
      integer :: inode
!------------------------------------------------------------------------------
!
!  ...check middle node type
      etype=NODES(Mdle)%type
      select case(etype)
      case('mdlp','mdlb','mdln','mdld')
      case default
        write(*,9999)Mdle,etype
 9999   format(' elem_face: inconsistent node type! Mdle,type = ',i7,2x,a4)
        stop
      endselect
!
!  ...determine face node
      call elem_nodes(Mdle, nodesl,void)
      inode = nvert(etype) + nedge(etype) + Iface
      mdlf  = nodesl(inode)
!
!
endsubroutine elem_face

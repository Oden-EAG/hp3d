!
subroutine initiate_order(Ntype, Norder)
!
   use node_types
   implicit none
!
   integer, intent(in)  :: Ntype
   integer, intent(out) :: Norder(19)
!
   integer,parameter,dimension(15) :: norder_prism =           &
      (/1,1,1, 1,1,1, 1,1,1, 1,1, 11,11,11, 11/)
   integer,parameter,dimension(19) :: norder_brick =           &
      (/1,1,1,1, 1,1,1,1, 1,1,1,1, 11,11,11,11,11,11, 111/)
   integer,parameter,dimension(11) :: norder_tetra =           &
      (/1,1,1,1,1,1, 1,1,1,1, 1/)
   integer,parameter,dimension(14) :: norder_pyram =           &
      (/1,1,1,1, 1,1,1,1, 11,1,1,1,1, 1/)
!
   select case(Ntype)
      case(MDLP,PRIS) ; Norder(1:15)=norder_prism(1:15)
      case(MDLB,BRIC) ; Norder(1:19)=norder_brick(1:19)
      case(MDLN,TETR) ; Norder(1:11)=norder_tetra(1:11)
      case(MDLD,PYRA) ; Norder(1:14)=norder_pyram(1:14)
   end select
!
end subroutine initiate_order

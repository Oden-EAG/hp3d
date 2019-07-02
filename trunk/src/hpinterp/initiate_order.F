      subroutine initiate_order(Type, Norder)
c
      implicit none
      character(len=4)     ,intent(in ) :: Type
      integer,dimension(19),intent(out) :: Norder
c
      integer,parameter,dimension(15) :: norder_prism =
     . (/1,1,1, 1,1,1, 1,1,1, 1,1, 11,11,11, 11/)
      integer,parameter,dimension(19) :: norder_brick =
     . (/1,1,1,1, 1,1,1,1, 1,1,1,1, 11,11,11,11,11,11, 111/)
      integer,parameter,dimension(11) :: norder_tetra =
     . (/1,1,1,1,1,1, 1,1,1,1, 1/)
      integer,parameter,dimension(14) :: norder_pyram = 
     . (/1,1,1,1, 1,1,1,1, 11,1,1,1,1, 1/)
c
      select case(Type)
      case('mdlp','pris') ; Norder(1:15)=norder_prism(1:15)
      case('mdlb','bric') ; Norder(1:19)=norder_brick(1:19)
      case('mdln','tetr') ; Norder(1:11)=norder_tetra(1:11)
      case('mdld','pyra') ; Norder(1:14)=norder_pyram(1:14)
      endselect
c
c      
      end subroutine

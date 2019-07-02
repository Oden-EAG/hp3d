!> Purpose : update the dirichle BC on the face
!!           H1 projection for ZnodH, L2 projection for ZnodE and ZnodV
!! @param[in]  Iflag        - a flag specifying which of the objects
!!                            5 pris, 6 hexa, 7 tetr, 8 pyra
!! @param[in]  No           - number of a specific object 
!! @param[in]  Etav         - reference coordinates of the element vertices
!! @param[in]  Type         - element (middle node) type
!! @param[in]  Icase        - the mid-face node case
!! @param[in]  Nedge_orient - edge orientation (not needed really)
!! @param[in]  Nface_orient - face orientation (not needed really)
!! @param[in]  Norder       - element order
!! @param[in]  Iface        - face number 
!! @param[in]  ZdofH        - H1 dof for the element (vertex values really)
!! @param[in]  ZdofE        - H(curl) dof for the element (edge values really)
!! 
!! @param[out] ZnodH        - H1 dof for the edge
!! @param[out] ZnodE        - H(curl) dof for the edge
!! @param[out] ZnodV        - H(div) dof for the edge
#include "implicit_none.h"
subroutine dhpface( &
     Iflag,No,Etav, Type,Icase, &
     Nedge_orient,Nface_orient,Norder,Iface, &
     ZdofH,ZnodH, &
     ZdofE,ZnodE, &
     ZnodV)
  use parameters
  use physics
  use element_data
  implicit none
  !
  ! ** Arguments
  !-----------------------------------------------------------------------
  integer,                                    intent(in)  :: Iflag,No
  integer,                                    intent(in)  :: Icase,Iface
  real*8,  dimension(3,8),                    intent(in)  :: Etav
  character(len=4),                           intent(in)  :: Type
  integer, dimension(12),                     intent(in)  :: Nedge_orient
  integer, dimension(6),                      intent(in)  :: Nface_orient
  integer, dimension(19),                     intent(in)  :: Norder

  VTYPE,   dimension(MAXEQNH,MAXbrickH),      intent(in)  :: ZdofH
  VTYPE,   dimension(MAXEQNE,MAXbrickE),      intent(in)  :: ZdofE
  VTYPE,   dimension(NRCOMS*NREQNH(Icase),*), intent(out) :: ZnodH
  VTYPE,   dimension(NRCOMS*NREQNE(Icase),*), intent(out) :: ZnodE
  VTYPE,   dimension(NRCOMS*NREQNV(Icase),*), intent(out) :: ZnodV
  !-----------------------------------------------------------------------
  call dhpfaceH( &
     Iflag,No,Etav, Type,Icase, &
     Nedge_orient,Nface_orient,Norder,Iface, &
     ZdofH,ZnodH, &
     ZnodV)
  !
  return

  call dhpfaceE( &
     Iflag,No,Etav, Type,Icase, &
     Nedge_orient,Nface_orient,Norder,Iface, &
     ZdofE,ZnodE, &
     ZnodV)
!
!
end subroutine dhpface

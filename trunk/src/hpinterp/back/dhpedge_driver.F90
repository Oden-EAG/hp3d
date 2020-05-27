!> Purpose : update the dirichle BC on the edge
!!           H1 projection for ZnodH, L2 projection for ZnodE
!! @param[in]  Iflag        - a flag specifying which of the objects
!!                            5 pris, 6 hexa, 7 tetr, 8 pyra
!! @param[in]  No           - number of a specific object
!! @param[in]  Etav         - reference coordinates of the element vertices
!! @param[in]  Type         - element (middle node) type
!! @param[in]  Icase        - the mid-edge node case
!! @param[in]  Nedge_orient - edge orientation (not needed really)
!! @param[in]  Nface_orient - face orientation (not needed really)
!! @param[in]  Norder       - element order
!! @param[in]  Iedg         - edge number
!! @param[in]  ZdofH        - H1 dof for the element (vertex values really)
!!
!! @param[out] ZnodH        - H1 dof for the edge
!! @param[out] ZnodE        - H(curl) dof for the edge

#include "typedefs.h"
subroutine dhpedge( &
     Iflag,No,Etav, Type,Icase, &
     Nedge_orient,Nface_orient,Norder,Iedg,  &
     ZdofH, &
     ZnodH, ZnodE)
  use parameters
  use physics
  use element_data
  implicit none
  ! ** Arguments
  !--------------------------------------------------------------------
  integer,                                    intent(in)  :: Iflag,No
  integer,                                    intent(in)  :: Icase,Iedg
  real*8,  dimension(3,8),                    intent(in)  :: Etav
  character(len=4),                           intent(in)  :: Type
  integer, dimension(12),                     intent(in)  :: Nedge_orient
  integer, dimension(6),                      intent(in)  :: Nface_orient
  integer, dimension(19),                     intent(in)  :: Norder
  VTYPE,   dimension(MAXEQNH,MAXbrickH),      intent(in)  :: ZdofH
  VTYPE,   dimension(NRCOMS*NREQNH(Icase),*), intent(out) :: ZnodH
  VTYPE,   dimension(NRCOMS*NREQNE(Icase),*), intent(out) :: ZnodE

  call dhpedgeH( &
       Iflag,No,Etav, Type,Icase, &
       Nedge_orient,Nface_orient,Norder,Iedg, &
       ZdofH,ZnodH)

  call dhpedgeE( &
       Iflag,No,Etav, Type,Icase, &
       Nedge_orient,Nface_orient,Norder,Iedg, &
       ZnodE)

end subroutine dhpedge


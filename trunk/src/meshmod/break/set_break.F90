!---------------------------------------------------------------------------------------------
!> Purpose : determine the infos necessary for node breaking
!!
!! @param[in ] Type_nod  - node type
!! @param[in ] Kref      - refinement kind
!! @param[in ] Nord      - order of approximation for the node
!! @param[in ] Nbc       - BC flags
!! @param[in ] Nbc       - BC flags
!! @param[in ] Subd      - subdomain of node
!! @param[out] Nrsons    - number of sons
!! @param[out] Type_sons - son type
!! @param[out] Norder    - order of approximation for sons
!! @param[out] Nbcond    - BC flags for sons
!! @param[out] Nsubd     - subdomain for sons
!!
!! rev@Dec 12
!---------------------------------------------------------------------------------------------
subroutine set_break(Type_nod,Kref,Nord,Nbc,Subd, Nrsons,Type_sons,Norder,Nbcond,Nsubd)
  implicit none
!
  character(len=4),                intent(in)  :: Type_nod
  integer,                         intent(in)  :: Kref, Nord, Nbc, Subd
!
  integer,          dimension(27), intent(out) :: Norder, Nbcond, Nsubd
  integer,                         intent(out) :: Nrsons
  character(len=4), dimension(27), intent(out) :: Type_sons
!
! initialize
  Norder = 0; Nbcond = 0; Type_sons(1:27) = 'none'
  select case (Type_nod)
! EDGE
  case('medg') ; call set_edge_break(Kref,Nord, Nrsons,Type_sons,Norder)
! TRIANGLE
  case('mdlt') ; call set_tria_break(Kref,Nord, Nrsons,Type_sons,Norder)
! QUAD
  case('mdlq') ; call set_quad_break(Kref,Nord, Nrsons,Type_sons,Norder)
! BRICK
  case('mdlb') ; call set_bric_break(Kref,Nord, Nrsons,Type_sons,Norder)
! TET
  case('mdln') ; call set_tetr_break(Kref,Nord, Nrsons,Type_sons,Norder)
! PRISM
  case('mdlp') ; call set_pris_break(Kref,Nord, Nrsons,Type_sons,Norder)
! PYRAMID
  case('mdld') ; call set_pyra_break(Kref,Nord, Nrsons,Type_sons,Norder)
  case default
     write(*,*) 'set_break: NOT SUPPORTED TYPE ', Type_nod
     stop
  endselect
!
! BC flags are inherited directly from the father node
  Nbcond(1:Nrsons)=Nbc
!
! inherit subdomain from middle node father
  Nsubd(1:27) = -1
  select case (Type_nod)
  case('mdlb','mdln','mdlp','mdld') ; Nsubd(1:Nrsons) = Subd
  endselect
!
endsubroutine set_break

!---------------------------------------------------------------------------------------------
!> @brief determine the infos necessary for node breaking
!!
!> @param[in ] Ntype      - node type
!> @param[in ] Kref       - refinement kind
!> @param[in ] Nord       - order of approximation for the node
!> @param[in ] Nbc        - BC flags
!> @param[in ] Subd       - subdomain of node
!> @param[out] Nrsons     - number of sons
!> @param[out] Ntype_sons - son type
!> @param[out] Norder     - order of approximation for sons
!> @param[out] Nbcond     - BC flags for sons
!> @param[out] Nsubd      - subdomain for sons
!!
!> @date Feb 2023
!---------------------------------------------------------------------------------------------
subroutine set_break(Ntype,Kref,Nord,Nbc,Subd, Nrsons,Ntype_sons,Norder,Nbcond,Nsubd)
  use node_types
  implicit none
  integer,                intent(in)  :: Ntype, Kref, Nord, Nbc, Subd
  integer,                intent(out) :: Nrsons
  integer, dimension(27), intent(out) :: Ntype_sons, Norder, Nbcond, Nsubd
!
! initialize
  Norder = 0; Nbcond = 0; Ntype_sons(1:27) = 0
  select case (Ntype)
! EDGE
  case(MEDG) ; call set_edge_break(Kref,Nord, Nrsons,Ntype_sons,Norder)
! TRIANGLE
  case(MDLT) ; call set_tria_break(Kref,Nord, Nrsons,Ntype_sons,Norder)
! QUAD
  case(MDLQ) ; call set_quad_break(Kref,Nord, Nrsons,Ntype_sons,Norder)
! BRICK
  case(MDLB) ; call set_bric_break(Kref,Nord, Nrsons,Ntype_sons,Norder)
! TET
  case(MDLN) ; call set_tetr_break(Kref,Nord, Nrsons,Ntype_sons,Norder)
! PRISM
  case(MDLP) ; call set_pris_break(Kref,Nord, Nrsons,Ntype_sons,Norder)
! PYRAMID
  case(MDLD) ; call set_pyra_break(Kref,Nord, Nrsons,Ntype_sons,Norder)
  case default
     write(*,*) 'set_break: NOT SUPPORTED TYPE ', S_Type(Ntype)
     stop
  endselect
!
! BC flags are inherited directly from the father node
  Nbcond(1:Nrsons) = Nbc
!
! inherit subdomain from middle node father
  Nsubd(1:27) = -1
  select case (Ntype)
    case(MDLB,MDLN,MDLP,MDLD) ; Nsubd(1:Nrsons) = Subd
  end select
!
endsubroutine set_break

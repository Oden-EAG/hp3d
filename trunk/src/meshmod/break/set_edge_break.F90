!----------------------------------------------------------------------------
!> @brief determines info for breaking an edge node
!!
!> @param[in ] Kref    - node refinement kind
!> @param[in ] Nord    - node order of approximation
!> @param[out] Nrsons  - number of sons
!> @param[out] Ntype   - sons' types
!> @param[out] Norder  - sons' orders of approximation
!!
!> @date Feb 2023
!----------------------------------------------------------------------------
subroutine set_edge_break(Kref,Nord, Nrsons,Ntype,Norder)
  use node_types
  implicit none
  integer, intent(in)  :: Kref,Nord
  integer, intent(out) :: Nrsons
  integer, intent(out) :: Ntype(27),Norder(27)
!----------------------------------------------------------------------------
!
! initialize
  Norder = 0; Ntype(1:27) = 0
!
! select refinement kind
  select case(Kref)
  case(1)
     Nrsons = 3
     Ntype (1:3) = (/ MEDG,MEDG,VERT /)
     Norder(1:3) = (/ Nord,Nord,1    /)
  case default
     Nrsons=0
  endselect
!
end subroutine set_edge_break


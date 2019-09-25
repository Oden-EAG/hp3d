!----------------------------------------------------------------------------
!> Purpose : determines info for breaking an edge node
!!
!> @param[in ] Kref    - node refinement kind
!> @param[in ] Nord    - node order of approximation
!> @param[out] Nrsons  - number of sons
!> @param[out] Type    - sons' types
!> @param[out] Norder  - sons' orders of approximation
!!
!> rev@Aug 2019
!----------------------------------------------------------------------------
subroutine set_edge_break(Kref,Nord, Nrsons,Type,Norder)
!
  implicit none
  integer,                         intent(in)  :: Kref,Nord
  integer,          dimension(27), intent(out) :: Norder
  integer,                         intent(out) :: Nrsons
  character(len=4), dimension(27), intent(out) :: Type
!----------------------------------------------------------------------------
!
! initialize
  Norder = 0; Type(1:27) = 'none'
!
! select refinement kind
  select case(Kref)
  case(1)
     Nrsons = 3
     Type(1:3) = (/'medg','medg','vert'/)
     Norder(1:3) = (/Nord,Nord,1/)
  case default
     Nrsons=0
  endselect
!
end subroutine set_edge_break


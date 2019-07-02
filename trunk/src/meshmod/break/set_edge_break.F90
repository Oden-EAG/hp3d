!----------------------------------------------------------------------------
!> Purpose : determines info for breaking an edge node
!!
!> @param[in ] Kref    - node refinement kind
!> @param[in ] Kfilter - node refinement filter
!> @param[in ] Nord    - node order of approximation
!> @param[out] Nrsons  - number of sons
!> @param[out] Type    - sons' types
!> @param[out] Norder  - sons' orders of approximation
!> @param[out] Nfilter - sons' refinement filters
!!
!> rev@Dec 12
!----------------------------------------------------------------------------
subroutine set_edge_break(Kref,Kfilter,Nord, Nrsons,Type,Norder,Nfilter)
!
!  
  implicit none
  integer,                         intent(in)  :: Kref, Kfilter, Nord
  integer,          dimension(27), intent(out) :: Norder, Nfilter
  integer,                         intent(out) :: Nrsons
  character(len=4), dimension(27), intent(out) :: Type
!----------------------------------------------------------------------------
!
! initialize
  Type(1:27)='none'
!
  Norder(1:27)=0 ; Nfilter(1:27)=0
!
! select refinement kind
  select case(Kref)
  case(1)
     Nrsons = 3
     Type(1:3) = (/'medg','medg','vert'/)
     Norder(1:3) = (/Nord,Nord,1/)
     Nfilter(1:2) = Kfilter
  case default
     Nrsons=0
  endselect
!
!
endsubroutine set_edge_break

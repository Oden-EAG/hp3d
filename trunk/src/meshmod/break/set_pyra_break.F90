!----------------------------------------------------------------------------
!> @brief determines info for breaking a pyramid's middle node
!!
!> @param[in ] Kref    - node refinement kind
!> @param[in ] Nord    - node order of approximation
!> @param[out] Nrsons  - number of sons
!> @param[out] Ntype   - sons' types
!> @param[out] Norder  - sons' orders of approximation
!!
!> @date Feb 2023
!----------------------------------------------------------------------------
subroutine set_pyra_break(Kref,Nord, Nrsons,Ntype,Norder)
  use node_types
  implicit none
  integer, intent(in)  :: Kref,Nord
  integer, intent(out) :: Nrsons
  integer, intent(out) :: Ntype(27),Norder(27)
  integer :: nordp
!----------------------------------------------------------------------------
!
  nordp=nord*10+nord
!
! initialize
  Norder = 0; Ntype(1:27) = 0
!
! select refinement kind
  select case(Kref)
  case(10)
     Nrsons = 7
!                 |    INTERIOR NODES     |    FACE NODES    |
     Ntype (1:7)=(/ MDLD,MDLP ,MDLP ,MDLP ,MDLQ ,MDLQ ,MDLQ  /)
     Norder(1:7)=(/ nord,nordp,nordp,nordp,nordp,nordp,nordp /)
  case default
     Nrsons=0
  end select
!
end subroutine set_pyra_break


!----------------------------------------------------------------------------
!> @brief determines info for breaking a quad node
!!
!> @param[in ] Kref    - node refinement kind
!> @param[in ] Nord    - node order of approximation
!> @param[out] Nrsons  - number of sons
!> @param[out] Ntype   - sons' types
!> @param[out] Norder  - sons' orders of approximation
!!
!> @date Feb 2023
!----------------------------------------------------------------------------
subroutine set_quad_break(Kref,Nord, Nrsons,Ntype,Norder)
  use node_types
  implicit none
  integer, intent(in)  :: Kref,Nord
  integer, intent(out) :: Nrsons
  integer, intent(out) :: Ntype(27),Norder(27)
  integer :: nordx,nordy
!----------------------------------------------------------------------------
!
! initialize
  Norder = 0; Ntype(1:27) = 0
!
  call decode(Nord, nordx,nordy)
!
! select refinement kind
  select case(Kref)
  case(01)
     Nrsons = 3
     Ntype (1:3) = (/ MDLQ,MDLQ,MEDG  /)
     Norder(1:3) = (/ nord,nord,nordx /)
  case(10)
     Nrsons = 3
     Ntype (1:3) = (/ MDLQ,MDLQ,MEDG  /)
     Norder(1:3) = (/ nord,nord,nordy /)
  case(11)
     Nrsons = 9
     Ntype (1:9) = (/ MDLQ,MDLQ,MDLQ,MDLQ,MEDG ,MEDG ,MEDG ,MEDG ,VERT /)
     Norder(1:9) = (/ nord,nord,nord,nord,nordy,nordx,nordy,nordx,1    /)
  case default
     Nrsons = 0
  end select
!
end subroutine set_quad_break


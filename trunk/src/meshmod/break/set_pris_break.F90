!----------------------------------------------------------------------------
!> @brief prepare the necessary information for breaking a prism middle node
!!
!> @param[in ] Kref    - node refinement kind
!> @param[in ] Nord    - node order of approximation
!> @param[out] Nrsons  - number of sons
!> @param[out] Ntype   - sons' types
!> @param[out] Norder  - sons' orders of approximation
!!
!> @date Feb 2023
!----------------------------------------------------------------------------
subroutine set_pris_break(Kref,Nord, Nrsons,Ntype,Norder)
  use node_types
  implicit none
  integer, intent(in)  :: Kref,Nord
  integer, intent(out) :: Nrsons
  integer, intent(out) :: Ntype(27),Norder(27)
!----------------------------------------------------------------------------
!
  integer :: nordh, nordz
  call decode(Nord, nordh,nordz)
!
  Norder = 0; Ntype(1:27) = 0
!
  select case(Kref)
  case(01)
     Nrsons       = 3
     Ntype(1:3)   = (/ MDLP,MDLP,MDLT  /)
     Norder(1:3)  = (/ nord,nord,nordh /)
  case(10)
     Nrsons = 7
     Ntype(1:4)   = MDLP
     Ntype(5:7)   = MDLQ
     Norder(1:7)  = nord
  case(11)
     Nrsons = 21
     Ntype(1:8)    = MDLP
     Ntype(9:12)   = MDLT
     Ntype(13:18)  = MDLQ
     Ntype(19:21)  = MEDG
     Norder(1:8)   = nord
     Norder(9:12)  = nordh
     Norder(13:18) = nord
     Norder(19:21) = nordh
  case default
     Nrsons = 0
  end select
!
end subroutine set_pris_break


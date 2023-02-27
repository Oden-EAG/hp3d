!----------------------------------------------------------------------------
!> @brief set the necessary infomation for the tetr break
!!
!> @param[in ] Kref    - node refinement kind
!> @param[in ] Nord    - node order of approximation
!> @param[out] Nrsons  - number of sons
!> @param[out] Ntype   - sons' types
!> @param[out] Norder  - sons' orders of approximation
!!
!> @date Feb 2023
!----------------------------------------------------------------------------
subroutine set_tetr_break(Kref,Nord, Nrsons,Ntype,Norder)
  use node_types
  implicit none
  integer, intent(in)  :: Kref,Nord
  integer, intent(out) :: Nrsons
  integer, intent(out) :: Ntype(27),Norder(27)
  integer :: nordp
!----------------------------------------------------------------------------
!
  Norder = 0; Ntype(1:27) = 0
!
  nordp = nord*10+nord
!
  select case(Kref)
  case(11,12,13)
     Nrsons       = 17
     Ntype(1:8)   = MDLN
     Ntype(9:16)  = MDLT
     Ntype(17)    = MEDG
     Norder(1:17) = nord
  case(24)
     Nrsons      = 7
     Ntype (1:7) = (/ MDLP ,MDLN,MDLN,MDLD,MDLQ ,MDLT,MDLT /)
     Norder(1:7) = (/ nordp,nord,nord,nord,nordp,nord,nord /)
  case(31,32,33,34)
     Nrsons = 3
     Ntype(1:3)  = (/ MDLN,MDLP ,MDLT /)
     Norder(1:3) = (/ nord,nordp,nord /)
  case default
     Nrsons = 0
  end select
!
end subroutine set_tetr_break


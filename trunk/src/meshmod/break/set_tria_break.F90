!----------------------------------------------------------------------------
!> @brief set the necessary infomation for the tria break
!!
!> @param[in ] Kref    - node refinement kind
!> @param[in ] Nord    - node order of approximation
!> @param[out] Nrsons  - number of sons
!> @param[out] Ntype   - sons' types
!> @param[out] Norder  - sons' orders of approximation
!!
!> @date Feb 2023
!----------------------------------------------------------------------------
subroutine set_tria_break(Kref,Nord, Nrsons,Ntype,Norder)
  use node_types
  implicit none
  integer, intent(in)  :: Kref,Nord
  integer, intent(out) :: Nrsons
  integer, intent(out) :: Ntype(27),Norder(27)
!----------------------------------------------------------------------------
!
  Norder = 0; Ntype(1:27) = 0
!
  select case(Kref)
  case(1)
     Nrsons = 7
     Ntype(1:7)  = (/ MDLT,MDLT,MDLT,MDLT, &
                      MEDG,MEDG,MEDG /)
     Norder(1:7) = Nord
  case(2,3,4)
     Nrsons = 3
     Ntype(1:3) = (/ MDLT,MDLQ,MEDG /)
     Norder(1) = Nord
     Norder(2) = Nord*10+Nord
     Norder(3) = Nord
  case default
     Nrsons = 0
  end select
!
end subroutine set_tria_break


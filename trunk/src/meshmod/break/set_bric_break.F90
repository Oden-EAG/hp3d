!----------------------------------------------------------------------------
!> @brief determines info for breaking a brick middle node
!!
!> @param[in ] Kref    - node refinement kind
!> @param[in ] Nord    - node order of approximation
!> @param[out] Nrsons  - number of sons
!> @param[out] Ntype   - sons' types
!> @param[out] Norder  - sons' orders of approximation
!!
!> @date Feb 2023
!----------------------------------------------------------------------------
subroutine set_bric_break(Kref,Nord, Nrsons,Ntype,Norder)
  use node_types
  implicit none
  integer, intent(in)  :: Kref,Nord
  integer, intent(out) :: Nrsons
  integer, intent(out) :: Ntype(27),Norder(27)
  integer :: nordxy, nordyz, nordxz, nordx, nordy, nordz
!----------------------------------------------------------------------------
!
  call decode(nord,   nordxy,nordz)
  call decode(nordxy, nordx, nordy)
!
  nordxz = nordx*10 + nordz
  nordyz = nordy*10 + nordz
!
  Norder = 0; Ntype(1:27) = 0
!
  select case (Kref)
  case (001,010,100)
     Nrsons     = 3
     Ntype(1:2) = MDLB
     Ntype(3)   = MDLQ
     select case (Kref)
     case(001)
        Norder(1:3) = (/nord,nord,nordxy/)
     case(010)
        Norder(1:3) = (/nord,nord,nordxz/)
     case(100)
        Norder(1:3) = (/nord,nord,nordyz/)
     end select
  case (011,101,110)
     Nrsons     = 9
     Ntype(1:4) = MDLB
     Ntype(5:8) = MDLQ
     Ntype(9)   = MEDG
     select case (Kref)
     case(011)
        Norder(1:9) = (/&
            nord,nord,nord,nord, &
            nordxz,nordxy,nordxz,nordxy,nordx/)
     case(101)
        Norder(1:9) = (/ &
            nord,nord,nord,nord, &
            nordyz,nordxy,nordyz,nordxy,nordy/)
     case(110)
        Norder(1:9) = (/ &
            nord,nord,nord,nord, &
            nordyz,nordxz,nordyz,nordxz,nordz/)
     end select
  case (111)
     Nrsons       = 27
     Ntype(1:8)   = MDLB
     Ntype(9:20)  = MDLQ
     Ntype(21:26) = MEDG
     Ntype(27)    = VERT
     ! should be as follows (according to hp book)
     ! Norder(1:27) = (/ &
     !      nord,nord,nord,nord,nord,nord,nord,nord, &
     !      nordyz,nordyz,nordyz,nordyz, &
     !      nordxz,nordxz,nordxz,nordxz, &
     !      nordxy,nordxy,nordxy,nordxy, &
     !      nordx,nordx,nordy,nordy,nordz,nordz,1/)
     !
     ! but is as follows (due to connectivity in files/ref/brick_111)
     Norder(1:27) = (/ &
          nord,nord,nord,nord,nord,nord,nord,nord, &
          nordyz,nordxz,nordyz,nordxz, &
          nordxy,nordxy,nordxy,nordxy, &
          nordyz,nordxz,nordyz,nordxz, &
          nordz,nordy,nordx,nordy,nordx,nordz,1/)
  case default
     Nrsons = 0
  end select
!
end subroutine set_bric_break


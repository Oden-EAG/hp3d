!----------------------------------------------------------------------------
!> Purpose : determines info for breaking a quad node
!!
!> @param[in ] Kref    - node refinement kind
!> @param[in ] Nord    - node order of approximation
!> @param[out] Nrsons  - number of sons
!> @param[out] Type    - sons' types
!> @param[out] Norder  - sons' orders of approximation
!!
!> rev@Aug 2019
!----------------------------------------------------------------------------
subroutine set_quad_break(Kref,Nord, Nrsons,Type,Norder)
!
  implicit none
  integer,                         intent(in)  :: Kref,Nord
  integer,          dimension(27), intent(out) :: Norder
  integer,                         intent(out) :: Nrsons
  character(len=4), dimension(27), intent(out) :: Type
  integer :: nordx, nordy
!----------------------------------------------------------------------------
!
! initialize
  Norder = 0; Type(1:27) = 'none'
!
  call decode(Nord, nordx,nordy)
!
! select refinement kind
  select case(Kref)
  case(01)
     Nrsons = 3
     Type(1:3) = (/'mdlq','mdlq','medg'/)
     Norder(1:3) = (/nord,nord,nordx/)
  case(10)
     Nrsons = 3
     Type(1:3) = (/'mdlq','mdlq','medg'/)
     Norder(1:3) = (/nord,nord,nordy/)
  case(11)
     Nrsons = 9
     Type(1:9) = (/ &
          'mdlq','mdlq','mdlq','mdlq', &
          'medg','medg','medg','medg', &
          'vert' &
          /)
     Norder(1:9)=(/nord,nord,nord,nord,nordy,nordx,nordy,nordx,1/)
  case default
     Nrsons = 0
  end select
!
end subroutine set_quad_break


!> Purpose : set the necessary infomation for the trig break
subroutine set_tria_break(Kref,Nord, Nrsons,Type,Norder)
  implicit none
  integer,                         intent(in)  :: Kref,Nord
  integer,          dimension(27), intent(out) :: Norder
  integer,                         intent(out) :: Nrsons
  character(len=4), dimension(27), intent(out) :: Type
!
  Norder = 0; Type(1:27) = 'none'
!
  select case(Kref)
  case(1)
     Nrsons = 7
     Type(1:7) = (/ &
          'mdlt','mdlt','mdlt','mdlt', &
          'medg','medg','medg' &
          /)
     Norder(1:7) = Nord
  case(2,3,4)
     Nrsons = 3
     Type(1:3) = (/'mdlt','mdlq','medg'/)
     Norder(1) = Nord
     Norder(2) = Nord*10+Nord
     Norder(3) = Nord
  case default
     Nrsons = 0
  end select
!
end subroutine set_tria_break


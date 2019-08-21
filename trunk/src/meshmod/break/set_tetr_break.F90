!> Purpose : set the necessary infomation for the tetr break
subroutine set_tetr_break(Kref,Nord, Nrsons,Type,Norder)
  implicit none
  integer,                         intent(in)  :: Kref,Nord
  integer,          dimension(27), intent(out) :: Norder
  integer,                         intent(out) :: Nrsons
  character(len=4), dimension(27), intent(out) :: Type
  integer :: nordp
!
  Norder = 0; Type(1:27) = 'none'
!
  nordp = nord*10+nord
!
  select case(Kref)
  case(11,12,13)
     Nrsons       = 17
     Type(1:8)    = 'mdln'
     Type(9:16)   = 'mdlt'
     Type(17)     = 'medg'
     Norder(1:17) = nord
  case(24)
     Nrsons    = 7
     Type(1:7) = (/ &
          'mdlp','mdln','mdln','mdld', &
          'mdlq','mdlt','mdlt' &
          /)
     Norder(1:7)  = (/nordp,nord,nord,nord, nordp,nord,nord/)
  case(31,32,33,34)
     Nrsons = 3
     Type(1:3) = (/'mdln','mdlp','mdlt'/)
     Norder(1:3) = (/nord,nordp,nord/)
  case default
     Nrsons = 0
  end select
!
end subroutine set_tetr_break


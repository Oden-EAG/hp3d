!> Purpose : set the necessary infomation for the bric break
subroutine set_tetr_break(&
     Kref, Kfilter, Nord, &
     Nrsons, Type, Norder, Nfilter)
  implicit none
  integer,                         intent(in)  :: Kref, Kfilter, Nord
  integer,          dimension(27), intent(out) :: Norder, Nfilter
  integer,                         intent(out) :: Nrsons
  character(len=4), dimension(27), intent(out) :: Type
  !!!! Add filter 
  integer :: nordp

  nordp = nord*10+nord
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

     if (Kfilter.ne.0) then
        Nfilter(1:7) = (/10,24,24,10,0,0,0/)
     endif

  case(31,32,33,34)
     Nrsons = 3
     Type(1:3) = (/'mdln','mdlp','mdlt'/)
     Norder(1:3) = (/nord,nordp,nord/)

     if (Kfilter.ne.0) then
        Nfilter(1:3) = (/Kfilter,01,0/)
     endif

  case default
     Nrsons = 0
  end select
end subroutine set_tetr_break

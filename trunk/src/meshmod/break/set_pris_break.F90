!> Purpose : prepare the necessary information for breaking a prism middle node
subroutine set_pris_break( &
           Kref, Kfilter, Nord, &
           Nrsons, Type, Norder, Nfilter)
  implicit none
  integer,                         intent(in)  :: Kref, Kfilter, Nord
  integer,          dimension(27), intent(out) :: Norder, Nfilter
  integer,                         intent(out) :: Nrsons
  character(len=4), dimension(27), intent(out) :: Type

  integer :: nordh, nordz
  call decode(Nord, nordh,nordz)
  select case(Kref)
  case(01)
     Nrsons       = 3
     Type(1:3)    = (/'mdlp','mdlp','mdlt'/)
     Norder(1:3)  = (/nord,nord,nordh/)
     Nfilter(1:2) = Kfilter
  case(10)
     Nrsons = 7
     Type(1:4)    = 'mdlp'
     Type(5:7)    = 'mdlq'
     Norder(1:7)  = nord
     Nfilter(1:4) = Kfilter
  case(11)
     Nrsons = 21
     Type(1:8)     = 'mdlp'
     Type(9:12)    = 'mdlt'
     Type(13:18)   = 'mdlq'
     Type(19:21)   = 'medg'
     Norder(1:8)   = nord
     Norder(9:12)  = nordh
     Norder(13:18) = nord 
     Norder(19:21) = nordh
     Nfilter(1:8)  = Kfilter
  case default
     Nrsons = 0
  end select
end subroutine set_pris_break

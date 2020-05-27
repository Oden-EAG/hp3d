!> Purpose : prepare the necessary information for breaking a brick middle node
subroutine set_bric_break(Kref,Nord, Nrsons,Type,Norder)
  implicit none
  integer,                         intent(in)  :: Kref,Nord
  integer,          dimension(27), intent(out) :: Norder
  integer,                         intent(out) :: Nrsons
  character(len=4), dimension(27), intent(out) :: Type
  integer :: nordxy, nordyz, nordxz, nordx, nordy, nordz
!
  call decode(nord,   nordxy,nordz)
  call decode(nordxy, nordx, nordy)
!
  nordxz = nordx*10 + nordz
  nordyz = nordy*10 + nordz
!
  Norder = 0; Type(1:27) = 'none'
!
  select case (Kref)
  case (001,010,100)
     Nrsons = 3
     Type(1:2) = 'mdlb'
     Type(3)   = 'mdlq'
     select case (Kref)
     case(001)
        Norder(1:3) = (/nord,nord,nordxy/)
     case(010)
        Norder(1:3) = (/nord,nord,nordxz/)
     case(100)
        Norder(1:3) = (/nord,nord,nordyz/)
     end select
  case (011,101,110)
     Nrsons      = 9
     Type(1:4)  = 'mdlb'
     Type(5:8)  = 'mdlq'
     Type(9)    = 'medg'
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
     Nrsons = 27
     Type(1:8)   = 'mdlb'
     Type(9:20)  = 'mdlq'
     Type(21:26) = 'medg'
     Type(27)    = 'vert'
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


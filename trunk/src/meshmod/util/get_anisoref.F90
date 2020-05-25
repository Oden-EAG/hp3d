subroutine get_anisoref(Mdle, Derr, Kref)
  use data_structure3D
  implicit none 
  integer,              intent(in)  :: Mdle
  real(8), dimension(3), intent(in) :: Derr
  integer,              intent(out) :: Kref
  real(8) :: derr_xy, derr_yz, derr_xz
  real(8), parameter :: ratio = 0.01d0

  derr_xy = sqrt(Derr(1)**2 + Derr(2)**2)
  derr_xz = sqrt(Derr(1)**2 + Derr(3)**2)
  derr_yz = sqrt(Derr(2)**2 + Derr(3)**2)

  Kref = 0
  select case(NODES(Mdle)%type)
  case('mdlp')
     if ((Derr(3)*ratio).gt.derr_xy) then
        Kref = Kref + 1
     end if
     if ((derr_xy*ratio).gt.Derr(3)) then
        Kref = Kref + 10
     end if
  case('mdlb')
     if ((Derr(3)*ratio).gt.derr_xy) then
        Kref = Kref + 1
     end if
     if ((Derr(2)*ratio).gt.derr_xz) then
        Kref = Kref + 10
     end if
     if ((Derr(3)*ratio).gt.derr_yz) then
        Kref = Kref + 100
     end if
     if (Kref.eq.0) then
        Kref = 111
     end if
  end select
  if (Kref.eq.0) then
     call get_isoref(Mdle, Kref)
  end if
end subroutine get_anisoref

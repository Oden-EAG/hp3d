!> Purpose : check initial mesh left oriented
!! @param[out] Ipass  - if all elements are right handed Ipass = 1
!!                          
subroutine check_left_oriented(Ipass)
  use data_structure3D
  implicit none
  integer, intent(out) :: Ipass
  integer :: iel
  
  Ipass = 1
  do iel=1, NRELIS
     if (.not.Is_right_handed(iel)) then
        write(*,*) 'Mdle = ', iel, ELEMS(iel)%type, '  LEFT ORIENTED'
        Ipass = 0
     end if
  end do

end subroutine check_left_oriented

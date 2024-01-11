!> @brief      check if the initial mesh is left oriented
!> @param[out] Ipass = 1, if all elements are right handed
!> @date       Feb 2023
subroutine check_left_oriented(Ipass)
   use data_structure3D
   implicit none
   integer, intent(out) :: Ipass
   integer :: iel
!
   Ipass = 1
   do iel=1, NRELIS
      if (.not.Is_right_handed(iel)) then
         write(*,*) 'Mdle = ', iel, S_Type(ELEMS(iel)%etype), '  LEFT ORIENTED'
         Ipass = 0
      end if
   end do
!
end subroutine check_left_oriented

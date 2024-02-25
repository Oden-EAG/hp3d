!-------------------------------------------------------------
!> @brief Returns factorial N!
!> @date Feb 2023
!-------------------------------------------------------------
integer function fact(N)
!
   implicit none
!
   integer, intent(in) :: N
!
   integer :: i
!
   if     (N.lt.0) then
      fact = 0
   elseif (N.eq.0) then
      fact = 1
   else
      fact = 1
      do i = 1,N
         fact = fact*i
      enddo
   endif
!
end function fact

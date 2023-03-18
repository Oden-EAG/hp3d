!-------------------------------------------------------------
!> @brief Returns factorial n!
!> @date Feb 2023
!-------------------------------------------------------------
integer function fact(n)
!
   implicit none
   integer :: n,i
!
   if (n.lt.0) then
      fact = 0
   elseif (n.eq.0) then
      fact = 1
   else
      fact = 1
      do i = 1,n
         fact = fact*i
      enddo
   endif
!
end function fact

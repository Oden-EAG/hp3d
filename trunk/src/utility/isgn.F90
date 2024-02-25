!----------------------------------------------------------------------
!
!   function name      - isgn
!
!----------------------------------------------------------------------
!> @brief Returns the sign of an integer value
!> @param[in] I - Integer value
!> @date  Feb 2023
!----------------------------------------------------------------------
integer function isgn(I)
!
   implicit none
!
   integer, intent(in) :: I
!
   if     (I.gt.0) then
      isgn= 1
   elseif (I.lt.0) then
      isgn=-1
   else
      isgn= 0
   endif
!
end function isgn

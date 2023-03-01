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
   integer I
!
   if (I.eq.0) isgn=  0; return
   if (I.gt.0) isgn=  1; return
   if (I.lt.0) isgn= -1; return
!
end function isgn

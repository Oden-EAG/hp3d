!----------------------------------------------------------------------
!
!   routine name       - locate_char
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine locates first instance of a
!                        string on a list of characters
!
!   arguments :
!     in:
!              Nel     - string to search for
!              List    - a list of strings
!              Nlist   - length of the list
!     out:
!              Number  - the string location
!
!----------------------------------------------------------------------
subroutine locate_char(Nel,List,Nlist, Number)
!
   implicit none
!
   integer         , intent(in)  :: Nlist
   character(len=*), intent(in)  :: Nel,List(Nlist)
   integer         , intent(out) :: Number
!
   integer :: i
!
   do i=1,Nlist
      if (List(i).eq.Nel) then
         Number=i
         return
      endif
   enddo
   Number=0
!
end subroutine locate_char

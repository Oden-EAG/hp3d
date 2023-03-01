c----------------------------------------------------------------------
c
c   routine name       - locate_char
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - routine locates first instance of a
c                        character on a list of characters
c
c   arguments :
c     in:
c              Nel     - character to search for
c              List    - a list of characters
c              Nlist   - length of the list
c     out:
c              Number  - the character location
c
c----------------------------------------------------------------------
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

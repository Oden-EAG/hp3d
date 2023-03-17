!----------------------------------------------------------------------
!
!   routine name       - locate
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine locates first instance of an integer
!                        on a list of integers
!
!   arguments :
!     in:
!              Nel     - an integer to be located
!              List    - a list of integers
!              Nlist   - length of the list
!     out:
!              Number  - the integer location
!
!----------------------------------------------------------------------
subroutine locate(Nel,List,Nlist, Number)
!
   implicit none
!
   integer, intent(in)  :: Nel
   integer, intent(in)  :: Nlist
   integer, intent(in)  :: List(Nlist)
   integer, intent(out) :: Number
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
end subroutine locate

!----------------------------------------------------------------------
!
!   routine name       - locate
!
!----------------------------------------------------------------------
!
!   computer           - machine independent
!
!   latest revision    - Jun 07
!
!   purpose            - routine locates an integer on a list
!
!   arguments :
!     in:
!              Nel     - an integer to be located
!              List    - a list of integers
!              Nlist   - length of the list
!     out:
!              Number  - the element location
!
!   required  routines -
!
!----------------------------------------------------------------------
!
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

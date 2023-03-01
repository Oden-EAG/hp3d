!----------------------------------------------------------------------
!
!   routine name       - decodg
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine decodes an integer into a sequence
!                        of digits
!
!   arguments :
!     in:
!               Nick   - a nickname to be decoded
!               Mod    - system number
!               N      - number of the digits stored in the nickname
!     out:
!               Narray - the decoded digits
!
!----------------------------------------------------------------------
subroutine decodg(Nick,Mod,N, Narray)
!
   implicit none
!
   integer, intent(in)  :: Nick,Mod,N
   integer, intent(out) :: Narray(N)
   integer :: i,nick1,nick2
!
   nick1 = Nick
   do i=1,N
      nick2 = nick1/Mod
      Narray(i) = nick1 - nick2*Mod
      nick1 = nick2
   enddo
!
end subroutine decodg

!----------------------------------------------------------------------
!
!   routine name       - decod
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine decodes an integer into a sequence
!                        of digits organized from left to right
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
subroutine decod(Nick,Mod,N, Narray)
!
   implicit none
!
   integer, intent(in)  :: Nick,Mod,N
   integer, intent(out) :: Narray(N)
!
!  ...local variables
   integer :: i,nick1,nick2
!
   nick1 = Nick
   do i=1,N-1
      nick2 = nick1/Mod
      Narray(N+1-i) = nick1 - nick2*Mod
      nick1 = nick2
   enddo
   Narray(1)=nick1
!
end subroutine decod
!
!
!----------------------------------------------------------------------
!
!   routine name       - decodLong
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine decodes a LONG integer into a sequence
!                        of digits organized from left to right
!
!   arguments :
!     in:
!               Nick   - a nickname to be decoded (LONG integer)
!               Mod    - system number
!               N      - number of the digits stored in the nickname
!     out:
!               Narray - the decoded digits
!
!----------------------------------------------------------------------
subroutine decodLong(Nick,Mod,N, Narray)
!
   implicit none
!
   integer(8), intent(in)  :: Nick
   integer,    intent(in)  :: Mod,N
   integer,    intent(out) :: Narray(N)
!
!..local variables
   integer    :: i
   integer(8) :: nick1,nick2,longmod
!
!..be careful with integer type - need to convert!
   longmod = INT(Mod,8)
   nick1 = Nick
   do i=1,N
      nick2 = nick1/longmod
!  ...be careful with integer type - need to convert!
      Narray(N+1-i) = INT(nick1 - nick2*longmod,4)
      nick1 = nick2
   enddo
!
end subroutine decodLong

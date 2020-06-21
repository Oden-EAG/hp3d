!-----------------------------------------------------------------------
!
!   routine name       - decode
!
!-----------------------------------------------------------------------
!
!   latest revision    - May 2020
!
!   purpose            - routine decodes a nickname consisting of
!                        2 digits
!
!   arguments :
!     in:       Nick   - a nickname
!     out:
!               J1     - the first  digit of the nickname
!               J2     - the second digit of the nickname
!
!-----------------------------------------------------------------------
subroutine decode(Nick, J1,J2)
   implicit none
   integer, intent(in)  :: Nick
   integer, intent(out) :: J1, J2
!
   J1 = Nick/10
   J2 = Nick - J1*10
!
end subroutine decode
!
!-----------------------------------------------------------------------
!
!   routine name       - ddecode
!
!-----------------------------------------------------------------------
!
!   latest revision    - May 2020
!
!   purpose            - routine decodes a nickname consisting of
!                        3 digits
!
!   arguments :
!     in:       Nick   - a nickname
!     out:
!               J1     - the first  digit of the nickname
!               J2     - the second digit of the nickname
!               J3     - the third  digit of the nickname
!
!-----------------------------------------------------------------------
subroutine ddecode(Nick, J1,J2,J3)
   implicit none
   integer, intent(in)  :: Nick
   integer, intent(out) :: J1, J2, J3
!
   integer :: nick2
!
   nick2 = Nick/10
   J3 = Nick - nick2*10
!
   J1 = nick2/10
   J2 = nick2 - J1*10
!
end subroutine ddecode

!----------------------------------------------------------------------
!
!   routine name       - decode2
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine decodes a nickname into two parts,
!                        using mod 100
!
!   arguments :
!     in:       Nick   - a nickname
!     out:
!               Nick1  - the first part of the nickname
!               Nick2  - the second two-digit part of the nickname
!
!----------------------------------------------------------------------
subroutine decode2(Nick, Nick1,Nick2)
!
   implicit none
!
   integer, intent(in)  :: Nick
   integer, intent(out) :: Nick1,Nick2
!
   Nick1 = Nick/100
   Nick2 = Nick - Nick1*100
!
end subroutine decode2

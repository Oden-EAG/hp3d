!----------------------------------------------------------------------
!
!   routine name       - encodg
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine encodes a sequence of digits
!                        into a nickname
!
!   arguments :
!     in:
!               Narray - digits to be encoded
!               Mod    - system number
!               N      - number of the digits stored in the nickname
!     out:
!               Nick   - the nickname
!
!----------------------------------------------------------------------
subroutine encodg(Narray,Mod,N, Nick)
!
   implicit none
!
   integer, intent(in)  :: Mod,N
   integer, intent(in)  :: Narray(N)
   integer, intent(out) :: Nick
!
   integer :: i
!
   Nick=Narray(N)
   do i=N-1,1,-1
      Nick = Nick*Mod + Narray(i)
   enddo
!
end subroutine encodg

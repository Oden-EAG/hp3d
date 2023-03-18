!----------------------------------------------------------------------
!
!   routine name       - encod
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine encodes a sequence of digits
!                        into a nickname from left to right
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
subroutine encod(Narray,Mod,N, Nick)
!
   implicit none
!
   integer, intent(in)  :: Mod,N
   integer, intent(in)  :: Narray(N)
   integer, intent(out) :: Nick
!
   integer :: i
!
   Nick=Narray(1)
   do i=2,N
      Nick = Nick*Mod + Narray(i)
   enddo
!
end subroutine encod
!
!
!----------------------------------------------------------------------
!
!   routine name       - encodLong
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine encodes a sequence of digits
!                        into a LONG nickname from left to right
!
!   arguments :
!     in:
!               Narray - digits to be encoded
!               Mod    - system number
!               N      - number of the digits stored in the nickname
!     out:
!               Nick   - the LONG (integer) nickname
!
!----------------------------------------------------------------------
subroutine encodLong(Narray,Mod,N, Nick)
!
   implicit none
!
   integer   , intent(in)  :: Mod,N
   integer   , intent(in)  :: Narray(N)
   integer(8), intent(out) :: Nick
!
!..local variables
   integer    :: i
   integer(8) :: longmod
!
!..be careful with integer type - need to convert!
   longmod = INT(Mod,8)
   Nick = INT(Narray(1),8)
   do i=2,N
!  ...be careful with integer type - need to convert!
      Nick = Nick*longmod + INT(Narray(i),8)
   enddo
!
end subroutine encodLong

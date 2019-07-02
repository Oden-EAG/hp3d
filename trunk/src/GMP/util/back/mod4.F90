!-----------------------------------------------------------------------
integer FUNCTION  mod4(i)
!-----------------------------------------------------------------------
!*** LATEST REVISION: Jul 09
!
!*** PURPOSE: function sets representative for equivalence class of 0 
!             equal to 4; hence mod4(n) belongs to set {1,2,3,4} instead 
!             of {0,1,2,3}, which is the case for function mod(n,3).
!-----------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------
  integer, intent(in) :: i
  mod4 = mod(i,4)
  if (mod4 .eq. 0) mod4 = 4
!
END FUNCTION mod4

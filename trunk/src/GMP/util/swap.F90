!-----------------------------------------------------------------------
!> @brief   routine swaps the values stored in two integers
!> @date    Oct 2024
!-----------------------------------------------------------------------
subroutine swap(I,J)
!
   implicit none
!
   integer, intent(inout) :: I,J
!
   integer                :: k
!
!-----------------------------------------------------------------------
!
   k=I ; I=J ; J=k
!
end subroutine swap

!-----------------------------------------------------------------------
!> @brief   routine computes the scalar product of two 3D vectors
!!
!> @param[in]  A    - 1st vector
!> @param[in]  B    - 2nd vector
!> @param[out] Prod - dot product
!!
!> @date    Oct 2024
!!
!> @note    the intrinsic Fortran function dot_product( , ) could be
!!          used instead
!-----------------------------------------------------------------------
subroutine scalar_product(A,B, Prod)
!
   implicit none
!
   real(8), intent(in)  :: A(3),B(3)
   real(8), intent(out) :: Prod
!
!-----------------------------------------------------------------------

   Prod = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)
!
end subroutine scalar_product

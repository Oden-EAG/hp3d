!-------------------------------------------------------------------------
!> @brief   routine computes Euclidean norm of a 3D vector
!
!> @param[in]  Vec   - 3D vector
!> @param[out] Rnorm - norm of the vector
!
!> @date    Oct 2024
!-------------------------------------------------------------------------
subroutine norm(Vec, Rnorm)
!
   implicit none
!
   real(8), intent(in ) :: Vec(3)
   real(8), intent(out) :: Rnorm
!
   real(8) :: s
!
!-------------------------------------------------------------------------
!
   call scalar_product(Vec(1:3),Vec(1:3), s) ; Rnorm=sqrt(s)
!
end subroutine norm

!-------------------------------------------------------------------------
!> @Purpose : routine computes Euclidean norm of a 3D vector
!
!> @param[in]  Vec   - 3D vector
!> @param[out] Rnorm - norm of the vector
!
!> @revision Feb 13
!-------------------------------------------------------------------------
subroutine norm(Vec, Rnorm)
!  
      implicit none
      real*8,dimension(3),intent(in ) :: Vec
      real*8,             intent(out) :: Rnorm
!
      real*8 :: s
!-------------------------------------------------------------------------
!
      call scalar_product(Vec(1:3),Vec(1:3), s) ; Rnorm=sqrt(s)
!
!
endsubroutine norm



!-----------------------------------------------------------------------
!> @Purpose : routine computes the scalar product of two 3D vectors
!
!> @param[in]  Vec1 - 1st vector
!> @param[in]  Vec2 - 2nd vector
!> @param[out] Prod - dot product
!
!> @revision Feb 13
!
!! Remark : the intrinsic Fortran fuction dot_product( , ) could be 
!!          used instead
!-----------------------------------------------------------------------
subroutine scalar_product(Vec1,Vec2, Prod)
!
      implicit none
      real*8,dimension(3),intent(in ) :: Vec1,Vec2
      real*8,             intent(out) :: Prod
!
      integer :: i
!-----------------------------------------------------------------------
!  
      Prod=0.d0
      do i=1,3 ; Prod=Prod+Vec1(i)*Vec2(i) ; enddo
!      
endsubroutine scalar_product

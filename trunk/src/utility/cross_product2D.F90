subroutine cross_product2D(Vec1,Vec2, Vec_result)
      implicit none
      real*8,dimension(2),intent(in ) :: Vec1,Vec2
      real*8,             intent(out) :: Vec_result
!      
      Vec_result=Vec1(1)*Vec2(2)-Vec1(2)*Vec2(1)
!      
endsubroutine cross_product2D

subroutine swap(I,J)
      implicit none
      integer, intent(inout) :: I,J
      integer                :: k
!      
      k=I ; I=J ; J=k
!      
endsubroutine swap

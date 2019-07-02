!--------------------------------------------------------------------      
subroutine normalize(Vec)
!
      implicit none
      real*8,dimension(3),intent(inout) :: Vec(3)
!
      real*8 :: s
!--------------------------------------------------------------------      
!
      call scalar_product(Vec(1:3),Vec(1:3), s) ; s=sqrt(s)
      Vec(1:3)=Vec(1:3)/s
!      
endsubroutine normalize
!
!
!
!--------------------------------------------------------------------      
subroutine normalize_real(Vec)
!
      implicit none
      real*8,dimension(3),intent(inout) :: Vec
!--------------------------------------------------------------------      
!
!  ...redirect to standard routine
      call normalize(Vec(1:3))
!      
endsubroutine normalize_real
!
!
!
!--------------------------------------------------------------------      
subroutine normalize_complex(Vec)
!
      implicit none
      complex*16,dimension(3),intent(inout) :: Vec
      real*8 :: rnorm
!--------------------------------------------------------------------      
!
!  ...use intrinsic function
      rnorm=sqrt( dot_product(Vec(1:3),Vec(1:3)) )
      Vec(1:3)=Vec(1:3)/rnorm
!      
endsubroutine normalize_complex



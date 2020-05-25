subroutine mixed_product(Vect1,Vect2,Vect3, Prod)
!
      implicit none
      real(8),dimension(3),intent(in ) :: Vect1,Vect2,Vect3
      real(8),             intent(out) :: Prod
!      
      real(8),dimension(3) :: vect
!
      call  cross_product(Vect1,Vect2, vect)
      call scalar_product(vect ,Vect3, Prod)
!
!
end subroutine mixed_product


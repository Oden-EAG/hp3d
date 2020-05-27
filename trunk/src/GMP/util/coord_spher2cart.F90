!--------------------------------------------------------------------------
!> @Purpose : convert spherical coordinates to Cartasian coordinates
!
!! Spherical coordinates:
!!   r     : radius
!!   theta : angle from positive z axis
!!   phi   : angle from positive x axis
!!
!> @param[in]  C_in  - spherical ( r, theta, phi )
!> @param[out] C_out - Cartesian ( x, y, z )
!
!> rev@Feb 13
!--------------------------------------------------------------------------
subroutine coord_spher2cart(C_in, C_out)
  !
  implicit none
  real(8),dimension(3),intent(in ) :: C_in
  real(8),dimension(3),intent(out) :: C_out
  call coord_sphere2cart(C_in, C_out)
end subroutine coord_spher2cart
!
!
subroutine coord_sphere2cart(C_in, C_out)
  implicit none
  real(8),dimension(3),intent(in ) :: C_in
  real(8),dimension(3),intent(out) :: C_out
  C_out(1) = C_in(1) * sin(C_in(2)) * cos(C_in(3))
  C_out(2) = C_in(1) * sin(C_in(2)) * sin(C_in(3))
  C_out(3) = C_in(1) * cos(C_in(2))
end subroutine coord_sphere2cart


!--------------------------------------------------------------------------
!> @Purpose : convert Cartesian coordinates to spherical coordinates
!
!! Spherical coordinates: 
!!   r     : radius
!!   theta : angle from positive z axis
!!   phi   : angle from positive x axis
!!
!> @param[in]  C_in  - Cartesian ( x, y, z )
!> @param[out] C_out - spherical ( r, theta, phi )
!
!> rev@Feb 13
!--------------------------------------------------------------------------
subroutine coord_cart2spher(C_in, C_out)
  !
  implicit none
  real*8,dimension(3),intent(in)  :: C_in
  real*8,dimension(3),intent(out) :: C_out
  !--------------------------------------------------------------------------
  !
  !  ...call old routine
  call coord_cart2sphere(C_in, C_out)
  !
end subroutine coord_cart2spher
!
!
!
subroutine coord_cart2sphere(C_in, C_out)
  implicit none
  real*8,dimension(3),intent(in)  :: C_in
  real*8,dimension(3),intent(out) :: C_out
  !    
  real*8, parameter :: eps=1.0e-15
  !--------------------------------------------------------------------------
  
  C_out(1) = sqrt(C_in(1)**2 + C_in(2)**2 + C_in(3)**2) 
  
  if (C_out(1) > eps) then
     C_out(2) = acos(C_in(3)/C_out(1))
     C_out(3) = atan2(C_in(2), C_in(1))
  else
     C_out(2:3) = 0.d0
  endif
  !
  !
end subroutine coord_cart2sphere

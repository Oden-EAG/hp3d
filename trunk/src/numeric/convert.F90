!> Purpose : convert unit vector from sphere to cartesian
!! @param Theta[in] - angle from positive z axis
!! @param Phi  [in] - angle from positive x axis
!! @param T    [out]- transform matrix
subroutine unit_cart2sphere(Theta, Phi, T)
  implicit none
  ! input arguments
  real(8), intent(in)                  :: Theta, Phi
  real(8), intent(out), dimension(3,3) :: T
  real(8) :: siTh, coTh, siPh, coPh
  
  siTh = sin(Theta); coTh = cos(Theta);
  siPh = sin(Phi);   coPh = cos(Phi);

  T(1,1)   =  siTh*coPh;      T(1,2)   =  siTh*siPh;       T(1,3)   =  coTh;
  T(2,1)   =  coTh*coPh;      T(2,2)   =  coTh*siPh;       T(2,3)   = -siTh;
  T(3,1)   = -siPh;           T(3,2)   =  coPh;            T(3,3)   =  0.d0;

end subroutine unit_cart2sphere

!> Purpose : convert gradient operator
!! @param R    [in] - radius
!! @param Theta[in] - angle from positive z axis
!! @param Phi  [in] - angle from positive x axis
!! @param T    [out]- transform matrix
subroutine grad_sphere2cart(R, Theta, Phi, T)
  implicit none
  ! input arguments
  real(8), intent(in)                  :: R, Theta, Phi
  real(8), intent(out), dimension(3,3) :: T
  real(8) :: siTh, coTh, siPh, coPh
  real(8), parameter :: epsilon = 1.0e-14, pi=4.*ATAN(1.D0)


  ! http://mathworld.wolfram.com/SphericalCoordinates.html
  if (Theta.lt.epsilon) then
     siTh = sin(epsilon)
  else
     siTh = sin(Theta)
  end if
  coTh = cos(Theta);
  siPh = sin(Phi);   
  coPh = cos(Phi);

  T(1,1)   =  coPh*siTh;    T(1,2)   =  coPh*coTh/R;     T(1,3)   = -siPh/(R*siTh)
  T(2,1)   =  siPh*siTh;    T(2,2)   =  siPh*coTh/R;     T(2,3)   =  coPh/(R*siTh)
  T(3,1)   =  coTh;         T(3,2)   =  -siTh/R;         T(3,3)   =  0.d0

end subroutine grad_sphere2cart


!> Purpose : convert displacement and its derivative from sphere to cart
!! @param S   [in]  - position vector in spherical coord
!! @param Zs  [in]  - displacement in spherical coord
!! @param Zds [in]  - derivative in spherical coord
!! @param Zc  [out] - displacement in cartesian coord
!! @param Zdc [out] - derivative in cartesian coord
#include "implicit_none.h"
subroutine disp_sphere2cart(S, Zs, Zds, Zc, Zdc)
  use parameters
  implicit none
  real(8),dimension(3),   intent(in)  :: S
  VTYPE,  dimension(3),   intent(in)  :: Zs
  VTYPE,  dimension(3,3), intent(in)  :: Zds
  VTYPE,  dimension(3),   intent(out) :: Zc
  VTYPE,  dimension(3,3), intent(out) :: Zdc
  
  integer                :: k1, k2 
  real(8),dimension(3,3) :: rot_u, rot_g
  VTYPE,  dimension(3,3) :: ztmp
  real(8)                :: coTh, siTh, coPh, siPh

  coTh = cos(S(2));
  siTh = sin(S(2));
  coPh = cos(S(3));
  siPh = sin(S(3));

  ! change of coordinates of spherical one to cartesian
  call unit_cart2sphere(S(2), S(3), rot_u)
  Zc = ZERO
  do k1=1,3
     Zc(k1) = Zc(k1) + dot_product(rot_u(1:3,k1),Zs(1:3))
  end do

  ! derivative is more tricky
  do k2=1,3
     do k1=1,3
        ztmp(k1,k2) = dot_product(rot_u(1:3,k1), Zds(1:3,k2))
     end do
  end do

  ! add its derivative r, theta and phi
  ztmp(1,2) = ztmp(1,2) + coTh*coPh*Zs(1) - siTh*coPh*Zs(2)
  ztmp(1,3) = ztmp(1,3) - siTh*siPh*Zs(1) - coTh*siPh*Zs(2) - coPh*Zs(3)
  
  ztmp(2,2) = ztmp(2,2) + coTh*siPh*Zs(1) - siTh*siPh*Zs(2)
  ztmp(2,3) = ztmp(2,3) + siTh*coPh*Zs(1) + coTh*coPh*Zs(2) - siPh*Zs(3)

  ztmp(3,2) = ztmp(3,2) -      siTh*Zs(1) -      coTh*Zs(2)

  call grad_sphere2cart(S(1), S(2), S(3), rot_g)

  do k1=1,3
     Zdc(1,k1) = dot_product(rot_g(1,1:3), ztmp(k1,1:3))
     Zdc(2,k1) = dot_product(rot_g(2,1:3), ztmp(k1,1:3))
     Zdc(3,k1) = dot_product(rot_g(3,1:3), ztmp(k1,1:3))
  end do
  
end subroutine disp_sphere2cart


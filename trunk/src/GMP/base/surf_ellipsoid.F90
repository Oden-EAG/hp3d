!------------------------------------------------------------------------------------
!> Purpose : ellipsoid parameterization
!!
!! @param[in]  X       - physical coordinates of a point
!! @param[in]  Center  - center
!! @param[in]  Rs1,2,3 - semiaxis
!! @param[out] Fval    - function value
!! @param[out] Dfdx    - function gradient
!!
!! @revision Nov 12
!------------------------------------------------------------------------------------
!
subroutine ellipsoid(X,Center,Rs1,Rs2,Rs3,  Fval,Dfdx)
!
      implicit none
      real(8), dimension(3), intent(in ) :: X,Center
      real(8)              , intent(in ) :: Rs1,Rs2,Rs3
      real(8)              , intent(out) :: Fval
      real(8), dimension(3), intent(out) :: Dfdx
!------------------------------------------------------------------------------------
!
      Fval =  ((X(1)-Center(1))/Rs1)**2 +        &
              ((X(2)-Center(2))/Rs2)**2 +        &
              ((X(3)-Center(3))/Rs3)**2 - 1.d0
!
      Dfdx(1) = (X(1)-Center(1))*2 /Rs1**2
      Dfdx(2) = (X(2)-Center(2))*2 /Rs2**2
      Dfdx(3) = (X(3)-Center(3))*2 /Rs3**2
!
!      
end subroutine ellipsoid

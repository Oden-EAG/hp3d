!----------------------------------------------------------------------
!> @Purpose : routine switches from 2D Cartesian to polar coordinates
!
!> @param[in]  X       - Cartesian coordinates of a point
!> @param[out] R,Theta - polar coordinates of the point
!
!> rev@Feb 13
!
!  REMARK : -pi <= Theta <= pi
!----------------------------------------------------------------------
subroutine coord_cart2polar(X, R,Theta)
!
      implicit none
!
      real(8), intent(in ) :: X(2)
      real(8), intent(out) :: R,Theta
!----------------------------------------------------------------------
!
!  ...call old routine
      call cart_to_polar(X, R,Theta)
!
!
end subroutine coord_cart2polar
!
!
!
!  ...for backward compatibility
subroutine cart_to_polar(Xp, R,Theta)
!
      implicit none
      real(8), intent(in ) :: Xp(2)
      real(8), intent(out) :: R,Theta
!
      real(8) :: x,y
      real(8), parameter :: PI = acos(-1.d0)
      real(8), parameter :: eps=1.d-13
      integer :: iprint
!----------------------------------------------------------------------
!
      iprint=0
!
      x=Xp(1) ; y=Xp(2)
!
      R=sqrt(x**2+y**2)
      if (R < eps) then
        Theta=0.d0
      else
        if (y == 0.d0) then
          if (x < 0.d0) then
            Theta=PI
          elseif (x > 0.d0) then
            Theta=0.d0
          endif
        elseif (y > 0.d0) then
          Theta= acos(x/R)
        else
          Theta=-acos(x/R)
        endif
      endif
!
      if (iprint.eq.1) then
        write(*,7001) x,y,R,Theta
 7001   format(' cart_to_polar: x,y,R,Theta = ',2f8.3,3x,2f8.3)
        call pause
      endif
!
!
end subroutine cart_to_polar

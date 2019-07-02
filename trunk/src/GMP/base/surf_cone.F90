!------------------------------------------------------------------------------------
!> Purpose : cone parameterization
!!
!! @param[in]  X      - physical coordinates of a point
!! @param[in]  Center - center
!! @param[in]  Dirvec - direction vector
!! @param[in]  Alpha  - aperture
!! @param[out] Fval   - function value
!! @param[out] Dfdx   - function gradient
!!
!! @revision Nov 12
!------------------------------------------------------------------------------------
!
subroutine cone(X,Center,Dirvec,Alpha,  Fval,Dfdx)
!
      implicit none
      real*8, dimension(3), intent(in ) :: X,Center,Dirvec
      real*8              , intent(in ) :: Alpha
      real*8              , intent(out) :: Fval
      real*8, dimension(3), intent(out) :: Dfdx
!------------------------------------------------------------------------------------
!  ...axis versor, relative position of X wrt cone vertex
      real*8, dimension(3) :: cvec,xvec
!  ...derivative of cone coordinate wrt Cartesian coordinates
      real*8, dimension(3) :: dzetdx,dr2dx
      real*8               :: s,rzet,r2,a
      integer              :: iprint
      real*8, parameter    :: eps=1.d-15
!------------------------------------------------------------------------------------
!
      iprint=0
      if (iprint.eq.1) then
        write(*,7001) X,Center,Dirvec,Alpha
 7001   format(' cone: X,Center,Dirvec,Alpha = ',3f8.3,2x,3f8.3,2x,3f8.3,2x,f8.3)
      endif
!
!  ...normalize the cone vector
      s = sqrt(Dirvec(1)**2+Dirvec(2)**2+Dirvec(3)**2)
      cvec(1:3) = Dirvec(1:3)/s
!
!  ...relative position vector
      xvec(1:3) = X(1:3) - Center(1:3)
!
!  ...derivatives of cone zeta coordinate wrt Cartesian coordinates
      rzet = xvec(1)*cvec(1)+xvec(2)*cvec(2)+xvec(3)*cvec(3)
      dzetdx(1:3) = cvec(1:3)
!
!  ...cone radial coordinate (squared) and its derivatives wrt Cartesian coordinates
      r2 = xvec(1)**2+xvec(2)**2+xvec(3)**2 - rzet**2
      dr2dx(1:3) = 2.d0*(xvec(1:3) - rzet*cvec(1:3))
!
!  ...printing
      if (iprint.eq.1) then
        write(*,7002) cvec,xvec,rzet,r2
 7002   format(' cone: cvec,xvec,rzet,r2 = ',3f8.3,2x,3f8.3,2x,2f8.3)
      endif

!  ...implicit cone function value
      a = tan(Alpha)**2
      Fval = r2 - a*rzet**2
!
!  ...gradient
      if (abs(rzet).le.eps) then
!  .....gradient not defined at vertex (return a large value)
        Dfdx(1:3) = 1.d20
      else
        Dfdx(1:3) = dr2dx(1:3) - 2.d0*a*rzet*dzetdx(1:3)
      endif
!
!  ...printing
      if (iprint.eq.1) then
        write(*,7003) Fval,Dfdx
 7003   format(' cone: Fval,Dfdx = ',e12.5,3x,3e12.5)
        call pause
      endif
!
!
endsubroutine cone



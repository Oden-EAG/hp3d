!------------------------------------------------------------------------------------
!> @brief superellipse parameterization
!!
!> @param[in]  X     - physical coordinates of a point
!> @param[in]  Point - a point on the axis
!> @param[in]  Axis1 - first semiaxis vector
!> @param[in]  Axis2 - second semiaxis vector
!> @param[in]  Rs1   - first semiaxis radius
!> @param[in]  Rs2   - second semiaxis radius
!> @param[in]  Pw1   - power of first coordinate
!> @param[in]  Pw2   - power of second coordinate
!> @param[out] Fval  - function value
!> @param[out] Dfdx  - function gradient
!!
!> @date Feb 18
!------------------------------------------------------------------------------------
subroutine superellipse(X,Point,Axis1,Axis2,Rs1,Rs2,Pw1,Pw2, Fval,Dfdx)
!
      implicit none
!
      real(8), dimension(3), intent(in ) :: X,Point,Axis1,Axis2
      real(8)              , intent(in ) :: Rs1,Rs2,Pw1,Pw2
      real(8)              , intent(out) :: Fval
      real(8), dimension(3), intent(out) :: Dfdx
!------------------------------------------------------------------------------------
!  ...superellipse axes unit vectors
      real(8), dimension(3)   :: unit1,unit2
!  ...relative position vector of X with respect to point of axis
      real(8), dimension(3)   :: vecV
!  ...projections onto unit axes
      real(8)                 :: xx,yy,sgnxx,sgnyy,s
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!------------------------------------------------------------------------------------
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7000) Point(1:3),Axis1(1:3),Axis2(1:3),Rs1,Rs2,Pw1,Pw2,X(1:3)
 7000   format(' superellipse: Point = ',3e12.5,/, &
               '               Axis1 = ',3e12.5,/, &
               '               Axis2 = ',3e12.5,/, &
               '               Rs1   = ', e12.5,/, &
               '               Rs2   = ', e12.5,/, &
               '               Pw1   = ', e12.5,/, &
               '               Pw2   = ', e12.5,/, &
               '               X     = ',3e12.5    )
      endif
#endif
!
!  ...superellipse has two orthogonal semiaxes, each with a radius and power
!     associated to it. In 2D it satisfies
!                   abs(xx/Rs1)**Pw1+abs(yy/Rs2)**Pw2=1
!     thereofore, in 2D when Rs1=Rs2 and Pw1=Pw2=p, it is an Lp ball, and in 3D if
!     Rs1=Rs2 and Pw1=Pw2=2, it is a cylinder.
!
!  ...evaluate axis unit vectors
      call norm(Axis1, s)
      unit1(1:3) = Axis1(1:3)/s
      call norm(Axis2, s)
      unit2(1:3) = Axis2(1:3)/s
!
!  ...evaluate the relative position vector of X
      vecV(1:3) = X(1:3) - Point(1:3)
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7001) unit1(1:3),unit2(1:3),vecV(1:3)
 7001   format(' superellipse: unit1 = ',3e12.5,' unit2 = ',3e12.5,' vecV = ',3e12.5)
      endif
#endif
!
!  ...project vecV to the relevant axes
      call scalar_product(vecV,unit1, xx)
      sgnxx = sign(1.d0,xx)
      call scalar_product(vecV,unit2, yy)
      sgnyy = sign(1.d0,yy)
!
!  ...function value and derivative
      Fval = (abs(xx)/Rs1)**Pw1+(abs(yy)/Rs2)**Pw2-1.d0
      Dfdx(1:3) = ((1.d0/Rs1)**Pw1)*Pw1*(abs(xx)**(Pw1-1))*unit1(1:3)*sgnxx &
                + ((1.d0/Rs2)**Pw2)*Pw2*(abs(yy)**(Pw2-1))*unit2(1:3)*sgnyy
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7004) Fval,Dfdx(1:3)
 7004   format(' superellipse: Fval,Dfdx = ',e12.5,2x,3e12.5)
        call pause
      endif
#endif
!
end subroutine superellipse

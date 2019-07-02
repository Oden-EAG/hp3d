!------------------------------------------------------------------------------------
!> Purpose : cylinder parameterization
!!
!! @param[in]  X     - physical coordinates of a point
!! @param[in]  Point - a point on the axis
!! @param[in]  Cvect - direction vector
!! @param[in]  Rad   - radius
!! @param[out] Fval  - function value
!! @param[out] Dfdx  - function gradient
!!
!! @revision Nov 12
!------------------------------------------------------------------------------------
!
subroutine cylinder(X,Point,Cvect,Rad, Fval,Dfdx)
!
      implicit none
      real*8, dimension(3), intent(in ) :: X,Point,Cvect
      real*8              , intent(in ) :: Rad
      real*8              , intent(out) :: Fval
      real*8, dimension(3), intent(out) :: Dfdx
!------------------------------------------------------------------------------------
!  ...cylinder axis unit vector
      real*8, dimension(3)   :: unit
!  ...relative position vector of X and its projection onto the plane normal to the
!     cylinder axis
      real*8, dimension(3)   :: xvec,xvecp
!  ...gradient of xvecp
      real*8, dimension(3,3) :: dxvecpdx
      integer                :: iprint,i
      real*8                 :: s,s1
!------------------------------------------------------------------------------------
!      
      iprint=0
!  ...printing      
      if (iprint.eq.1) then
        write(*,7000) Point(1:3),Cvect(1:3),Rad,X(1:3)
 7000   format(' cylinder: Point = ',3e12.5,/, &
               '           Cvect = ',3e12.5,/, &
               '           Rad   = ', e12.5,/, &
               '           X     = ',3e12.5    )
      endif
!
!  ...evaluate cylinder unit vector
      call norm(Cvect, s)
      unit(1:3) = Cvect(1:3)/s
!
!  ...evaluate the relative position vector of X
      xvec(1:3) = X(1:3) - Point(1:3)
!
!  ...printing      
      if (iprint.eq.1) then
        write(*,7001) unit(1:3),xvec(1:3)
 7001   format(' cylinder: unit = ',3e12.5,' xvec = ',3e12.5)
      endif
!
!  ...project the relative position vector onto the plane normal to the
!     cylinder axis
      call scalar_product(xvec,unit, s)
      xvecp(1:3) = xvec(1:3) - s*unit(1:3)
!
!  ...printing      
      if (iprint.eq.1) then
        write(*,7002) xvecp(1:3)
 7002   format('cylinder: xvecp = ',3e12.5)
      endif
!
!  ...compute the gradient of the projection vector
      do i=1,3
        dxvecpdx(1:3,i) = -unit(i)*unit(1:3)
        dxvecpdx(i,i) = dxvecpdx(i,i) + 1
      enddo
!
      call scalar_product(xvecp,xvecp, s)
      if (iprint.eq.1) then
        call norm(xvecp, s1)
        write(*,7003) s1
 7003   format(' cylinder: |xvecp| = ',e12.5)
      endif
      Fval = s - Rad**2
      do i=1,3
        call scalar_product(xvecp,dxvecpdx(1:3,i), s)
        Dfdx(i) = 2.d0*s
      enddo
      if (iprint.eq.1) then
        write(*,7004) Fval,Dfdx(1:3)
 7004   format(' cylinder: Fval,Dfdx = ',e12.5,2x,3e12.5)
        call pause
      endif
!
!
endsubroutine cylinder

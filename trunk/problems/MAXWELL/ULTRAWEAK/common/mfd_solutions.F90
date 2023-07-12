!
!----------------------------------------------------------------------
!
!     routine name      - maxwell_solution
!
!----------------------------------------------------------------------
!
!     latest revision:  - June 18
!
!     purpose:          - compute all relevant quantities of the exact
!                         solutions for the Maxwell problem
!
!     arguments:
!
!     in:
!             Xp        - a point in physical space
!     out:
!              p        - value of the solution pressure
!          Gradp        - corresponding first derivatives
!         Grad2p        - corresponding second derivatives
!
!----------------------------------------------------------------------
!
   subroutine mfd_solutions(Xp, p,Gradp,Grad2p)
!
      use data_structure3D
      use commonParam
      use parameters, only : ZERO, ZONE, ZIMG
!
      implicit none
!
      real(8),    intent(in)  :: Xp(3)
      complex(8), intent(out) :: p
      complex(8), intent(out) :: Gradp(3)
      complex(8), intent(out) :: Grad2p(3,3)
!
!  ...intermediate variables
      real(8) :: rl, w0, p0, z_R, rk, pi_mod
      real(8) :: theta_x, theta_y, theta_z
      real(8) :: rad_x, rad_y, rad_z
      real(8) :: sinx, siny, cosx, cosy, sinz, cosz
      real(8) :: x, y, z
      real(8) :: Rx(3,3), Ry(3,3), Rz(3,3), Rxy(3,3), Rxyz(3,3)
!
!  ...functions and their derivatives
      complex(8) :: zf_x, zf_y, zf_z, cn
      complex(8) :: dzf_x,dzf_y, dzf_z
      complex(8) :: ddzf_x, ddzf_y, ddzf_z
!
      complex(8) :: p_x, p_xx, p_xy, p_xz
      complex(8) :: p_y, p_yx, p_yy, p_yz
      complex(8) :: p_z, p_zx, p_zy, p_zz
!
      real(8) :: r_x, r_xx, r_xy, r_xz
      real(8) :: r_y, r_yx, r_yy, r_yz
      real(8) :: r_z, r_zx, r_zy, r_zz
!
      real(8) :: dxdxs, dxdys, dxdzs
      real(8) :: dydxs, dydys, dydzs
      real(8) :: dzdxs, dzdys, dzdzs
!
      real(8) :: x1, x2, x3, xshift, yshift, zshift
      real(8) :: r, alpha, zm
!
!TODO: Fix this
      integer :: iprob
!
!---------------------------------------------------------------------------------------
!
      x1 = Xp(1); x2 = Xp(2); x3 = Xp(3)
!
      iprob = Isol
!
      select case (iprob)
!
!  ...polynomial of order NPX*NPY*NPZ
      case(0)
!
         p = ZERO
         Gradp = ZERO
         Grad2p = ZERO
!
!..sin solution
      case(1)
!
         pi_mod = OMEGA
         p =dsin(x1*pi_mod)*dsin(x2*pi_mod)*dsin(x3*pi_mod)*(1.0D0,1.0D0)
!
!     ...1st order derivatives
         Gradp(1) = pi_mod*dcos(x1*pi_mod)*dsin(x2*pi_mod)*dsin(x3*pi_mod)*(1.0D0,1.0D0)
         Gradp(2) = pi_mod*dcos(x2*pi_mod)*dsin(x1*pi_mod)*dsin(x3*pi_mod)*(1.0D0,1.0D0)
         Gradp(3) = pi_mod*dcos(x3*pi_mod)*dsin(x1*pi_mod)*dsin(x2*pi_mod)*(1.0D0,1.0D0)
!
!     ...2nd derivative (3,3) matrix - Hessian
         Grad2p(1,1) = -pi_mod**2*dsin(x1*pi_mod)*dsin(x2*pi_mod)*dsin(x3*pi_mod)*(1.0D0,1.0D0)
         Grad2p(1,2) =  pi_mod**2*dcos(x1*pi_mod)*dcos(x2*pi_mod)*dsin(x3*pi_mod)*(1.0D0,1.0D0)
         Grad2p(1,3) =  pi_mod**2*dcos(x1*pi_mod)*dcos(x3*pi_mod)*dsin(x2*pi_mod)*(1.0D0,1.0D0)
         Grad2p(2,1) =  pi_mod**2*dcos(x1*pi_mod)*dcos(x2*pi_mod)*dsin(x3*pi_mod)*(1.0D0,1.0D0)
         Grad2p(2,2) = -pi_mod**2*dsin(x1*pi_mod)*dsin(x2*pi_mod)*dsin(x3*pi_mod)*(1.0D0,1.0D0)
         Grad2p(2,3) =  pi_mod**2*dcos(x2*pi_mod)*dcos(x3*pi_mod)*dsin(x1*pi_mod)*(1.0D0,1.0D0)
         Grad2p(3,1) =  pi_mod**2*dcos(x1*pi_mod)*dcos(x3*pi_mod)*dsin(x2*pi_mod)*(1.0D0,1.0D0)
         Grad2p(3,2) =  pi_mod**2*dcos(x2*pi_mod)*dcos(x3*pi_mod)*dsin(x1*pi_mod)*(1.0D0,1.0D0)
         Grad2p(3,3) = -pi_mod**2*dsin(x1*pi_mod)*dsin(x2*pi_mod)*dsin(x3*pi_mod)*(1.0D0,1.0D0)
!
!  ...polynomial solution vanishing on the boundary
      case(2)
!
         cn=1.D2*ZONE;     !  cn=1.D2*(ZONE+ZIMG);
!     ...displacement
         p = -cn*x1*x2*x3*(x1-1.0D0)*(x2-1.0D0)*(x3-1.0D0)
!
!     ...1st order derivatives
         Gradp(1) =-cn*x2*x3*(x1*2.0D0-1.0D0)*(x2-1.0D0)*(x3-1.0D0)
         Gradp(2) =-cn*x1*x3*(x2*2.0D0-1.0D0)*(x1-1.0D0)*(x3-1.0D0)
         Gradp(3) =-cn*x1*x2*(x3*2.0D0-1.0D0)*(x1-1.0D0)*(x2-1.0D0)
!     ...2nd derivative (3,3) matrix - hessian
         Grad2p(1,1) =cn*x2*x3*(x2-1.0D0)*(x3-1.0D0)*(-2.0D0)
         Grad2p(1,2) =-cn*x3*(x1*2.0D0-1.0D0)*(x2*2.0D0-1.0D0)*(x3-1.0D0)
         Grad2p(1,3) =-cn*x2*(x1*2.0D0-1.0D0)*(x3*2.0D0-1.0D0)*(x2-1.0D0)
         Grad2p(2,1) =-cn*x3*(x1*2.0D0-1.0D0)*(x2*2.0D0-1.0D0)*(x3-1.0D0)
         Grad2p(2,2) =cn*x1*x3*(x1-1.0D0)*(x3-1.0D0)*(-2.0D0)
         Grad2p(2,3) =-cn*x1*(x2*2.0D0-1.0D0)*(x3*2.0D0-1.0D0)*(x1-1.0D0)
         Grad2p(3,1) =-cn*x2*(x1*2.0D0-1.0D0)*(x3*2.0D0-1.0D0)*(x2-1.0D0)
         Grad2p(3,2) =-cn*x1*(x2*2.0D0-1.0D0)*(x3*2.0D0-1.0D0)*(x1-1.0D0)
         Grad2p(3,3) =cn*x1*x2*(x1-1.0D0)*(x2-1.0D0)*(-2.0D0)
!
!  ...plane wave
      case(3)
!
         zf_x = exp(-ZIMG*OMEGA*x1)
         zf_y = exp(-ZIMG*OMEGA*x2)
         zf_z = exp(-ZIMG*OMEGA*x3)
!
!     ...first order derivatives
         dzf_x = -ZIMG*OMEGA*zf_x
         dzf_y = -ZIMG*OMEGA*zf_y
         dzf_z = -ZIMG*OMEGA*zf_z
!
!  ...second order derivatives
         ddzf_x = -OMEGA**2*zf_x
         ddzf_y = -OMEGA**2*zf_y
         ddzf_z = -OMEGA**2*zf_z
!
         p = zf_x*zf_y*zf_z

!  ...1st order derivatives
         Gradp(1) = dzf_x*zf_y*zf_z
         Gradp(2) = zf_x*dzf_y*zf_z
         Gradp(3) = zf_x*zf_y*dzf_z
!
         Grad2p(1,1) = ddzf_x*zf_y*zf_z
         Grad2p(1,2) = dzf_x*dzf_y*zf_z
         Grad2p(1,3) = dzf_x*zf_y*dzf_z
         Grad2p(2,1) = dzf_x*dzf_y*zf_z
         Grad2p(2,2) = zf_x*ddzf_y*zf_z
         Grad2p(2,3) = zf_x*dzf_y*dzf_z
         Grad2p(3,1) = dzf_x*zf_y*dzf_z
         Grad2p(3,2) = zf_x*dzf_y*dzf_z
         Grad2p(3,3) = zf_x*zf_y*ddzf_z
!
!  ...point source
      case(4)
!
!     ...shift origin
         xshift = x1 + 0.1d0; yshift = x2 + 0.1d0; zshift = x3 + 0.1d0
         r = dsqrt(xshift**2 + yshift**2 + zshift**2)
         r_x = xshift/r ; r_y = yshift/r ;   r_z = zshift/r
!
         r_xx = 1.d0/r*(1.d0-r_x**2)
         r_xy = -r_x/r*r_y
         r_xz = -r_x/r*r_z
         r_yx = r_xy
         r_yy = 1.d0/r*(1.d0-r_y**2)
         r_yz = -r_y/r*r_z
         r_zx = r_xz
         r_zy = r_yz
         r_zz = 1.d0/r*(1.d0-r_z**2)
!
         alpha = -1.d0/(4.d0*PI)
!
         p = alpha*exp(ZIMG*OMEGA*r)/r
!
         Gradp(1)    = -r_x*(-ZIMG*OMEGA+1.d0/r) * p
         Gradp(2)    = -r_y*(-ZIMG*OMEGA+1.d0/r) * p
         Gradp(3)    = -r_z*(-ZIMG*OMEGA+1.d0/r) * p
!
         Grad2p(1,1) = (r_x**2/r**2-r_xx*(-ZIMG*OMEGA+1.D0/r))*p   &
                     -  r_x*(-ZIMG*OMEGA+1.d0/r)*Gradp(1)
         Grad2p(1,2) = (r_x*r_y/r**2-r_xy*(-ZIMG*OMEGA+1.D0/r))*p   &
                     -  r_x*(-ZIMG*OMEGA+1.d0/r)*Gradp(2)
         Grad2p(1,3) = (r_x*r_z/r**2-r_xz*(-ZIMG*OMEGA+1.D0/r))*p   &
                     -  r_x*(-ZIMG*OMEGA+1.d0/r)*Gradp(3)
         Grad2p(2,1) = Grad2p(1,2)
         Grad2p(2,2) = (r_y**2/r**2-r_yy*(-ZIMG*OMEGA+1.D0/r))*p   &
                     -  r_y*(-ZIMG*OMEGA+1.d0/r)*Gradp(2)
         Grad2p(2,3) = (r_y*r_z/r**2-r_yz*(-ZIMG*OMEGA+1.D0/r))*p   &
                     -  r_y*(-ZIMG*OMEGA+1.d0/r)*Gradp(3)
         Grad2p(3,1) = Grad2p(1,3)
         Grad2p(3,2) = Grad2p(2,3)
         Grad2p(3,3) = (r_z**2/r**2-r_zz*(-ZIMG*OMEGA+1.D0/r))*p   &
                     -  r_z*(-ZIMG*OMEGA+1.d0/r)*Gradp(3)
!
!  ...Gaussian beam
      case(5)
!
!     ...wavenumber
         rk = OMEGA
!     ...shift the origin
         xshift = x1 + 0.02d0; yshift = x2 + 0.02d0; zshift = x3 + 0.02d0
         theta_x = 35.d0; theta_y=-45.d0;
         rad_x = theta_x*pi/180.d0;  rad_y = theta_y*pi/180.d0;
!
         sinx=sin(rad_x); cosx = cos(rad_x);
         siny=sin(rad_y); cosy = cos(rad_y);
!
         Rx(1,1:3) = (/1.d0, 0.d0,  0.d0/)
         Rx(2,1:3) = (/0.d0, cosx, -sinx/)
         Rx(3,1:3) = (/0.d0, sinx,  cosx/)
!
         Ry(1,1:3) = (/cosy , 0.0d0, siny/)
         Ry(2,1:3) = (/0.0d0, 1.0d0, 0.0d0/)
         Ry(3,1:3) = (/-siny, 0.0d0, cosy/)

!     ...initialize rotation matrix R
         Rxy = 0.d0
!
         call DGEMM('N','N',3,3,3,1.0d0,Rx,3,Ry,3,0.0d0,Rxy,3)
!
!     ...new coordinates after rotation
         x = Rxy(1,1)*xshift+Rxy(1,2)*yshift+Rxy(1,3)*zshift
         y = Rxy(2,1)*xshift+Rxy(2,2)*yshift+Rxy(2,3)*zshift
         z = Rxy(3,1)*xshift+Rxy(3,2)*yshift+Rxy(3,3)*zshift
!
!     ...change of variable derivatives
         dxdxs = Rxy(1,1) ; dxdys = Rxy(1,2) ; dxdzs = Rxy(1,3)
         dydxs = Rxy(2,1) ; dydys = Rxy(2,2) ; dydzs = Rxy(2,3)
         dzdxs = Rxy(3,1) ; dzdys = Rxy(3,2) ; dzdzs = Rxy(3,3)
!
!     ...beam waist radius
         w0 = 0.05d0
!
         r = dsqrt(x**2+y**2)
!     ...first derivatives (only dependence on x and y)
         r_x = x/r ; r_y = y/r
!     ...second derivatives (only dependence on x and y)
         r_xx = (1.d0-r_x**2)/r ; r_xy = -r_y*r_x/r
         r_yx = r_xy;             r_yy =  (1.d0-r_y**2)/r
!
!     ...pressure
!
         p0 = 1.0d0*(2.0d0/pi/w0**2)**0.25;
         p  = p0*exp(-ZIMG*OMEGA*z)*exp(-r**2/w0**2)
!
!     ...derivatives with respect to x,y,z
         p_x = (-2.d0/w0**2)* x*p
         p_y = (-2.d0/w0**2)* y*p
         p_z = (-ZIMG*OMEGA)* p
!     ...second derivatives with respect to x,y,z
         p_xx = (-2.d0/w0**2)*(p+x*p_x)
         p_xy = (-2.d0/w0**2)* x*p_y
         p_xz = (-2.d0/w0**2)* x*p_z
         p_yx = p_xy
         p_yy = (-2.d0/w0**2)*(p+y*p_y)
         p_yz = (-2.d0/w0**2)* y*p_z
         p_zx = p_xz
         p_zy = p_yz
         p_zz = (-ZIMG*OMEGA)* p_z
!
!     ...derivatives with respect to xshift,yshift,zshift
         Gradp(1) = p_x*dxdxs + p_y*dydxs + p_z*dzdxs
         Gradp(2) = p_x*dxdys + p_y*dydys + p_z*dzdys
         Gradp(3) = p_x*dxdzs + p_y*dydzs + p_z*dzdzs
!     ...second derivatives  with respect to xshift,yshift,zshift
         Grad2p(1,1) = (p_xx*dxdxs + p_yx*dydxs + p_zx*dzdxs)*dxdxs    &
                     + (p_xy*dxdxs + p_yy*dydxs + p_zy*dzdxs)*dydxs    &
                     + (p_xz*dxdxs + p_yz*dydxs + p_zz*dzdxs)*dzdxs

         Grad2p(1,2) = (p_xx*dxdxs + p_yx*dydxs + p_zx*dzdxs)*dxdys    &
                     + (p_xy*dxdxs + p_yy*dydxs + p_zy*dzdxs)*dydys    &
                     + (p_xz*dxdxs + p_yz*dydxs + p_zz*dzdxs)*dzdys

         Grad2p(1,3) = (p_xx*dxdxs + p_yx*dydxs + p_zx*dzdxs)*dxdzs    &
                     + (p_xy*dxdxs + p_yy*dydxs + p_zy*dzdxs)*dydzs    &
                     + (p_xz*dxdxs + p_yz*dydxs + p_zz*dzdxs)*dzdzs

         Grad2p(2,1) = (p_xx*dxdys + p_yx*dydys + p_zx*dzdys)*dxdxs    &
                     + (p_xy*dxdys + p_yy*dydys + p_zy*dzdys)*dydxs    &
                     + (p_xz*dxdys + p_yz*dydys + p_zz*dzdys)*dzdxs

         Grad2p(2,2) = (p_xx*dxdys + p_yx*dydys + p_zx*dzdys)*dxdys    &
                     + (p_xy*dxdys + p_yy*dydys + p_zy*dzdys)*dydys    &
                     + (p_xz*dxdys + p_yz*dydys + p_zz*dzdys)*dzdys

         Grad2p(2,3) = (p_xx*dxdys + p_yx*dydys + p_zx*dzdys)*dxdzs    &
                     + (p_xy*dxdys + p_yy*dydys + p_zy*dzdys)*dydzs    &
                     + (p_xz*dxdys + p_yz*dydys + p_zz*dzdys)*dzdzs

         Grad2p(3,1) = (p_xx*dxdzs + p_yx*dydzs + p_zx*dzdzs)*dxdxs    &
                     + (p_xy*dxdzs + p_yy*dydzs + p_zy*dzdzs)*dydxs    &
                     + (p_xz*dxdzs + p_yz*dydzs + p_zz*dzdzs)*dzdxs

         Grad2p(3,2) = (p_xx*dxdzs + p_yx*dydzs + p_zx*dzdzs)*dxdys    &
                     + (p_xy*dxdzs + p_yy*dydzs + p_zy*dzdzs)*dydys    &
                     + (p_xz*dxdzs + p_yz*dydzs + p_zz*dzdzs)*dzdys

         Grad2p(3,3) = (p_xx*dxdzs + p_yx*dydzs + p_zx*dzdzs)*dxdzs    &
                     + (p_xy*dxdzs + p_yy*dydzs + p_zy*dzdzs)*dydzs    &
                     + (p_xz*dxdzs + p_yz*dydzs + p_zz*dzdzs)*dzdzs
!
   end select
!
!
   end subroutine mfd_solutions


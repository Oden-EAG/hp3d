!-----------------------------------------------------------------------------------
!> Purpose : compute all relevant quantities of the exact solutions for the
!            linear elasticty problems
!!
!! @param[in]  X        - a point in physical space
!! @param[out] p        - value of the solution displacement
!! @param[out] Gradp    - corresponding first derivatives
!! @param[out] Grad2p   - corresponding second derivatives
!-----------------------------------------------------------------------------------
!
!----------------------------------------------------------------------
!                                                                     
!     routine name      - acoustics_solution
!                                                                     
!---------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - July 17
!                                                                     
!     purpose:          - compute all relevant quantities of the exact
!                         solutions for the linear acoustics problems
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
   subroutine acoustics_solution(Xp, p,Gradp,Grad2p)
!   
   use data_structure3D
   use common_prob_data, only : PI, OMEGA,IEXACT_PROB
   use parameters, only : ZERO, ZONE, ZIMG
!
   implicit none  
!
   real*8,     dimension(3),   intent(in)  :: Xp
   real*8 ::   x1,x2,x3, xshift, yshift,zshift,r,alpha,r_x, r_y, r_z
   real*8 ::   r_xx, r_xy, r_xz, r_yx, r_yy, r_yz, r_zx, r_zy, r_zz 
   complex*16,                 intent(out) :: p
   complex*16, dimension(3),   intent(out) :: Gradp  ! 1st derivative - gradient
   complex*16, dimension(3,3), intent(out) :: Grad2p ! 2nd derivative - Hessian
   complex*16  cn, zf_x,zf_y,zf_z,dzf_x,dzf_y, dzf_z, ddzf_x, ddzf_y, ddzf_z
   real*8 :: nn, x10, x20, x30
!..required for the Gaussian beam   
!..parameters
   real*8 :: rl, w0, p0, z_R,rk, theta_x, theta_y,rad_x,rad_y, pi_mod
!..rotation variables and matrices   
   real*8 :: sinx,siny, cosx, cosy
   real*8 :: x, y, z
   real*8 :: dxdxs, dxdys, dxdzs, dydxs, dydys, dydzs, dzdxs, dzdys, dzdzs 
   real*8, dimension(3, 3) :: Rx, Ry, Rxy   
!..functions and their derivatives
   real*8     :: w,w_z,w_zz
   real*8     :: f,f_z,f_zz
   real*8     :: A,A_x,A_xx,A_xy,A_xz,A_y,A_yx,A_yy,A_yz,A_z,A_zx,A_zy,A_zz
   real*8     :: g,g_x,g_xx,g_xy,g_xz,g_y,g_yx,g_yy,g_yz,g_z,g_zx,g_zy,g_zz
   real*8     :: phi, phi_z, phi_zz
   real*8     :: D, D_z, D_zz
   real*8     :: C,C_x,C_xx,C_xy,C_xz,C_y,C_yx,C_yy,C_yz,C_z,C_zx,C_zy,C_zz
   real*8     :: B,B_x,B_xx,B_xy,B_xz,B_y,B_yx,B_yy,B_yz,B_z,B_zx,B_zy,B_zz
   complex*16 :: h,h_x,h_xx,h_xy,h_xz,h_y,h_yx,h_yy,h_yz,h_z,h_zx,h_zy,h_zz 
   complex*16 ::   p_x,p_xx,p_xy,p_xz,p_y,p_yx,p_yy,p_yz,p_z,p_zx,p_zy,p_zz 
! 
!---------------------------------------------------------------------------------------
!
   x1 = Xp(1); x2 = Xp(2); x3 = Xp(3)


   select case (IEXACT_PROB)
!
!..polynomial in x dimension
   case(0)  

      p = x1**2 
      Gradp(1) = 2.d0*x1
      Gradp(2) = ZERO
      Gradp(3) = ZERO
      Grad2p(1,1) = 2.d0
      Grad2p(1,2) = ZERO
      Grad2p(1,3) = ZERO
      Grad2p(2,1) = ZERO
      Grad2p(2,2) = ZERO
      Grad2p(2,3) = ZERO
      Grad2p(3,1) = ZERO
      Grad2p(3,2) = ZERO
      Grad2p(3,3) = ZERO

!
!..sin solution vanishing on the boundary
   case(1)   
!
      pi_mod = OMEGA
      p =dsin(x1*pi_mod)*dsin(x2*pi_mod)*dsin(x3*pi_mod)*(1.0D0,1.0D0)
!   
!  ...1st order derivatives
      Gradp(1) = pi_mod*dcos(x1*pi_mod)*dsin(x2*pi_mod)*dsin(x3*pi_mod)*(1.0D0,1.0D0)
      Gradp(2) = pi_mod*dcos(x2*pi_mod)*dsin(x1*pi_mod)*dsin(x3*pi_mod)*(1.0D0,1.0D0)
      Gradp(3) = pi_mod*dcos(x3*pi_mod)*dsin(x1*pi_mod)*dsin(x2*pi_mod)*(1.0D0,1.0D0)
! 
!  ...2nd derivative (3,3) matrix - Hessian
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
!..polynomial solution vanishing on the boundary
   case(2)
!
      cn=1.D2*ZONE;     !  cn=1.D2*(ZONE+ZIMG);
!  ...displacement
      p = -cn*x1*x2*x3*(x1-1.0D0)*(x2-1.0D0)*(x3-1.0D0)
!  
!  ...1st order derivatives
      Gradp(1) =-cn*x2*x3*(x1*2.0D0-1.0D0)*(x2-1.0D0)*(x3-1.0D0)
      Gradp(2) =-cn*x1*x3*(x2*2.0D0-1.0D0)*(x1-1.0D0)*(x3-1.0D0)
      Gradp(3) =-cn*x1*x2*(x3*2.0D0-1.0D0)*(x1-1.0D0)*(x2-1.0D0)
!  ...2nd derivative (3,3) matrix - hessian
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

!..plane wave
   case(3)
      zf_x = exp(-ZIMG*OMEGA*x1)
      zf_y = exp(-ZIMG*OMEGA*x2)
      zf_z = exp(-ZIMG*OMEGA*x3)
!      
!  ..first order derivatives      
      dzf_x = -ZIMG*OMEGA*zf_x
      dzf_y = -ZIMG*OMEGA*zf_y
      dzf_z = -ZIMG*OMEGA*zf_z
!
! ...second order derivatives      
      ddzf_x = -OMEGA**2*zf_x
      ddzf_y = -OMEGA**2*zf_y
      ddzf_z = -OMEGA**2*zf_z
!
      p = zf_x*zf_y*zf_z

! ...1st order derivatives
      Gradp(1) = dzf_x*zf_y*zf_z
      Gradp(2) = zf_x*dzf_y*zf_z
      Gradp(3) = zf_x*zf_y*dzf_z
! 
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

!..point source 
   case(4)
!
!  ...shift origin
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
!..Gaussian beam
   case(5)
!
!  ...wavenumber
      rk = OMEGA
!  ...shift the origin   
      xshift = x1 + 0.05d0; yshift = x2 + 0.05d0; zshift = x3 + 0.05d0
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
! 
!  ...initialize rotation matrix R
      Rxy = 0.d0
!   
      call DGEMM('N','N',3,3,3,1.0d0,Rx,3,Ry,3,0.0d0,Rxy,3)
!   
!  ...new coordinates after rotation
      x = Rxy(1,1)*xshift+Rxy(1,2)*yshift+Rxy(1,3)*zshift
      y = Rxy(2,1)*xshift+Rxy(2,2)*yshift+Rxy(2,3)*zshift
      z = Rxy(3,1)*xshift+Rxy(3,2)*yshift+Rxy(3,3)*zshift

!  ...change of variable derivatives
      dxdxs = Rxy(1,1) ; dxdys = Rxy(1,2) ; dxdzs = Rxy(1,3)
      dydxs = Rxy(2,1) ; dydys = Rxy(2,2) ; dydzs = Rxy(2,3)
      dzdxs = Rxy(3,1) ; dzdys = Rxy(3,2) ; dzdzs = Rxy(3,3)

!  ...wavelength
      rl = 2.d0*pi/rk
! 
!  ...beam waist radius
      w0 = 0.09d0
!  ...z_R
      z_R = pi*w0**2/rl   
!
!  ...functions and their derivatives
!  ...function w
      w  = w0*dsqrt(1.d0 + (z/z_R)**2)
!  ...first derivative (only dependence on z)
      w_z = (w0**2/z_R**2)  *( z / w )
!  ...second derivative (only dependence on z)
      w_zz = (w0**2/z_R**2) * (w - z*w_z) / w**2  
!
!  ...function f
      f = w0/w
!  ...first derivative (only dependence on z)
      f_z = - w0 * (w_z / w**2)
!  ...second derivative (only dependence on z)
      f_zz = - w0 * (w_zz*w - 2.d0*w_z**2)/w**3 
!
!
!  ...function r (radial distance)
      r = dsqrt(x**2+y**2)
!  ...first derivatives (only dependence on x and y)
      r_x = x/r ; r_y = y/r
!  ...second derivatives (only dependence on x and y)
      r_xx = (1.d0-r_x**2)/r ; r_xy = -r_y*r_x/r
      r_yx = r_xy;             r_yy =  (1.d0-r_y**2)/r
!
!
!  ...function A
      A = r**2/w**2
!  ...first derivatives 
      A_x =  2.d0*r*r_x/w**2  
      A_y =  2.d0*r*r_y/w**2  
      A_z = -2.d0*r**2*w_z/w**3  
!  ...second derivatives 
      A_xx =  2.d0 * (r_x**2 + r* r_xx)/w**2
      A_xy =  2.d0 * (r_y*r_x + r* r_xy)/w**2
      A_xz = -4.d0 * r*r_x*w_z/w**3
      A_yx =  A_xy
      A_yy =  2.d0 * (r_y**2 + r* r_yy)/w**2
      A_yz = -4.d0 * r*r_y*w_z/w**3
      A_zx = A_xz
      A_zy = A_yz
      A_zz = -2.d0*r**2 *(w_zz*w-3.d0*w_z**2)/w**4
!
!
!  ...function g
      g = exp(-A)
!  ...first derivatives   
      g_x = -A_x*g ;  g_y = -A_y*g ;  g_z = -A_z*g ;
!  ...second derivatives
      g_xx =-A_xx*g-A_x*g_x ; g_xy =-A_xy*g-A_x*g_y ; g_xz =-A_xz*g-A_x*g_z
      g_yx = g_xy           ; g_yy =-A_yy*g-A_y*g_y ; g_yz =-A_yz*g-A_y*g_z
      g_zx = g_xz           ; g_zy = g_yz           ; g_zz =-A_zz*g-A_z*g_z
!
!
!  ...function D
      D = z + z_R**2/z
!  ...first derivative (only dependence on z)
      D_z = 1.d0 - z_R**2/z**2
!  ...second derivative (only dependence on z)
      D_zz = 2.d0*z_R**2/z**3
!
!   
!  ...function phi
      phi = atan(z/z_R)
!  ...first derivative (only dependence on z)
      phi_z = z_R/(z_R**2+z**2)
!  ...second derivative (only dependence on z)
      phi_zz = -2.d0*z_R*z/(z_R**2+z**2)**2
!
!
!  ...function C
      C = r**2/D
!  ...first derivatives
      C_x = 2.0*r*r_x/D ; C_y = 2.d0*r*r_y/D; C_z = -r**2*D_z/D**2
!  ...second derivatives
      C_xx = 2.d0 *(r_x**2 +  r*r_xx)/D
      C_xy = 2.d0 *(r_x*r_y + r*r_xy)/D
      C_xz = -2.d0 *r* r_x * D_z / D**2
      C_yx = C_xy
      C_yy = 2.d0 *(r_y**2 + r*r_yy)/D
      C_yz = -2.d0 *r* r_y * D_z / D**2
      C_zx = C_xz
      C_zy = C_yz
      C_zz = -r**2 *(D_zz*D - 2.d0*D_z**2)/D**3
!
!  ...function B
      B = rk*z+rk/2.d0*C-phi
!  ...first derivatives
      B_x = rk/2.d0*C_x
      B_y = rk/2.d0*C_y
      B_z = rk + rk/2.d0*C_z - phi_z

!  ...second derivatives
      B_xx = rk/2.d0*C_xx ; B_xy = rk/2.d0*C_xy ; B_xz = rk/2.d0*C_xz
      B_yx = B_xy         ; B_yy = rk/2.d0*C_yy ; B_yz = rk/2.d0*C_yz
      B_zx = B_xz         ; B_zy = B_yz         ; B_zz =  rk/2.d0*C_zz - phi_zz
   
   
!  ...function h
      h = exp(-ZIMG*B)
!  ...first derivatives
      h_x = -ZIMG*B_x*h ; h_y = -ZIMG*B_y*h ; h_z = -ZIMG*B_z*h
!  ...second derivatives
      h_xx = -ZIMG*(B_xx*h + B_x*h_x)
      h_xy = -ZIMG*(B_xy*h + B_x*h_y)
      h_xz = -ZIMG*(B_xz*h + B_x*h_z)
      h_yx = h_xy
      h_yy = -ZIMG*(B_yy*h + B_y*h_y)
      h_yz = -ZIMG*(B_yz*h + B_y*h_z)
      h_zx = h_xz
      h_zy = h_yz
      h_zz = -ZIMG*(B_zz*h + B_z*h_z)
!
!
!  ...pressure
      p0 = 1.0d0
      p = p0*f*g*h   
!  ...derivatives with respect to x,y,z
      p_x = p0*f*(g_x*h+g*h_x)
      p_y = p0*f*(g_y*h+g*h_y)
      p_z = p0*(f_z*g*h + f*g_z*h + f*g*h_z)
!  ...second derivatives with respect to x,y,z
      p_xx = p0*f*(g_xx*h+2.d0*g_x*h_x+g*h_xx)
      p_xy = p0*f*(g_xy*h+g_x*h_y + g_y*h_x +g*h_xy)
      p_xz = p0*f_z*(g_x*h+g*h_x)+p0*f*(g_xz*h+g_x*h_z+g_z*h_x+g*h_xz)
      p_yx = p_xy
      p_yy = p0*f*(g_yy*h+2.d0*g_y*h_y+g*h_yy)
      p_yz = p0*f_z*(g_y*h+g*h_y)+p0*f*(g_yz*h+g_y*h_z+g_z*h_y+g*h_yz)
      p_zx = p_xz
      p_zy = p_yz
      p_zz = p0*(  f_zz*g*h  + f_z*g_z*h  + f_z*g*h_z  &    
                 + f_z*g_z*h + f*g_zz*h   + f*g_z*h_z  & 
                 + f_z*g*h_z + f*g_z*h_z  + f*g*h_zz  )

!  ...derivatives with respect to xshift,yshift,zshift
      Gradp(1) = p_x*dxdxs + p_y*dydxs + p_z*dzdxs
      Gradp(2) = p_x*dxdys + p_y*dydys + p_z*dzdys
      Gradp(3) = p_x*dxdzs + p_y*dydzs + p_z*dzdzs
!  ...second derivatives  with respect to xshift,yshift,zshift
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
!..Gaussian beam simplified
   case(6)
!
!  ...wavenumber
      rk = OMEGA
!  ...shift the origin   
      xshift = x1 + 0.05d0; yshift = x2 + 0.05d0; zshift = x3 + 0.05d0
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

!  ...initialize rotation matrix R
      Rxy = 0.d0
!   
      call DGEMM('N','N',3,3,3,1.0d0,Rx,3,Ry,3,0.0d0,Rxy,3)
!
!  ...new coordinates after rotation
      x = Rxy(1,1)*xshift+Rxy(1,2)*yshift+Rxy(1,3)*zshift
      y = Rxy(2,1)*xshift+Rxy(2,2)*yshift+Rxy(2,3)*zshift
      z = Rxy(3,1)*xshift+Rxy(3,2)*yshift+Rxy(3,3)*zshift

!  ...change of variable derivatives
      dxdxs = Rxy(1,1) ; dxdys = Rxy(1,2) ; dxdzs = Rxy(1,3)
      dydxs = Rxy(2,1) ; dydys = Rxy(2,2) ; dydzs = Rxy(2,3)
      dzdxs = Rxy(3,1) ; dzdys = Rxy(3,2) ; dzdzs = Rxy(3,3)
! 
!  ...beam waist radius
      w0 = 0.05d0
!
      r = dsqrt(x**2+y**2)
!  ...first derivatives (only dependence on x and y)
      r_x = x/r ; r_y = y/r
!  ...second derivatives (only dependence on x and y)
      r_xx = (1.d0-r_x**2)/r ; r_xy = -r_y*r_x/r
      r_yx = r_xy;             r_yy =  (1.d0-r_y**2)/r
!
!  ...pressure
! 
      p0 = 1.0d0*(2.0d0/pi/w0**2)**0.25;
      p  = p0*exp(-ZIMG*OMEGA*z)*exp(-r**2/w0**2)
!   
!  ...derivatives with respect to x,y,z
      p_x = (-2.d0/w0**2)* x*p
      p_y = (-2.d0/w0**2)* y*p
      p_z = (-ZIMG*OMEGA)* p
!  ...second derivatives with respect to x,y,z
      p_xx = (-2.d0/w0**2)*(p+x*p_x)
      p_xy = (-2.d0/w0**2)* x*p_y
      p_xz = (-2.d0/w0**2)* x*p_z
      p_yx = p_xy
      p_yy = (-2.d0/w0**2)*(p+y*p_y)
      p_yz = (-2.d0/w0**2)* y*p_z
      p_zx = p_xz
      p_zy = p_yz
      p_zz = (-ZIMG*OMEGA)* p_z

!  ...derivatives with respect to xshift,yshift,zshift
      Gradp(1) = p_x*dxdxs + p_y*dydxs + p_z*dzdxs
      Gradp(2) = p_x*dxdys + p_y*dydys + p_z*dzdys
      Gradp(3) = p_x*dxdzs + p_y*dydzs + p_z*dzdzs
!  ...second derivatives  with respect to xshift,yshift,zshift
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
   end subroutine acoustics_solution





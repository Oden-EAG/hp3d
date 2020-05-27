      subroutine hankel_source_field_E(
     .  X, ZRad, Theta, Phi,
     .  Z_int, ZD_int)
#include "syscom.blk"
c
c     ... Internal field on spherical coord
      complex*16, dimension(3)            :: Z_int
      complex*16, dimension(3,3)          :: ZD_int
c
      real*8, parameter :: eps = 1.0e-12
c
      z_i = cmplx(0.d0, 1.d0)
      z_0 = cmplx(0.d0, 0.d0)
c
      if (abs(ZRad).lt.eps) then
        ZRad = cmplx(eps, 0.d0)
      endif
c
c     ...used in Z_int which is influenced by Zref
      z_rx   = ZRad*X
c
c     ...spherical bessel function and its derivatvies
c     ...first kind
ccc      zh_rx  = -z_i*exp(z_i*z_rx)/z_rx
ccc      zhd_rx = exp(zh_rx*z_i)/zh_rx + (exp(zh_rx*z_i)*z_i)/zh_rx**2
c     ...second kind... this is outgoing wave
      zh_rx  =  z_i*exp(-z_i*z_rx)/z_rx
      zhd_rx = (z_1/exp(z_rx*z_i))/z_rx -
     .  ((z_1/exp(z_rx*z_i))*z_i)/z_rx**2
c
      Z_int(1:3) = z_0; ZD_int(1:3,1:3) = z_0
c
      Z_int(1)    = zh_rx
      ZD_int(1,1) = zhd_rx
c
c
      end subroutine
c
c

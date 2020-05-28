      subroutine mie_near_field_cart( Ipec,
     .  N, Zref, X,
     .  Xp,
     .  ZE_inc, ZE_int, ZE_sca,
     .  ZC_inc, ZC_int, ZC_sca)
#include "syscom.blk"
      real*8,     dimension(3) :: Xp
      complex*16, dimension(3) :: ZE_inc, ZE_int, ZE_sca
      complex*16, dimension(3) :: ZC_inc, ZC_int, ZC_sca
c
      complex*16, dimension(3)   :: Z_s_inc,  Z_s_int,  Z_s_sca
      complex*16, dimension(3,3) :: ZD_s_inc, ZD_s_int, ZD_s_sca
c
      complex*16, dimension(3)   :: Z_c_inc,  Z_c_int,  Z_c_sca
      complex*16, dimension(3,3) :: ZD_c_inc, ZD_c_int, ZD_c_sca
c
      complex*16, parameter :: z_i = cmplx(0.d0,1.d0)
      complex*16, parameter :: z_0 = cmplx(0.d0,0.d0)
      real*8 :: sp(3)
c
      call coord_cart2sphere(Xp, sp)
      call mie_near_field_E(N, Zref,
     .  X, sp(1), sp(2), sp(3),
     .  Z_s_inc,  Z_s_int,  Z_s_sca,
     .  ZD_s_inc, ZD_s_int, ZD_s_sca)
c
      if (Ipec.eq.1) then
        call mie_pec_near_field_E(N,
     .    X, sp(1), sp(2), sp(3),
     .    Z_s_sca, ZD_s_sca)
        Z_s_int = z_0; ZD_s_int = z_0;
      end if
c
      call disp_sphere2cart(sp, Z_s_inc, ZD_s_inc, Z_c_inc, ZD_c_inc)
      call disp_sphere2cart(sp, Z_s_int, ZD_s_int, Z_c_int, ZD_c_int)
      call disp_sphere2cart(sp, Z_s_sca, ZD_s_sca, Z_c_sca, ZD_c_sca)
c
      ZE_inc = Z_c_inc
      ZE_int = Z_c_int
      ZE_sca = Z_c_sca
c
      do k2=1,3
        do k1=1,3
          if (isnan(real(ZD_c_inc(k1,k2)))) then
            ZD_c_inc(k1,k2) = z_0
          end if
          if (isnan(real(ZD_c_int(k1,k2)))) then
            ZD_c_int(k1,k2) = z_0
          end if
          if (isnan(real(ZD_c_sca(k1,k2)))) then
            ZD_c_sca(k1,k2) = z_0
          end if
        end do
      end do
c
      call compute_curl(ZD_c_inc, ZC_inc)
      call compute_curl(ZD_c_int, ZC_int)
      call compute_curl(ZD_c_sca, ZC_sca)
c     ...if u use mie_sph_back following is not necessary..
c     ...I probably did some sign mistake in derivative, cannot find it out yet.
c     ...Or we can see whether mie is wrong or my solution is wrong later
      ZC_inc = -ZC_inc
      ZC_int = -ZC_int
      ZC_sca = -ZC_sca
c
      end subroutine

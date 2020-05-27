c> Purpose : find out the truncation n for mie series.
c     @param[in] X - k*a (a-r of the sphere, k-wave number)
c------------------------------------------------------
      integer function mie_n_max(X)
#include "syscom.blk"
c
      mie_n_max = nint( X+4.d0*X**(1.d0/3.d0)+2.d0 )
      end function mie_n_max
c------------------------------------------------------
c> Purpose : Computes a matrix of Mie coefficients, a,b,c,d
c     from order n=1 to N
c     @param[in]  N    - nmax to compute coefficients
c     @param[in]  Zref - complex refractive index
c     @param[in]  X    - x = k0*a, where k0 = wave number
c     @param[out] za,zb,zc,zd - mie coefficients
c------------------------------------------------------
      subroutine mie_coeff(
     .  Anzat,
     .  N, Zref,X,
     .  Za,Zb,Zc,Zd)
#include "syscom.blk"
c
      integer, parameter :: idelta = 2
c
      dimension Za(N), Zb(N), Zc(N), Zd(N)
c
      character(len=3) :: Anzat
c     local variables
      dimension
     .  zs1   (0:N+idelta),    zd1(0:N+idelta),
     .  zs2   (0:N+idelta),    zd2(0:N+idelta),
     .  zs1_zr(0:N+idelta), zd1_zr(0:N+idelta),
     .  zs2_zr(0:N+idelta), zd2_zr(0:N+idelta),
     .  zs1_h (0:N+idelta),  zd1_h(0:N+idelta)
c
      dimension ztemp(2,2)
c
      iprint = 0
c
      z_i = cmplx(0.d0, 1.d0)
      z_0 = cmplx(0.d0, 0.d0)
c
      z = cmplx(X, 0.d0)
      zref2 = Zref*Zref
      zr = Zref*z
c
c     spherical bessel functions
      call csphjy(N,  z, N, zs1,    zd1,    zs2,    zd2)
      call csphjy(N, zr, N, zs1_zr, zd1_zr, zs2_zr, zd2_zr)
c
c     hankel function 1st kind
      select case(Anzat)
      case('neg');
        zs1_h(0:N+idelta) = zs1(0:N+idelta) + z_i*zs2(0:N+idelta)
        zd1_h(0:N+idelta) = zd1(0:N+idelta) + z_i*zd2(0:N+idelta)
      case('pos');
        zs1_h(0:N+idelta) = zs1(0:N+idelta) - z_i*zs2(0:N+idelta)
        zd1_h(0:N+idelta) = zd1(0:N+idelta) - z_i*zd2(0:N+idelta)
      end select
c
c
      do i=1, N
c
        ztemp(1,1) = zs1_zr(i)*(z* zs1   (i-1) - i*zs1   (i))
        ztemp(1,2) = zs1(i)   *(zr*zs1_zr(i-1) - i*zs1_zr(i))
        ztemp(2,1) = zs1_zr(i)*(z* zs1_h (i-1) - i*zs1_h (i))
        ztemp(2,2) = zs1_h(i) *(zr*zs1_zr(i-1) - i*zs1_zr(i))
c
        Za(i) = ( zref2*ztemp(1,1) - ztemp(1,2) )/
     .          ( zref2*ztemp(2,1) - ztemp(2,2) )
        Zb(i) = (       ztemp(1,1) - ztemp(1,2) )/
     .          (       ztemp(2,1) - ztemp(2,2) )
c
        ztemp(1,1) = zs1  (i)*(z* zs1_h (i-1) - i*zs1_h (i))
        ztemp(1,2) = zs1_h(i)*(z* zs1   (i-1) - i*zs1   (i))
c
        Zc(i) = ( ztemp(1,1) - ztemp(1,2) )/
     .          ( ztemp(2,1) - ztemp(2,2) )
        Zd(i) = ( Zref*ztemp(1,1) - Zref*ztemp(1,2) )/
     .          ( zref2*ztemp(2,1) - ztemp(2,2) )
c
      enddo
c
      if (iprint.eq.1) then
        do i=1, N
          write(*,6000) i, Za(i), Zb(i), Zc(i), Zd(i)
 6000     format(i3, ' an,bn,cn,dn = ', 4(2e12.9,2x))
        enddo
      endif
c
c
      end subroutine mie_coeff
c------------------------------------------------------
c> Computation of Mie Efficiencies for given
cc complex refractive-index ratio m=m'+im"
cc and size parameter x=k0*a, where k0= wave number in ambient
cc medium, a=sphere radius, using complex Mie Coefficients
cc an and bn for n=1 to nmax,
cc s. Bohren and Huffman (1983) BEWI:TDD122, p. 103,119-122,477.
c------------------------------------------------------
      subroutine mie( Anzat,
     .  N, Zref, X,
     .  Qext, Qsca, Qabs, Qb)
#include "syscom.blk"
      character(len=3) :: Anzat
      integer, parameter :: idelta = 3
      complex*16 :: zqb
      complex*16, dimension(N+idelta) :: za, zb, zc, zd
c
      z_i = cmplx(0.d0, 1.d0)
      z_0 = cmplx(0.d0, 0.d0)
c
      Qext = 0.d0; Qsca = 0.d0; Qabs = 0.d0; Qb = 0.d0;
      zqb  = z_0
      if (X.eq.0.d0) return
c
      call mie_coeff(Anzat, N,Zref,X,  za,zb,zc,zd)
c
      do i=1, nmax
        Qsca = Qsca + (2.d0*i+1.d0)*
     .                ( abs(za(i))**2 + abs(zb(i))**2 )
        Qext = Qext + (2.d0*i+1.d0)*
     .                ( real(za(i)+zb(i)) )
        zqb  = zqb  + (2.d0*i+1.d0)*((-1)**i)*
     .                ( za(i) - zb(i) )
      enddo

      Qsca = Qsca*2.d0/(x**2)
      Qext = Qext*2.d0/(x**2)
      Qabs = Qext - Qsca
      Qb   = (zqb*conjg(zqb))/(x**2)
c
      end subroutine
c------------------------------------------------------
c     Purpose : compute angular functions for Mie
c     N     (in) - n_max for Mie series
c     U     (in) - cos(theta)
c     P, DP (out)- phi function and its derivative
c     T, DT (out)- tau function and its derivative
c------------------------------------------------------
      subroutine mie_phitau(N, U, P, T, DP, DT)
#include "syscom.blk"
      integer, parameter :: idelta = 3
      real*8, dimension(N+idelta) :: P, T, DP, DT
      integer :: i
      real*8 :: tmp1, tmp2
c
      P(1) = 1.d0;        T(1) = U;
      P(2) = 3.d0*U;      T(2) = 6.d0*(U**2) - 3.d0

      DP(1)= 0.d0;        DT(1)= 1
      DP(2)= 3.d0;        DT(2)= 12.d0*U

      do i=3, N
        tmp1 = (2.d0*i-1.d0)/(i-1.d0)*P(i-1)*U;
        tmp2 = i/(i-1.d0)*P(i-2);
        P(i) = tmp1 - tmp2;
c
        tmp1 = i*U*P(i);
        tmp2 = (i+1.d0)*P(i-1);
        T(i) = tmp1 - tmp2;
c
        tmp1 = (2.d0*i-1.d0)/(i-1.d0)*(DP(i-1)*U+P(i-1))
        tmp2 = i*DP(i-2)/(i-1.d0);
        DP(i)= tmp1 - tmp2;
c
        tmp1 = i*(U*DP(i)+P(i));
        tmp2 = (i+1.d0)*DP(i-1);
        DT(i)= tmp1 - tmp2;
c
      enddo
c
      end subroutine
c------------------------------------------------------
c
c------------------------------------------------------
      subroutine mie_s12( Anzat,
     .  N, Zref, X, U,
     .  ZS1, ZS2)
#include "syscom.blk"
c
      character(len=3) :: Anzat
c
      integer, parameter :: idelta = 3
      real*8,     dimension(N+idelta) :: p, t, dp, dt
      complex*16, dimension(N+idelta) :: za, zb, zc, zd
c
      z_i = cmplx(0.d0, 1.d0)
      z_0 = cmplx(0.d0, 0.d0)
c
      call mie_coeff(Anzat, N,Zref,X,  za,zb,zc,zd)
      call mie_phitau(N, U, p, t, dp, dt);
c
      ZS1 = z_0
      ZS2 = z_0
      do i=1, N
        d = (2.d0*i+1.d0)/i/(i+1.d0)
ccc        write(*,*) d, p(i), t(i)
        ZS1 = ZS1 + d*(za(i)*p(i) + zb(i)*t(i))
        ZS2 = ZS2 + d*(za(i)*t(i) + zb(i)*p(i))
      enddo
c
      end subroutine
c
c------------------------------------------------------
c
c------------------------------------------------------
      subroutine mie_far_field(Anzat,
     .  N, Zref, X,
     .  Rk, Rad, Theta, Phi,
     .  ZST, ZSP)
#include "syscom.blk"
      character(len=3) :: Anzat
c
      z_i = cmplx(0.d0, 1.d0)
      z_0 = cmplx(0.d0, 0.d0)
c
      call mie_s12(Anzat, N, Zref, X, dcos(Theta), zs1, zs2)
c
      ZST =
     .  exp(z_i*Rk*Rad)
     .  /(-z_i*Rk*Rad)
     .  *dcos(Phi)*zs2
c
      ZSP =
     .  exp(z_i*Rk*Rad)
     .  /(z_i*Rk*Rad)
     .  *dsin(Phi)*zs1
c
      end subroutine
c
c
c------------------------------------------------------
c
c------------------------------------------------------
      subroutine Mie_mn(
     .  N, Theta, Phi,
     .  Pn, Tn, PDn, TDn,
     .  Z, ZD,
     .  Zj_a, Zj_b, Zjd_a, Zjd_b,
     .  Zm, Zn, ZDm, ZDn)
#include "syscom.blk"
      integer                    :: N
      real*8                     :: Theta, Phi, Pn, Tn
      complex*16                 :: Z, ZD, Zj_a, Zj_b, Zjd_a, Zjd_b
      complex*16, dimension(3)   :: Zm, Zn
      complex*16, dimension(3,3) :: ZDm, ZDn
c

      z_i = cmplx(0.d0, 1.d0)
      z_0 = cmplx(0.d0, 0.d0)
c
      coTh = dcos(Theta)
      siTh = dsin(Theta)
      coPh = dcos(Phi)
      siPh = dsin(Phi)
c
      Zm(1:3)      = z_0
      Zn(1:3)      = z_0
      ZDm(1:3,1:3) = z_0
      ZDm(1:3,1:3) = z_0
c
      zjd   = ( Z*Zj_b - N*Zj_a )
      zjd2  = ( Zj_b + Z*Zjd_b - N*Zjd_a )
c
      Zm(2)    =  coPh*Pn*Zj_a
      Zm(3)    = -siPh*Tn*Zj_a
c
      Zn(1)    =  N*(N+1)*coPh*siTh*Pn*Zj_a/Z
      Zn(2)    =  coPh*Tn*zjd/Z
      Zn(3)    = -siPh*Pn*zjd/Z
c
      ZDm(2,1) =  coPh*Pn*Zjd_a*ZD
      ZDm(2,2) =  coPh*PDn*(-siTh)*Zj_a
      ZDm(2,3) = -siPh*Pn*Zj_a
c
      ZDm(3,1) = -siPh*Tn*zjd_a*ZD
      ZDm(3,2) = -siPh*TDn*(-siTh)*Zj_a
      ZDm(3,3) = -coPh*Tn*Zj_a
c
c
      ZDn(1,1) = N*(N+1)*coPh*siTh*Pn*(Zjd_a*Z-Zj_a)*ZD/(Z**2)
      ZDn(1,2) = N*(N+1)*coPh*(coTh*Pn-siTh**2*PDn)*Zj_a/Z
      ZDn(1,3) = N*(N+1)*(-siPh)*siTh*Pn*Zj_a/Z
c
      ZDn(2,1) = coPh*Tn*(zjd2*Z-zjd)*ZD/(Z**2)
      ZDn(2,2) = coPh*TDn*(-siTh)*zjd/Z
      ZDn(2,3) = -siPh*Tn*zjd/Z
c
      ZDn(3,1) = -siPh*Pn*(zjd2*Z-zjd)*ZD/(Z**2)
      ZDn(3,2) = -siPh*PDn*(-siTh)*zjd/Z
      ZDn(3,3) = -coPh*Pn*zjd/Z
c
      end subroutine
c
c------------------------------------------------------
c
c------------------------------------------------------
      subroutine mie_near_field_E( Anzat,
     .  N, Zref, X,
     .  Rad, Theta, Phi,
     .  Z_inc,  Z_int,  Z_sca,
     .  ZD_inc, ZD_int, ZD_sca)
#include "syscom.blk"
c
      character(len=3) :: Anzat
c     ... Internal field on spherical coord
      complex*16, dimension(3)            :: Z_inc,  Z_int,  Z_sca
      complex*16, dimension(3,3)          :: ZD_inc, ZD_int, ZD_sca
c
      integer, parameter                  :: idelta = 3
c
      real*8,     dimension(N+idelta)     :: p, t, dp, dt
      complex*16, dimension(N+idelta)     :: zca, zcb, zcc, zcd
c
      complex*16, dimension(0:N+idelta-1) ::
     .  zjd_rzx, zj_rzx, zjd_rx, zj_rx,
     .  zyd_rzx, zy_rzx, zyd_rx, zy_rx,
     .  zhd_rx,  zh_rx
c
      complex*16, dimension(3)            :: zm, zn
      complex*16, dimension(3,3)          :: zdm, zdn
      real*8, parameter :: epsilon = 1.0e-15, pi=4.*ATAN(1.D0)
c
      z_i = cmplx(0.d0, 1.d0)
      z_0 = cmplx(0.d0, 0.d0)
c
      if (abs(Rad).lt.eps) then
        Rad = eps
      endif

      select case(Anzat)
      case('pos');
        Theta = pi - Theta
      end select
c
c     ...used in Z_int which is influenced by Zref
      z_rzx   = Rad*Zref*X
      z_rx    = Rad*X
c
c     ...used in Z_inc, Z_sca in free space
      zd_rzx  = Zref*X
      zd_rx   = X
c
c     ...mie coefficients and scattering pattern functions
      call mie_coeff  (Anzat, N,Zref,X,  zca,zcb,zcc,zcd)
      call mie_phitau (N, dcos(Theta), p, t, dp, dt);
c
c     ...spherical bessel function and its derivatvies
      call csphjy(N, z_rzx, N, zj_rzx, zjd_rzx, zy_rzx, zyd_rzx)
      call csphjy(N, z_rx,  N, zj_rx,  zjd_rx,  zy_rx,  zyd_rx )
c
c
ccc   ...use hankel first kind
      select case(Anzat)
      case('neg');
        zh_rx (:) = zj_rx (:) + z_i*zy_rx (:)
        zhd_rx(:) = zjd_rx(:) + z_i*zyd_rx(:)
      case('pos');
        zh_rx (:) = zj_rx (:) - z_i*zy_rx (:)
        zhd_rx(:) = zjd_rx(:) - z_i*zyd_rx(:)
      end select
c
ccc   ...initialization with zero
      Z_inc(1:3)      = z_0
      Z_int(1:3)      = z_0
      Z_sca(1:3)      = z_0
c
      ZD_inc(1:3,1:3) = z_0
      ZD_int(1:3,1:3) = z_0
      ZD_sca(1:3,1:3) = z_0
c
      do i=1, N
c
        z_n = (z_i**(i))*(2*i+1.d0)/(i*(i+1.d0))
ccc   ...Variables
ccc   zca,zcb,zcc,zcd = coefficients
ccc   z    = spherical coordinate r direction
ccc   zd   = derivative of z to the direction of r
ccc   zj   = spherical bessel value for the z
ccc   zjd  = derivative of (z*zj(z))
ccc   zjd2 = double derivative of (z*zj(z)) to phi and theta direction
ccc          this is not to r direction
ccc
ccc
ccc   ...Incident wave
        z     =   z_rx
        zd    =   zd_rx
        zj_a  =   zj_rx(i);        zj_b  =   zj_rx(i-1)
        zjd_a =   zjd_rx(i);       zjd_b =   zjd_rx(i-1)
ccc
        call Mie_mn(
     .    i, Theta, Phi,
     .    p(i), t(i), dp(i), dt(i),
     .    z, zd,
     .    zj_a, zj_b, zjd_a, zjd_b,
     .    zm, zn, zdm, zdn)
c
        Z_inc(1:3) = Z_inc(1:3)
     .    + z_n * ( zm(1:3) - z_i*zn(1:3) )
        ZD_inc(1:3,1:3) = ZD_inc(1:3,1:3)
     .    + z_n*( zdm(1:3,1:3) - z_i*zdn(1:3,1:3) )
ccc
ccc   ....Internal wave
        z     =   z_rzx
        zd    =   zd_rzx
        zj_a  =   zj_rzx(i);        zj_b  =   zj_rzx(i-1)
        zjd_a =   zjd_rzx(i);       zjd_b =   zjd_rzx(i-1)
ccc
        call Mie_mn(
     .    i, Theta, Phi,
     .    p(i), t(i), dp(i), dt(i),
     .    z, zd,
     .    zj_a, zj_b, zjd_a, zjd_b,
     .    zm, zn, zdm, zdn)
c
        Z_int(1:3) = Z_int(1:3)
     .    + z_n * ( zcc(i)*zm(1:3) - z_i*zcd(i)*zn(1:3) )
        ZD_int(1:3,1:3) = ZD_int(1:3,1:3)
     .    + z_n*( zcc(i)*zdm(1:3,1:3) - z_i*zcd(i)*zdn(1:3,1:3) )
ccc
ccc   ....Scattering wave
        z     =   z_rx
        zd    =   zd_rx
        zj_a  =   zh_rx(i);        zj_b  =   zh_rx(i-1)
        zjd_a =   zhd_rx(i);       zjd_b =   zhd_rx(i-1)
c
        call Mie_mn(
     .    i, Theta, Phi,
     .    p(i), t(i), dp(i), dt(i),
     .    z, zd,
     .    zj_a, zj_b, zjd_a, zjd_b,
     .    zm, zn, zdm, zdn)
c
        Z_sca(1:3) = Z_sca(1:3)
     .    + z_n * (-zcb(i)*zm(1:3) + z_i*zca(i)*zn(1:3) )
        ZD_sca(1:3,1:3) = ZD_sca(1:3,1:3)
     .    + z_n*(-zcb(i)*zdm(1:3,1:3) + z_i*zca(i)*zdn(1:3,1:3) )
c
      enddo
c
      end subroutine
c
c

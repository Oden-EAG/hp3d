c> Purpose : find out the truncation n for mie series.
cc @param[in] X - k*a (a-r of the sphere, k-wave number)
c------------------------------------------------------
      integer function mie_n_max(X)
#include "syscom.blk"
c
      mie_n_max = nint( X+4.d0*X**(1.d0/3.d0)+2.d0 )+1
      end function mie_n_max
c------------------------------------------------------
c> Purpose : Computes a matrix of Mie coefficients, a,b,c, and d
cc of orders n=1 to nmax, complex refractive index Zref
cc and size parameter x=k0*a, where k0= wave number
cc in the ambient medium, a=sphere radius;
cc p. 100, 477 in Bohren and Huffman (1983) BEWI:TDD122
cc C. Mätzler, June 2002
c------------------------------------------------------      
      subroutine mie_coeff(
     .  N, Zref,X,
     .  Za,Zb,Zc,Zd) 
#include "syscom.blk"
c      
      integer, parameter :: idelta = 1
c      
      dimension Za(N), Zb(N), Zc(N), Zd(N) 
c
c     local variables
      dimension
     .  zs1   (0:N+idelta),    zd1(0:N+idelta),
     .  zs2   (0:N+idelta),    zd2(0:N+idelta),
     .  zs1_zr(0:N+idelta), zd1_zr(0:N+idelta),
     .  zs1_h (0:N+idelta),  zd1_h(0:N+idelta)
c
      dimension
     .  s1    (0:N+idelta),     d1(0:N+idelta),
     .  s2    (0:N+idelta),     d2(0:N+idelta),
     .  s1_r  (0:N+idelta),   d1_r(0:N+idelta)
c
      dimension ztemp(2,2)
c
      z = cmplx(X, 0.d0)
      zref2 = Zref*Zref
      zr = Zref*z
c
      call sphj(N, real(z),  N, s1, d1)
      call sphy(N, real(z),  N, s2, d2)
      call sphj(N, real(zr), N, s1_r, d1_r)
c
      zs1    = s1;      zd1    = d1
      zs2    = s2;      zd2    = d2
      zs1_zr = s1_r;    zd1_zr = d1_r
c
      zs1_h(0:N+idelta) = zs1(0:N+idelta) +
     .                    cmplx(0.d0, 1.d0)*zs2(0:N+idelta)
      zd1_h(0:N+idelta) = zd1(0:N+idelta) +
     .                    cmplx(0.d0, 1.d0)*zd2(0:N+idelta)
c
      do i=1, N
c
        ztemp(1,1) = zs1_zr(i)*(zs1(i)+z*zd1(i))
        ztemp(1,2) = zs1(i)*(zs1_zr(i)+zr*zd1_zr(i))
        ztemp(2,1) = zs1_zr(i)*(zs1_h(i)+z*zd1_h(i))
        ztemp(2,2) = zs1_h(i)*(zs1_zr(i)+zr*zd1_zr(i))
c
        Za(i) = ( zref2*ztemp(1,1) - ztemp(1,2) )/
     .          ( zref2*ztemp(2,1) - ztemp(2,2) )
        
        Zb(i) = ( ztemp(1,1) - ztemp(1,2) )/
     .          ( ztemp(2,1) - ztemp(2,2) )
c
        ztemp(1,1) = zs1(i)*(zs1_h(i)+z*zd1_h(i))
        ztemp(1,2) = zs1_h(i)*(zs1(i)+z*zd1(i))
c
        Zc(i) = ( ztemp(1,1) - ztemp(1,2) )/
     .          ( ztemp(2,1) - ztemp(2,2) )
        Zd(i) = ( Zref*ztemp(1,1) - Zref*ztemp(1,2) )/
     .          ( zref2*ztemp(2,1) - ztemp(2,2) )
c
      enddo
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
      subroutine mie(
     .  N, Zref, X,
     .  Qext, Qsca, Qabs, Qb)
#include "syscom.blk"
      integer, parameter :: idelta = 2
      complex*16 :: zqb
      complex*16, dimension(N+idelta) :: za, zb, zc, zd
c      
      Qext = 0.d0; Qsca = 0.d0; Qabs = 0.d0; Qb = 0.d0;
      zqb  = cmplx(0.d0, 0.d0)
      if (X.eq.0.d0) return
c      
      call mie_coeff(N,Zref,X,  za,zb,zc,zd)
c
      write(*,*) zqb
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
c
c------------------------------------------------------            
      subroutine mie_phitau(N, U, P, T)
#include "syscom.blk"
      integer, parameter :: idelta = 2
      real*8, dimension(N+idelta) :: P, T
c
ccc      write(*,*) 'U,N ', U, N
      P(1) = 1.d0;        T(1) = U;
      P(2) = 3.d0*U;      T(2) = 3.d0*dcos(2.d0*acos(U));
      do i=3, N
        tmp1 = (2.d0*i-1.d0)/(i-1.d0)*P(i-1)*U;
        tmp2 = i/(i-1.d0)*P(i-2);
        P(i) = tmp1 - tmp2;
c
        tmp1 = i*U*P(i);
        tmp2 = (i+1.d0)*P(i-1);
        T(i) = tmp1 - tmp2;
      enddo
c
      end subroutine
c------------------------------------------------------            
c
c------------------------------------------------------      
      subroutine mie_s12(
     .  N, Zref, X, U,
     .  ZS1, ZS2)
#include "syscom.blk"
c
      integer, parameter :: idelta = 2
      real*8,     dimension(N+idelta) :: p, t
      complex*16, dimension(N+idelta) :: za, zb, zc, zd
c      
      call mie_coeff(N,Zref,X,  za,zb,zc,zd)
      call mie_phitau(N, U, p, t);
c
      ZS1 = cmplx(0.d0, 0.d0)
      ZS2 = cmplx(0.d0, 0.d0)
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
      subroutine mie_far_field(
     .  N, Zref, X, 
     .  Rk, Rad, Theta, Phi,
     .  ZST, ZSP)
#include "syscom.blk"
c
      call mie_s12(N, Zref, X, dcos(Theta), zs1, zs2)
c      
      ZST =
     .  exp(cmplx(0.d0, 1.d0)*Rk*Rad)
     .  /(-cmplx(0.d0, 1.d0)*Rk*Rad)
     .  *dcos(Phi)*zs2
c
      ZSP =
     .  exp(cmplx(0.d0, 1.d0)*Rk*Rad)
     .  /(cmplx(0.d0, 1.d0)*Rk*Rad)
     .  *dsin(Phi)*zs1
c     
      end subroutine
c
c
c------------------------------------------------------            
c
c------------------------------------------------------      
      subroutine mie_near_field(
     .  N, Zref, X,  
     .  Rad, Theta, Phi,
     .  Z_int, Z_sca)
#include "syscom.blk"
c      
c     ... Internal field on spherical coord
      complex*16, dimension(3)            :: Z_int, Z_sca
      
      integer, parameter                  :: idelta = 2
c      
      real*8,     dimension(N+idelta)     :: p, t
      complex*16, dimension(N+idelta)     :: za, zb, zc, zd
c
      complex*16, dimension(0:N+idelta-1) :: zjd, zj
      real*8,     dimension(0:N+idelta-1) :: rjd, rj
c
      complex*16, dimension(3)            :: zm, zn
      real*8, parameter :: eps = 1.0e-14
c     
      z   = Rad*Zref*X
c
      if (abs(z).lt.eps) then
        z = z + eps
      endif
c      
      call mie_coeff(N,Zref,X,  za,zb,zc,zd)
      call mie_phitau(N, dcos(Theta), p, t);
c
      call sphj(N, real(z), N, rj, rjd)
      zj = rj;    zjd = rjd      
c
      Z_int(1:3) = cmplx(0.d0, 0.d0)
      Z_sca(1:3) = cmplx(0.d0, 0.d0)
c
      do i=1, N
c
        zm(1) =  cmplx(0.d0, 0.d0)
        zm(2) =  dcos(Phi)*p(i)*zj(i)
        zm(3) = -dsin(Phi)*t(i)*zj(i)
c
        zdd   = (zj(i)+z*zjd(i))
c        
        zn(1) =  i*(i+1)*dcos(Phi)*dsin(Theta)*p(i)*zj(i)/z
        zn(2) =  dcos(Phi)*t(i)*zdd/z
        zn(3) = -dsin(Phi)*p(i)*zdd/z
c
        Z_n = (cmplx(0.d0, 1.d0)**i)*(2*i+1.d0)/(i*(i+1.d0))
        
        Z_int(1:3) = Z_int(1:3) 
     .    + Z_n * ( zc(i)*zm(1:3) - cmplx(0.d0, 1.d0)*zd(i)*zn(1:3) )

        Z_sca(1:3) = Z_sca(1:3)
     .    + Z_n * (-zb(i)*zm(1:3) + cmplx(0.d0, 1.d0)*za(i)*zn(1:3) )
c        
      enddo
c     
      end subroutine
c
c

c--------------------------------------------------------------------
c
c     routine name      - shape3Q
c
c--------------------------------------------------------------------
c
c     latest revision:  - Feb 14
c
c     purpose:          - routine returns values of 3D L2
c                         shape functions
c
c     arguments:
c
c     in:
c
c             Type      - element type
c             Xi        - master element coordinates
c             Norder    - polynomial order
c
c     out:
c             NrdofQ    - number of dof
c             ShapQ     - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape3Q_1(Type,Xi,Norder, NrdofQ,ShapQ)
c
      use parameters
#include "syscom.blk"
c
      character(len=4) :: Type
      dimension Norder(19)
      dimension Xi(3),ShapQ(MAXbrickQ)
c
c
      select case(Type)
      case('bric','mdlb')
        call shape3Q_bric(Xi,Norder(19), NrdofQ,ShapQ)
      case('tetr','mdln')
        call shape3Q_tetr(Xi,Norder(11), NrdofQ,ShapQ)
      case('pris','mdlp')
        call shape3Q_pris(Xi,Norder(15), NrdofQ,ShapQ)
      case('pyra','mdld')
        call shape3Q_pyra(Xi,Norder(14), NrdofQ,ShapQ)
      case default
        write(*,*)'shape3Q: Type = ', Type
        stop 1
      end select
c
      end subroutine shape3Q_1
c
c--------------------------------------------------------------------
c
c     routine name      - shape3Q_bric
c
c--------------------------------------------------------------------
c
c     latest revision:  - Dec 13
c
c     purpose:          - routine returns values of brick L2
c                         shape functions
c
c     arguments:
c
c     in:
c             Xi        - master bric coordinates from (0,1)^3
c             Nord      - polynomial order
c
c     out:
c             NrdofQ    - number of dof
c             ShapQ     - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape3Q_bric(Xi,Nord, NrdofQ,ShapQ)
c
      use parameters
#include "syscom.blk"
c
      dimension Xi(3),ShapQ(MAXbrickQ)
c
c  ...work space for shape1Q
      dimension shapQ1(MAXP),shapQ2(MAXP),shapQ3(MAXP)
c
      iprint=0
c
      call decode(Nord, norda,nord3)
      call decode(norda, nord1,nord2)
      call shape1Q(Xi(1),nord1, nrdofQ1,shapQ1)
      call shape1Q(Xi(2),nord2, nrdofQ2,shapQ2)
      call shape1Q(Xi(3),nord3, nrdofQ3,shapQ3)
c
      k=0
c
      do l=1,nrdofQ3
        do j=1,nrdofQ2
          do i=1,nrdofQ1
            k=k+1
            ShapQ(k) = shapQ1(i)*shapQ2(j)*shapQ3(l)
          enddo
        enddo
      enddo
c
      NrdofQ=k
c
      if (iprint.eq.1) then
        write(*,7001) Xi(1:3),Nord
 7001   format('shap3Q_bric: Xi = ',3f8.3,' Nord = ',i3)
        do k=1,NrdofQ
          write(*,7002) k,ShapQ(k)
 7002     format('k = ',i3,' ShapQ, = ',e12.5)
        enddo
        call pause
      endif
c
      end subroutine shape3Q_bric
c
c--------------------------------------------------------------------
c
c     routine name      - shape3Q_tetr
c
c--------------------------------------------------------------------
c
c     latest revision:  - Nov 13
c
c     purpose:          - routine returns values of tetrahedron L2
c                         shape functions
c
c     arguments:
c
c     in:
c             Xi        - master tetrahedron coordinates
c             Nord      - polynomial order
c
c     out:
c             NrdofQ    - number of dof
c             ShapQ     - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape3Q_tetr(Xi,Nord, NrdofQ,ShapQ)
c
      use parameters
#include "syscom.blk"
c
      dimension Xi(3),ShapQ(MAXtetraQ)
c
c  ...barycentric coordinates
      dimension vshap(0:3)
c
c  ...work space for routine Legendre
      dimension p1(0:MAXP-1),
     .          p2(0:MAXP+1,0:MAXP),
     .          p3(0:MAXP+1,0:MAXP),
     .          void(0:MAXP-1)
c
      double precision z
      integer alpha1,alpha2
c
      iprint=0
c
c
c  ...barycentric coordinates = vertex shape functions
      vshap(0) = 1.d0-Xi(1)-Xi(2)-Xi(3)
      vshap(1) = Xi(1)
      vshap(2) = Xi(2)
      vshap(3) = Xi(3)
c
c  ...determine shifted scaled Legendre polynomials in lambda_1
      x = vshap(1); t = vshap(0) + vshap(1)
      call LegendreNEW(x,t,Nord-1, p1)
c
c  ...determine shifted scaled Jacobi polynomials in lambda_2
      y = vshap(2); t = vshap(0) + vshap(1) + vshap(2)
      maxalpha = 2*Nord-1
      call Jacobi(y,t,maxalpha,Nord-1, p2)
c
c  ...determine shifted Jacobi polynomials in lambda_3
      z = vshap(3); t = 1.d0
      call Jacobi(z,t,maxalpha,Nord-1, p3)
c
      k=0
c
c  ...total order
      do np=0,Nord-1
c
c  .....order in lambda_3
        do l=0,np
c
c  .......order in lambda_2
          do j=0,np-l
c
c  .........order in lambda_1
            i=np-l-j
c
c  .........Jacobi weights in y and z
c            alpha1 = 2*i+1
c            alpha2 = 2*(i+j)+1
            k=k+1
            ShapQ(k) = p1(i)*p2(i,j)*p3(i+j,l)
          enddo
        enddo
      enddo
c
      NrdofQ=k
c
      if (iprint.eq.1) then
        write(*,7001) Xi(1:3),Nord
 7001   format('shap3Q_tetr: Xi = ',3f8.3,' Nord = ',i3)
        do k=1,NrdofQ
          write(*,7002) k,ShapQ(k)
 7002     format('k = ',i3,' ShapQ, = ',e12.5)
        enddo
        call pause
      endif
c
      end subroutine shape3Q_tetr
c
c--------------------------------------------------------------------
c
c     routine name      - shape3Q_pris
c
c--------------------------------------------------------------------
c
c     latest revision:  - Dec 13
c
c     purpose:          - routine returns values of prism L2
c                         shape functions
c
c     arguments:
c
c     in:
c             Xi        - master prism coordinates
c             Nord      - polynomial order
c
c     out:
c             NrdofQ    - number of dof
c             ShapQ     - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape3Q_pris(Xi,Nord, NrdofQ,ShapQ)
c
      use parameters
#include "syscom.blk"
c
      dimension Xi(3),ShapQ(MAXprismQ)
c
c  ...work space for routines shape2Q_tria, shape1Q
      dimension shapQt(MAXtriaQ),shapQ1(MAXP)
c
      iprint=0
c
      call decode(Nord, nord1,nord2)
      call shape2Q_tria(Xi   ,nord1, nrdofQt,shapQt)
      call shape1Q(     Xi(3),nord2, nrdofQ1,shapQ1)
c
      k=0
c
      do j=1,nrdofQ1
        do i=1,nrdofQt
          k=k+1
          ShapQ(k) = shapQt(i)*shapQ1(j)
        enddo
      enddo
c
      NrdofQ=k
c
      if (iprint.eq.1) then
        write(*,7001) Xi(1:3),Nord
 7001   format('shap3Q_pris: Xi = ',3f8.3,' Nord = ',i3)
        do k=1,NrdofQ
          write(*,7002) k,ShapQ(k)
 7002     format('k = ',i3,' ShapQ, = ',e12.5)
        enddo
        call pause
      endif
c
      end subroutine shape3Q_pris
c
c--------------------------------------------------------------------
c
c     routine name      - shape3Q_pyra
c
c--------------------------------------------------------------------
c
c     latest revision:  - Feb 14
c
c     purpose:          - routine returns values of pyramid L2
c                         shape functions
c
c     arguments:
c
c     in:
c             Xip       - master pyramid coordinates
c             Nord      - polynomial order
c
c     out:
c             NrdofQ    - number of dof
c             ShapQ     - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape3Q_pyra(Xip,Nord, NrdofQ,ShapQ)
c
      use parameters
#include "syscom.blk"
c
      dimension Xip(3),ShapQ(MAXpyramQ)
c
      double precision zeta,z
c
c  ...work space for Legendre and Jacobi
      dimension p1(0:MAXP),p2(0:MAXP),p3(0:MAXP+1,0:MAXP)
c
c  ...Jacobi weight
      integer a
c
c
      iprint=0
c
      xi = Xip(1); eta = Xip(2); zeta = Xip(3)
      x = xi/(1.d0-zeta)
      y = eta/(1.d0-zeta)
c
      k=0
c
c  ...just bubbles
      call LegendreNEW(   x,1.d0,Nord-1, p1)
      call LegendreNEW(   y,1.d0,Nord-1, p2)
      maxalpha = 2*Nord-1
      call Jacobi(zeta,1.d0,maxalpha,Nord, p3)
      do m=0,nord
        do j=0,nord-1
          do i=0,nord-1
            n = max(i+1,j+1,m)+1
c            alpha = 2*max(i+1,j+1)-1
            a = max(i+1,j+1)
            k=k+1
c
c  .........watch for the pullback factor
            ShapQ(k) = p1(i)*p2(j)*p3(a,m)*(1-zeta)**(n-m-2)
          enddo
        enddo
      enddo
c
c
      NrdofQ=k
c
      if (iprint.eq.1) then
        write(*,7001) Xip(1:3),Nord
 7001   format('shap3Q_pyra: Xi = ',3f8.3,' Nord = ',i3)
        do k=1,NrdofQ
          write(*,7002) k,ShapQ(k)
 7002     format('k = ',i3,' ShapQ, = ',e12.5)
        enddo
        call pause
      endif
c
      end subroutine shape3Q_pyra

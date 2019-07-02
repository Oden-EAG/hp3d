c--------------------------------------------------------------------
c
c     routine name      - shape2Q_1
c
c--------------------------------------------------------------------
c
c     latest revision:  - Feb 14
c
c     purpose:          - routine returns values of a 2D element
c                         L2 shape functions
c
c     arguments:
c
c     in:
c             Type      - element type
c             Xi        - master quad coordinates
c             Norder    - polynomial order
c
c     out:
c             NrdofQ    - number of dof
c             ShapQ     - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape2Q_1(Type,Xi,Norder, NrdofQ,ShapQ)
c
      use parameters
#include "syscom.blk"
c
      character(len=4) :: Type
      dimension Xi(2),Norder(5),
     .          ShapQ(MAXP**2)
c
      select case(Type)
      case('tria','mdlt')
        call shape2Q_tria(Xi,Norder(4), NrdofQ,ShapQ)
      case('quad','mdlq')
        call shape2Q_quad(Xi,Norder(5), NrdofQ,ShapQ)
      case default
        write(*,*)'shape2E: Type = ', Type
        stop 1
      end select
c
      end subroutine shape2Q_1
c
c--------------------------------------------------------------------
c
c     routine name      - shape2Q_quad
c
c--------------------------------------------------------------------
c
c     latest revision:  - Feb 14
c
c     purpose:          - routine returns values of quad L2
c                         shape functions
c
c     arguments:
c
c     in:
c             Xi        - master quad coordinates from (0,1)x(0,1)
c             Nord      - polynomial order
c
c     out:
c             NrdofQ    - number of dof
c             ShapQ     - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape2Q_quad(Xi,Nord, NrdofQ,ShapQ)
c
      use parameters
#include "syscom.blk"
c
      dimension Xi(2),ShapQ(MAXP**2)
c
c  ...work space for routine shape1Q
      dimension shapQ1(MAXP),shapQ2(MAXP)
c
      iprint=0
c
      call decode(Nord, nord1,nord2)
      call shape1Q(Xi(1),nord1, nrdofQ1,shapQ1)
      call shape1Q(Xi(2),nord2, nrdofQ2,shapQ2)
c
      k=0
c
      do j=1,nrdofQ2
        do i=1,nrdofQ1
          k=k+1
          ShapQ(k) = shapQ1(i)*shapQ2(j)
        enddo
      enddo
c
      NrdofQ=k
c
      if (iprint.eq.1) then
        write(*,7001) Xi(1:2),Nord
 7001   format('shape2Q_quad: Xi = ',2f8.3,' Nord = ',i2)
        do k=1,NrdofQ
          write(*,7002) k,ShapQ(k)
 7002     format('k = ',i3,' ShapQ = ',e12.5)
        enddo
        call pause
      endif
c
      end subroutine shape2Q_quad
c
c--------------------------------------------------------------------
c
c     routine name      - shape2Q_tria
c
c--------------------------------------------------------------------
c
c     latest revision:  - Feb 14
c
c     purpose:          - routine returns values of quad L2
c                         shape functions
c
c     arguments:
c
c     in:
c             Xi        - master quad coordinates from (0,1)x(0,1)
c             Nord      - polynomial order
c
c     out:
c             NrdofQ    - number of dof
c             ShapQ     - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape2Q_tria(Xi,Nord, NrdofQ,ShapQ)
c
      use parameters
#include "syscom.blk"
c
      dimension Xi(2),ShapQ(MAXtriaQ)
c
c  ...work space for routine Legendreand Jacobi
      dimension p1(0:MAXP-1),p2(0:MAXP+1,0:MAXP)
c
c  ...Jacobi weight
      integer alpha
c
      iprint=0
c
c  ...determine shifted Legendre polynomials in the first variable
      x = Xi(1); t = 1.d0 - Xi(2)
      call LegendreNEW(x,t,Nord-1, p1)
c
c  ...determine shifted Jacobi polynomials in the second variable
      x = Xi(2); t = 1.d0
      maxalpha = 2*Nord-1
      call Jacobi(x,t,maxalpha,Nord-1, p2)
c
c  ...shape function counter
      k=0
c
c  ...watch for the structure of the loop to enforce hierarchy
c     in polynomial order, np = total order
      do np=0,Nord-1
c
c  .....order in y
        do j=0,np
c
c  .......order in x
          i=np-j
c
c  .......weight in y
c          alpha = 2*i+1
          k=k+1
          ShapQ(k) = p1(i)*p2(i,j)
        enddo
c
c  ...end of loop through order of polynomials
      enddo
      NrdofQ=k
c
c
      if (iprint.eq.1) then
        write(*,7001) Xi(1:2),Nord
 7001   format('shape2Q_tria: Xi = ',2f8.3,' Nord = ',i2)
        do k=1,NrdofQ
          write(*,7002) k,ShapQ(k)
 7002     format('k = ',i3,' ShapQ = ',e12.5)
        enddo
        call pause
      endif
c
      end subroutine shape2Q_tria
c--------------------------------------------------------------------
c
c     routine name      - shape2Q_orient
c
c--------------------------------------------------------------------
c
c     latest revision:  - May 14
c
c     purpose:          - routine returns values of a 2D element
c                         L2 shape functions accounting for orientation
c
c     arguments:
c
c     in:
c             Type      - element type
c             Xi        - master quad coordinates
c             Norder    - polynomial order of element nodes
c             Norient   - element orientation
c
c     out:
c             NrdofQ    - number of dof
c             ShapQ     - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape2Q_orient(Type,Xi,Norder,Norient, NrdofQ,ShapQ)
c
      use parameters
      use element_data
#include "syscom.blk"
c
      character(len=4) :: Type
      dimension Xi(2),Norder(5),
     .          ShapQ(MAXP**2)
c
c  ...1D L2 shape functions
      dimension shapQ1(MAXP),shapQ2(MAXP)
c
c  ...translating quad orientation into product of 1D orientations
      dimension norient221(1:3,0:7)
      data norient221/0,0,0, 1,1,0, 0,1,1, 1,0,1,
     .                1,0,0, 0,1,0, 1,1,1, 0,0,1/
c
      iprint=0
c
      select case(Type)
      case('tria','mdlt')
        call  shape2Q_tria_orient(Xi,Norder(4),Norient, NrdofQ,ShapQ)
      case('quad','mdlq','rect')
        call decode(Norder(5), nord1,nord2)
        nor1 = norient221(1,Norient)
        call shape1Q_orient(Xi(1),nord1,nor1, nrdofQ1,shapQ1)
        nor2 = norient221(2,Norient)
        call shape1Q_orient(Xi(2),nord2,nor2, nrdofQ2,shapQ2)
        k=0
        do j=1,nrdofQ2
          do i=1,nrdofQ1
            k=k+1
            ShapQ(k) = shapQ1(i)*shapQ2(j)
          enddo
        enddo
        NrdofQ = k
      case default
        write(*,*)'shape2Q_orient: Type = ', Type
        stop 1
      end select
      if (iprint.eq.1) then
        write(*,7001) Type,Xi(1:2),Norient
 7001   format('shape2Q_orient: Type = ',a4,' Xi = ',2f8.3,
     .         ' Norient = ',i2)
        write(*,7002) Norder(1:nedge(Type)+1)
 7002   format('                Norder = ',4i2,i3)
        write(*,7003)
 7003   format('shape2Q_orient: ShapQ = ')
        write(*,7010) ShapQ(1:NrdofQ)
 7010   format(10(10e12.5,/))
        call pause
      endif
c
      end subroutine shape2Q_orient
c
c--------------------------------------------------------------------
c
c     routine name      - shape2Q_tria_orient
c
c--------------------------------------------------------------------
c
c     latest revision:  - Feb 14
c
c     purpose:          - routine returns values of triangle L2
c                         shape functions accounting for orientation
c
c     arguments:
c
c     in:
c             Xi        - master quad coordinates from (0,1)x(0,1)
c             Nord      - polynomial order
c             Norient   - orientation of the triangle
c
c     out:
c             NrdofQ    - number of dof
c             ShapQ     - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape2Q_tria_orient(Xi,Nord,Norient, NrdofQ,ShapQ)
c
      use parameters
#include "syscom.blk"
c
      dimension Xi(2),ShapQ(MAXmdltQ)
c
c  ...face coordinates
      dimension xif(2)
c
      iprint=0
c
c  ...switch to face coordinates
      x = Xi(1); y=Xi(2)
      select case(Norient)
      case(0)
        xif(1) = x; xif(2) = y; det = 1.d0
      case(1)
        xif(1) = y; xif(2) = 1.d0-x-y; det = 1.d0
      case(2)
        xif(1) = 1.d0-x-y; xif(2) = x; det = 1.d0
      case(3)
        xif(1) = y; xif(2) = x; det = -1.d0
      case(4)
        xif(1) = 1.d0-x-y; xif(2) = y; det = -1.d0
      case(5)
        xif(1) = x; xif(2) = 1.d0-x-y; det = -1.d0
      end select
c
      call shape2Q_tria(xif,Nord, NrdofQ,ShapQ)
c
c  ...Piola transform
      ShapQ(1:NrdofQ) = ShapQ(1:NrdofQ)/det
c
c
      if (iprint.eq.1) then
        write(*,7001) Xi(1:2),Nord,Norient
 7001   format('shape2Q_tria_orient: Xi = ',2f8.3,' Nord = ',i2,
     .         ' Norient = ',i2)
        do k=1,NrdofQ
          write(*,7002) k,ShapQ(k)
 7002     format('k = ',i3,' ShapQ = ',e12.5)
        enddo
        call pause
      endif
c
      end subroutine shape2Q_tria_orient



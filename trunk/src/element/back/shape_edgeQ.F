c--------------------------------------------------------------------
c
c     routine name      - shape_edgeQ
c
c--------------------------------------------------------------------
c
c     latest revision:  - Feb 14
c
c     purpose:          - routine returns extensions of edge L2 shape
c                         functions using triangle affine coordinates
c
c     arguments:
c
c     in:
c             Vshap     - edge vertex shape functions (affine
c                         coordinates)
c             Nord      - polynomial order
c
c     out:
c             NrdofQ    - number of shape functions
c             ShapQ     - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape_edgeQ(Vshap,Nord, NrdofQ,ShapQ)
c
      use parameters
#include "syscom.blk"
c
      dimension Vshap(0:1), ShapQ(MAXP)
c
c  ...work space for routine Legendre
      dimension P(0:MAXP-1)
c
      iprint=0
c
c  ...use shifted scaled Legendre polynomials
      x =  Vshap(1); t = Vshap(0) + Vshap(1)
      call LegendreNEW(x,t,Nord-1, P)
      NrdofQ = Nord
c
      ShapQ(1:NrdofQ) = P(0:NrdofQ-1)
c
c
      if (iprint.eq.1) then
        write(*,7001) Nord, Vshap(0:1)
 7001   format('shape_edgeQ: Nord = ',i2,' Vshap = ',2F8.3)
        do k=1,NrdofQ
          write(*,7002) k,ShapQ(k)
 7002     format(' k = ',i2,' ShapQ = ',e12.5)
        enddo
        call pause
      endif
c
      end subroutine shape_edgeQ
c
c--------------------------------------------------------------------
c
c     routine name      - shape_edgeQ_orient
c
c--------------------------------------------------------------------
c
c     latest revision:  - Feb 14
c
c     purpose:          - routine returns extensions of edge L2 shape
c                         functions using triangle affine coordinates
c                         accounting for orientation
c
c     arguments:
c
c     in:
c             Vshap     - edge vertex shape functions (affine
c                         coordinates)
c             Nord      - polynomial order
c             Norient   - edge orientation
c
c     out:
c             NrdofQ    - number of shape functions
c             ShapQ     - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape_edgeQ_orient(Vshap,Nord,Norient, NrdofQ,ShapQ)
c
      use parameters
#include "syscom.blk"
c
      dimension Vshap(0:1), ShapQ(MAXP)
      dimension vshap_loc(0:1)
c
      iprint=0
c
c  ...account for orientation
      select case(Norient)
      case(0)
        n0=0; n1=1
      case(1)
        n0=1; n1=0
      end select
      vshap_loc(0) = Vshap(n0); vshap_loc(1) = Vshap(n1)
      call shape_edgeQ(vshap_loc,Nord, NrdofQ,ShapQ)
c
c
      if (iprint.eq.1) then
        write(*,7001) Nord, Norient,Vshap
 7001   format('shape_edgeQ_orient: Nord = ',i2,' Norient = ',i1,
     .         ' Vshap = ',2F8.3)
        do k=1,NrdofQ
          write(*,7002) k,ShapQ(k)
 7002     format(' k = ',i2,' ShapQ = ',e25.15)
        enddo
        call pause
      endif
c
c
      end subroutine shape_edgeQ_orient



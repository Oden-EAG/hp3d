c--------------------------------------------------------------------
c
c     routine name      - shape1Q_1
c
c--------------------------------------------------------------------
c
c     latest revision:  - Feb 14
c
c     purpose:          - routine returns values of 1D L2 shape
c                         functions
c
c     arguments:
c
c     in:
c             Xi        - master element coordinate from [0,1]
c             Nord      - polynomial order
c
c     out:
c             NrdofQ    - number of shape functions
c             ShapQ     - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape1Q_1(Xi,Nord, NrdofQ,ShapQ)
c
      use parameters
#include "syscom.blk"
c
      dimension ShapQ(MAXP)
c
      iprint=0
c
c  ...number of L2 shape functions
      NrdofQ = Nord
c
c  ...no scaling
      t = 1.d0
      call LegendreNEW(Xi,t,Nord-1, ShapQ)
c
c
      if (iprint.eq.1) then
        write(*,7001) Nord, Xi
 7001   format('shape1Q: Nord = ',i2,' Xi = ',F8.3)
        do k=1,NrdofQ
          write(*,7002) k,ShapQ(k)
 7002     format(' k = ',i2,' ShapQ = ',e25.15)
        enddo
        call pause
      endif
c
c
      end subroutine shape1Q_1
c
c--------------------------------------------------------------------
c
c     routine name      - shape1Q_orient
c
c--------------------------------------------------------------------
c
c     latest revision:  - Feb 14
c
c     purpose:          - routine returns values of 1D L2 shape
c                         functions accounting for orientation
c
c     arguments:
c
c     in:
c             Xi        - master element coordinate from [0,1]
c             Nord      - polynomial order
c             Norient   - orientation
c
c     out:
c             NrdofQ    - number of shape functions
c             ShapQ     - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape1Q_orient(Xi,Nord,Norient, NrdofQ,ShapQ)
c
      use parameters
#include "syscom.blk"
c
      dimension ShapQ(MAXP)
c
      iprint=0
c
      select case(Norient)
      case(0)
        call shape1Q(Xi,Nord, NrdofQ,ShapQ)
      case(1)
        xi1 = 1.d0 - Xi
        call shape1Q(xi1,Nord, NrdofQ,ShapQ)
        ShapQ(1:NrdofQ) = - ShapQ(1:NrdofQ)
      end select
c
c
      if (iprint.eq.1) then
        write(*,7001) Nord, Xi
 7001   format('shape1Q: Nord = ',i2,' Xi = ',F8.3)
        do k=1,NrdofQ
          write(*,7002) k,ShapQ(k)
 7002     format(' k = ',i2,' ShapQ = ',e25.15)
        enddo
        call pause
      endif
c
c
      end subroutine shape1Q_orient


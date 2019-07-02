c--------------------------------------------------------------------
c
c     routine name      - shape1H
c
c--------------------------------------------------------------------
c
c     latest revision:  - Feb 14
c
c     purpose:          - routine returns values of 1D H1 shape
c                         functions and their derivatives
c
c     arguments:
c
c     in:
c             Xi        - master edge coordinate (0,1)
c             Nord      - polynomial order
c
c     out:
c             NrdofH    - number of dof
c             ShapH     - values of the shape functions at the point
c             DshapH    - derivatives of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape1H_1(Xi,Nord, NrdofH,ShapH,DshapH)
c
      use parameters
#include "syscom.blk"
c
      dimension ShapH(MAXP+1),DshapH(MAXP+1)
c
      iprint=0
c
      NrdofH = Nord+1
c
c  ...vertex shape functions (1D affine coordinates)
      ShapH(1) = 1.d0-Xi; DshapH(1) = -1.d0
      ShapH(2) =      Xi; DshapH(2) =  1.d0
      if (Nord.eq.1) go to 999
      call shape1Hb(Xi,Nord, Nvoid,ShapH(3:MAXP+1),DshapH(3:MAXP+1))
c
 999  continue
      if (iprint.eq.1) then
        write(*,7001) Xi,Nord
 7001   format('shape1H: Xi,Nord = ',f8.3,2x,i2)
        do k=1,NrdofH
          write(*,7002) k,ShapH(k),DshapH(k)
 7002     format(' k = ',i2,' ShapH,DshapH = ',2e12.5)
        enddo
        call pause
      endif
c
c
      end subroutine shape1H_1
c
c--------------------------------------------------------------------
c
c     routine name      - shape1Hb
c
c--------------------------------------------------------------------
c
c     latest revision:  - Feb 14
c
c     purpose:          - routine returns values of 1D H1 bubble
c                         shape functions and their derivatives
c
c     arguments:
c
c     in:
c             Xi        - master edge coordinate (0,1)
c             Nord      - polynomial order
c
c     out:
c             NrdofH    - number of dof
c             ShapH     - values of the shape functions at the point
c             DshapH    - derivatives of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape1Hb(Xi,Nord, NrdofH,ShapH,DshapH)
c
      use parameters
#include "syscom.blk"
c
      dimension ShapH(MAXP-1),DshapH(MAXP-1)
c
      iprint=0
c
      NrdofH = Nord-1
c
c  ...no scaling
      t = 1.d0
      call ILegendre(Xi,t,Nord,0, ShapH,DshapH,void)
c
 999  continue
      if (iprint.eq.1) then
        write(*,7001) Xi,Nord
 7001   format('shape1Hb: Xi,Nord = ',f8.3,2x,i2)
        do k=1,NrdofH
          write(*,7002) k,ShapH(k),DshapH(k)
 7002     format(' k = ',i2,' ShapH,DshapH = ',2e12.5)
        enddo
        call pause
      endif
c
c
      end subroutine shape1Hb
c
c--------------------------------------------------------------------
c
c     routine name      - shape1Hb_orient
c
c--------------------------------------------------------------------
c
c     latest revision:  - Feb 14
c
c     purpose:          - routine returns values of 1D H1 bubble
c                         shape functions and their derivatives
c                         accounting for orientation
c
c     arguments:
c
c     in:
c             Xi        - master edge coordinate (0,1)
c             Nord      - polynomial order
c             Norient   - orientation
c
c     out:
c             NrdofH    - number of dof
c             ShapH     - values of the shape functions at the point
c             DshapH    - derivatives of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape1Hb_orient(Xi,Nord,Norient, NrdofH,ShapH,DshapH)
c
      use parameters
#include "syscom.blk"
c
      dimension ShapH(MAXP-1),DshapH(MAXP-1)
c
      iprint=0
c
      select case(Norient)
      case(0)
        call shape1Hb(Xi,Nord, NrdofH,ShapH,DshapH)
      case(1)
        xi1 = 1.d0 - Xi
        call shape1Hb(xi1,Nord, NrdofH,ShapH,DshapH)
        DshapH(1:NrdofH) = - DshapH(1:NrdofH)
      end select
c
 999  continue
      if (iprint.eq.1) then
        write(*,7001) Xi,Nord,Norient
 7001   format('shape1Hb_orient: Xi,Nord,Norient = ',f8.3,2x,2i2)
        do k=1,NrdofH
          write(*,7002) k,ShapH(k),DshapH(k)
 7002     format(' k = ',i2,' ShapH,DshapH = ',2e12.5)
        enddo
        call pause
      endif
c
c
      end subroutine shape1Hb_orient

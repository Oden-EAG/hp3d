c------------------------------------------------------------------------
c
c   routine name       - recta
c
c------------------------------------------------------------------------
c
c   latest revision    - Jun 07
c
c   purpose            - routine evaluates physical coordinates
c                        and its derivatives wrt to reference
c                        coordinates for a point in a GMP rectangle
c                        USING A BILINEAR APPROXIMATION
c
c   arguments :
c     in:
c               No     - the rectangle number
c               Eta    - reference coordinates of a point
c     out:
c               X      - physical coordinates of the point
c               Dxdeta - derivatives of the physical coordinates wrt
c                        to the reference coordinates
c
c  REMARK, Nov 12 : do not remove! needed for graphics
c------------------------------------------------------------------------
c
      subroutine recta_linear(No,Eta, X,Dxdeta)
c
      use GMP
#include "syscom.blk"
c
      dimension Eta(2),X(3),Dxdeta(NDIM,2)
c
c  ...vertex coordinates, work space
      dimension xv(NDIM,4),aux(NDIM,3)
c
c------------------------------------------------------------------------
c       
      iprint=0
c
c  ...evaluate vertex coordinates
      do i=1,4
        np = RECTANGLES(No)%VertNo(i)
        call pointr(np, xv(1:NDIM,i))
      enddo
c
c  ...evaluate coefficients defining the bilinear map
      do j=1,NDIM
        aux(j,1) = xv(j,2)-xv(j,1)
        aux(j,2) = xv(j,4)-xv(j,1)
        aux(j,3) = xv(j,1)-xv(j,2)+xv(j,3)-xv(j,4)
      enddo
c
c  ...evaluate physical coordinates and their derivatives
      do j=1,NDIM
        Dxdeta(j,1) = aux(j,1) + Eta(2)*aux(j,3)
        Dxdeta(j,2) = aux(j,2) + Eta(1)*aux(j,3)
        X(j) = xv(j,1) + Eta(1)*aux(j,1) + Eta(2)*aux(j,2)
     .       + Eta(1)*Eta(2)*aux(j,3)
      enddo
c
c
      end









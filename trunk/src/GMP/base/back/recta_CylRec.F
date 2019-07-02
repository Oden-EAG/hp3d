c//////////////////////////////////////////////////////////////////////
c P. Gatto, Nov 12 : this is a LEGACY routine. Do not use it, unless
c                    you revise it!
c//////////////////////////////////////////////////////////////////////
c
c
c-----------------------------------------------------------------------
c
c   routine name       - recta_CylRec
c
c-----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Jan 06
c
c   purpose            - routine evaluates physical coordinates
c                        and its derivatives wrt to reference
c                        coordinates (theta,z) for a point in a quadrant
c                        of a cylinder
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
c-----------------------------------------------------------------------
c
      subroutine recta_CylRec(No,Eta, X,Dxdeta)
c
      use GMP
#include "syscom.blk"
c
      dimension Eta(2),X(3),Dxdeta(3,2)
c
c  ...vertex coordinates, work space
      dimension xv(3,4),aux(3,3)
c
c  ...center coordinates
      dimension center(3)
c
c  ...transformation matrix from local to global coordinates
      dimension aij(3,3)
c
c  ...local coordinates and their derivatives wrt reference
c     coordinates
      dimension xprim(3),dxprdxi(3,2)
c
c-----------------------------------------------------------------------
c
      iprint = 0
c
      if ((RECTANGLES(No)%Type.ne.'CylRec').or.(NDIM.ne.3)) then
        write(*,7001)
 7001   format('recta_CylRec: WRONG CALL TO recta_CylRec')
        stop 1
      endif
c
c  ...evaluate vertex coordinates
      do i=1,4
        np = RECTANGLES(No)%VertNo(i)
        call pointr(np, xv(1:NDIM,i))
      enddo
c
c  ...evaluate height and radius
      height=0.d0
      rad = 0.d0
      do i=1,3
        aux(i,1) = xv(i,2)-xv(i,1)
        aux(i,2) = xv(i,4)-xv(i,1)
        rad = rad + (aux(i,1))**2
        height = height + (aux(i,2))**2
      enddo
      height = sqrt(height)
      rad = sqrt(rad/2.d0)
c
c  ...evaluate coordinates of the center
      aux(1,3) = aux(2,1)*aux(3,2) - aux(3,1)*aux(2,2)
      aux(2,3) = aux(3,1)*aux(1,2) - aux(1,1)*aux(3,2)
      aux(3,3) = aux(1,1)*aux(2,2) - aux(2,1)*aux(1,2)
      s = 0.d0
      do i=1,3
        s = s + aux(i,3)**2
      enddo
      s = sqrt(s)
      sr = sqrt(2.d0)
      do i=1,3
        center(i) = (xv(i,1)+xv(i,2))/2.d0 - aux(i,3)/s*rad/sr
      enddo
c
c  ...evaluate the transformation matrix from local to global
c     coordinates
      do i=1,3
        aij(i,1) = (xv(i,1) - center(i))/rad
        aij(i,2) = (xv(i,2) - center(i))/rad
        aij(i,3) = (xv(i,3) - xv(i,2))/height
      enddo
c
c  ...evaluate coordinates and their derivatives in the auxilary
c     system of coordinates (see manual for explanation)
      pihalf = dacos(0.d0)
      theta = pihalf*Eta(1)
      costh = cos(theta)
      sinth = sin(theta)
      xprim(1) = rad*costh
      xprim(2) = rad*sinth
      xprim(3) = height*Eta(2)
      dxprdxi(1,1) = -rad*pihalf*sinth
      dxprdxi(2,1) =  rad*pihalf*costh
      dxprdxi(3,1) =  0.d0
      dxprdxi(1,2) =  0.d0
      dxprdxi(2,2) =  0.d0
      dxprdxi(3,2) =  height
c
c  ...evaluate physical coordinates and their derivatives wrt the
c     reference coordinates
      do i=1,3
        s  = 0.d0
        s1 = 0.d0
        s2 = 0.d0
        do j=1,3
          s  = s  + aij(i,j)*xprim(j)
          s1 = s1 + aij(i,j)*dxprdxi(j,1)
          s2 = s2 + aij(i,j)*dxprdxi(j,2)
        enddo
        X(i) = s + center(i)
        Dxdeta(i,1) = s1
        Dxdeta(i,2) = s2
      enddo
c
      return
      end









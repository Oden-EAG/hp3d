c-----------------------------------------------------------------------
c
c   routine name       - trian_QtCirTri
c
c-----------------------------------------------------------------------
c
c   latest revision    - Jan 06
c
c   purpose            - routine evaluates physical coordinates
c                        and its derivatives wrt to reference
c                        coordinates for a point in a circular
c                        (planar) quadrant of a circle
c                        (old routine: use with caution...)
c
c   arguments :
c     in:
c               No     - a GMP triangle number
c               Eta    - reference coordinates of a point
c                        in the triangle
c     out:
c               X      - physical coordinates of the point
c               Dxdeta - derivatives of the physical coordinates wrt
c                        to the parameters
c
c
c-----------------------------------------------------------------------
c
      subroutine trian_QtCirTri(No,Eta, X,Dxdeta)
c
      use GMP
#include "syscom.blk"
c
      dimension Eta(2),X(NDIM),Dxdeta(NDIM,2)
c
c  ...vertex coordinates, edge curves numbers
      dimension xv(3,3),ncurv(3)
c
c  ...transformatiom matrix
      dimension aij(1:NDIM,3)
c
c  ...auxiliary quantities
      dimension drdeta(2),center(3),dthdeta(2),xprim(3),dxprdeta(3,2)
c
      data one /1.d0/
      iprint=0
c
      if (TRIANGLES(No)%Type.ne.'QtCirTr') then
        write(*,7001) TRIANGLES(No)%Type
 7001   format('tria_QtCirTr: WRONG TRIANGLE TYPE = ',a10)
        stop 1
      endif
c
      if (iprint.eq.1) then
        write(*,7002) No,Eta
 7002   format('trian_QtCirTr: No,Eta = ',i4,2x,2f8.3)
      endif
c
c  ...get edge curves and vertex coordinates
      do i=1,3
        ncurv(i) = TRIANGLES(No)%EdgeNo(i)
        np = TRIANGLES(No)%VertNo(i)
        call pointr(np, xv(1:NDIM,i))
      enddo
c
c  ...check compatibility
      do i=1,3
        nc = abs(ncurv(i))
        select case(i)
        case(1,3)
        if (CURVES(nc)%Type.ne.'Seglin') then
          write(*,7004) No,i,nc,CURVES(nc)%Type
 7004       format('trian_QtCirTr: INCOMPATIBLE DEFINITIONS: No,i,ic=',
     .              i4,2x,i2,i4,' type = ',a10)
          stop 1
        endif
        case(2)
        if (CURVES(nc)%Type.ne.'QuaCir') then
          write(*,7004) No,i,nc,CURVES(nc)%Type
          stop 1
        endif
        end select
      enddo
c
c-----------------------------------------------------------------------
c
c  ...evaluate radius and center of the circle
      rad=0.d0
      do i=1,NDIM
        rad = rad + (xv(i,2)-xv(i,1))**2
        center(i) = xv(i,1)
      enddo
      rad = sqrt(rad)
c
c  ...evaluate the transformation matrix
      do i=1,NDIM
        aij(i,1) = (xv(i,2) - center(i))/rad
        aij(i,2) = (xv(i,3) - center(i))/rad
      enddo
      if (NDIM.eq.3) then
        aij(1,3) = aij(2,1)*aij(3,2) - aij(3,1)*aij(2,2)
        aij(2,3) = aij(3,1)*aij(1,2) - aij(1,1)*aij(3,2)
        aij(3,3) = aij(1,1)*aij(2,2) - aij(2,1)*aij(1,2)
      endif
c
c  ...evaluate coordinates and their derivatives in the auxilary
c     system of coordinates (see manual for explanation)
      pihalf = dacos(0.d0)
      aa = Eta(1) + Eta(2)
      if (aa.lt..001d0) aa=.001d0
      aa2 = aa*aa
      theta = pihalf*Eta(2)/aa
      dthdeta(1) = -pihalf*Eta(2)/aa2
      dthdeta(2) =  pihalf*Eta(1)/aa2
      r = rad*aa
      drdeta(1) = rad
      drdeta(2) = rad
c
      costh = cos(theta)
      sinth = sin(theta)
      xprim(1) = r*costh
      xprim(2) = r*sinth
      xprim(3) = 0.d0
      do j=1,2
        dxprdeta(1,j) = costh*drdeta(j) - r*sinth*dthdeta(j)
        dxprdeta(2,j) = sinth*drdeta(j) + r*costh*dthdeta(j)
        dxprdeta(3,j) = 0.d0
      enddo
c
c  ...evaluate physical coordinates and their derivatives wrt the
c     parameters
      do i=1,NDIM
        s  = 0.d0
        s1 = 0.d0
        s2 = 0.d0
        do j=1,NDIM
          s  = s  + aij(i,j)*xprim(j)
          s1 = s1 + aij(i,j)*dxprdeta(j,1)
          s2 = s2 + aij(i,j)*dxprdeta(j,2)
        enddo
        X(i) = s + center(i)
        Dxdeta(i,1) = s1
        Dxdeta(i,2) = s2
      enddo
c
      end



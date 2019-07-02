c-----------------------------------------------------------------------
c
c   routine name       - trian_PaSphTri
c
c-----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Jan 06
c
c   purpose            - routine evaluates physical coordinates
c                        and its derivatives wrt to reference
c                        coordinates for a point in a triangular
c                        part of a sphere (old routine, use with
c                        caution...)
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
      subroutine trian_PaSphTri(No,Eta, X,Dxdeta)
c
      use GMP
      use control
#include "syscom.blk"
c
      dimension Eta(2),X(3),Dxdeta(3,2)
c
c  ...vertex coordinates, edge curves
      dimension xv(3,3),ncurv(3)
c
c  ...transformation matrix
      dimension aij(3,3)
c
c  ...centers for the three circular segments
      dimension centers(3,3)
c
c  ...auxiliary quantities
      dimension dthdeta(2),xprim(3),dpsdeta(2),dxprdeta(3,2),aux(3)
c
c-----------------------------------------------------------------------
c
      iprint=0
c
      if (TRIANGLES(No)%Type.ne.'PaSphTri') then
        write(*,7000) TRIANGLES(No)%Type
 7000   format('tria_PaSphTri: WRONG TRIANGLE TYPE = ',a10)
        stop 1
      endif
c
      if (iprint.eq.1) then
        write(*,7001) No,Eta
 7001   format('trian_PaSphTri: No,Eta = ',i4,2x,2f8.3)
      endif
c
c  ...get the edges and then the vertices
      do i=1,3
        ncurv(i) = TRIANGLES(No)%EdgeNo(i)
        np = TRIANGLES(No)%VertNo(i)
        call pointr(np, xv(1:NDIM,i))
      enddo
c
c  ...check compatibility
      do i=1,3
        nc = iabs(ncurv(i))
        select case(i)
        case(1,3)
        if (CURVES(nc)%Type.ne.'SegCir') then
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
        centers(1:3,i) = CURVES(iabs(nc))%Rdata(1:3)
        end select
      enddo
      d = 0.d0
      do i=1,3
        d = d + (centers(j,1) - centers(j,3))**2
      enddo
      d = sqrt(d)
      if (d.gt.GEOM_TOL) then
        write(*,7007) centers
 7007   format('trian_QtCirTr: INCOMPATIBLE centers = ',3(3f8.3,2x))
        stop 1
      endif
c
c-----------------------------------------------------------------------
c
c  ...evaluate radius and beta
      radc=0.d0
      rads=0.d0
      skpr=0.d0
      do i=1,3
        aij(i,1) = xv(i,2) - centers(i,2)
        aij(i,2) = xv(i,3) - centers(i,2)
        aij(i,3) = xv(i,1) - centers(i,1)
        aux(i)  = xv(i,2) - centers(i,1)
        radc = radc + (aij(i,1))**2
        rads = rads + (aij(i,3))**2
        skpr = skpr + (aij(i,3))*aux(i)
      enddo
      radc = sqrt(radc)
      rads = sqrt(rads)
      cobeta = skpr/(rads*rads)
c
      if (iprint.eq.1) then
        do j=1,3
          write(*,7005) j,xv(j,1:3)
 7005     format('trian_PaSphTri: j,xv(j,*) = ',i2,2x,3f8.3)
        enddo
        write(*,7006) centers(1:3,1:2), radc,rads
 7006   format('trian_PaSphTri: centers, radc,rads = ',
     .          2(3f8.3,2x),2f8.3)
      endif
c
c  ...evaluate the transformation matrix
      do i=1,3
        aij(i,1) = (aij(i,1))/radc
        aij(i,2) = (aij(i,2))/radc
        aij(i,3) = (aij(i,3))/rads
      enddo
c
c  ...evaluate coordinates and their derivatives in the auxilary
c     system of coordinates (see manual for explanation)
      pihalf = dacos(0.d0)
      beta   = pihalf - dacos(cobeta)
      aa = Eta(1) + Eta(2)
      if (aa.lt..001d0) aa=.001d0
      aa2 = aa*aa
      theta = pihalf*Eta(1)/aa
      dthdeta(1) =  pihalf*Eta(2)/aa2
      dthdeta(2) = -pihalf*Eta(1)/aa2
      psi   = (beta-pihalf)*(Eta(1)+Eta(2))+pihalf
      dpsdeta(1) = beta - pihalf
      dpsdeta(2) = beta - pihalf
c
      costh = cos(theta)
      sinth = sin(theta)
      cosps = cos(psi)
      sinps = sin(psi)
      xprim(1) = rads*cosps*sinth
      xprim(2) = rads*cosps*costh
      xprim(3) = rads*sinps
      do j=1,2
        dxprdeta(1,j) =  rads*cosps*costh*dthdeta(j)
     .                  -rads*sinps*sinth*dpsdeta(j)
        dxprdeta(2,j) = -rads*cosps*sinth*dthdeta(j)
     .                  -rads*sinps*costh*dpsdeta(j)
        dxprdeta(3,j) =  rads*cosps*dpsdeta(j)
      enddo
c
c  ...evaluate physical coordinates and their derivatives wrt the
c     parameters
      do i=1,3
        s  = 0.d0
        s1 = 0.d0
        s2 = 0.d0
        do j=1,3
          s  = s  + aij(i,j)*xprim(j)
          s1 = s1 + aij(i,j)*dxprdeta(j,1)
          s2 = s2 + aij(i,j)*dxprdeta(j,2)
        enddo
        X(i) = s + centers(i,1)
        Dxdeta(i,1) = s1
        Dxdeta(i,2) = s2
      enddo
c
c
      end



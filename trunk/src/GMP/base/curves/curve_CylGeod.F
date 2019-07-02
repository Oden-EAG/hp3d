c-----------------------------------------------------------------------
c
c   routine name       - curve_CylGeod
c
c-----------------------------------------------------------------------
c
c   latest revision    - Jan 06
c
c   purpose            - routine evaluates physical coordinates
c                        for a culinder-geodesics-curve parametrization,
c                        and their  derivative wrt to reference
c                        coordinate
c
c   arguments :
c     in:
c               No     - the curve number
c               Eta     - reference coordinate  (between 0 and 1)
c     out:
c               X      - physical coordinates of the point
c               Dxdeta  - derivatives of the physical coordinates wrt
c                        reference coordinate
c
c-----------------------------------------------------------------------
c
      subroutine curve_CylGeod(No,Eta, X,Dxdeta)
c
      use GMP
#include "syscom.blk"
      common /cCylGeod/ iprint
c
      dimension X(3),Dxdeta(3)
c
c  ...cylinder origin, transformation matrix, local coordinates
c     and corresponding derivatives
      dimension xcenter(3),transf(3,3),xp(3),dxpdeta(3)
c
ccc      iprint=0
      if (iprint.eq.1) then
        write(*,7000) No,Eta
 7000   format('curve_CylGeod: DEBUGGING FOR No,Eta = ',i4,2x,f8.3)
      endif
c
c  ...check the type
      if (CURVES(No)%Type.ne.'CylGeod') then
        write(*,7001) No, CURVES(No)%Type
 7001   format('curve_CylGeod: WRONG CURVE TYPE, No, Type = ',
     .         i4,a10)
        stop 1
      endif
c
c  ...get the curve data
      theta2 = CURVES(No)%Rdata(1)
      rz1    = CURVES(No)%Rdata(2)
      rz2    = CURVES(No)%Rdata(3)
      l=3
      do i=1,3
        do j=1,3
          l=l+1
          transf(j,i) = CURVES(No)%Rdata(l)
        enddo
      enddo
c
c  ...get the cylinder origin
      ns = CURVES(No)%Idata(1)
      if (SURFACES(ns)%Type.ne.'Cylinder') then
        write(*,7002) ns,SURFACES(ns)%Type
 7002   format('curve_CylGeod: WRONG SURFACE TYPE, No, Type = ',
     .         i4,a10)
        stop 1
      endif
      xcenter(1:3) = SURFACES(ns)%Rdata(1:3)
      radius = SURFACES(ns)%Rdata(7)
c
c  ...compute the physical coordinates and derivatives
      theta = theta2*Eta
      xp(1) = radius*cos(theta)
      xp(2) = radius*sin(theta)
      xp(3) = rz1 + (rz2-rz1)*Eta
      dxpdeta(1) = - xp(2)*theta2
      dxpdeta(2) =   xp(1)*theta2
      dxpdeta(3) =   rz2-rz1
c
c  ...transform to the global system of coordinates
      call transform(21,transf,X,xp)
      X = X + xcenter
      call transform(21,transf,Dxdeta,dxpdeta)
c
      if (iprint.eq.1) then
        write(*,7003) No,X,Dxdeta
 7003   format('curve_CylGeod: No,X,Dxdeta = ', i4,2x,2(3e12.5,2x))
        call pause
      endif
c
c
      return
      end

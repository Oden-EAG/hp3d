c-----------------------------------------------------------------------
c
c   latest revision    - Jan 08
c
c   purpose            - routine evaluates physical coordinates
c                        and its derivatives wrt to reference
c                        coordinates for a point on a curve resulting
c                        form intersection of two surfaces
c
c   arguments :
c     in:
c              No      - curve number
c              I1,I2   - indices of surface numbers
c              Eta     - reference coordinate
c     out:
c              X       - physical coordinates of the point
c              Dxdeta  - derivatives of the physical coordinates wrt
c                        to the parameters
c              Ierror  = 0 if no error occurred
c                      = 1 otherwise
c
c-----------------------------------------------------------------------
c
      subroutine curve_2SurfsCur(No,I1,I2,Eta, X,Dxdeta,Ierror)
c
      use control
      use GMP
#include "syscom.blk"
      common /cmnewt/ iprint_mnewt
      common /ccone/ iprint_cone
c
      dimension X(3),Dxdeta(3)
c
      dimension nsurf(4)
      dimension xv(3,2),dxv(3),fgrad(3),f3grad(3),f4grad(3),
     .          xs(3),sfact(4),void(4),aij(3,3),aux(3)
c
c  ...local variables:
c     nsurf     - surfaces defining the curve
c     xv        - endpoints coordinates
c     dxv       - vector connecting the endpoints
c     fgrad,f3grad,f4grad - surface gradients
c     xs        - a starting point for NR iterations
c     sfact     - renormalization factor to adjust the definition
c                 of surfaces cutting of the segment
c     aij,aux   - matrices to determine derivatives
c
c
ccc      if (No.eq.39) then
ccc        iprint=1
ccc      else
        iprint=0
ccc      endif
      Ierror=0
 5    continue
      iprint_mnewt=iprint
      iprint_cone=iprint
      if (iprint.eq.1) then
        write(*,7001) No
 7001   format('curve_2SurfsCur: No = ',i5)
      endif
c
c  ...get endpoints for the curve
      do i=1,2
        np = CURVES(No)%EndPoNo(i)
        call pointr(np, xv(1:3,i))
      enddo
      dxv(1:3) = xv(1:3,2) - xv(1:3,1)
c
c  ...determine the linear interpolant of vertex values used
c     as a starting point for NR iterations, and the mid-point
c     used to adjust surface orientations
      do j=1,NDIM
        xs(j) = Eta*xv(j,2) + (1.d0-Eta)*xv(j,1)
      enddo
c
c  ...get surface numbers
      ns1 = CURVES(No)%Idata(I1)
      ns2 = CURVES(No)%Idata(I2)
      nsurf(1)=ns1; nsurf(2)=ns2
c
c  ...create two temporary planes
      do i=1,2
        NRSURFS=NRSURFS+1
        SURFACES(NRSURFS)%Type = 'VecPt'
        allocate(SURFACES(NRSURFS)%Rdata(6))
        SURFACES(NRSURFS)%Rdata(1:3) = xv(1:3,i)
        call surf(ns1,xv(1:3,i), fval,f3grad)
        call surf(ns2,xv(1:3,i), fval,f4grad)
        call cross_product(f3grad,f4grad, aux)
        call norm(aux,s)
        if (s.lt.GEOM_TOL) then
          write(*,*) 'curve_2SurfsCur: SINGULAR POINT'
          stop 1
        endif
        SURFACES(NRSURFS)%Rdata(4:6) = aux(1:3)/s
        nsurf(2+i)=NRSURFS
      enddo
c
      if (iprint.eq.1) then
        write(*,7002) xv, xs, Eta,nsurf
 7002   format('curve_2SurfsCur: ENDPOINS        = ',2(3e12.5,2x),
     .       /,'                 STARTING POINT  = ',3e12.5,
     .       /,'                 Eta             = ',f8.3,
     .       /,'                 SURFACE NUMBERS = ',4i5)
      endif
c
c  ...check consistency of data and determine the renormalization
c     factors
      do i=1,2
        call surf(nsurf(2+i),xv(1:3,i), fval,fgrad)
        s1 = dot_product(dxv,fgrad)
        s2 = dot_product(fgrad,fgrad)
        s2 = 1.d0/sqrt(s2)
c
c  .....after the renormalization the gradient of both surface
c       at the endpoints is a unit vector with direction
c       compatible with the direction of the line segment
c
c       sign(x,y) = abs(x) , y >= 0 ; -abs(x) , y < 0
c
        sfact(i) = sign(s2,s1)
        if (fval.gt.GEOM_TOL) then
          write(*,7003) No
 7003     format('curve_2SurfsCur: 2 SURFACE CURVE = ',i4)
          write(*,7004) i,CURVES(No)%EndPoNo(i),2+i,nsurf(2+i),
     ,                  fval
 7004     format('curve_2SurfsCur: INCOMPATIBILITY ! POINT ',i2,'(',i4,
     .           ')',' DOES NOT LIE ON SURFACE ',i2,'(',i3,')',
     .           ' VALUE = ',e12.5)
          call pause
        endif
      enddo
c
c  ...use NR method to determine the point
      call mnewt(2,nsurf,Eta,void,xs,sfact, X)
c
c  ...compute derivatives
      do i=1,2
        call surf(nsurf(i),X, fval,fgrad)
        if (fval.gt.GEOM_TOL) then
          write(*,7005) i,fval
 7005     format('curve_2SurfsCur: INCONSISTENCY i,fval = ',i2,2x,e12.5)
          stop 1
        endif
        aij(i,1:3) = fgrad(1:3)
        aux(i) = 0.d0
      enddo
c
      call surf(nsurf(3),X, fval3,f3grad)
      call surf(nsurf(4),X, fval4,f4grad)
      fval3 = fval3*sfact(1); f3grad(1:3) = f3grad(1:3)*sfact(1)
      fval4 = fval4*sfact(2); f4grad(1:3) = f4grad(1:3)*sfact(2)
c
      aij(3,1:3) = (1.d0-Eta)*f3grad(1:3)+Eta*f4grad(1:3)
      aux(3) = fval3 - fval4
c
c  ...solve for derivatives
      call saruss(aij,aux, Dxdeta,ifl)
      if (ifl.ne.0) then
        iprint=1
        go to 5
      endif
      if (iprint.eq.1) then
        write(*,7006) X(1:3)
 7006   format('curve_2SurfsCur: X      = ',3e12.5)
        write(*,7007) Dxdeta(1:NDIM)
 7007   format('                 Dxdeta = ',3e12.5)
        call pause
        call print_GMP
      endif
c
c  ...delete the two temporary planes
      do i=1,2
        deallocate(SURFACES(NRSURFS)%Rdata)
        NRSURFS=NRSURFS-1
      enddo
c
c
      end subroutine curve_2SurfsCur

c----------------------------------------------------------------------
c
c   routine name       - curve_ImpCir
c
c----------------------------------------------------------------------
c
c   latest revision    - Jan 06
c
c   purpose            - implicit curve parametrization routine
c
c   arguments :
c     in:
c               No     - the curve number
c               Eta    - reference coordinate  (between 0 and 1)
c     out:
c               X      - physical coordinates of the point
c               Dxdeta - derivatives of the physical coordinates wrt
c                        reference coordinate
c
c----------------------------------------------------------------------
c
      subroutine curve_ImpCir(No,Eta, X,Dxdeta)
c
      use GMP
      use control
#include "syscom.blk"
c
      dimension X(3),Dxdeta(3)
c
      dimension nsurf(6)
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
      iprint=0
 5    continue
      nsurf(1:6)=0 ; void(1:4)=0.d0 ; sfact(1:4)=0.d0
c
      if ((CURVES(No)%Type.ne.'ImpCir').or.(NDIM.ne.3)) then
        write(*,7001)
 7001   format('curve_ImpCir: WRONG CALL')
        stop 1
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
      nsurf(1:4) = CURVES(No)%Idata(1:4)
      if (iprint.eq.1) then
        write(*,7002) xv, xs, Eta,nsurf
 7002   format('curve_ImpCir: ENDPOINS        = ',2(3e12.5,2x),
     .       /,'              STARTING POINT  = ',3e12.5,
     .       /,'              Eta             = ',f8.3,
     .       /,'              SURFACE NUMBERS = ',4i5)
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
        sfact(i) = sign(s2,s1)
        if (fval.gt.GEOM_TOL) then
          write(*,7003) No
 7003     format('curve_ImpCir: IMPLICIT CURVE = ',i4)
          write(*,7004) i,CURVES(No)%EndPoNo(i),2+i,nsurf(2+i),
     ,                  fval
 7004     format('curve_ImpCir:: INCOMPATIBILITY ! POINT ',i2,'(',i4,
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
 7005     format('curve_IMPCiR: INCONSISTENCY i,fval = ',i2,2x,e12.5)
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
 7006   format('curve_IMPCir: X      = ',3e12.5)
        write(*,7007) Dxdeta(1:NDIM)
 7007   format('              Dxdeta = ',3e12.5)
        call pause
      endif
c
c
      end

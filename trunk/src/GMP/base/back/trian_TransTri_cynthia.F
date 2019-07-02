c----------------------------------------------------------------------
c
c   routine name       - trian_TransTri
c
c----------------------------------------------------------------------
c
c   latest revision    - Jan 06
c
c   purpose            - routine evaluates physical coordinates
c                        and its derivatives wrt to reference
c                        triangle parameters for a point in
c                        a triangle using transfinite interpolation
c
c   arguments :
c     in:
c               No     - the triangle number
c               Eta    - coordinates of a point in the reference
c                        triangle
c     out:
c               X      - physical coordinates of the point
c               Dxdeta - derivatives of the physical coordinates wrt
c                        to the reference triangle coordinates
c
c----------------------------------------------------------------------
c
      subroutine trian_TransTri(No,Eta, X,Dxdeta)
c
      use GMP
#include "syscom.blk"
c
      dimension Eta(2),X(NDIM),Dxdeta(NDIM,3)
c
      dimension noc(3),norientc(3),xv(NDIM,3),xc(NDIM),dxc(NDIM),
     .          rlambda(3),drlambda(2,3)
c
c  ...explanation of local variables
c     noc      - edge curve numbers
c     norientc - edge orientations
c     xv       - vertex coordinates
c     xc       - coordinates of a point on a curve
c     dxc      - derivatives of coordinates of a curve point wrt the
c                curve coordinate
c     rlambda,drlambda - affine coordinates and their derivatives
c
c----------------------------------------------------------------------
c
      data eps /1e-13/
c
c  ...set the derivatives of the affine coordinates
      data drlambda /-1.d0,-1.d0,1.d0,0.d0,0.d0,1.d0/
c
      imod(j,mod) = j-(j-1)/mod*mod
c
      iprint=0
c
      if (TRIANGLES(No)%Type.ne.'TransTri') then
        write(*,7000) TRIANGLES(No)%Type
 7000   format('tria_TransTri: WRONG TRIANGLE TYPE = ',a10)
        stop 1
      endif
c
c  ...get the edges and vertex coordinates
      norientc(1:3)=0
      do i=1,3
        noc(i) = iabs(TRIANGLES(No)%EdgeNo(i))
        if (TRIANGLES(No)%EdgeNo(i).lt.0) norientc(i)=1
        np = TRIANGLES(No)%VertNo(i)
        call pointr(np, xv(1:NDIM,i))
      enddo
c
      if (iprint.eq.1 )then
        write(*,7001) Eta
 7001   format('trian_TransTri:  Eta = ',3e12.5)
        write(*,7004)
 7004   format('trian_TransTri: vertex coordinates = ')
        write(*,7002) ((xv(j,i),j=1,NDIM),i=1,3)
 7002   format(3(2e12.5,2x))
      endif
c
c  ...calculate the affine coordinates
      rlambda(1) = 1.d0 - Eta(1) - Eta(2)
      rlambda(2) = Eta(1)
      rlambda(3) = Eta(2)
c
c  ...initiate
      X(1:NDIM) = 0.d0
      Dxdeta(1:NDIM,1:2) = 0.d0
c
c  ...calculate the linear, vertex contribution
      do i=1,3
        do j=1,NDIM
c
c  .......coordinates
          X(j) = X(j) + xv(j,i)*rlambda(i)
c
c  .......derivatives
          Dxdeta(j,1) = Dxdeta(j,1) + xv(j,i)*drlambda(1,i)
          Dxdeta(j,2) = Dxdeta(j,2) + xv(j,i)*drlambda(2,i)
        enddo
      enddo
      if (iprint.eq.1) then
        write(*,*) 'trian_TransTri: AFTER LINEAR INTERPOLATION'
        do j=1,NDIM
          write(*,7003) j,X(j),(Dxdeta(j,i),i=1,2)
 7003     format('j,X(j),Dxdeta(j,*) = ',i2,e12.5,2x,2e12.5)
        enddo
        call pause
      endif
c
c  ...exit if the point is in the singular region
c
c  ...loop through the edges
      do i=1,3
        i1 = imod(i+1,3)
c
c  .....exit if in a singular region
        if ((rlambda(i).lt.eps).or.(rlambda(i1).lt.eps)) cycle
c
c  .....first singular contribution
        xic = rlambda(i1)
        call curve_local(noc(i),norientc(i),xic, xc,dxc)
        xc(1:NDIM) = xc(1:NDIM)
     .             - (xv(1:NDIM,i)*(1.d0-xic) + xv(1:NDIM,i1)*xic)
        dxc(1:NDIM) = dxc(1:NDIM)
     .              - (xv(1:NDIM,i1) - xv(1:NDIM,i))
        xc(1:NDIM) = 0.5d0*xc(1:NDIM)
        dxc(1:NDIM) = 0.5d0*dxc(1:NDIM)
c
c  .....accumulate
        X(1:NDIM) = X(1:NDIM)
     .            + xc(1:NDIM)*rlambda(i)/(1.d0-rlambda(i1))
        do j=1,2
          Dxdeta(1:NDIM,j) = Dxdeta(1:NDIM,j)
     .    + dxc(1:NDIM)*drlambda(j,i1)*rlambda(i)/(1.d0-rlambda(i1))
     .    + xc(1:NDIM)*(drlambda(j,i)*(1.d0-rlambda(i1))
     .      + rlambda(i)*drlambda(j,i1))/(1.d0-rlambda(i1))**2
        enddo
c
c  .....second singular contribution
        xic = 1.d0 - rlambda(i)
        call curve_local(noc(i),norientc(i),xic, xc,dxc)
        xc(1:NDIM) = xc(1:NDIM)
     .             - (xv(1:NDIM,i)*(1.d0-xic) + xv(1:NDIM,i1)*xic)
        dxc(1:NDIM) = dxc(1:NDIM)
     .              - (xv(1:NDIM,i1) - xv(1:NDIM,i))
        xc(1:NDIM) = 0.5d0*xc(1:NDIM)
        dxc(1:NDIM) = 0.5d0*dxc(1:NDIM)
c
c  .....accumulate
        X(1:NDIM) = X(1:NDIM)
     .            + xc(1:NDIM)*rlambda(i1)/(1.d0-rlambda(i))
        do j=1,2
          Dxdeta(1:NDIM,j) = Dxdeta(1:NDIM,j)
     .    - dxc(1:NDIM)*drlambda(j,i)*rlambda(i1)/(1.d0-rlambda(i))
     .    + xc(1:NDIM)*(drlambda(j,i1)*(1.d0-rlambda(i))
     .      + rlambda(i1)*drlambda(j,i))/(1.d0-rlambda(i))**2
        enddo
        if (iprint.eq.1) then
          write(*,7005) i
 7005     format('trian_TransTri: AFTER ',i1,'TH SIDE  INTERPOLATION')
          do j=1,NDIM
            write(*,7003) j,X(j),(Dxdeta(j,ii),ii=1,2)
          enddo
          call pause
        endif
c
c  ...end of loop through edges
      enddo
c
c
      return
      end
c

c----------------------------------------------------------------------
c
c   routine name       - curve_QuaCir
c
c----------------------------------------------------------------------
c
c   latest revision    - Jan 06
c
c   purpose            - routine defines a parametrization for
c                        a quarter of circle
c
c   arguments :
c     in:
c               No     - the curve number
c               Eta     - reference coordinate  (between 0 and 1)
c     out:
c               X      - physical coordinates of the point
c               Dxdeta - derivatives of the physical coordinates wrt
c                        reference coordinate
c
c----------------------------------------------------------------------
c
      subroutine curve_QuaCir(No,Eta, X,Dxdeta)
c
      use GMP
#include "syscom.blk"
c
      dimension X(3),Dxdeta(3)
c
c  ...coordinates of the endpoints
      dimension xv(3,2)
c
c  ...center coordinates
      dimension center(3)
c
c  ...transformation matrix
      dimension aij(3,3)
c
c  ...local coordinates and their derivatives
      dimension xprim(3),dxprdeta(3)
c
      iprint=0
c
      if (CURVES(No)%Type.ne.'QuaCir') then
        write(*,7001)
 7001   format(' curve_QuaCir: WRONG CALL')
        stop
      endif
c
      if (iprint.eq.1) then
        write(*,*)'No,Eta = ',No,Eta
ccc        write(*,7002) No, Eta
ccc 7002   format(' curve_QuaCir: No,Eta = ',i4,2x,f8.3)
        call pause
      endif
c
c  ...get the endpoints coordinates
      do i=1,2
        np=CURVES(No)%EndPoNo(i)
        call pointr(np, xv(1:3,i))
      enddo
c
c  ...get coordinates of the center
      center(1:3)=CURVES(No)%Rdata(1:3)
c
c  ...evaluate the radius
      rad = 0.d0
      do i=1,3
        rad = rad + (xv(i,1) - center(i))**2
      enddo
      rad = sqrt(rad)
c
c  ...evaluate the transformation matrix
      do i=1,3
        aij(i,1) = (xv(i,1) - center(i))/rad
        aij(i,2) = (xv(i,2) - center(i))/rad
      enddo
      call cross_product(aij(1:3,1),aij(1:3,2), aij(1:3,3))
c
c  ...printing
      if (iprint.eq.1) then
        do i=1,3
          write(*,7005) i,aij(i,1:3)
 7005     format(' i,a(i,:) = ',i1,2x,3(e12.5,2x))
        enddo
      endif
c
c  ...evaluate coordinates and their derivatives in the auxiliary
c     system of coordinates (see manual for explanation)
      pihalf = dacos(0.d0)
      alpha = pihalf*Eta
      xprim(1) = rad*dcos(alpha)
      xprim(2) = rad*dsin(alpha)
      xprim(3) = 0.d0
      dxprdeta(1) = -xprim(2)*pihalf
      dxprdeta(2) =  xprim(1)*pihalf
      dxprdeta(3) =  0.d0
c
c  ...printing
      if (iprint.eq.1) then
        write(*,7006) xprim(1:3)
 7006   format(' xprim    = ',3(e12.5,2x))
        write(*,7007) dxprdeta(1:3)
 7007   format(' dxprdeta = ',3(e12.5,2x))
      endif
c
c  ...evaluate physical coordinates and their derivatives wrt the
c     reference coordinate
      do i=1,3
        s =0.d0
        s1=0.d0
        do j=1,3
          s  = s  + aij(i,j) * xprim(j)
          s1 = s1 + aij(i,j) * dxprdeta(j)
        enddo
        X(i) = s + center(i)
        Dxdeta(i) = s1
      enddo
c
c  ...printing
      if (iprint.eq.1) then
        write(*,7003) X(1:3)
        write(*,7004) Dxdeta(1:3)
 7003   format(' X      = ',3(e12.5,2x))
 7004   format(' Dxdeta = ',3(e12.5,2x))
      endif
c
c
      end subroutine curve_QuaCir

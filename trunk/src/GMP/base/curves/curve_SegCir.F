c//////////////////////////////////////////////////////////////////////
c P. Gatto, Nov 12 : this is a LEGACY routine. Do not use it, unless
c                    you revise it!
c//////////////////////////////////////////////////////////////////////
c
c
c----------------------------------------------------------------------
c
c   routine name       - curve_SegCir
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Jan 06
c
c   purpose            - routine defines a parametrization for
c                        a segment of a circle
c
c   arguments :
c     in:
c               No     - the curve number
c               Eta    - reference coordinate  (between 0 and 1)
c     out:
c               X      - physical coordinates of the point
c               Dxdeta - derivatives of the physical coordinates wrt
c                        the reference coordinate
c
c   required  routines -
c
c----------------------------------------------------------------------
c
      subroutine curve_SegCir(No,Eta, X,Dxdeta)
c
      use GMP
      use control
#include "syscom.blk"
ccc      common /common_curve_SegCir/ iprint
c
      dimension X(NDIM),Dxdeta(NDIM)
c
c  ...coordinates of the endpoints
      dimension xv(NDIM,2)
c
c  ...center coordinates
      dimension center(NDIM)
c
c  ...vectors connecting the center with endpoints
      dimension radv1(3),radv2(3)
c
c  ...local coordinates and their derivatives
      dimension xprim(3),dxprdeta(3)
c
c  ...transformation matrix
      dimension aij(3,3)
c
      iprint=0
c
c      if (CURVES(No)%Type.ne.'SegCir') then
c        write(*,7001)
c 7001   format('curve_SegCir: WRONG CALL')
c        stop 1
c      endif
c
      if (iprint.eq.1) then
        write(*,7002) No,Eta
 7002   format('curve_SegCir: No,Eta = ',i4,2x,f8.3)
      endif
c
c  ...get the endpoints coordinates
      do i=1,2
        np = CURVES(No)%EndPoNo(i)
        call pointr(np, xv(1:NDIM,i))
      enddo
c
c  ...get coordinates of the center
      nsurface =  CURVES(No)%Idata(1)
      center(1:NDIM) = SURFACES(nsurface)%Rdata(1:NDIM)
      if (iprint.eq.1) then
        write(*,7006) xv(1:NDIM,1),xv(1:NDIM,2), center(1:NDIM)
 7006   format('curve_SegCir: xv,center = ',9f8.3)
      endif
c
c  ...evaluate radius and segment angle beta
      radv1(3) = 0.d0; radv2(3) = 0.d0
      rad1 = 0.d0; rad2 = 0.d0; skpr = 0.d0
      do i=1,NDIM
        radv1(i) = xv(i,1) - center(i)
        radv2(i) = xv(i,2) - center(i)
        rad1 = rad1 + (radv1(i))**2              !*** accumulate for radius to 1st vertex
        rad2 = rad2 + (radv2(i))**2              !*** accumulate for radius to 2nd vertex
        skpr = skpr + radv1(i)*radv2(i)          !*** accumulate for scalar product
      enddo
c
c  ...check compatibility
      rad1 = sqrt(rad1); rad2 = sqrt(rad2)
      if (abs(rad1 - rad2) .gt. GEOM_TOL) then
        write(*,7003) rad1,rad2
 7003   format('curve_SegCir: INCOMPATIBLE DATA, rad1,rad2 = ',2e12.5)
        stop 1
      endif
      radv1(1:NDIM) = radv1(1:NDIM)/rad1         !*** normalize 1st vector
      radv2(1:NDIM) = radv2(1:NDIM)/rad2         !*** normalize 2nd vector
      cobeta = skpr/(rad1*rad2)
      beta = dacos(cobeta)                       !*** get cosine of angle
      if (iprint.eq.1) then
        write(*,*) 'curve_SegCir: beta = ',beta
      endif
c
c  ...evaluate the transformation matrix
      call cross_product(radv1,radv2, aij(1:3,3))
      call norm(aij(1:3,3),rn)
      aij(1:3,3) = aij(1:3,3)/rn
      aij(1:3,1) = radv1(1:3)
      call cross_product(aij(1:3,3),aij(1:3,1), aij(1:3,2))
c
c  ...evaluate coordinates and their derivatives in the local
c     system of coordinates (see manual for explanation)
      alpha = beta*Eta
      xprim(1) = rad1*cos(alpha)
      xprim(2) = rad1*sin(alpha)
      xprim(3) = 0.d0
      dxprdeta(1) = -xprim(2)*beta
      dxprdeta(2) =  xprim(1)*beta
      dxprdeta(3) =  0.d0
      if (iprint.eq.1) then
        write(*,7007) xprim(1:NDIM), dxprdeta(1:NDIM)
 7007   format('curve_SegCir: xprim,dxprdeta = ',6f8.3)
      endif
c
c  ...transform to the global coordinates
      do i=1,NDIM
        s =0.d0
        s1=0.d0
        do j=1,NDIM
          s  = s  + aij(i,j)*xprim(j)
          s1 = s1 + aij(i,j)*dxprdeta(j)
        enddo
        X(i) = s + center(i)
        Dxdeta(i) = s1
      enddo
      if (iprint.eq.1) then
        write(*,7004) X(1:NDIM)
 7004   format('curve_SegCir: X = ',3e12.5)
        write(*,7005) Dxdeta(1:NDIM)
 7005   format('curve_SegCir: Dxdeta = ',3e12.5)
        call pause
      endif
c
c
      end subroutine curve_SegCir



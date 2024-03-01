!----------------------------------------------------------------------
!
!   routine name       - curve_QuaCir
!
!----------------------------------------------------------------------
!
!   latest revision    - Mar 2023
!
!   purpose            - routine defines a parametrization for
!                        a quarter of circle
!
!   arguments :
!     in:
!               No     - the curve number
!               Eta    - reference coordinate (between 0 and 1)
!     out:
!               X      - physical coordinates of the point
!               Dxdeta - derivatives of the physical coordinates wrt
!                        reference coordinate
!
!----------------------------------------------------------------------
!
      subroutine curve_QuaCir(No,Eta, X,Dxdeta)
!
      use GMP
!
      implicit none
!
      integer :: No
      real(8) :: Eta,X(3),Dxdeta(3)
!
!  ...coordinates of the endpoints
      real(8) :: xv(3,2)
!
!  ...center coordinates
      real(8) :: center(3)
!
!  ...transformation matrix
      real(8) :: aij(3,3)
!
!  ...local coordinates and their derivatives
      real(8) :: xprim(3),dxprdeta(3)
!
      integer :: i,j,np
      real(8) :: alpha,pihalf,rad,s,s1
!
      integer :: iprint
      iprint=0
!
      if (CURVES(No)%Type.ne.'QuaCir') then
        write(*,7001)
 7001   format(' curve_QuaCir: WRONG CALL')
        stop
      endif
!
      if (iprint.eq.1) then
        write(*,*)'No,Eta = ',No,Eta
!!!        write(*,7002) No, Eta
!!! 7002   format(' curve_QuaCir: No,Eta = ',i4,2x,f8.3)
        call pause
      endif
!
!  ...get the endpoints coordinates
      do i=1,2
        np=CURVES(No)%EndPoNo(i)
        call pointr(np, xv(1:3,i))
      enddo
!
!  ...get coordinates of the center
      center(1:3)=CURVES(No)%Rdata(1:3)
!
!  ...evaluate the radius
      rad = 0.d0
      do i=1,3
        rad = rad + (xv(i,1) - center(i))**2
      enddo
      rad = sqrt(rad)
!
!  ...evaluate the transformation matrix
      do i=1,3
        aij(i,1) = (xv(i,1) - center(i))/rad
        aij(i,2) = (xv(i,2) - center(i))/rad
      enddo
      call cross_product(aij(1:3,1),aij(1:3,2), aij(1:3,3))
!
!  ...printing
      if (iprint.eq.1) then
        do i=1,3
          write(*,7005) i,aij(i,1:3)
 7005     format(' i,a(i,:) = ',i1,2x,3(e12.5,2x))
        enddo
      endif
!
!  ...evaluate coordinates and their derivatives in the auxiliary
!     system of coordinates (see manual for explanation)
      pihalf = dacos(0.d0)
      alpha = pihalf*Eta
      xprim(1) = rad*dcos(alpha)
      xprim(2) = rad*dsin(alpha)
      xprim(3) = 0.d0
      dxprdeta(1) = -xprim(2)*pihalf
      dxprdeta(2) =  xprim(1)*pihalf
      dxprdeta(3) =  0.d0
!
!  ...printing
      if (iprint.eq.1) then
        write(*,7006) xprim(1:3)
 7006   format(' xprim    = ',3(e12.5,2x))
        write(*,7007) dxprdeta(1:3)
 7007   format(' dxprdeta = ',3(e12.5,2x))
      endif
!
!  ...evaluate physical coordinates and their derivatives wrt the
!     reference coordinate
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
!
!  ...printing
      if (iprint.eq.1) then
        write(*,7003) X(1:3)
        write(*,7004) Dxdeta(1:3)
 7003   format(' X      = ',3(e12.5,2x))
 7004   format(' Dxdeta = ',3(e12.5,2x))
      endif
!
!
      end subroutine curve_QuaCir

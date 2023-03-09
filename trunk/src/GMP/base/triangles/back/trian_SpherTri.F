c-----------------------------------------------------------------------
c
c   routine name       - trian_SpherTri
c
c-----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - routine evaluates physical coordinates
c                        and its derivatives wrt to reference
c                        coordinates for a point on a spherical
c                        octant (old technology...)
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
c-----------------------------------------------------------------------
c
      subroutine trian_SpherTri(No,Eta, X,Dxdeta)
c
      use GMP
c
      implicit none
c
      integer :: No
      real(8) :: Eta(2),X(NDIM),Dxdeta(NDIM,2)
c
c  ...vertex coordinates and edge curves
      real(8) :: xv(3,3)
      integer :: ncurv(3)
c
c  ...transformation matrix
      real(8) ::  aij(3,3)
c
c  ...auxiliary quantities
      real(8) :: xprim(3),dxprdxi(3,2),dthdxi(2),dpsdxi(2),
     .           aux(3,3),center(3)
c
      real(8) :: aa,aa2,cosps,costh,psi,rad,s,s1,s2,sinps,sinth,theta
      integer :: i,j,np
c
      real(8) :: pihalf
c
      integer :: iprint=0
c
      if ((TRIANGLES(No)%Type.ne.'SpherTri'.or.(NDIM.ne.3))) then
        write(*,7001) TRIANGLES(No)%Type
 7001   format('tria_SpherTri: WRONG TRIANGLE TYPE = ',a10)
        stop 1
      endif
c
      if (iprint.eq.1) then
        write(*,7002) No,Eta
 7002   format('trian_SpherTri: No,Eta = ',i4,2x,2f8.3)
      endif
c
c  ...get the edges and vertex coordinates
      do i=1,3
        ncurv(i) = TRIANGLES(No)%EdgeNo(i)
        np = TRIANGLES(No)%VertNo(i)
        call pointr(np, xv(1:NDIM,i))
      enddo
      if (iprint.eq.1) then
        write(*,7003)
 7003   format('trian_SpherTri: VERTEX COORDINATES = ')
        do j=1,NDIM
          write(*,7004) xv(j,1:3)
 7004     format(3f8.3)
        enddo
      endif
c
c  ...evaluate radius and center coordinates
      rad=0.d0
      do i=1,3
        aux(i,1) = xv(i,2) - xv(i,1)
        aux(i,2) = xv(i,3) - xv(i,1)
        rad = rad + (aux(i,1))**2
      enddo
      rad = sqrt(rad/2.d0)
      call cross_product(aux(1:3,1),aux(1:3,2), aux(1:3,3))
      s = sqrt(aux(1,3)**2+aux(2,3)**2+aux(3,3)**2)
      do i=1,3
        center(i) = (xv(i,1)+xv(i,2)+xv(i,3))/3.d0
     .            - aux(i,3)/s*rad/sqrt(3.d0)
      enddo
c
      if (iprint.eq.1) then
        write(*,7005) center,rad
 7005   format('trian_SpherTri: center,rad = ',3f8.3,2x,f8.3)
      endif
c
c  ...evaluate the transformation matrix
      do i=1,3
        aij(i,1) = (xv(i,3) - center(i))/rad
        aij(i,2) = (xv(i,2) - center(i))/rad
        aij(i,3) = (center(i) - xv(i,1))/rad
      enddo
c
c  ...evaluate coordinates and their derivatives in the auxiliary
c     system of coordinates (see manual for explanation)
      pihalf = dacos(0.d0)
      aa = Eta(1) + Eta(2)
      if (aa.lt..001d0) aa=.001d0
      aa2 = aa*aa
      theta = pihalf*Eta(1)/aa
      dthdxi(1) =  pihalf*Eta(2)/aa2
      dthdxi(2) = -pihalf*Eta(1)/aa2
      psi   = pihalf*(Eta(1)+Eta(2)-1.d0)
      dpsdxi(1) = pihalf
      dpsdxi(2) = pihalf
c
      costh = cos(theta)
      sinth = sin(theta)
      cosps = cos(psi)
      sinps = sin(psi)
      xprim(1) = rad*cosps*costh
      xprim(2) = rad*cosps*sinth
      xprim(3) = rad*sinps
      do j=1,2
        dxprdxi(1,j) = -rad*cosps*sinth*dthdxi(j)
     .                 -rad*sinps*costh*dpsdxi(j)
        dxprdxi(2,j) =  rad*cosps*costh*dthdxi(j)
     .                 -rad*sinps*sinth*dpsdxi(j)
        dxprdxi(3,j) =  rad*cosps*dpsdxi(j)
      enddo
c
c  ...evaluate physical coordinates and their derivatives wrt
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
c
      end subroutine trian_SpherTri

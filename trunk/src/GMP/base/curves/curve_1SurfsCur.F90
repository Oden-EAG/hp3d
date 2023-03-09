!----------------------------------------------------------------------------
!> Purpose : parametrization of a curve lying on a single algebraic surface
!!
!! @param[in ] No     - the curve number
!! @param[in ] Eta    - reference coordinate  (between 0 and 1)
!! @param[out] X      - physical coordinates of the point
!! @param[out] Dxdeta - derivatives of the physical coordinates
!!
!! @revision Nov 12
!----------------------------------------------------------------------------
subroutine curve_1SurfsCur(Nc,Eta, R,Dr)
!
      use GMP
      implicit none
      integer             ,intent(in ) :: Nc
      real(8)             ,intent(in ) :: Eta
      real(8),dimension(3),intent(out) :: R
      real(8),dimension(3),intent(out) :: Dr
      !
      real(8),dimension(3,2) :: v
      integer                :: nv,ns,i
!----------------------------------------------------------------------------
!
      ns=CURVES(nc)%Idata(1)
!
      select case (SURFACES(ns)%Type)
!
      case('VecPt','ThreePt')
        do i=1,2
          nv=CURVES(nc)%EndPoNo(i) ; v(1:3,i)=POINTS(nv)%Rdata(1:3)
        enddo
        do i=1,3
          R(i) = (1.d0 - Eta)*v(i,1) + Eta*v(i,2)
        enddo
        Dr(1:3) = v(1:3,2) - v(1:3,1)
      case('PPwCC')    ; call diag_segment(nc, eta, r, Dr)
      case('Sphere')   ; call circular_segment(nc, eta, r, Dr)
      case('Cylinder') ; call cylinder_geodesic(nc, eta, r, Dr)
      case('Cone')     ; call cone_geodesic(nc, eta, r, Dr)
!
      endselect
!
!
end subroutine curve_1SurfsCur
!
!
!
!----------------------------------------------------------------------
!> Purpose : parametrization for a curve lying on a yz plane
!            parametrized with cylindrical coordinates
!!
!! @param[in]  No     - curve number
!! @param[in]  Eta    - reference coordinates of a point
!! @param[out] X      - physical coordinates of the point
!! @param[out] Dxdeta - derivatives wrt reference coordinates
!!
!! @date Mar 2023
!----------------------------------------------------------------------
subroutine diag_segment(Nc,Eta, X,dX_dEta)
!
      use control
      use GMP
!
      implicit none
!
      integer :: Nc
      real(8) :: Eta,X(3),dX_dEta(3)
!
!  ...x coordinate, min and max radius
!      x0,rmin,rmax
!
!  ...endpoints, point in the parametric domain,
      real(8) :: eta_v(2,2),etap(2)
!
      integer :: i,iv,ivar,ns,nv
      real(8) :: dr,dr_deta,dtheta_deta
      real(8) :: r,r2,rmax,rmin,theta,theta1,x0
!
      real(8) :: pi
      integer :: iprint=0
!----------------------------------------------------------------------
!
      pi = acos(-1.d0)
!
      if (iprint.eq.1) then
        write(*,7000)Nc,Eta
 7000   format(' diag_segment: Nc,Eta = ',i5,2x,e12.5)
      endif
!
!  ...check surface
      ns = CURVES(Nc)%Idata(1)
      if (SURFACES(ns)%Type.ne.'PPwCC') then
        write(*,*) '---------------------------------'
        write(*,*) 'diag_segment:'
        write(*,*) 'ERROR: inconsistent surface type.'
        write(*,*) '---------------------------------'
        stop
      end if
!
!  ...get the data for the plane
      x0       = SURFACES(ns)%Rdata(1)
      rmin     = SURFACES(ns)%Rdata(2)
      rmax     = SURFACES(ns)%Rdata(3)
      dr = rmax - rmin
      if (iprint.eq.1) then
        write(*,7001) x0,rmin,rmax,dr
 7001   format('diag_segment: x0,rmin,rmax,dr = ',4f8.3)
      endif
!
!  ...check radii consistency
      if ((rmin.le.GEOM_TOL).or.(rmax-rmin.le.GEOM_TOL)) then
        write(*,*) 'diag_segment: rmin,rmax = ',rmin,rmax
        stop 1
      endif
!
!  ...determine point coordinates in the parametric plane
      do iv=1,2
        nv = CURVES(Nc)%EndPoNo(iv)
        if (abs(POINTS(nv)%Rdata(1)-x0).gt.GEOM_TOL) then
          write(*,*) 'diag_segment: Nc,iv,nv = ',Nc,iv,nv
          stop 1
        endif
        r2 = POINTS(nv)%Rdata(2)**2 + POINTS(nv)%Rdata(3)**2
        r = sqrt(r2)
        if ((r.lt.rmin).or.(r.gt.rmax)) then
          write(*,*) 'diag_segment: rmin,rmax,r = ',rmin,rmax,r
          stop 1
        endif
        eta_v(2,iv) = (r - rmin)/dr
        select case(iv)
        case(1)
          theta = datan2(POINTS(nv)%Rdata(3),POINTS(nv)%Rdata(2))
        case(2)
          theta1 = datan2(POINTS(nv)%Rdata(3),POINTS(nv)%Rdata(2))
          if (theta1-theta.gt.pi) then
            theta1 = theta1 - 2.d0*pi
          elseif (theta-theta1.gt.pi) then
            theta1 = theta1 + 2.d0*pi
          endif
          theta = theta1
        end select
        eta_v(1,iv) = (pi/2.d0 - theta)/2.d0/pi
        if (iprint.eq.1) then
          write(*,7002) iv,eta_v(1:2,iv)
 7002     format('diag_segment: iv, eta(1:2,iv) = ',i2,2x,2f8.3)
        endif
      enddo
!
!  ...plane reference coordinates of the point
      etap(1:2) = eta_v(1:2,1)*(1.d0-Eta) + eta_v(1:2,2)*Eta
      r = rmin + etap(2)*dr
      theta = 2.d0*pi*etap(1)
!
!  ...point coordinates
      X(1) = x0
      X(2) = r*sin(theta)
      X(3) = r*cos(theta)
!
!  ...derivatives
      dX_dEta(1) = 0.d0
      dr_deta = (eta_v(2,2)-eta_v(2,1))*dr
      dtheta_deta = (eta_v(1,2)-eta_v(1,1))*2.d0*pi
      dX_dEta(2) = dr_deta*sin(theta) + r*cos(theta)*dtheta_deta
      dX_dEta(3) = dr_deta*cos(theta) - r*sin(theta)*dtheta_deta
!
      if (iprint.eq.1) then
        write(*,7020)
 7020   format(' diag_segment: X,dX_dEta = ')
        do ivar=1,3
          write(*,7030) ivar,X(ivar),dX_dEta(ivar)
 7030     format(i3,3x,f8.3,3x,f8.3)
        enddo
        call pause
      endif
!
!
end subroutine diag_segment
!
!!
!----------------------------------------------------------------------
!> Purpose : parametrization for a geodesics on a sphere
!!
!! @param[in]  No     - curve number
!! @param[in]  Eta    - reference coordinates of a point
!! @param[out] X      - physical coordinates of the point
!! @param[out] Dxdeta - derivatives wrt reference coordinates
!!
!! @date Mar 2023
!----------------------------------------------------------------------
subroutine circular_segment(Nc,Eta, X,dX_dEta)
!
      use control
      use GMP
!
      implicit none
!
      integer :: Nc
      real(8) :: Eta,X(3),dX_dEta(3)
      real(8) :: cen(3),r1(3),r2(3),ax(3),dX_dalpha(3)
!
      integer :: i,ns,nv1,nv2
      real(8) :: alpha,dalpha_dEta,rad,rnorm1,rnorm2,sp,theta
!
      integer :: iprint=0
!----------------------------------------------------------------------
      if (iprint.eq.1) then
        write(*,7000)Nc,Eta
 7000   format(' circular_segment: Nc,Eta = ',i5,2x,e12.5)
      endif
!
!  ...check surface
      ns = CURVES(Nc)%Idata(1)
      if (SURFACES(ns)%Type.ne.'Sphere') then
        write(*,*) '---------------------------------'
        write(*,*) 'circular_segment:'
        write(*,*) 'ERROR: inconsistent surface type.'
        write(*,*) '---------------------------------'
        stop
      end if
!
!  ...get sphere data
      cen(1:3) = SURFACES(ns)%Rdata(1:3)
      rad      = SURFACES(ns)%Rdata(4)
      nv1 = CURVES(Nc)%EndPoNo(1) ; nv2 = CURVES(Nc)%EndPoNo(2)
      r1(1:3) = POINTS(nv1)%Rdata(1:3) - cen(1:3)
      r2(1:3) = POINTS(nv2)%Rdata(1:3) - cen(1:3)
!
!  ...check radii consistency
      call norm(r1(1:3), rnorm1)
      if(abs(rnorm1-rad).gt.GEOM_TOL) then
        write(*,*) 'circular_segment: inconsistent radii'
        call pause
      endif
      call norm(r2(1:3), rnorm2)
      if(abs(rnorm2-rad).gt.GEOM_TOL) then
        write(*,*) 'circular_segment: inconsistent radii'
        call pause
      endif
!
!  ...set up rotation
      call cross_product(r1(1:3),r2(1:3), ax(1:3))
      call scalar_product(r1(1:3),r2(1:3), sp)
      theta = acos(sp/(rnorm1*rnorm2))
      alpha = Eta*theta  ; dalpha_dEta = theta
      if (iprint.eq.1) then
        write(*,7001)theta,alpha
 7001   format(' circular_segment: theta,alpha = ',2(e12.5,2x))
      endif
!
      X = POINTS(nv1)%Rdata(1:3)
      call rotation_about_axis(ax,cen,alpha, X,dX_dalpha)
      dX_dEta = dX_dalpha(1:3)*dalpha_dEta
      if(iprint.eq.1) then
        write(*,7002)X,dX_dEta
 7002   format(' circular_segment: X,dX_dEta = ',4(e12.5,2x))
      endif
!
!
end subroutine circular_segment
!
!
!
!----------------------------------------------------------------------
!> Purpose : parametrization for a geodesics on a cone
!!
!! @param[in]  No     - curve number
!! @param[in]  Eta    - reference coordinates of a point
!! @param[out] X      - physical coordinates of the point
!! @param[out] Dxdeta - derivatives wrt reference coordinates
!!
!! @date Mar 2023
!----------------------------------------------------------------------
subroutine cylinder_geodesic(No,Eta, X,Dxdeta)
!
      use control
      use GMP
      use element_data
!
      implicit none
!
      integer :: No
      real(8) :: Eta,X(3),Dxdeta(3)
!
!  ...cylinder center and direction vector
      real(8) :: center(3),direction(3)
!
!  ...relative position vectors for curve endpoints in the global
!     cartesian system
      real(8) :: xrelv(3,2)
!
!  ...endpoint coordinates in the cylinder Cartesian system and in the
!     cylinder reference coordinates
      real(8) :: xrelsv(3,2),xparv(2,2)
!
!  ...transformation matrix from global Cartesian to cylinder Cartesian
!     coordinates; rows of the transformation matrix are unit
!     vectors of the cylinder system computed in the global system
      real(8) :: transf(3,3)
!
!  ...point in the parameter space, derivatives
      real(8) :: xpar(2),dxpardeta(2)
!
!  ...local cartesian coordinates, derivatives
      real(8) :: xrels(3),dxrelsdeta(3)
!
!  ...work space
      real(8) :: void(3)
!
      integer :: i,iv,ivar,j,np,ns
      real(8) :: fval,rad,r,s
!
      integer :: iprint=0
!-----------------------------------------------------------------------
!
      if (iprint.eq.1) then
        write(*,*) '-----------------------------------'
        write(*,7001) No,Eta
 7001   format('cylinder_geodesic: No,Eta = ',i4,2x,e12.5)
      endif
!
      ns = CURVES(No)%Idata(1)
      if (SURFACES(ns)%Type.ne.'Cylinder') then
        write(*,*) 'cylinder_geodesic: INCONSISTENT SURFACE TYPE'
        stop
      endif
!
!  ...cylinder parameters in human-readable form
      center(   1:3)=SURFACES(ns)%Rdata(1:3)
      direction(1:3)=SURFACES(ns)%Rdata(4:6)
      rad           =SURFACES(ns)%Rdata(7)
!
!  ...shift the originally specified center to a new location such
!     that the new z coordinate of the first vertex will be zero
      call normalize(direction(1:3))
      np = CURVES(No)%EndPoNo(1)
      void(1:3) = POINTS(np)%Rdata(1:3) - center(1:3)
      call scalar_product(void,direction, s)
      center(1:3) = center(1:3) + s*direction(1:3)
!
!  ...compute the relative position vector for vertices
      do iv=1,2
        np = CURVES(No)%EndPoNo(iv)
!
!  .....check if on the surface
        call surf(ns,POINTS(np)%Rdata(1:3), fval,void)
        if (abs(fval).gt.GEOM_TOL) then
          write(*,7011) No,ns,iv,np,fval
 7011     format(' cylinder_geodesic: No,ns,iv,np,fval = ', &
                                      i5,2i3,i6,e12.5)
          call pause
          call print_GMP
        endif
        xrelv(1:3,iv) = POINTS(np)%Rdata(1:3) - center(1:3)
      enddo
!
!  ...compute the transformation matrix from the global Cartesian
!     to the cylinder Cartesian system
!
!  ...z axis coincides with the cylinder axis
      transf(3,1:3) = direction(1:3)
!
!  ...the cylinder x axis goes through the first vertex
      call norm(xrelv(1:3,1), s)
      transf(1,1:3) = xrelv(1:3,1)/s
!
!  ...y axis
      call cross_product(transf(3,1:3),transf(1,1:3), transf(2,1:3))
!
!  ...printing
      if (iprint.eq.1) then
        do i=1,3
          write(*,9999)i,transf(i,1:3)
9999      format(' i,transf(i,1:3) = ',i1,2x,3(e12.5,2x))
        enddo
      endif
!
!  ...compute endpoint coordinates in the cylinder Cartesian system
      do iv=1,2
        do i=1,3
          s = 0.d0
          do j=1,3
            s = s + transf(i,j)*xrelv(j,iv)
          enddo
          xrelsv(i,iv) = s
        enddo
      enddo
!
!  ...printing
      if (iprint.eq.1) then
        do i=1,2
          write(*,9998)i,xrelsv(1:3,i)
9998      format(' i,xrelsv(1:3,i) = ',i1,2x,3(e12.5,2x))
       enddo
      endif
!
!  ...compute the vertex coordinates in the cylinder parametric space
      do iv=1,2
        call cart_to_polar(xrelsv(1:2,iv), r,xparv(1,iv))
        if (abs(rad-r).gt.GEOM_TOL) then
          write(*,7012) iv,rad,r
 7012     format('cylinder_geodesic: iv,rad,r = ',i2,2e12.5)
          call pause
        endif
        xparv(2,iv) = xrelsv(3,iv)
      enddo
!
!  ...printing
      if (iprint.eq.1) then
        do iv=1,2
          write(*,7031) iv,xparv(1:2,iv)
 7031     format(' i,xparv(1:2,i)  = ',i1,2x,2(e12.5,2x))
        enddo
      endif
!
!  ...geodesics on the cylinder is simply a straight line in the
!     (\theta,z) parameter space
!cc      xpar(1:2) = xparv(1:2,1)*(1.d0-Eta) + xparv(1:2,2)*Eta
      xpar(     1:2) = xparv(1:2,2)*Eta
      dxpardeta(1:2) = xparv(1:2,2)
!
!  ...local cartesian system
      xrels(1) = rad*cos(xpar(1))
      dxrelsdeta(1) = -rad*sin(xpar(1))*dxpardeta(1)
      xrels(2) = rad*sin(xpar(1))
      dxrelsdeta(2) =  rad*cos(xpar(1))*dxpardeta(1)
      xrels(3) = xpar(2)
      dxrelsdeta(3) = dxpardeta(2)
!
!  ...printing
      if (iprint.eq.1) then
        write(*,9997)xpar(1:2)
 9997   format(' xpar(1:2)  = ',2(e12.5,2x))
        write(*,9996)xrels(1:3)
 9996   format(' xrels(1:3) = ',3(e12.5,2x))
      endif
!
!  ...transform to the global Cartesian system
      X(1:3)=0.d0 ; Dxdeta(1:3)=0.d0
      do i=1,3
        do j=1,3
          X(i)      = X(i)      + transf(j,i)*xrels(j)
          Dxdeta(i) = Dxdeta(i) + transf(j,i)*dxrelsdeta(j)
        enddo
      enddo
      X(1:3)= X(1:3) + center(1:3)
!
      if (iprint.eq.1) then
        write(*,*) 'cylinder_geodesic: X,Dxdeta = '
        do ivar=1,3
          write(*,7035) X(ivar),Dxdeta(ivar)
 7035     format(e12.5,3x,2e12.5)
        enddo
        call pause
      endif
!
!
end subroutine cylinder_geodesic
!
!
!----------------------------------------------------------------------
!!!      subroutine check_cyl_geodesic(Nc, Iflag)
!!!c
!!!      use control
!!!      use GMP
!!!      use element_data
!!!#include "syscom.blk"
!!!c
!!!      dimension center(1:3),axis(3),x(3),dxdeta(3),void(3)
!!!      dimension xrel(3)
!!!c
!!!      nsub=100 ; Iflag=0
!!!c
!!!      ns=CURVES(Nc)%Idata(1)
!!!      center(1:3)=SURFACES(ns)%Rdata(1:3)
!!!      axis(  1:3)=SURFACES(ns)%Rdata(4:6)
!!!      call normalize(axis(1:3))
!!!c
!!!c  ...shift center
!!!      np=CURVES(Nc)%EndPoNo(1)
!!!      void(1:3) = POINTS(np)%Rdata(1:3) - center(1:3)
!!!      call scalar_product(void,axis, s)
!!!      center(1:3) = center(1:3) + s*axis(1:3)
!!!c
!!!c  ...xrel is prep to axis and vector from center through
!!!c     first endpoint
!!!      void(1:3) = POINTS(np)%Rdata(1:3) - center(1:3)
!!!      call cross_product(axis,void,xrel)
!!!      call normalize(xrel)
!!!c
!!!c  ...compute values at second endpoint
!!!      np=CURVES(Nc)%EndPoNo(2)
!!!      void(1:3) = POINTS(np)%Rdata(1:3) - center(1:3)
!!!      call scalar_product(void,axis, x_max)
!!!      call scalar_product(void,xrel, y_max)
!!!c
!!!      do i=0,nsub
!!!        eta=float(i)*1.d0/nsub
!!!        call cylinder_geodesic(Nc,eta, x,dxdeta)
!!!c
!!!c  .....compute cartesian coordinates of projection
!!!        x(1:3)=x(1:3)-center(1:3)
!!!        call scalar_product(axis,x, s_x)
!!!        call scalar_product(xrel,x, s_y)
!!!        s_x=s_x/x_max ; s_y=s_y/y_max
!!!ccc        write(*,*)s_y/s_x
!!!        write(*,*)'s_x,s_y = ',s_x,s_y
!!!ccc        s=s/s_max
!!!ccc        if (abs(s-eta).gt.1.d-10) then
!!!ccc          write(*,*)'s,eta = ',s,eta
!!!ccc        endif
!!!c
!!!c  .....second check
!!!        call scalar_product(axis,dxdeta, s_new)
!!!        if (i.eq.0) then
!!!          s_old=s_new ; cycle
!!!        endif
!!!        if (abs(s_new-s_old).gt.1.d-10) then
!!!          Iflag=1
!!!          write(*,7111)eta,s_old,s_new
!!!7111      format(' eta,s_old,s_new = ',f4.2,2x,2(e12.5,2x))
!!!        endif
!!!        s_old=s_new
!!!      enddo
!!!      if (Iflag.eq.0) then
!!!        write(*,*)'PASS!'
!!!      else
!!!        write(*,*)'FAIL!'
!!!      endif
!!!c
!!!c
!!!      end
!
!
!----------------------------------------------------------------------
!> Purpose : parametrization for a geodesics on a cone
!!
!! @param[in]  No     - curve number
!! @param[in]  Eta    - reference coordinates of a point
!! @param[out] X      - physical coordinates of the point
!! @param[out] Dxdeta - derivatives wrt reference coordinates
!!
!! @date Mar 2023
!----------------------------------------------------------------------
subroutine cone_geodesic(No,Eta, X,Dxdeta)
!
      use control
      use GMP
!
      implicit none
!
      integer :: No
      real(8) :: Eta,X(3),Dxdeta(3)
!
!  ...cone vertex and direction vector
      real(8) :: vertex(3),direction(3)
!
!  ...relative position vectors for curve endpoints in the global
!     cartesian system
      real(8) :: xrelv(3,2)
!
!  ...endpoint coordinates in the cone Cartesian system and in the
!     cone reference coordinates
      real(8) :: xrelsv(3,2),xparv(2,2)
!
!  ...transformation matrix from global Cartesian to cone Cartesian
!     coordinates; rows of the transformation matrix are unit
!     vectors of the cone system computed in the global system
      real(8) :: transf(3,3)
!
!  ...point in the parameter space, derivatives
      real(8) :: xpar(2),dxpardeta(2)
!
!  ...local cartesian coordinates, derivatives
      real(8) :: xrels(3),dxrelsdeta(3)
!
!  ...work space
      real(8) :: void(3)
!
      integer :: i,iv,ivar,j,np,ns
      real(8) :: ap,fval,r,s
!
      integer :: iprint=0
!-----------------------------------------------------------------------
!
      if (iprint.eq.1) then
        write(*,*) '-----------------------------------'
        write(*,7001) No,Eta
 7001   format('cone_geodesic: No,Eta = ',i4,2x,e12.5)
      endif
!
      ns = CURVES(No)%Idata(1)
      if (SURFACES(ns)%Type.ne.'Cone') then
        write(*,*) 'cone_geodesic: INCONSISTENT SURFACE TYPE'
        stop
      endif
!
!  ...get cone Rdata
      vertex(1:3) = SURFACES(ns)%Rdata(1:3)
      direction(1:3) = SURFACES(ns)%Rdata(4:6)
      call normalize(direction(1:3))
      ap = SURFACES(ns)%Rdata(7)
!
!  ...compute the relative position vector for vertices
      do iv=1,2
        np = CURVES(No)%EndPoNo(iv)
!
!  .....check if on the surface
        call surf(ns,POINTS(np)%Rdata(1:3), fval,void)
        if (abs(fval).gt.GEOM_TOL) then
          write(*,7011) No,ns,iv,np,fval
 7011     format(' cone_geodesic: No,ns,iv,np,fval = ', &
                                      i5,2i3,i6,e12.5)
          call pause
          call print_GMP
        endif
        xrelv(1:3,iv) = POINTS(np)%Rdata(1:3) - vertex(1:3)
        call scalar_product(xrelv(1:3,iv),direction(1:3), s)
        if (s.lt.0.d0) then
          write(*,*)'cone_geodesic: flip direction of axis'
          write(*,*)'for cone = ',ns
          write(*,*)'Then pray and try again!'
          stop
        endif
      enddo
!
!  ...compute the transformation matrix from the global Cartesian
!     to the cone Cartesian system
!
!  ...z axis coincides with the cone axis
      transf(3,1:3) = direction(1:3)
!
!  ...cone (x,z) plane goes through "average" of endpoints (to avoid vertex)
      void(1:3) = xrelv(1:3,1) + xrelv(1:3,2)
      call scalar_product(direction(1:3),void(1:3), s)
      void(1:3) = void(1:3) - s*direction(1:3)
!  ...check
      call scalar_product(void(1:3),direction(1:3), s)
      if (abs(s).gt.GEOM_TOL) then
        write(*,*)'geodesic_cone: wrong cone system of coordinates!'
        stop
      endif
      call norm(void(1:3), s)
      transf(1,1:3) = void(1:3)/s
!
!  ...y axis
      call cross_product(transf(3,1:3),transf(1,1:3), transf(2,1:3))
!
!  ...compute endpoint coordinates in the cylinder Cartesian system
      do iv=1,2
        do i=1,3
          s = 0.d0
          do j=1,3
            s = s + transf(i,j)*xrelv(j,iv)
          enddo
          xrelsv(i,iv) = s
        enddo
      enddo
!
!  ...compute the vertex coordinates in the cone parametric space
      do iv=1,2
        call cart_to_polar(xrelsv(1:2,iv), r,xparv(1,iv))
        xparv(2,iv) = xrelsv(3,iv)*tan(ap)
      enddo
      if (iprint.eq.1) then
        write(*,*) 'cone_geodesic: RELATIVE VERTEX COORDINATES ='
        do iv=1,2
          write(*,7031) iv,xrelsv(1:3,iv),xparv(1:2,iv)
 7031     format('VERTEX ',i1,' COORDINATES = ',3e12.5, &
                 ' THETA,Z = ',2e12.5)
        enddo
      endif
!
!  ...geodesics on the cone is simply a straight line in the
!     (\theta,r) parameter space
      xpar(1:2) = xparv(1:2,1)*(1.d0-Eta) + xparv(1:2,2)*Eta
      dxpardeta(1:2) = xparv(1:2,2) - xparv(1:2,1)
!
!  ...local cartesian system
      xrels(1) = xpar(2)*cos(xpar(1))
      dxrelsdeta(1) = dxpardeta(2)*cos(xpar(1))              - &
                           xpar(2)*sin(xpar(1))*dxpardeta(1)
      xrels(2) = xpar(2)*sin(xpar(1))
      dxrelsdeta(2) = dxpardeta(2)*sin(xpar(1))              + &
                           xpar(2)*cos(xpar(1))*dxpardeta(1)
      xrels(3) = xpar(2)/tan(ap)
      dxrelsdeta(3) = dxpardeta(2)/tan(ap)
!
!  ...transform to the global Cartesian system
      X(1:3) = 0.d0; Dxdeta(1:3) = 0.d0
      do i=1,3
        do j=1,3
          X(i) = X(i) + transf(j,i)*xrels(j)
          Dxdeta(i) = Dxdeta(i) + transf(j,i)*dxrelsdeta(j)
        enddo
        X(i) = X(i) + vertex(i)
      enddo
      if (iprint.eq.1) then
        write(*,*) 'cone_geodesic: X,Dxdeta = '
        do ivar=1,3
          write(*,7035) X(ivar),Dxdeta(ivar)
 7035     format(e12.5,3x,2e12.5)
        enddo
        call pause
      endif
!
!
end subroutine cone_geodesic

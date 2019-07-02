!-------------------------------------------------------------------------------
!> Purpose : parametric transfinite interpolation for a triangle
!!
!! @param[in]  No     - a GMP triangle number
!! @param[in]  Eta    - reference coordinates of a point
!! @param[out] X      - physical coordinates of the point
!! @param[out] Dxdeta - derivatives of the physical coordinates
!!
!! @ revision Mar 11
!-------------------------------------------------------------------------------
!
subroutine trian_PTITri(No,Eta, X,Dxdeta)
!
      use GMP
!-------------------------------------------------------------------------------
      implicit none
      integer,              intent(in)  :: No
      real*8,dimension(2),  intent(in)  :: Eta
      real*8,dimension(3),  intent(out) :: X
      real*8,dimension(3,2),intent(out) :: Dxdeta
!-------------------------------------------------------------------------------
      integer :: ns, iprint
!-------------------------------------------------------------------------------
!
      iprint=0
!
!  ...get surface number
      ns = TRIANGLES(No)%Idata(1)
      if (iprint.eq.1) then
        write(*,7001) No,Eta
 7001   format('trian_PTITri: No,Eta = ',i4,2x,2f8.3)
        write(*,7002) ns,SURFACES(ns)%Type
 7002   format('              ns,Type = ',i3,2x,a10)
      endif
!
!  ...select surface type
      select case(SURFACES(ns)%Type)
      case('Sphere'         ) ; call spherical_triangle(  No,Eta, X,Dxdeta)
      case('Cylinder'       ) ; call cylindrical_triangle(No,Eta, X,Dxdeta)
      case('Cone'           ) ; call conic_triangle(      No,Eta, X,Dxdeta)
      case('VecPt','ThreePt') ; call trian_TransTri(      No,Eta, X,Dxdeta)
      case('PPwCC'          ) ; call PPwCC_triangle(      No,Eta, X,Dxdeta)
      case default
        write(*,7020) ns,SURFACES(ns)%Type
 7020   format('trian_PTITri: surface type not supported, ns,Type = ', &
                i5,2x,a10 )
        stop
      end select
!
!
end subroutine trian_PTITri
!
!
!
!----------------------------------------------------------------------
!> Purpose : parametric transfinite interpolation for a triangle on a
!!           sphere
!!
!! @param[in]  No        - triangle number
!! @param[in]  Eta       - reference coordinates of a point
!! @param[out] X         - physical coordinates of the point
!! @param[out] Dxdeta    - derivatives wrt reference coordinates
!!
!! @revision Mar 11
!----------------------------------------------------------------------
!
subroutine spherical_triangle(No,Eta, X,Dxdeta)
!
      use control
      use GMP
      use element_data
#include "syscom.blk"
      common /cspherical_triangle/ iprint
!
      dimension Eta(2), X(3),Dxdeta(3,2)
!
!  ...relative position vectors wrt sphere origin
      dimension xrels(3),xpar(2),dxrels(3,2),dxpardeta(2,2)
!
!  ...relative position vectors for triangle vertices in the global
!     cartesian system
      dimension xrelv(3,3)
!
!  ...vertex coordinates in the sphere Cartesian system and in the
!     sphere reference coordinates
      dimension xrelsv(3,3),xparv(2,3)
!
!  ...transformation matrix from global Cartesian to sphere Cartesian
!     coordinates; rows of the transformation matrix are unit
!     vectors of the sphere system computed in the global system
      dimension transf(3,3)
!
!  ...linear shape functions in the parametric space (psi,theta)
      dimension shapH(3),dshapH(2,3)
!
!  ...point on an edge
      dimension dsedeta(2),                                  &
                xe(3),dxedeta(3),xerels(3),dxerelsdeta(3),   &
                xepar(2),dxepardeta(2)
!
!  ...workspace
      dimension void(3),dblend(2)
!
!-----------------------------------------------------------------------
!
      iprint=0
      if (iprint.eq.1) then
        write(*,*) '-----------------------------------'
        write(*,7001) No,Eta(1:2)
 7001   format('spherical triangle: No,Eta = ',i4,2x,2e12.5)
      endif
!
!  ...check surface type
      ns = TRIANGLES(No)%Idata(1)
      if (SURFACES(ns)%Type.ne.'Sphere') then
        write(*,*) 'spherical_triangle: INCONSISTENT SURFACE TYPE'
        stop
      endif
      rad = SURFACES(ns)%Rdata(4)
!
!----------------------------------------------------------------------
!  ...compute the relative position vector for vertices
!
!  ...loop through vertices
      do iv=1,3
!  .....get global point number
        np = TRIANGLES(No)%VertNo(iv)
!
!  .....check if point is on the surface
        call surf(ns,POINTS(np)%Rdata(1:3), fval,void)
        if (abs(fval) .gt. GEOM_TOL) then
          write(*,7011) No,ns,iv,np,fval
 7011     format('spherical_triangle: No,ns,iv,np,fval = ',  &
                                      i5,2i3,i6,2x,e12.5)
          call norm(POINTS(np)%Rdata(1:3), r_aux)
          write(*,9000) SURFACES(ns)%Rdata(4),r_aux
 9000     format('spherical_triangle: actual rad, computed rad = ',   &
                                                       e12.5,2x e12.5)
          call pause
          call print_GMP
        endif
!
        xrelv(1:3,iv) = POINTS(np)%Rdata(1:3) - SURFACES(ns)%Rdata(1:3)
!  ...end of loop through vertices
      enddo
!
!----------------------------------------------------------------------
!  ...compute the transformation matrix from the global Cartesian
!     system to the sphere Cartesian system
!
!  ...the sphere x axis goes through the center of the triangle
      void(1:3)=xrelv(1:3,1)+xrelv(1:3,2)+xrelv(1:3,3)
      call norm(void(1:3), s)
      transf(1,1:3) = void(1:3)/s
!
!  ...second vertex lies on the sphere xy plane
      call cross_product(xrelv(1:3,3),transf(1,1:3), void)
      call norm(void, s)
      transf(2,1:3) = void(1:3)/s
      call cross_product(transf(1,1:3),transf(2,1:3), transf(3,1:3))
!
!----------------------------------------------------------------------
!  ...compute vertex coordinates in the sphere Cartesian system
      do iv = 1, 3
        do i = 1, 3
          s = 0.d0
          do j = 1, 3
            s = s + transf(i,j)*xrelv(j,iv)
          enddo
          xrelsv(i,iv) = s
        enddo
      enddo
!
!----------------------------------------------------------------------
!  ...compute the vertex coordinates in the sphere parametric space
      do iv=1,3
!  .....project on horizontal plane to determine theta
        call cart_to_polar(xrelsv(1:2,iv), rsinpsi,xparv(2,iv))
!  .....project on vertical plane to determine psi
        call cart_to_polar((/xrelsv(3,iv),rsinpsi/), r,xparv(1,iv))
        if (abs(rad-r) .gt. GEOM_TOL) then
          write(*,7012) iv,rad,r
 7012     format('spherical_triangle: iv,rad,r = ',i2,2e12.5)
          call pause
        endif
      enddo
!
!  ...printing statement
      if (iprint .eq. 1) then
        write(*,*) 'spherical_triangle: RELATIVE VERTEX COORDINATES = '
        do iv=1,3
          write(*,7031) iv,xrelsv(1:3,iv),xparv(1:2,iv)
 7031     format('VERTEX ',i1,' COORDINATES = ',3e12.5,  &
                 ' PSI,THETA = ',2e12.5)
        enddo
      endif
!-----------------------------------------------------------------------
!  ...evaluate linear shape functions
      shapH(1) = 1.d0 - Eta(1) - Eta(2); dshapH(1:2,1) = -1.d0
      shapH(2) = Eta(1); dshapH(1,2) = 1.d0; dshapH(2,2) = 0.d0
      shapH(3) = Eta(2); dshapH(1,3) = 0.d0; dshapH(2,3) = 1.d0
!
!-----------------------------------------------------------------------
!     L I N E A R    I N T E R P O L A T I O N
!-----------------------------------------------------------------------
      xpar(1:2) = 0.d0; dxpardeta(1:2,1:2) = 0.d0
      do iv=1,3
        xpar(1:2) = xpar(1:2) + xparv(1:2,iv)*shapH(iv)
        do j=1,2
          dxpardeta(1:2,j) = dxpardeta(1:2,j)  &
                           + xparv(1:2,iv)*dshapH(j,iv)
        enddo
      enddo
!
!  ...printing statement
      if (iprint .eq. 1) then
        write(*,*) 'spherical_triangle: AFTER VERT xpar,dxpardeta = '
        do ivar=1,2
          write(*,7035) xpar(ivar),dxpardeta(ivar,1:2)
        enddo
      endif
!
!----------------------------------------------------------------------
!     E D G E    B U B B L E S
!----------------------------------------------------------------------
!  ...loop over edges
      do ie=1,3
!
!  .....get the curve number
        nc = TRIANGLES(No)%EdgeNo(ie); norient = 0
        if (nc.lt.0) then
          nc = -nc; norient = 1
        endif
        if (iprint .eq. 1) then
          write(*,7003) ie, nc, CURVES(nc)%Type
 7003     format('spherical_triangle: ie,nc,Type = ',i2,i5,2x,a5)
        endif
!
!  .....get the edge vertices specifying the local edge orientation
        iv1=TRIAN_EDGE_TO_VERT(1,ie) ; iv2=TRIAN_EDGE_TO_VERT(2,ie)
!
!  .....project Eta onto the edge ie
        call proj_t2e(iv1,iv2,shapH,dshapH, se,dsedeta)
!  .....cycle if start or end point
        if ((se.lt.GEOM_TOL).or.(se.gt.1.d0-GEOM_TOL)) cycle
!
!  .....compute the edge coordinate and its derivative
        select case(norient)
        case(0); s = se; dsdse = 1.d0
        case(1); s = 1.d0 - se; dsdse = -1.d0
        end select
        call curve(nc,s, xe,dxedeta)
        dxedeta(1:3) = dxedeta(1:3)*dsdse
!
!  .....transform to the sphere Cartesian system
        xerels(1:3) = 0.d0; dxerelsdeta(1:3) = 0.d0
        do j = 1,3
          xerels(1:3) = xerels(1:3)  &
                      + transf(1:3,j)*(xe(j) - SURFACES(ns)%Rdata(j))
          dxerelsdeta(1:3) = dxerelsdeta(1:3)  &
                           + transf(1:3,j)*dxedeta(j)
        enddo
!
!  .....transform to the parametric space
!  .....xepar(1) = psi  ;  xepar(2) = theta  ;  -pi <= psi,theta <= pi
!       (usual spherical polars, NO Gilardi's style!!!)
        call cart_to_polar(xerels(1:2), rsinpsi,xepar(2))
        call cart_to_polar((/xerels(3),rsinpsi/), r,xepar(1))
!
!  .....geometry consistency check and printing statements
        if (iprint.eq.1) then
          write(*,7032) s,xerels(1:3),xepar(1:2)
 7032     format('s = ',e12.5,' COORDINATES = ',3e12.5, &
                 ' PSI,THETA = ',2e12.5)
        endif
        if (abs(rad-r).gt.GEOM_TOL) then
          write(*,7004) rad,r
 7004     format('spherical_triangle: rad,r = ',2e12.5)
!c          call pause
        endif
        if (abs(sin(xepar(1))).lt.GEOM_TOL) then
          write(*,*) 'spherical triangle: POINT ON A POLE'
          stop
        endif
!
!  .....compute dxepardeta(1:2)
        dxepardeta(1) = -dxerelsdeta(3)/(rad*sin(xepar(1)))
        s1 = dxerelsdeta(1) - xerels(3)*dxepardeta(1)*cos(xepar(2))
        s2 = dxerelsdeta(2) - xerels(3)*dxepardeta(1)*sin(xepar(2))
        ile=0
        if (abs(xerels(2)).gt.GEOM_TOL) then
          ile=ile+1
          s1 = -s1/xerels(2); dxepardeta(2) = s1
        endif
        if (abs(xerels(1)).gt.GEOM_TOL) then
          ile = ile + 1
          s2 =  s2/xerels(1); dxepardeta(2) = s2
        endif
        if (iprint.eq.1) then
          write(*,7027) dxepardeta(1),s1,s2
 7027     format('spherical_triangle: dxepardeta(1),s1,s2 = ',3e12.5)
        endif
        if (ile.eq.2) then
          if (abs(s1-s2).gt.GEOM_TOL) then
            write(*,7005) s1,s2
 7005       format('spherical_triangle: s1,s2 = ',2e12.5)
!c            call pause
          endif
        endif
        if (iprint.eq.1) then
          write(*,*) 'ie,xepar,dxepardeta = ',ie
          do ivar=1,2
            write(*,7038) xepar(ivar),dxepardeta(ivar)
 7038       format(e12.5,2x,2e12.5)
          enddo
        endif
!
!  .....compute the bubble
        xepar(1:2) = xepar(1:2)  &
                   - (xparv(1:2,iv1)*(1.d0-se)+xparv(1:2,iv2)*se)
        dxepardeta(1:2) = dxepardeta(1:2)  &
                   - (xparv(1:2,iv2) - xparv(1:2,iv1))
!
!  .....compute the kernel
        b = (1.d0-se)*se; db = 1.d0-2.d0*se
        dxepardeta(1:2) = dxepardeta(1:2)/b   &
                        - xepar(1:2)*db/b**2
        xepar(1:2) = xepar(1:2)/b
!
!  .....add the edge contribution
        blend = shapH(iv1)*shapH(iv2)
        dblend(1:2) = dshapH(1:2,iv1)*shapH(iv2)  &
                    + shapH(iv1)*dshapH(1:2,iv2)
        xpar(1:2) = xpar(1:2) + xepar(1:2)*blend
        do j = 1,2
          dxpardeta(1:2,j) = dxpardeta(1:2,j)                   &
                           + dxepardeta(1:2)*dsedeta(j)*blend   &
                           + xepar(1:2)*dblend(j)
        enddo
        if (iprint .eq. 1) then
          write(*,*) 'ie,xpar,dxpardeta = ',ie
          do ivar=1,2
            write(*,7038) xpar(ivar),dxpardeta(ivar,1:2)
          enddo
        endif
!
!  ...loop over edges
      enddo
!
!  ...printing
      if (iprint.eq.1) then
        write(*,*) 'spherical_triangle: FINAL xpar,dxpardeta = '
        do ivar=1,2
          write(*,7035) xpar(ivar),dxpardeta(ivar,1:2)
        enddo
      endif
!
!-----------------------------------------------------------------------
!
!  ...transform to sphere Cartesian coordinates
      xrels(1) = rad*sin(xpar(1))*cos(xpar(2))
      dxrels(1,1:2) = rad*cos(xpar(1))*dxpardeta(1,1:2)*cos(xpar(2))  &
                    - rad*sin(xpar(1))*sin(xpar(2))*dxpardeta(2,1:2)
      xrels(2) = rad*sin(xpar(1))*sin(xpar(2))
      dxrels(2,1:2) = rad*cos(xpar(1))*dxpardeta(1,1:2)*sin(xpar(2))  &
                    + rad*sin(xpar(1))*cos(xpar(2))*dxpardeta(2,1:2)
      xrels(3) = rad*cos(xpar(1))
      dxrels(3,1:2) = -rad*sin(xpar(1))*dxpardeta(1,1:2)
!
!  ...transform to global Cartesian coordinates
      X(1:3) = 0.d0; Dxdeta(1:3,1:2) = 0.d0
      do i=1,3
        do j=1,3
          X(i) = X(i) + transf(j,i)*xrels(j)
          Dxdeta(i,1:2) = Dxdeta(i,1:2) + transf(j,i)*dxrels(j,1:2)
        enddo
        X(i) = X(i) + SURFACES(ns)%Rdata(i)
      enddo
!
!  ...verification....
      call norm(X, radnew)
      if (abs(radnew-rad).gt.GEOM_TOL) then
        write(*,8001) rad,radnew
 8001   format('spherical_triangle: rad,radnew = ',2e12.5)
        call pause
      endif
      do ieta=1,2
        call scalar_product(X,Dxdeta(1:3,ieta), s)
        if (abs(s).gt.GEOM_TOL) then
          write(*,8002) ieta,s
 8002     format('spherical_triangle: ieta,s = ',i1,e12.5)
          call pause
        endif
      enddo
!
!  ...printing statement
      if (iprint.eq.1) then
        write(*,*) 'spherical_triangle: X,Dxdeta = '
        do ivar=1,3
          write(*,7035) X(ivar),Dxdeta(ivar,1:2)
 7035     format(e12.5,3x,3e12.5)
        enddo
        call pause
      endif
!
end subroutine spherical_triangle
!
!
!
!----------------------------------------------------------------------
!> Purpose : parametric transfinite interpolation for a triangle
!!           on a cylinder
!!
!! @param[in]  No     - triangle number
!! @param[in]  Eta    - reference coordinates of a point
!! @param[out] X      - physical coordinates of the point
!! @param[out] Dxdeta - derivatives wrt reference coordinates
!!
!! @revision Mar 11
!----------------------------------------------------------------------
!
subroutine cylindrical_triangle(No,Eta, X,Dxdeta)
!
      use control
      use GMP
      use element_data
#include "syscom.blk"
      common /ccylindrical_triangle/ iprint
!
      dimension Eta(2), X(3),Dxdeta(3,2)
!
!  ...cylinder center and direction vector
      dimension center(3), direction(3)
!
!  ...relative position vectors wrt center of the cylinder, cylinder
!     parameters (\theta,z), derivatives
      dimension xrels(3),xpar(2),dxrels(3,2),dxpardeta(2,2)
!
!  ...relative position vectors for triangle vertices in the global
!     cartesian system
      dimension xrelv(3,3)
!
!  ...vertex coordinates in the cylinder Cartesian system and in the
!     cylinder reference coordinates
      dimension xrelsv(3,3),xparv(2,3)
!
!  ...transformation matrix from global Cartesian to cylinder Cartesian
!     coordinates; rows of the transformation matrix are unit
!     vectors of the cylinder system computed in the global system
      dimension transf(3,3)
!
!  ...linear shape functions in the parametric space (psi,theta)
      dimension shapH(3),dshapH(2,3)
!
!  ...point on an edge
      dimension dsedeta(2),                                 &
                xe(3),dxedeta(3),xerels(3),dxerelsdeta(3),  &
                xepar(2),dxepardeta(2)
!
!  ...work space
      dimension void(3),dblend(2)
!
!-----------------------------------------------------------------------
!
      iprint=0
      if (iprint.eq.1) then
        write(*,*) '-----------------------------------'
        write(*,7001) No,Eta(1:2)
 7001   format('cylindrical triangle: No,Eta = ',i4,2x,2e12.5)
      endif
!
      ns = TRIANGLES(No)%Idata(1)
      if (SURFACES(ns)%Type.ne.'Cylinder') then
        write(*,*) 'cylindrical_triangle: INCONSISTENT SURFACE TYPE'
        stop
      endif
      center(1:3) = SURFACES(ns)%Rdata(1:3)
!
!  ...shift the originally specified center to a new location such
!     that the new z coordinate of the first vertex will be zero
      direction(1:3) = SURFACES(ns)%Rdata(4:6)
      call normalize(direction)
      np = TRIANGLES(No)%VertNo(1)
      void(1:3) = POINTS(np)%Rdata(1:3) - center(1:3)
      call scalar_product(void,direction, s)
      center(1:3) = center(1:3) + s*direction(1:3)
      rad = SURFACES(ns)%Rdata(7)
!
!  ...compute the relative position vector for vertices
      do iv=1,3
        np = TRIANGLES(No)%VertNo(iv)
!
!  .....check if on the surface
        call surf(ns,POINTS(np)%Rdata(1:3), fval,void)
        if (abs(fval).gt.GEOM_TOL) then
          write(*,7011) No,ns,iv,np,fval
 7011     format('cylindrical_triangle: No,ns,iv,np,fval = ',  &
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
!  ...compute vertex coordinates in the cylinder Cartesian system
      do iv=1,3
        do i=1,3
          s = 0.d0
          do j=1,3
            s = s + transf(i,j)*xrelv(j,iv)
          enddo
          xrelsv(i,iv) = s
        enddo
      enddo
!
!  ...compute the vertex coordinates in the cylinder parametric space
      do iv=1,3
        call cart_to_polar(xrelsv(1:2,iv), r,xparv(1,iv))
        if (abs(rad-r).gt.GEOM_TOL) then
          write(*,7012) iv,rad,r
 7012     format('cylindrical_triangle: iv,rad,r = ',i2,2e12.5)
          call pause
        endif
        xparv(2,iv) = xrelsv(3,iv)
      enddo
      if (iprint.eq.1) then
        write(*,*) 'cylindrical_triangle: RELATIVE VERTEX COORDINATES ='
        do iv=1,3
          write(*,7031) iv,xrelsv(1:3,iv),xparv(1:2,iv)
 7031     format('VERTEX ',i1,' COORDINATES = ',3e12.5,    &
                 ' THETA,Z = ',2e12.5)
        enddo
      endif
!
!  ...evaluate linear shape functions
      shapH(1) = 1.d0-Eta(1)-Eta(2); dshapH(1:2,1) = -1.d0
      shapH(2) = Eta(1); dshapH(1,2) = 1.d0; dshapH(2,2) = 0.d0
      shapH(3) = Eta(2); dshapH(1,3) = 0.d0; dshapH(2,3) = 1.d0
!
!-----------------------------------------------------------------------
!
!  ...linear interpolant
      xpar(1:2) = 0.d0; dxpardeta(1:2,1:2) = 0.d0
      do iv=1,3
        xpar(1:2) = xpar(1:2) + xparv(1:2,iv)*shapH(iv)
        do j=1,2
          dxpardeta(1:2,j) = dxpardeta(1:2,j)    &
                           + xparv(1:2,iv)*dshapH(j,iv)
        enddo
      enddo
      if (iprint.eq.1) then
        write(*,*) 'cylindrical_triangle: AFTER VERT xpar,dxpardeta ='
        do ivar=1,2
          write(*,7035) xpar(ivar),dxpardeta(ivar,1:2)
        enddo
      endif
!
!  ...add edge bubbles
      do ie=1,3
!
!  .....get the curve number
        nc = TRIANGLES(No)%EdgeNo(ie); norient=0
        if (nc.lt.0) then
          nc = -nc; norient=1
        endif
        if (iprint.eq.1) then
          write(*,7003) ie,nc,CURVES(nc)%Type
 7003     format('cylindrical_triangle: ie,nc,Type = ',i2,i5,2x,a5)
        endif
!
!  .....get the edge vertices specifying the local edge orientation
        iv1=TRIAN_EDGE_TO_VERT(1,ie) ; iv2=TRIAN_EDGE_TO_VERT(2,ie)
!
!  .....project s onto the edge
        call proj_t2e(iv1,iv2,shapH,dshapH, se,dsedeta)
        if ((se.lt.GEOM_TOL).or.(se.gt.1.d0-GEOM_TOL)) cycle
!
!  .....compute the edge parametrization
        select case(norient)
        case(0); s = se; dsdse = 1.d0
        case(1); s = 1.d0-se; dsdse = -1.d0
        end select
        call curve(nc,s, xe,dxedeta)
        dxedeta(1:3) = dxedeta(1:3)*dsdse
!
!  .....transform to the cylinder Cartesian system
        xerels(1:3) = 0.d0; dxerelsdeta(1:3) = 0.d0
        do j=1,3
          xerels(1:3) = xerels(1:3)   &
                      + transf(1:3,j)*(xe(j)- center(j))
          dxerelsdeta(1:3) = dxerelsdeta(1:3)  &
                           + transf(1:3,j)*dxedeta(j)
        enddo
!
!  .....transform to the parametric space
        call cart_to_polar(xerels(1:2), r,xepar(1))
        xepar(2) = xerels(3)
        if (iprint.eq.1) then
          write(*,7032) s,xerels(1:3),xepar(1:2)
 7032     format('s = ',e12.5,' COORDINATES = ',3e12.5,  &
                 ' THETA,Z = ',2e12.5)
        endif
        if (abs(rad-r).gt.GEOM_TOL) then
          write(*,7004) rad,r
 7004     format('cylindrical_triangle: rad,r = ',2e12.5)
          call pause
        endif
!
!  .....this part is specific for cylindrical coordinates
!cc        write(*,*)'dxerelsdeta(1)*xerels(2) = ',dxerelsdeta(1)*xerels(2)
!cc        write(*,*)'dxerelsdeta(2)*xerels(1) = ',dxerelsdeta(2)*xerels(1)
        dxepardeta(2) = dxerelsdeta(3)
        s1 = dxerelsdeta(1)
        s2 = dxerelsdeta(2)
        ile=0
        if (abs(xerels(2)).gt.GEOM_TOL) then
          ile=ile+1
          s1 = -s1/xerels(2); dxepardeta(1) = s1
        endif
        if (abs(xerels(1)).gt.GEOM_TOL) then
          ile=ile+1
          s2 =  s2/xerels(1); dxepardeta(1) = s2
        endif
        if (iprint.eq.1) then
          write(*,7027) s1,s2,dxepardeta(2)
 7027     format('cylindrical_triangle: s1,s2,dxepardeta(2) = ',3e12.5)
        endif
        if (ile.eq.2) then
          if (abs(s1-s2).gt.GEOM_TOL) then
            write(*,7005) s1,s2
 7005       format('cylindrical_triangle: s1,s2 = ',2e12.5)
            call pause
          endif
        endif
        if (iprint.eq.1) then
          write(*,*) 'ie,xepar,dxepardeta = ',ie
          do ivar=1,2
            write(*,7038) xepar(ivar),dxepardeta(ivar)
 7038       format(e12.5,2x,2e12.5)
          enddo
        endif
!
!  .....compute the bubble
        xepar(1:2) = xepar(1:2)  &
                   - (xparv(1:2,iv1)*(1.d0-se)+xparv(1:2,iv2)*se)
        dxepardeta(1:2) = dxepardeta(1:2)  &
                   - (xparv(1:2,iv2) - xparv(1:2,iv1))
!
!  .....compute the kernel
        b = (1.d0-se)*se; db = 1.d0-2.d0*se
        dxepardeta(1:2) = dxepardeta(1:2)/b &
                        - xepar(1:2)*db/b**2
        xepar(1:2) = xepar(1:2)/b
!
!  .....add the edge contribution
        blend = shapH(iv1)*shapH(iv2)
        dblend(1:2) = dshapH(1:2,iv1)*shapH(iv2)  &
                    + shapH(iv1)*dshapH(1:2,iv2)
        xpar(1:2) = xpar(1:2) + xepar(1:2)*blend
        do j=1,2
          dxpardeta(1:2,j) = dxpardeta(1:2,j)     &
                           + dxepardeta(1:2)*dsedeta(j)*blend  &
                           + xepar(1:2)*dblend(j)
        enddo
        if (iprint.eq.1) then
          write(*,*) 'ie,xpar,dxpardeta = ',ie
          do ivar=1,2
            write(*,7038) xpar(ivar),dxpardeta(ivar,1:2)
          enddo
        endif
      enddo
      if (iprint.eq.1) then
        write(*,*) 'cylindrical_triangle: FINAL xpar,dxpardeta = '
        do ivar=1,2
          write(*,7035) xpar(ivar),dxpardeta(ivar,1:2)
        enddo
      endif
!
!-----------------------------------------------------------------------
!
!  ...transform to sphere Cartesian coordinates
      xrels(1) = rad*cos(xpar(1))
      dxrels(1,1:2) = -rad*sin(xpar(1))*dxpardeta(1,1:2)
      xrels(2) = rad*sin(xpar(1))
      dxrels(2,1:2) =  rad*cos(xpar(1))*dxpardeta(1,1:2)
      xrels(3) = xpar(2)
      dxrels(3,1:2) = dxpardeta(2,1:2)
!
!  ...transform to global Cartesian coordinates
      X(1:3) = 0.d0; Dxdeta(1:3,1:2) = 0.d0
      do i=1,3
        do j=1,3
          X(i) = X(i) + transf(j,i)*xrels(j)
          Dxdeta(i,1:2) = Dxdeta(i,1:2) + transf(j,i)*dxrels(j,1:2)
        enddo
        X(i) = X(i) + center(i)
      enddo
      if (iprint.eq.1) then
        write(*,*) 'cylindrical_triangle: X,Dxdeta = '
        do ivar=1,3
          write(*,7035) X(ivar),Dxdeta(ivar,1:2)
 7035     format(e12.5,3x,3e12.5)
        enddo
        call pause
      endif
!
!
end subroutine cylindrical_triangle
!
!
!
!----------------------------------------------------------------------
!> Purpose : parametric transfinite interpolation for a triangle on
!!           a cone
!!
!! @param[in]  No     - triangle number
!! @param[in]  Eta    - reference coordinates of a point
!! @param[out] X      - physical coordinates of the point
!! @param[out] Dxdeta - derivatives wrt reference coordinates
!!
!! @revision Mar 11
!----------------------------------------------------------------------
!
subroutine conic_triangle(No,Eta, X,Dxdeta)
!
      use control
      use GMP
      use element_data
#include "syscom.blk"
      common /cconic_triangle/ iprint
!
      dimension Eta(2), X(3),Dxdeta(3,2)
!
!  ...cone axis vector normalized
      dimension cvec(3)
!
!  ...relative position vectors wrt cone origin
      dimension xrelc(3),xpar(2),dxrelc(3,2),dxpardeta(2,2)
!
!  ...relative position vectors for triangle vertices in the global
!     cartesian system
      dimension xrelv(3,3)
!
!  ...vertex coordinates in the cone Cartesian system and in the
!     cone reference coordinates
      dimension xrelcv(3,3),xparv(2,3)
!
!  ...transformation matrix from global Cartesian to cone Cartesian
!     coordinates; rows of the transformation matrix are unit
!     vectors of the cone system computed in the global system
      dimension transf(3,3)
!
!  ...linear shape functions in the parametric space
      dimension shapH(3),dshapH(2,3)
!
!  ...derivative of cone coordinate wrt Cartesian coordinates
      dimension dzetdx(3),dr2dx(3)
!
!  ...point on an edge
      dimension dsedeta(2),  &
                xe(3),dxedeta(3),xerelc(3),dxerelcdeta(3),  &
                xepar(2),dxepardeta(2)
!
!  ...work space
      dimension void(3),dblend(2)
!
!-----------------------------------------------------------------------
!
!cc      iprint=0
      if (iprint.eq.1) then
        write(*,*) '-----------------------------------'
        write(*,7001) No,Eta(1:2)
 7001   format('conic triangle: No,Eta = ',i4,2x,2e12.5)
      endif
!
      ns = TRIANGLES(No)%Idata(1)
      if (SURFACES(ns)%Type.ne.'Cone') then
        write(*,*) 'conic_triangle: INCONSISTENT SURFACE TYPE'
        stop
      endif
!
!  ...compute the cone axis unit vector
      call norm(SURFACES(ns)%Rdata(4:6), rnorm)
      cvec(1:3) = SURFACES(ns)%Rdata(4:6)/rnorm
!
!  ...compute the relative position vector for vertices
      do iv=1,3
        np = TRIANGLES(No)%VertNo(iv)
!
!  .....check if on the surface
        call surf(ns,POINTS(np)%Rdata(1:3), fval,void)
        if (abs(fval).gt.GEOM_TOL) then
          write(*,7011) No,ns,iv,np,fval
 7011     format('conic_triangle: No,ns,iv,np,fval = ',i5,2i3,i6,e12.5)
          call pause
          call print_GMP
        endif
        xrelv(1:3,iv) = POINTS(np)%Rdata(1:3) - SURFACES(ns)%Rdata(1:3)
      enddo
!
!  ...compute the transformation matrix from the global Cartesian
!     to the cone Cartesian system
!
!  ...the cone z axis coincides with the cone axis
      transf(3,1:3) = cvec(1:3)
!
!  ...the cone y coordinate of the first vertex is zero
      call scalar_product(xrelv(1:3,1),cvec, s)
      transf(1,1:3) = xrelv(1:3,1) - s*cvec(1:3)
      call normalize(transf(1,1:3))
      call cross_product(transf(3,1:3),transf(1,1:3), transf(2,1:3))
!
!  ...compute vertex coordinates in the cone Cartesian system
      do iv=1,3
        do i=1,3
          s = 0.d0
          do j=1,3
            s = s + transf(i,j)*xrelv(j,iv)
          enddo
          xrelcv(i,iv) = s
        enddo
      enddo
!
!  ...compute the vertex coordinates in the cone parametric space
      c = tan(SURFACES(ns)%Rdata(7))
      do iv=1,3
        call cart_to_polar(xrelcv(1:2,iv), r,xparv(1,iv))
        if (abs(c*xrelcv(3,iv)-r).gt.GEOM_TOL) then
          write(*,7012) iv,c*xrelcv(3,iv),r
 7012     format('conic_triangle: iv,c*xrelcv(3,iv),r = ',i2,2e12.5)
          call pause
        endif
        xparv(2,iv) = xrelcv(3,iv)
      enddo
      if (iprint.eq.1) then
        write(*,*) 'conic_triangle: RELATIVE VERTEX COORDINATES = '
        do iv=1,3
          write(*,7031) iv,xrelcv(1:3,iv),xparv(1:2,iv)
 7031     format('VERTEX ',i1,' COORDINATES = ',3e12.5,  &
                 ' THETA,Z = ',2e12.5)
        enddo
      endif
!
!  ...evaluate linear shape functions
      shapH(1) = 1.d0-Eta(1)-Eta(2); dshapH(1:2,1) = -1.d0
      shapH(2) = Eta(1); dshapH(1,2) = 1.d0; dshapH(2,2) = 0.d0
      shapH(3) = Eta(2); dshapH(1,3) = 0.d0; dshapH(2,3) = 1.d0
!
!-----------------------------------------------------------------------
!
!  ...linear interpolant
      xpar(1:2) = 0.d0; dxpardeta(1:2,1:2) = 0.d0
      do iv=1,3
        xpar(1:2) = xpar(1:2) + xparv(1:2,iv)*shapH(iv)
        do j=1,2
          dxpardeta(1:2,j) = dxpardeta(1:2,j)   &
                           + xparv(1:2,iv)*dshapH(j,iv)
        enddo
      enddo
      if (iprint.eq.1) then
        write(*,*) 'conic_triangle: AFTER VERT xpar,dxpardeta = '
        do ivar=1,2
          write(*,7035) xpar(ivar),dxpardeta(ivar,1:2)
        enddo
      endif
!
!  ...add edge bubbles
      do ie=1,3
!
!  .....get the curve number
        nc = TRIANGLES(No)%EdgeNo(ie); norient=0
        if (nc.lt.0) then
          nc = -nc; norient=1
        endif
!
!  .....get the edge vertices specifying the local edge orientation
        iv1=TRIAN_EDGE_TO_VERT(1,ie) ; iv2=TRIAN_EDGE_TO_VERT(2,ie)
!
!  .....project s onto the edge
        call proj_t2e(iv1,iv2,shapH,dshapH, se,dsedeta)
        if ((se.lt.GEOM_TOL).or.(se.gt.1.d0-GEOM_TOL)) cycle
        if (iprint.eq.1) then
          write(*,7003) ie,nc,CURVES(nc)%Type
 7003     format('conic_triangle: ie,nc,Type = ',i2,i5,2x,a5)
        endif
!
!  .....compute the edge parametrization
        select case(norient)
        case(0); s = se; dsdse = 1.d0
        case(1); s = 1.d0-se; dsdse = -1.d0
        end select
        call curve(nc,s, xe,dxedeta)
        dxedeta(1:3) = dxedeta(1:3)*dsdse
!
!  .....transform to the cone Cartesian system
        xerelc(1:3) = 0.d0; dxerelcdeta(1:3) = 0.d0
        do j=1,3
          xerelc(1:3) = xerelc(1:3)  &
                      + transf(1:3,j)*(xe(j)- SURFACES(ns)%Rdata(j))
          dxerelcdeta(1:3) = dxerelcdeta(1:3)  &
                           + transf(1:3,j)*dxedeta(j)
        enddo
!
!  .....transform to the parametric space
        call cart_to_polar(xerelc(1:2), r,xepar(1))
        if (iprint.eq.1) then
          write(*,7032) s,xerelc(1:3),xepar(1:2)
 7032     format('s = ',e12.5,' COORDINATES = ',3e12.5,  &
                 ' THETA,Z = ',2e12.5)
        endif
        if (abs(c*xerelc(3)-r).gt.GEOM_TOL) then
          write(*,7004) c*xerelc(3),r
 7004     format('conic_triangle: c*xerelc(3),r = ',2e12.5)
          call pause
        endif
        xepar(2) = xerelc(3)
        dxepardeta(2) = dxerelcdeta(3)
        s1 = dxerelcdeta(1) - xerelc(1)*dxepardeta(2)/xepar(2)
        s2 = dxerelcdeta(2) - xerelc(2)*dxepardeta(2)/xepar(2)
        ile=0
        if (abs(xerelc(2)).gt.GEOM_TOL) then
          ile=ile+1
          s1 = -s1/xerelc(2); dxepardeta(1) = s1
        endif
        if (abs(xerelc(1)).gt.GEOM_TOL) then
          ile=ile+1
          s2 =  s2/xerelc(1); dxepardeta(1) = s2
        endif
        if (ile.eq.2) then
          if (abs(s1-s2).gt.GEOM_TOL) then
            write(*,7005) s1,s2
 7005       format('conic_triangle: s1,s2 = ',2e12.5)
            call pause
          endif
        endif
!
!  .....compute the bubble
        xepar(1:2) = xepar(1:2)  &
                   - (xparv(1:2,iv1)*(1.d0-se)+xparv(1:2,iv2)*se)
        dxepardeta(1:2) = dxepardeta(1:2)  &
                   - (xparv(1:2,iv2) - xparv(1:2,iv1))
!
!  .....compute the kernel
        b = (1.d0-se)*se; db = 1.d0-2.d0*se
        dxepardeta(1:2) = dxepardeta(1:2)/b &
                           - xepar(1:2)*db/b**2
        xepar(1:2) = xepar(1:2)/b
!
!  .....add the edge contribution
        blend = shapH(iv1)*shapH(iv2)
        dblend(1:2) = dshapH(1:2,iv1)*shapH(iv2)  &
                    + shapH(iv1)*dshapH(1:2,iv2)
        xpar(1:2) = xpar(1:2) + xepar(1:2)*blend
        do j=1,2
          dxpardeta(1:2,j) = dxpardeta(1:2,j)                   &
                           + dxepardeta(1:2)*dsedeta(j)*blend   &
                           + xepar(1:2)*dblend(j)
        enddo
      enddo

!-----------------------------------------------------------------------
!
!  ...transform to cone Cartesian coordinates
      xrelc(1) = c*xpar(2)*cos(xpar(1))
      dxrelc(1,1:2) = c*dxpardeta(2,1:2)*cos(xpar(1)) &
                    - c*xpar(2)*sin(xpar(1))*dxpardeta(1,1:2)
      xrelc(2) = c*xpar(2)*sin(xpar(1))
      dxrelc(2,1:2) = c*dxpardeta(2,1:2)*sin(xpar(1)) &
                    + c*xpar(2)*cos(xpar(1))*dxpardeta(1,1:2)
      xrelc(3) = xpar(2)
      dxrelc(3,1:2) = dxpardeta(2,1:2)
!
!  ...transform to global Cartesian coordinates
      X(1:3) = 0.d0; Dxdeta(1:3,1:2) = 0.d0
      do i=1,3
        do j=1,3
          X(i) = X(i) + transf(j,i)*xrelc(j)
          Dxdeta(i,1:2) = Dxdeta(i,1:2) + transf(j,i)*dxrelc(j,1:2)
        enddo
        X(i) = X(i) + SURFACES(ns)%Rdata(i)
      enddo
      if (iprint.eq.1) then
        write(*,*) 'conic_triangle: X,Dxdeta = '
        do ivar=1,3
          write(*,7035) X(ivar),Dxdeta(ivar,1:2)
 7035     format(e12.5,3x,3e12.5)
        enddo
        call pause
      endif
!
!
end subroutine conic_triangle




!
!----------------------------------------------------------------------
!> Purpose : provide for a triangle on a
!!           a plane parametrized with Cylindrical coordinates
!!           This is a ``dirty''routine, as it works only
!!           for a very specific class of triangles
!!
!! @param[in]  No        - triangle number
!! @param[in]  Eta       - reference coordinates of a point
!! @param[out] X         - physical coordinates of the point
!! @param[out] Dxdeta    - derivatives wrt reference coordinates
!!
!! @revision Mar 11
!----------------------------------------------------------------------
!
subroutine PPwCC_triangle(No,Eta, X,Dxdeta)
!
      use control
      use GMP
      use element_data
#include "syscom.blk"
!
      dimension Eta(2), X(3),Dxdeta(3,2)
!
!  ...x coordinate, min and max radius
!      x0,rmin,rmax
!
!  ...endpoints, point in the parametric domain,
      dimension eta_v(2,3), etap(2)
!
!  ...derivatives
      dimension dr_deta(2), dtheta_deta(2)
!
      pi = acos(-1.d0)
      iprint=0
!
!  ...check surface
      ns = TRIANGLES(No)%Idata(1)
      if (SURFACES(ns)%Type.ne.'PPwCC') then
        write(*,*) '---------------------------------'
        write(*,*) 'PPwCC_triangle:'
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
 7001   format('PPwCC_triangle: x0,rmin,rmax,dr = ',4f8.3)
      endif
!
!  ...check radii consistency
      if ((rmin.le.GEOM_TOL).or.(rmax-rmin.le.GEOM_TOL)) then
        write(*,*) 'PPwCC_triangle: rmin,rmax = ',rmin,rmax
        stop 1
      endif
!
!  ...determine point coordinates in the parametric plane
      do iv=1,3
        nv = TRIANGLES(No)%VertNo(iv)
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
        case(2,3)
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
 7002     format('PPwCC_triangle: iv, eta(1:2,iv) = ',i2,2x,2f8.3)
        endif
      enddo
!
!  ...plane reference coordinates of the point
      etap(1:2) = eta_v(1:2,1)*(1.d0-Eta(1)-Eta(2)) &
                + eta_v(1:2,2)*Eta(1) + eta_v(1:2,3)*Eta(2)
      r = rmin + etap(2)*dr
      theta = 2.d0*pi*etap(1)
!
!  ...point coordinates
      X(1) = x0
      X(2) = r*sin(theta)
      X(3) = r*cos(theta)
!
!  ...derivatives
      Dxdeta(1,1:2) = 0.d0
      dr_deta(1) = (eta_v(2,2)-eta_v(2,1))*dr
      dr_deta(2) = (eta_v(2,3)-eta_v(2,1))*dr
      dtheta_deta(1) = (eta_v(1,2)-eta_v(1,1))*2.d0*pi
      dtheta_deta(2) = (eta_v(1,3)-eta_v(1,1))*2.d0*pi
      Dxdeta(2,1:2) = dr_deta(1:2)*sin(theta) + r*cos(theta)*dtheta_deta(1:2)
      Dxdeta(3,1:2) = dr_deta(1:2)*cos(theta) - r*sin(theta)*dtheta_deta(1:2)
!
!  ...printing statement
      if (iprint.eq.1) then
        write(*,*) 'PPwCC_triangle: X,Dxdeta = '
        do ivar=1,3
          write(*,7035) X(ivar),Dxdeta(ivar,1:2)
 7035     format(e12.5,3x,3e12.5)
        enddo
        call pause
      endif
!
end subroutine PPwCC_triangle
!


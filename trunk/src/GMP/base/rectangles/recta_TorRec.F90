!-----------------------------------------------------------------------
!
!   routine name       - recta_TorRec
!
!-----------------------------------------------------------------------
!
!   latest revision    - Dec 2015
!
!   purpose            - routine evaluates physical coordinates
!                        and its derivatives wrt to reference
!                        coordinates for a point on the image
!                        of a linear rectagle through a global system
!                        of coordinates: 
!                        x = Rmin*sin(\phi),
!                        y = (Rmaj+ Rmin*cos(\phi))cos(\theta),
!                        z = (Rmaj+ Rmin*cos(\phi))sin(\theta)
!                        and their  derivative wrt to reference
!                        coordinates.
!                        The torus centre lies in the origin, revolving
!                        about axis x. 
!                        Major radius is Rmaj, minor radius is Rmin.
!                        Idata(1:2) gives the point numbers for centers
!                        of arcs P1P2 and P3P4, respectively.
!
!   arguments :
!     in:
!               No     - a GMP rectagle number
!               Eta    - reference coordinates of a point
!                        in the rectangle
!     out:
!               X      - physical coordinates of the point
!               Dxdeta - derivatives of the physical coordinates wrt
!                        to the parameters
!
!
!-----------------------------------------------------------------------
!
   subroutine recta_TorRec(No,Eta, X,Dxdeta)
!
      use control
      use GMP          , only : RECTANGLES,POINTS,NDIM
      use node_types   , only : QUAD
      implicit none
!----------------------------------------------------------------------
      integer,                 intent(in)  :: No
      real(8), dimension(2),   intent(in)  :: Eta
      real(8), dimension(3),   intent(out) :: X(3)
      real(8), dimension(3,2), intent(out) :: Dxdeta
!----------------------------------------------------------------------
!  ...vertex shape functions
      real(8), dimension(4)   :: vshape
      real(8), dimension(2,4) :: dvshape
!  ...toroidal coordinates
      real(8)               :: phi,theta,phip,thetap,phitmp
      real(8), dimension(2) :: dphideta,dthetadeta
!----------------------------------------------------------------------
!     misc.
      integer :: iprint,iv,np,i
      real(8) :: pi,twopi,raux,px,costheta,sintheta,cosphi,sinphi
      real(8) :: theta12,theta34,rmaj,rmin
      real(8), dimension(3) :: c12,c34,oc,cv
!----------------------------------------------------------------------
!
      select case(No)
      case(6)
        iprint=0
      case default
        iprint=0
      end select
!
      if ((RECTANGLES(No)%Type.ne.'TorRec'.or.(NDIM.ne.3))) then
        write(*,7001) RECTANGLES(No)%Type
 7001   format('recta_TorRec: WRONG RECTANGLE TYPE = ',a10)
        stop 1
      endif
!
      if (iprint.eq.1) then
        write(*,7002) No,Eta
 7002   format('recta_TorRec: No,Eta = ',i4,2x,2f8.3)
      endif
!
!  ...initiate
      X(1:3) = 0.d0; Dxdeta(1:3,1:2) = 0.d0
      phi = 0.d0; dphideta(1:2) = 0.d0
      theta = 0.d0; dthetadeta(1:2) = 0.d0
!
      pi = acos(-1.d0)
      twopi = pi*2.d0
!
!  ...save coordinates of arcs' centers. 
!     c12 is the center of arc P1P2, c34 is the center of arc P3P4
      c12 = RECTANGLES(No)%Rdata(1:3)
      c34 = RECTANGLES(No)%Rdata(4:6)
!  ...check that coordinate x of both centers is zero
      if (abs(c12(1)).gt.GEOM_TOL .or. abs(c34(1)).gt.GEOM_TOL) then
        write(*,*) 'recta_TorRec: coord x of both arc centers must be 0'
        write(*,*) '              c12(1), c34(1) =', c12(1), c34(1)
        stop 1
      endif
!
!  ...determine major radius
      rmaj = NORM2(c12)
      raux = NORM2(c34)
      if (abs(rmaj - raux).gt.GEOM_TOL) then
        write(*,7003) rmaj,raux
 7003   format('recta_TorRec: INCONSISTENCY IN MAJOR RADIUS = ',2(f15.10))
        stop 1
      else 
        rmaj = 0.5*(rmaj+raux)
      endif
!
!  ...determine minor radius
!     take radius from c12 to p1
      np = RECTANGLES(No)%VertNo(1)
      rmin = NORM2(POINTS(np)%Rdata(1:3)-c12)
!     compare current value against distance from c12 to p2
      np = RECTANGLES(No)%VertNo(2)
      raux = NORM2(POINTS(np)%Rdata(1:3)-c12)
!     check consistency
      if (abs(rmin - raux).gt.GEOM_TOL) then
        write(*,7004) rmaj,raux
 7004   format('recta_TorRec: INCONSISTENCY IN MINOR RADIUS = ',2(f15.10))
        stop 1
      else 
        rmin = 0.5*(rmin+raux)
      endif
!     compare current value against distance from c34 to p3
      np = RECTANGLES(No)%VertNo(3)
      raux = NORM2(POINTS(np)%Rdata(1:3)-c34)
      if (abs(rmin - raux).gt.GEOM_TOL) then
        write(*,7004) rmin,raux
        stop 1
      else 
        rmin = 0.5*(rmin+raux)
      endif
!     compare current value against distance from c34 to p4
      np = RECTANGLES(No)%VertNo(4)
      raux = NORM2(POINTS(np)%Rdata(1:3)-c34)
      if (abs(rmin - raux).gt.GEOM_TOL) then
        write(*,7004) rmin,raux
        stop 1
      else 
        rmin = 0.5*(rmin+raux)
      endif
! 
!  ...determine theta angles of both arcs
      if (abs(c12(2)).gt.GEOM_TOL) then
        theta12 = ATAN(c12(3)/c12(2))
      else
        theta12 = pi/2.d0 * SIGN(1.d0,c12(3))
      endif 
      if (c12(2).lt. -GEOM_TOL) then
        theta12 = theta12 + pi
      endif
      if (abs(c34(2)).gt.GEOM_TOL) then
        theta34 = ATAN(c34(3)/c34(2))
      else
        theta34 = pi/2.d0 * SIGN(1.d0,c34(3))
      endif 
      if (c34(2).lt. -GEOM_TOL) then
        theta34 = theta34 + pi
      endif
      if (theta34.le.theta12) then 
        theta12 = theta12 - twopi
      endif
      ! theta12 = 0.d0; theta34 = pi/2.d0
!
!  ...vertex shape functions
      call vshape2(QUAD,Eta, vshape,dvshape)
!
!  ...compute the cylindrical coordinates of the endpoints
      do iv=1,4
        np = RECTANGLES(No)%VertNo(iv)
        px = POINTS(np)%Rdata(1)
        select case(iv)
!    ...select appropriate vector origin--arc center
        case(1,2)
          oc = c12
        case(3,4)
          oc = c34
        end select
!    ...store vector arc center--point iv
        cv = POINTS(np)%Rdata(1:3)-oc
!    ...evaluate dot product between major radial unit vector and 
!       minor radial unit vector.
        call dot_product(oc/rmaj,cv/rmin,raux)
!    ...Get phi. If px is negative, phi must be corrected
        if (px.ge.0.d0) then
          phitmp = ACOS(raux)
        else
          phitmp = twopi - ACOS(raux)
        endif
!
!  .....set theta and, if necessary, adjust phi
        select case(iv)
        case(1)
          if (iprint.eq.2) write(*,*) 'oc,cv,raux=',oc,cv,raux
          thetap = theta12
          phip = phitmp
        case(2)
          thetap = theta12
!      ...if phi2 is less than phi1, adjust by adding 2pi
          if (phitmp.lt.phip) then
            phip = phitmp + twopi
          else
            phip = phitmp
          endif
        case(3)
          thetap = theta34
!      ...if phi3 is twopi away from phi2, adjust by equaling phi3 to phi2
          if (abs(abs(phip-phitmp)-twopi).lt.GEOM_TOL) then
            phip = phip
          else
            phip = phitmp
          endif
        case(4)
          thetap = theta34
!      ...if phi4 is more than phi1, adjust by subtracting 2pi
          if (phitmp.gt.phip) then
            phip = phitmp - twopi
          else
            phip = phitmp
          endif
        end select
!    ...evaluate phi and theta according to vertex shape functions
        phi   = phi   + phip*vshape(iv)
        theta = theta + thetap*vshape(iv)
!    ...evaluate derivatives as well (w.r.t. eta)
        dphideta(1:2)   = dphideta(1:2)   + phip*dvshape(1:2,iv)
        dthetadeta(1:2) = dthetadeta(1:2) + thetap*dvshape(1:2,iv)

        if (iprint.eq.1) then
          write(*,*) 'recta_TorRec: iv, phip,thetap=',iv,phip,thetap
        endif

      enddo

      if (iprint.eq.1) then
        write(*,*) 'recta_TorRec: shape functions grad'
        write(*,*) '              dvshape(1) =',dvshape(:,1)
        write(*,*) '              dvshape(2) =',dvshape(:,2)
        write(*,*) '              dvshape(3) =',dvshape(:,3)
        write(*,*) '              dvshape(4) =',dvshape(:,4)
        write(*,*) 'recta_TorRec: completed'
        write(*,*) '              phi        =',phi
        write(*,*) '              dphideta   =',dphideta
        write(*,*) '              theta      =',phi
        write(*,*) '              dthetadeta =',dthetadeta
      endif
!
!  ...find Cartesian coordinates
      costheta = COS(theta); sintheta = SIN(theta)
      cosphi = COS(phi); sinphi = SIN(phi)
      X(1) = rmin*sinphi
      X(2) = (rmaj+rmin*cosphi)*costheta
      X(3) = (rmaj+rmin*cosphi)*sintheta
      Dxdeta(1,1:2) =  rmin*cosphi*dphideta(1:2)
      Dxdeta(2,1:2) = -rmin*costheta*sinphi*dphideta(1:2)   &
                      -X(3)*dthetadeta(1:2)
      Dxdeta(3,1:2) = -rmin*sintheta*sinphi*dphideta(1:2)   &
                      +X(2)*dthetadeta(1:2)
      if (iprint.eq.1) then
        write(*,*) 'theta = ', theta
        write(*,*) 'phi = ', phi
        write(*,7005) X
        write(*,7006) (Dxdeta(i,1:2),i=1,3)
 7005   format('recta_TorRec: ',/,'X      = ',3f8.3)
 7006   format('Dxdeta = ',2f8.3,/,  &
               '         ',2f8.3,/,  &
               '         ',2f8.3)
        write(*,*) 'c12 = ', c12
        write(*,*) 'c34 = ', c34
        write(*,*) 'rmaj, rmin = ',rmaj, rmin
        write(*,*) 'theta12, theta34 = ', theta12, theta34
        ! call pause
      endif
!
   end subroutine recta_TorRec


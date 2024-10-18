!-----------------------------------------------------------------------
!
!   routine name       - recta_CylRec
!
!-----------------------------------------------------------------------
!
!   latest revision    - Dec 2015
!
!   purpose            - routine evaluates physical coordinates
!                        and its derivatives wrt to reference
!                        coordinates for a point on the image
!                        of a linear rectangle through a global system
!                        of coordinates: x,y=rcos(\theta),z=rsin(\theta)
!                        and their  derivative wrt to reference
!                        coordinates
!
!   arguments :
!     in:
!               No     - a GMP rectangle number
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
   subroutine recta_CylRec(No,Eta, X,Dxdeta)
!
      use control
      use GMP          , only : RECTANGLES,POINTS,NDIM
      use node_types   , only : QUAD
      implicit none
!----------------------------------------------------------------------
      integer,                 intent(in)  :: No
      real(8), dimension(3),   intent(in)  :: Eta(2)
      real(8), dimension(3),   intent(out) :: X(3)
      real(8), dimension(3,3), intent(out) :: Dxdeta(3,3)
!----------------------------------------------------------------------
!  ...vertex shape functions
      real(8), dimension(4)   :: vshape
      real(8), dimension(3,4) :: dvshape
!  ...cylindrical coordinates
      real(8)               :: r,theta,rp,thetap,thetaTmp
      real(8), dimension(3) :: drdeta,dthetadeta
!----------------------------------------------------------------------
!     misc.
      integer :: iv,np
      real(8) :: pi,twopi,costhet,sinthet
!
#if HP3D_DEBUG
      integer :: i,iprint
      iprint=0
#endif
!----------------------------------------------------------------------
!
#if HP3D_DEBUG
      select case(No)
      case(1)
        iprint=0
      case default
        iprint=0
      end select
#endif
!
      if ((RECTANGLES(No)%Type.ne.'CylRec'.or.(NDIM.ne.3))) then
        write(*,7001) RECTANGLES(No)%Type
 7001   format('recta_CylRec: WRONG RECTANGLE TYPE = ',a10)
        stop 1
      endif
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7002) No,Eta
 7002   format('recta_CylRec: No,Eta = ',i4,2x,2f8.3)
      endif
#endif
!
!  ...initiate
      X(1:3) = 0.d0; Dxdeta(1:3,1:3) = 0.d0
      r = 0.d0; drdeta(1:3) = 0.d0
      theta = 0.d0; dthetadeta(1:3) = 0.d0
!
      pi = acos(-1.d0)
      twopi = pi*2.d0
!
!  ...vertex shape functions
      call vshape2(QUAD,Eta, vshape,dvshape)
!
!  ...compute the cylindrical coordinates of the endpoints
      do iv=1,4
        np = RECTANGLES(No)%VertNo(iv)
        !  x-coord
        X(1) = X(1) + POINTS(np)%Rdata(1)*vshape(iv)
        Dxdeta(1,1:3) = Dxdeta(1,1:3)  &
                      + POINTS(np)%Rdata(1)*dvshape(1:3,iv)
        !  y,z-coords
        call coord_cart2polar(POINTS(np)%Rdata(2:3), rp,thetaTmp)
!
!  .....adjust the angle, if necessary
        select case(iv)
        case(1)
          thetap = thetaTmp
        case default
          if (thetaTmp-thetap.gt.pi) then
            thetap = thetaTmp - twopi
          elseif (thetaTmp-thetap.lt.-pi) then
            thetap = thetaTmp + twopi
          else
            thetap = thetaTmp
          endif
        end select
        ! write(*,*) 'thetap = ', thetap
        r     = r     + rp*vshape(iv)
        theta = theta + thetap*vshape(iv)
        drdeta(1:3)     = drdeta(1:3)     + rp*dvshape(1:3,iv)
        dthetadeta(1:3) = dthetadeta(1:3) + thetap*dvshape(1:3,iv)
      enddo
!
      costhet = cos(theta); sinthet = sin(theta)
      X(2) = r*costhet
      X(3) = r*sinthet
      Dxdeta(2,1:3) = drdeta(1:3)*costhet - r*sinthet*dthetadeta(1:3)
      Dxdeta(3,1:3) = drdeta(1:3)*sinthet + r*costhet*dthetadeta(1:3)
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) 'theta = ', theta
        write(*,*) 'r = ', r
        write(*,7003) X
        write(*,7004) (Dxdeta(i,1:3),i=1,3)
 7003   format('recta_CylRec: ',/,'X      = ',3f8.3)
 7004   format('Dxdeta = ',3f8.3,/,  &
               '         ',3f8.3,/,  &
               '         ',3f8.3)
        call pause
      endif
#endif
!
!
   end subroutine recta_CylRec


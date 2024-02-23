c-----------------------------------------------------------------------
c
c   routine name       - curve_CylCoord
c
c-----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - routine evaluates physical coordinates
c                        for a point on the image of a straight line
c                        segment through a global cylicrical system
c                        of coordinates: x,y=rcos(\theta),z=rsin(\theta)
c                        and their  derivative wrt to reference
c                        coordinate
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
c-----------------------------------------------------------------------
c
      subroutine curve_CylCoord(No,Eta, X,Dxdeta)
c
      use GMP
c
      implicit none
c
      integer :: No
      real(8) :: Eta,X(3),Dxdeta(3)
c
c  ...cylindrical coordinates of the endpoints of the curve
      real(8) :: xp(3,2)
c
      integer :: iv,np
      real(8) :: costhet,dr_deta,dthet_deta,pi,r,sinthet,theta,twopi
c
      integer :: iprint
      iprint=0
c
      if (iprint.eq.1) then
        write(*,7000) No,Eta
 7000   format('curve_CylCoord: DEBUGGING FOR No,Eta = ',i4,2x,f8.3)
      endif
c
      pi = acos(-1.d0)
      twopi = pi*2.d0
c
c  ...check the type
      if (CURVES(No)%Type.ne.'CylCoord') then
        write(*,7001) No, CURVES(No)%Type
 7001   format('curve_CylCoord: WRONG CURVE TYPE, No, Type = ',
     .         i4,a10)
        stop 1
      endif
c
c  ...compute the cylindrical coordinates of the endpoints
      do iv=1,2
        np = CURVES(No)%EndPoNo(iv)
        xp(1,iv) = POINTS(np)%Rdata(1)
        call coord_cart2polar(POINTS(np)%Rdata(2:3), r,theta)
        xp(2,iv) = r
c
c  .....adjust the angle, if necessary
        select case(iv)
        case(2)
          if (theta-xp(3,1).gt.pi) theta = theta - twopi
          if (theta-xp(3,1).lt.-pi) theta = theta + twopi
        end select
        xp(3,iv) = theta
      enddo
c
c  ...interpolate
      X(1)  = (1.d0-Eta)*xp(1,1) + Eta*xp(1,2)
      r     = (1.d0-Eta)*xp(2,1) + Eta*xp(2,2)
      theta = (1.d0-Eta)*xp(3,1) + Eta*xp(3,2)
      dr_deta    = xp(2,2) - xp(2,1)
      dthet_deta = xp(3,2) - xp(3,1)
      costhet = cos(theta); sinthet = sin(theta)
      X(2) = r*costhet
      X(3) = r*sinthet
      Dxdeta(1) = xp(1,2) - xp(1,1)
      Dxdeta(2) = dr_deta*costhet - r*sinthet*dthet_deta
      Dxdeta(3) = dr_deta*sinthet + r*costhet*dthet_deta
c
      if (iprint.eq.1) then
        write(*,7003) No,X,Dxdeta
 7003   format('curve_CylCoord: No,X,Dxdeta = ', i4,2x,2(3e12.5,2x))
        call pause
      endif
c
c
      end subroutine curve_CylCoord


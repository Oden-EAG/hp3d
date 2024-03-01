!-----------------------------------------------------------------------
!
!   routine name       - curve_CylCoord
!
!-----------------------------------------------------------------------
!
!   latest revision    - Mar 2023
!
!   purpose            - routine evaluates physical coordinates
!                        for a point on the image of a straight line
!                        segment through a global cylicrical system
!                        of coordinates: x,y=rcos(\theta),z=rsin(\theta)
!                        and their  derivative wrt to reference
!                        coordinate
!
!   arguments :
!     in:
!               No     - the curve number
!               Eta    - reference coordinate  (between 0 and 1)
!     out:
!               X      - physical coordinates of the point
!               Dxdeta - derivatives of the physical coordinates wrt
!                        reference coordinate
!
!-----------------------------------------------------------------------
!
      subroutine curve_CylCoord(No,Eta, X,Dxdeta)
!
      use GMP
!
      implicit none
!
      integer :: No
      real(8) :: Eta,X(3),Dxdeta(3)
!
!  ...cylindrical coordinates of the endpoints of the curve
      real(8) :: xp(3,2)
!
      integer :: iv,np
      real(8) :: costhet,dr_deta,dthet_deta,pi,r,sinthet,theta,twopi
!
      integer :: iprint
      iprint=0
!
      if (iprint.eq.1) then
        write(*,7000) No,Eta
 7000   format('curve_CylCoord: DEBUGGING FOR No,Eta = ',i4,2x,f8.3)
      endif
!
      pi = acos(-1.d0)
      twopi = pi*2.d0
!
!  ...check the type
      if (CURVES(No)%Type.ne.'CylCoord') then
        write(*,7001) No, CURVES(No)%Type
 7001   format('curve_CylCoord: WRONG CURVE TYPE, No, Type = ', &
               i4,a10)
        stop 1
      endif
!
!  ...compute the cylindrical coordinates of the endpoints
      do iv=1,2
        np = CURVES(No)%EndPoNo(iv)
        xp(1,iv) = POINTS(np)%Rdata(1)
        call coord_cart2polar(POINTS(np)%Rdata(2:3), r,theta)
        xp(2,iv) = r
!
!  .....adjust the angle, if necessary
        select case(iv)
        case(2)
          if (theta-xp(3,1).gt.pi) theta = theta - twopi
          if (theta-xp(3,1).lt.-pi) theta = theta + twopi
        end select
        xp(3,iv) = theta
      enddo
!
!  ...interpolate
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
!
      if (iprint.eq.1) then
        write(*,7003) No,X,Dxdeta
 7003   format('curve_CylCoord: No,X,Dxdeta = ', i4,2x,2(3e12.5,2x))
        call pause
      endif
!
!
      end subroutine curve_CylCoord


c-----------------------------------------------------------------------
c
c   routine name       - trian_CylTri
c
c-----------------------------------------------------------------------
c
c   latest revision    - Dec 15
c
c   purpose            - routine evaluates physical coordinates
c                        and its derivatives wrt to reference
c                        coordinates for a point on the image
c                        of a linear triangle through a global system
c                        of coordinates: x,y=rcos(\theta),z=rsin(\theta)
c                        and their  derivative wrt to reference
c                        coordinates
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
c
c-----------------------------------------------------------------------
c
      subroutine trian_CylTri(No,Eta, X,Dxdeta)
c
      use GMP
#include "syscom.blk"
c
      dimension Eta(2),X(NDIM),Dxdeta(NDIM,2)
c
c  ...cylindrical coordinates of the endpoints of the triangle
      dimension xp(3,3)
c
      iprint=0
c
      if ((TRIANGLES(No)%Type.ne.'CylTri'.or.(NDIM.ne.3))) then
        write(*,7001) TRIANGLES(No)%Type
 7001   format('tria_CylTri: WRONG TRIANGLE TYPE = ',a10)
        stop 1
      endif
c
      if (iprint.eq.1) then
        write(*,7002) No,Eta
 7002   format('trian_CylTri: No,Eta = ',i4,2x,2f8.3)
      endif
c
      pi = acos(-1.d0)
      twopi = pi*2.d0
c
c  ...compute the cylindrical coordinates of the endpoints
      do iv=1,3
        np = TRIANGLES(No)%VertNo(iv)
        xp(1,iv) = POINTS(np)%Rdata(1)
        call coord_cart2polar(POINTS(np)%Rdata(2:3), r,theta)
        xp(2,iv) = r
c
c  .....adjust the angle, if necessary
        select case(iv)
        case(2,3)
          if (theta-xp(3,1).gt.pi) theta = theta - twopi
          if (theta-xp(3,1).lt.-pi) theta = theta + twopi
        end select
        xp(3,iv) = theta
      enddo
c
c  ...interpolate
      X(1)  = (1.d0-Eta(1)-Eta(2))*xp(1,1)
     .      + Eta(1)*xp(1,2)+ Eta(2)*xp(1,3)
      r     = (1.d0-Eta(1)-Eta(2))*xp(2,1)
     .      + Eta(1)*xp(2,2)+ Eta(2)*xp(2,3)
      theta = (1.d0-Eta(1)-Eta(2))*xp(3,1)
     .      + Eta(1)*xp(3,2)+ Eta(2)*xp(3,3)
      dr_deta1    = xp(2,2) - xp(2,1)
      dr_deta2    = xp(2,3) - xp(2,1)
      dthet_deta1 = xp(3,2) - xp(3,1)
      dthet_deta2 = xp(3,3) - xp(3,1)
      costhet = cos(theta); sinthet = sin(theta)
      X(2) = r*costhet
      X(3) = r*sinthet
      Dxdeta(1,1) = xp(1,2) - xp(1,1)
      Dxdeta(1,2) = xp(1,3) - xp(1,1)
      Dxdeta(2,1) = dr_deta1*costhet - r*sinthet*dthet_deta1
      Dxdeta(2,2) = dr_deta2*costhet - r*sinthet*dthet_deta2
      Dxdeta(3,1) = dr_deta1*sinthet + r*costhet*dthet_deta1
      Dxdeta(3,2) = dr_deta2*sinthet + r*costhet*dthet_deta2
c
c
      end subroutine trian_CylTri


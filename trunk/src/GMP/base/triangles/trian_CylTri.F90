!-----------------------------------------------------------------------
!
!   routine name       - trian_CylTri
!
!-----------------------------------------------------------------------
!
!   latest revision    - Mar 2023
!
!   purpose            - routine evaluates physical coordinates
!                        and its derivatives wrt to reference
!                        coordinates for a point on the image
!                        of a linear triangle through a global system
!                        of coordinates: x,y=rcos(\theta),z=rsin(\theta)
!                        and their  derivative wrt to reference
!                        coordinates
!
!   arguments :
!     in:
!               No     - a GMP triangle number
!               Eta    - reference coordinates of a point
!                        in the triangle
!     out:
!               X      - physical coordinates of the point
!               Dxdeta - derivatives of the physical coordinates wrt
!                        to the parameters
!
!-----------------------------------------------------------------------
!
   subroutine trian_CylTri(No,Eta, X,Dxdeta)
!
      use GMP
!
      implicit none
!
      integer :: No
      real(8) :: Eta(2),X(NDIM),Dxdeta(NDIM,2)
!
!  ...cylindrical coordinates of the endpoints of the triangle
      real(8) :: xp(3,3)
!
      integer :: iv,np
      real(8) :: costhet,dr_deta1,dr_deta2,dthet_deta1,dthet_deta2
      real(8) :: pi,twopi,r,sinthet,theta
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
      if ((TRIANGLES(No)%Type.ne.'CylTri'.or.(NDIM.ne.3))) then
        write(*,7001) TRIANGLES(No)%Type
 7001   format('tria_CylTri: WRONG TRIANGLE TYPE = ',a10)
        stop 1
      endif
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7002) No,Eta
 7002   format('trian_CylTri: No,Eta = ',i4,2x,2f8.3)
      endif
#endif
!
      pi = acos(-1.d0)
      twopi = pi*2.d0
!
!  ...compute the cylindrical coordinates of the endpoints
      do iv=1,3
        np = TRIANGLES(No)%VertNo(iv)
        xp(1,iv) = POINTS(np)%Rdata(1)
        call coord_cart2polar(POINTS(np)%Rdata(2:3), r,theta)
        xp(2,iv) = r
!
!  .....adjust the angle, if necessary
        select case(iv)
        case(2,3)
          if (theta-xp(3,1).gt.pi) theta = theta - twopi
          if (theta-xp(3,1).lt.-pi) theta = theta + twopi
        end select
        xp(3,iv) = theta
      enddo
!
!  ...interpolate
      X(1)  = (1.d0-Eta(1)-Eta(2))*xp(1,1) &
            + Eta(1)*xp(1,2)+ Eta(2)*xp(1,3)
      r     = (1.d0-Eta(1)-Eta(2))*xp(2,1) &
            + Eta(1)*xp(2,2)+ Eta(2)*xp(2,3)
      theta = (1.d0-Eta(1)-Eta(2))*xp(3,1) &
            + Eta(1)*xp(3,2)+ Eta(2)*xp(3,3)
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
!
!
   end subroutine trian_CylTri

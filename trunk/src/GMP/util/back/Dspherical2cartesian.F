!*****************************************************************************
      SUBROUTINE Dspherical2cartesian(rho, theta, phi, r, Dr)
!*****************************************************************************
      !*** LATEST REVISION: Set 08
      !
      !*** PURPOSE: routine evaluates the cartesian coordinates of a point
      !             from spherical coordinates; it also evaluates derivatives
      !             of cartesian coordinates wrt spherical ones.

!*****************************************************************************
      IMPLICIT NONE
!*****************************************************************************
      !*** dummy arguments
      double precision, intent(in) :: rho
      double precision, intent(in) :: theta
      double precision, intent(in) :: phi
      double precision, dimension(3), intent(out) :: r
      double precision, dimension(3,3), intent(out) :: Dr
!*****************************************************************************
      !
      r(1) = rho*sin(phi)*cos(theta)
      r(2) = rho*sin(phi)*sin(theta)
      r(3) = rho*cos(phi)
      !
      Dr(1,1) = sin(phi)*cos(theta)
      Dr(1,2) = -rho*sin(phi)*sin(theta)
      Dr(1,3) = rho*cos(phi)*cos(theta)
      Dr(2,1) = sin(phi)*sin(theta)
      Dr(2,2) = rho*sin(phi)*cos(theta)
      Dr(2,3) = rho*cos(phi)*sin(theta)
      Dr(3,1) = cos(phi)
      Dr(3,2) = 0.d0
      Dr(3,3) = -rho*sin(phi)
      !
      END SUBROUTINE Dspherical2cartesian








!************************************************************************
      SUBROUTINE spherical2cartesian(rho, theta, phi, r)
!************************************************************************
      !*** LATEST REVISION: Sep 08
      !
      !*** PURPOSE: given the spherical coordinates of a point routine
      !             evaluates the cartesian coodinates.
      !
      !*** REQUIRED PROCEDURES: none
!************************************************************************
      IMPLICIT NONE
!************************************************************************
      !*** DUMMY ARGUMENTS
      double precision, intent(in) :: rho
      double precision, intent(in) :: theta
      double precision, intent(in) :: phi
      double precision, dimension(3), intent(out) :: r
!************************************************************************
      !
      r(1) = rho*sin(phi)*cos(theta)
      r(2) = rho*sin(phi)*sin(theta)
      r(3) = rho*cos(phi)
      !
      END SUBROUTINE spherical2cartesian








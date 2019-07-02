!****************************************************************************
      SUBROUTINE cylindrical2cartesian(rho, theta, z, r)
!****************************************************************************
      !*** LATEST REVISION: Sep 08
      !
      !*** PURPOSE: routine evaluate the cartsian coordinates of a point, 
      !             given the cylindrical coordinates
      !
      !*** REQUIRED PROCEDURES: none
!****************************************************************************
      IMPLICIT NONE
!****************************************************************************
      double precision, intent(in) :: rho
      double precision, intent(in) :: theta
      double precision, intent(in) :: z
      double precision, dimension(3), intent(out) :: r  
!****************************************************************************
      r(1) = rho*cos(theta)
      r(2) = rho*sin(theta)
      r(3) = z
      !
      END SUBROUTINE cylindrical2cartesian

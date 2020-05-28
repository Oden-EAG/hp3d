!***********************************************************************************************************
      SUBROUTINE Dcylindrical2cartesian(rho, theta, z, r, Dr)
!***********************************************************************************************************
      !*** LATEST REVISION:     Aug 08
      !
      !*** PURPOSE:             given the cylindrical coordinates of a point, routine return the cartesian
      !                         coordinates and the derivatives of cartesian coordinates wrt cylindrical
      !
      !*** REQUIRED PROCEDURES: none
!***********************************************************************************************************
      IMPLICIT NONE
!***********************************************************************************************************
      !*** DUMMY ARGUMENTS
      double precision, intent(in) :: rho, theta, z            !*** cylindrical coordinates
      double precision, dimension(3), intent(out) :: r         !*** cartesian coordinates
      double precision, dimension(3,3), intent(out) :: Dr      !*** derivatives
!***********************************************************************************************************
      !
      r(1) = rho*cos(theta)
      r(2) = rho*sin(theta)
      r(3) = z
      !
      Dr(1,1) = cos(theta)
      Dr(1,2) = - rho*sin(theta)
      Dr(1,3) = 0.d0
      Dr(2,1) = sin(theta)
      Dr(2,2) = rho*cos(theta)
      Dr(2,3) = 0.d0
      Dr(3,1) = 0.d0
      Dr(3,2) = 0.d0
      Dr(3,3) = 1.d0
      !
      END SUBROUTINE Dcylindrical2cartesian

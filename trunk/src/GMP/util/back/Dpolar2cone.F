!**************************************************************************************************************
      SUBROUTINE Dpolar2cone(alpha, rho, theta, r, Dr)
!**************************************************************************************************************
      !*** LATEST REVISION:     Aug 08
      !
      !*** PURPOSE:             subroutine maps a circular sector of amplitude beta to a cone of half-aperture
      !                         alpha, vertex at the origin, and axis oriented as the z-axis. Derivatives wrt
      !                         rho and theta are also returned.
      !
      !*** REQUIRED PROCEDURES: none
      !
      !*** REMARKS:             angle alpha must belong to (0, pi); the relationship between alpha and beta is:
      !                         beta = 2pi*sin(alpha).
!**************************************************************************************************************
      IMPLICIT NONE
!**************************************************************************************************************
      !*** DUMMY ARGUMENTS
      double precision, intent(in) :: alpha                    !*** half-aperture for cone
      double precision, intent(in) :: rho, theta               !*** polar cordinates for a point
      double precision, dimension(3), intent(out) :: r         !*** cartesian coordinates for image point
      double precision, dimension(3,2), intent(out) :: Dr      !*** derivatives wrt rho and theta
      !
      !*** LOCAL VARIABLES
      double precision :: beta                                 !*** amplitude of circular sector
      double precision :: pi                                   !*** pi
!**************************************************************************************************************
      !
      pi = acos(-1.d0)
      !
      beta = pi*sin(alpha)                                     !*** get HALF amplitude of circular sector
      !
      if ((theta .le. -beta) .or. (theta .gt. beta)) then
        write(*,*)'--------------------------------'
        write(*,*)'Dpolar2cone:'
        write(*,*)'ERROR: invalid theta coordinate.'
        write(*,*)'--------------------------------'
        stop
      end if
      !
      r(1) = rho*sin(alpha)*cos(theta/sin(alpha))              !*** compute (x,y,z)
      r(2) = rho*sin(alpha)*sin(theta/sin(alpha))
      r(3) = rho*cos(alpha)
      !
      Dr(1,1) = sin(alpha)*cos(theta/sin(alpha))               !*** compute d(x,y,z)/d(rho,theta)
      Dr(2,1) = sin(alpha)*sin(theta/sin(alpha))
      Dr(3,1) = cos(alpha)
      Dr(1,2) = -rho*sin(theta/sin(alpha))
      Dr(2,2) = rho*cos(theta/sin(alpha))
      Dr(3,2) = 0.d0
      !
      END SUBROUTINE Dpolar2cone

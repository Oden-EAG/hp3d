!**************************************************************************************************************
      SUBROUTINE polar2cone(alpha, rho, theta, r)
!**************************************************************************************************************
      !*** LATEST REVISION:     Aug 08
      !
      !*** PURPOSE:             subroutine maps a circular sector of amplitude beta to a cone of half-aperture
      !                         alpha, vertex at the origin, and axis oriented as the z-axis.
      !
      !*** REQUIRED PROCEDURES: none
      !
      !*** REMARKS:             angle alpha must belong to (0, pi); the relationship between alpha and beta is:
      !                         beta = 2pi*sin(alpha).
!**************************************************************************************************************
      IMPLICIT NONE
!**************************************************************************************************************
      !*** DUMMY ARGUMENTS
      double precision, intent(in) :: alpha               !*** half-aperture for cone
      double precision, intent(in) :: rho, theta          !*** polar cordinates for a point
      double precision, dimension(3), intent(out) :: r    !*** cartesian coordinates for image point
      !
      !*** LOCAL VARIABLES
      double precision :: beta                            !*** amplitude of circular sector
      double precision :: pi                              !*** pi
!**************************************************************************************************************
      !
      pi = acos(-1.d0)
      !
      beta = pi*sin(alpha)                                     !*** get HALF amplitude of circular sector
      !
      if ((theta .le. -beta) .OR. (theta .gt. beta)) then
        write(*,*)'--------------------------------'
        write(*,*)'polar2cone:'
        write(*,*)'ERROR: invalid theta coordinate.'
        write(*,*)'--------------------------------'
        stop
      end if
      !
      r(1) = rho*sin(alpha)*cos(theta/sin(alpha))
      r(2) = rho*sin(alpha)*sin(theta/sin(alpha))
      r(3) = rho*cos(alpha)
      !
      END SUBROUTINE polar2cone


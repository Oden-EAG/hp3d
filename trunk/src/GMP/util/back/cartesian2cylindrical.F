!************************************************************************************************************************      
      SUBROUTINE cartesian2cylindrical(r, rho ,theta, z)
!************************************************************************************************************************
      !*** LATEST REVISION:     Aug 08
      !
      !*** PURPOSE:             given a point with cartesian coordinates (r(1), r(2), r(3)), routines returns cylindrical
      !                         coordinates (rho, theta, z) for the point.
      !
      !*** REQUIRED PROCEDURES: none
      !
      !*** REMARKS:             angle theta belongs to (-pi, pi].
!************************************************************************************************************************
      IMPLICIT NONE
!************************************************************************************************************************
      !*** DUMMY ARGUMENTS
      double precision, dimension(3), intent(in) ::  r      !*** cartesian coordinates for a point
      double precision, intent(out) :: rho, theta, z        !*** cylindrical coordinate for a point
      !
      !*** LOCAL VARIABLES
      double precision :: pi                                !*** pi
      integer :: iverb                                      !*** verbosity flag: 0, 1
      !
      !*** PARAMETERS
      double precision, parameter ::  eps = 1.d-13          !*** geometrical tolerance
!************************************************************************************************************************
      !
      iverb = 0
      !
      rho = sqrt(r(1)**2 + r(2)**2)                                   !*** compute rho
      !
      if (iverb .ge. 1) then
        write(*,1000) rho
 1000   format (1X, 'cartesian2cylindrical: rho = ', F7.3)
      end if
      !
      if (rho .lt. eps) then                                          !*** compute theta
        theta = 0.d0
        write(*,*) 'cartesian2cylindrical: WARNING:
     . theta is undetermined.'
      else
        if (r(2) .eq. 0.d0) then
          if (r(1) .lt. 0.d0) then
            pi = acos(-1.d0)
            theta = pi
          else if (r(1) .gt. 0.d0) then
            theta = 0.d0
          end if
      !
        elseif (r(2) .gt. 0.d0) then
          theta = acos(r(1)/rho)
      !
        else
          theta = - acos(r(1)/rho)
        endif
      endif
      !
      if (iverb .ge. 1) then
        write(*,1001) theta
 1001   format (1X, 'cartesian2cylindrical: theta = ', F5.3)
      end if
      !
      z = r(3)                                                       !*** compute z
      !
      if (iverb .ge. 1) then
        write(*,1002) z
 1002   format (1X, 'cartesian2cylindrical: z = ', F7.3)
      end if
      !
      END SUBROUTINE cartesian2cylindrical

!**************************************************************************************************************************
      SUBROUTINE cone2polar(alpha, r, rho, theta)
!**************************************************************************************************************************
      !*** LATEST REVISION: Sep 08
      !
      !*** PURPOSE: let us consider a cone of half-aperture alpha, vertex at the origin, and axis coinciding
      !             with the z-axis; when we unfold such a surface on a plane we obatain a cicular sector of amplitude
      !             2pi/sin(alpha). Given a point on the cone, routine returns the polar coordinates of the image
      !             point when the cylinder is unfolded on the plane. For mapping the circular sector onto the cone
      !             cf. subroutine POLAR2CONE.
      !
      !*** REQUIRED PROCEDURES: none
!**************************************************************************************************************************
      IMPLICIT NONE
!**************************************************************************************************************************
      !*** DUMMY ARGUMENTS
      double precision, intent(in) :: alpha               !*** half aperture for cone
      double precision, dimension(3), intent(in) :: r     !*** cartesian coordinates in 3D space
      double precision, intent(out) :: rho, theta         !*** polar coordinates in 2D space
      !
      !*** LOCAL VARIABLES
      double precision :: pi                              !*** pi
      double precision :: x                               !*** auxiliary variable
      integer :: iverb                                    !*** verbosity flag 0,1
      !
      !*** PARAMETERS
      double precision, parameter :: eps = 1.d-13         !*** geometry tolerance
!**************************************************************************************************************************
      iverb = 0
      !
      rho = sqrt(r(1)**2 + r(2)**2 + r(3)**2)                  !*** compute rho
      !
      if (iverb .ge. 1) then
        write(*,1000) rho
 1000   format (1X, 'cone2polar: rho = ', F7.3)
      end if
      !
      if (rho .lt. eps) then                                   !*** compute theta
        theta = 0.d0
        write(*,*) 'cone2polar: WARNING:
     . theta is undetermined.'
      else
        if (r(2) .eq. 0.d0) then                               !*** IF on x-axis
          if (r(1) .lt. 0.d0) then
            pi = acos(-1.d0)
            theta = pi
          else if (r(1) .gt. 0.d0) then
            theta = 0.d0
          end if
      !
        elseif (r(2) .gt. 0.d0) then                           !*** ELSE IF strictly above x-axis
          theta = acos(r(1)/sqrt(r(1)**2 + r(2)**2))
      !
        else                                                   !*** ELSE strictly below x-axis
          theta = - acos(r(1)/sqrt(r(1)**2 + r(2)**2))
        endif
      endif
      !
      theta = theta*sin(alpha)                                 !*** rescale theta to interval (-pi*sin(alpha), pi*sin(alpha)]
      !
      if (iverb .ge. 1) then
        write(*,1001) theta
 1001   format (1X, 'cone2polar: theta = ', F5.3)
      end if
      !
      END SUBROUTINE cone2polar

!****************************************************************************************************************
      SUBROUTINE cartesian2spherical(r, rho, theta, phi)
!****************************************************************************************************************
      !*** LATEST REVISION: Sep 08
      !
      !*** PURPOSE: given a point with cartesian coordinates (r(1), r(2), r(3)), routines returns spherical
      !             coordiantes (rho, theta, phi) for the point; angle theta belongs to (-pi, pi] while angle
      !             phi belongs to [0, pi).
      !
      !*** REQUIRED PROCEDURES: none
!****************************************************************************************************************
      IMPLICIT NONE
!*****************************************************************************************************************
      !*** dummy arguments
      double precision, dimension(3), intent(in) :: r           !*** cartesian coordinates of a point
      double precision, intent(out) :: rho                      !*** radial distance
      double precision, intent(out) :: theta                    !*** longitude from xz plane [0,2pi)
      double precision, intent(out) :: phi                      !*** colatitude from z axis [0,pi)
      !
      !*** local variables
      double precision :: pi                                    !*** pi
      double precision :: rho_aux                               !*** auxiliary variable
      integer :: iverb                                          !*** verbosity flag: 0, 1, 2
      !
      !*** constants
      double precision, parameter :: eps = 1.d-13               !*** geometrical tolerance    
!*****************************************************************************************************************
      iverb = 0
      !
      rho = sqrt(r(1)**2 + r(2)**2 + r(3)**2)
      !
      if (iverb .ge. 1) then
        write(*,1000) rho
 1000   format(1X, 'cartesian2spherical: rho = ', F7.3)
      end if
      !
      if (rho .lt. eps) then                                         !*** IF less than geometrical tolerance
        theta = 0.d0
        phi = 0.d0
        write(*,*) 'cartesian2spherical: WARNING:
     . theta and phi are undetermined.'
      else                                                           !*** ELSE more than geometrical tolerance
        phi = acos(r(3)/rho)                                         !*** get phi
        rho_aux = sqrt(r(1)**2 + r(2)**2)                            !*** get rho_aux
      !
        if (iverb .ge. 2) then
          write(*,1003) rho_aux
 1003     format(1X, 'cartesian2spherical: rho_aux = ', F7.3)
        end if
      !
        if (rho_aux .lt. eps) then                                   !***** IF less than geometrical tolerance
          theta = 0.d0
          write(*,*) 'cartesian2spherical: WARNING:
     . theta is undetermined.'
        else if (r(2) .eq. 0.d0) then                                !***** ELSE IF on x-axis
          if (r(1) .lt. 0.d0) then
            pi = acos(-1.d0)                                         !*** get pi
            theta = pi
          else if (r(1) .gt. 0.d0) then
            theta = 0.d0
          end if
        else if (r(2) .gt. 0.d0) then                                !***** ELSE IF strictly above x-axis 
          theta = acos(r(1)/rho_aux)
        else                                                         !***** ELSE strictly below the x-axis
          theta = - acos(r(1)/rho_aux)              
        end if                                                       !***** END IF
      end if                                                         !*** END IF
      !
      if (iverb .ge. 1) then
        write(*,1001) theta
 1001   format(1X, 'cartesian2spherical: theta = ', F5.3)
        write(*,1002) phi
 1002   format(1X, 'cartesian2spherical: phi = ', F5.3)
      end if
      !
      END SUBROUTINE cartesian2spherical

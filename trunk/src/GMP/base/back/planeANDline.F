!*****************************************************************************************************
      SUBROUTINE planeANDline(n_plane, x_plane, n_line, x_line, r)
!*****************************************************************************************************
      !*** LATEST REVISION: Sep 08
      !
      !*** PURPOSE: routine evaluates the intersection of a plane and a line.
      !
      !*** REQUIRED PROCEDURES: scalar_product
!*****************************************************************************************************
      IMPLICIT NONE
!*****************************************************************************************************      
      !*** DUMMY ARGUMENTS
      double precision, dimension(3), intent(in) :: n_plane        !*** normal to plane
      double precision, dimension(3), intent(in) :: x_plane        !*** a point on the plane
      double precision, dimension(3), intent(in) :: n_line         !*** direction of line        
      double precision, dimension(3), intent(in) :: x_line         !*** a point on the line
      double precision, dimension(3), intent(out) :: r             !*** coordinates of point of intersection
      !
      !*** LOCAL VARIABLES
      double precision :: prod1, prod2, prod3
      double precision :: t
      double precision, parameter :: eps = 1.d-13                  !*** geometry tolerance
      integer :: i
!*****************************************************************************************************
      !
      call scalar_product(n_plane, n_line, prod1)
      !
      if (abs(prod1) .le. eps) then
        write(*,*) 'planeANDline: line is parallel to plane'
        stop
      end if
      !     
      call scalar_product(x_plane, n_plane, prod2)
      call scalar_product(x_line, n_plane, prod3) 
      !
      t = (prod2 - prod3)/prod1  
      !
      do i = 1, 3
        r(i) = x_line(i) + t*n_line(i)
      end do
      !                  
      END SUBROUTINE planeANDline

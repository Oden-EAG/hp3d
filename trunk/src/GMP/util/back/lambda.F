c---------------------------------------------------------------------
c
      double precision FUNCTION lambda(i, eta)
c
c---------------------------------------------------------------------
c     LATEST REVISION: Nov 08
c
c     PURPOSE:         routines evaluate baricentrical coordinates
c                      with respect to a 2D simplex (triangle)
c
c     IN:              i - baricentrical coordinate number
c                      eta - cartesian coordinates
c
c     OUT:             lambda - baricentrical coordinate
c---------------------------------------------------------------------
      IMPLICIT NONE
!*********************************************************************
      !*** dummmy arguments
      integer, intent(in) :: i
      double precision, dimension(2), intent(in) :: eta
!*********************************************************************

      select case (i)
        case (1)
	lambda = 1.d0 - eta(1) - eta(2)

        case (2)
	lambda = eta(1)

	case (3)
	lambda = eta(2)

	case default
	lambda = 0.d0
	write(*,*)'lambda: unknown baricentrical coordinate'
	stop

      end select

      END FUNCTION lambda

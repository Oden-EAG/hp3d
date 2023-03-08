c-----------------------------------------------------------------------
c
c   routine name       - trobs
c
c-----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - routine transforms cartesian coordinates
c                        of a point into observer's coordinates; the
c                        first two correspond to the projection plane,
c                        the third to the direction of projection;
c                        the third coordinate decides about the order
c                        in which objects are displayed
c
c   arguments :
c     in:
c              Xcart   - 3D cartesian coordinates of a point
c     out:
c              Xobs    - 2D observer's coordinates of a point
c
c-----------------------------------------------------------------------
c
      subroutine trobs(Xcart, Xobs)
c
      use graphmod, only: RMTR
c
      implicit none
c
      real(8), intent(in)  :: Xcart(3)
      real(8), intent(out) :: Xobs(3)
c
      real(8) :: s
      integer :: i,j
c
      do i=1,3
        s = 0.d0
        do j=1,3
          s = s + RMTR(i,j)*Xcart(j)
        enddo
        Xobs(i) = s
      enddo
c
      end subroutine trobs

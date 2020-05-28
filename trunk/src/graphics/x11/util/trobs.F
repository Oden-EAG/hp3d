c-----------------------------------------------------------------------
c
c   routine name       - trobs
c
c-----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Sep 07
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
c   required  routines - none
c
c-----------------------------------------------------------------------
c
      subroutine trobs(Xcart, Xobs)
      use graphmod
c
#include "syscom.blk"
c
      dimension Xcart(3),Xobs(3)
c
      do i=1,3
        s = 0.d0
        do j=1,3
          s = s + RMTR(i,j)*Xcart(j)
        enddo
        Xobs(i) = s
      enddo
c
c
      end subroutine trobs

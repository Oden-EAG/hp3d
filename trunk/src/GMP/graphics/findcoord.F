#if HP3D_USE_X11

c----------------------------------------------------------------------
c
c   routine name       - findcoord
c
c   latest revision    - Mar 2023
c
c   purpose            - given a triangle or rectangle number and a point
c                        both  projected onto a plane, find the correspon-
c                        ding master element coordinates of the point
c                        wrt to the figure
c
c   arguments :
c     in:
c                No    - triangle or rectangle number
c                Lab   = 1 triangle
c                      = 2 rectangle
c                Xp    - coordinates of a projected point
c
c     out:
c                Eta   - master element coordinates of the point
c                X     - physical coordinates of the point (before
c                        the projection)
c                Nfl   = 0, if the procedure has converged
c                        1  otherwise
c
c----------------------------------------------------------------------
c
      subroutine findcoord(No,Lab,Xp, Eta,X,Nfl)
c
      use control
      use GMP
c
      implicit none
c
      integer :: No,Lab,Nfl
      real(8) :: Xp(2),Eta(2),X(3)
c
c  ...derivatives of the physical coordinates wrt reference coordinates
      real(8) :: dxdeta(3,2)
c
c  ...work space for NR iterations
      real(8) :: xy(2),aa(2,2),bb(2),deta(2)
c
      real(8) :: d
      integer :: iter
c
      integer :: iprint
      iprint=0
c
c-----------------------------------------------------------------------
c
      select case(Lab)
      case(1)
c
        write(*,*) 'findcoord: UNFINISHED'
        stop 1
c
c-----------------------------------------------------------------------
c
      case(2)
c
c  .....start with the midpoint
        Eta(1:2) = 0.d0
c
c  .....Newton-Raphson iterations
        do iter=1,10
c
          call recta_linear(No,Eta,X,Dxdeta)
          call trobs(X, xy)
c
c  .......residual
          bb(1:2) = Xp(1:2) - xy(1:2)
          d = sqrt(bb(1)**2+bb(2)**2)
          if (d.lt.GEOM_TOL) then
            Nfl = 1
            return
          endif
c
c  .......tangent matrix
          call trobs(Dxdeta(1:3,1), aa(1:2,1))
          call trobs(Dxdeta(1:3,2), aa(1:2,2))
c
c  .......solve for the increment
          call gausse(aa,2,bb, deta,2)
c
c  .......update
          Eta(1:2) = Eta(1:2) + deta(1:2)
        enddo
c
        Nfl = 0
        Eta(1:2) = 0.d0
        X(1:3) = 0.d0
c
      end select
c
      end subroutine findcoord

#endif

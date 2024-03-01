#if HP3D_USE_X11

!----------------------------------------------------------------------
!
!   routine name       - findcoord
!
!   latest revision    - Mar 2023
!
!   purpose            - given a triangle or rectangle number and a point
!                        both  projected onto a plane, find the correspon-
!                        ding master element coordinates of the point
!                        wrt to the figure
!
!   arguments :
!     in:
!                No    - triangle or rectangle number
!                Lab   = 1 triangle
!                      = 2 rectangle
!                Xp    - coordinates of a projected point
!
!     out:
!                Eta   - master element coordinates of the point
!                X     - physical coordinates of the point (before
!                        the projection)
!                Nfl   = 0, if the procedure has converged
!                        1  otherwise
!
!----------------------------------------------------------------------
!
   subroutine findcoord(No,Lab,Xp, Eta,X,Nfl)
!
      use control
      use GMP
!
      implicit none
!
      integer :: No,Lab,Nfl
      real(8) :: Xp(2),Eta(2),X(3)
!
!  ...derivatives of the physical coordinates wrt reference coordinates
      real(8) :: dxdeta(3,2)
!
!  ...work space for NR iterations
      real(8) :: xy(2),aa(2,2),bb(2),deta(2),aux(3,2)
!
      real(8) :: d
      integer :: iter
!
      integer :: iprint
      iprint=0
!
!-----------------------------------------------------------------------
!
      select case(Lab)
      case(1)
!
        write(*,*) 'findcoord: UNFINISHED'
        stop 1
!
!-----------------------------------------------------------------------
!
      case(2)
!
!  .....start with the midpoint
        Eta(1:2) = 0.d0
!
!  .....Newton-Raphson iterations
        do iter=1,10
!
          call recta_linear(No,Eta,X,Dxdeta)
          call trobs(X, aux(1:3,1))
          xy(1:2) = aux(1:2,1)
!
!  .......residual
          bb(1:2) = Xp(1:2) - xy(1:2)
          d = sqrt(bb(1)**2+bb(2)**2)
          if (d.lt.GEOM_TOL) then
            Nfl = 1
            return
          endif
!
!  .......tangent matrix
          call trobs(Dxdeta(1:3,1), aux(1:3,1))
          call trobs(Dxdeta(1:3,2), aux(1:3,2))
          aa(1:2,1) = aux(1:2,1)
          aa(1:2,2) = aux(1:2,2)
!
!  .......solve for the increment
          call gausse(aa,2,bb, deta,2)
!
!  .......update
          Eta(1:2) = Eta(1:2) + deta(1:2)
        enddo
!
        Nfl = 0
        Eta(1:2) = 0.d0
        X(1:3) = 0.d0
!
      end select
!
   end subroutine findcoord

#endif

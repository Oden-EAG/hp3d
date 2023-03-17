c----------------------------------------------------------------------
c
c   routine name       - shapeQb
c
c----------------------------------------------------------------------
c
c   latest revision    - Oct 12
c
c   purpose            - routine evaluates L2 shape functions for master
c                        brick
c
c   arguments :
c     in:
c         Xi           - master element coordinates
c         Norder       - order of approximation for the nodes (in the H1 sense)
c     out:
c         Nrdof        - number of dof
c         Vshap        - values of shape functions
c
c   Note: Inspired in shapeHb but simplified to adapt it to L2 case
c----------------------------------------------------------------------
c
      subroutine shapeQb(Xi,Norder,Nrdof,ShapQ)
c
      use parameters

      implicit none

c  ...subroutine arguments
      double precision,intent(in) :: Xi(3)
      integer, intent(in) :: Norder(19)
      integer, intent(out) :: Nrdof
      double precision, intent(inout) ::  ShapQ(MAXbrickQ)

c  ...local variables
      double precision :: shapQ1(MAXP),shapQ2(MAXP),shapQ3(MAXP)
      integer :: nordh,nord1,nord2,nord3
      integer :: nrdofQ1,nrdofQ2,nrdofQ3
      integer :: i,j,k,kk

c
c
c----------------------------------------------------------------------
c
c  ...there are only middle node shape functions for L2

      call decode(Norder(19), nordh, nord3)
      call decode(nordh, nord1,nord2)

      call shape1Q(Xi(1),nord1,nrdofQ1,shapQ1)
      call shape1Q(Xi(2),nord2,nrdofQ2,shapQ2)
      call shape1Q(Xi(3),nord3,nrdofQ3,shapQ3)

      k=0
      do kk=1,nrdofQ3
        do j=1,nrdofQ2
          do i=1,nrdofQ1
            k=k+1
            ShapQ(k) = shapQ1(i)*shapQ2(j)*shapQ3(kk)
          enddo
        enddo
      enddo
c
      Nrdof = k
c
c
      end
c

c----------------------------------------------------------------------
c
c     routines to be implemented here (not done yet):

c                        - shapeQp (prism)
c                        - shapeQn (tetr)
c                        - shapeQd (pyra)
c
c----------------------------------------------------------------------
c
c   latest revision    - Oct 12
c
c   purpose            - routine evaluates L2 shape functions for master
c                        XXXXX
c
c   arguments :
c     in:
c         Xi           - master element coordinates
c         Norder       - order of approximation for the nodes
c     out:
c         Nrdof        - number of dof
c         Vshap        - values of shape functions
c
c   Note: Inspired in shapeHx but simplify to adapted it for L2 case
c----------------------------------------------------------------------



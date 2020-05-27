c----------------------------------------------------------------------
c----------------------------------------------------------------------
c DEPRECATED INTERFACE:
c previously used to fork between old and new shape functions package.
c
c Use the new shape functions, located in:
c     shape_1/ContExactSequence.F
c     shape_1/Segment.F
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
c   routine name       - shape1H
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 14
c
c   purpose            - evaluate 1D H1 shape functions
c
c   arguments :
c     in:
c        Xi            - master element coordinates
c        Nord          - polynomial order of approximation (H1 sense)
c     out:
c        NrdofH        - number of shape functions
c        ShapH         - values of shape functions
c        GradH         - values of derivatives of the shape functions
c                        wrt to master element coordinates
c
c----------------------------------------------------------------------
c
      subroutine shape1H(Xi,Nord, NrdofH,ShapH,GradH)
c
      use parameters , only : MAXP
c
      implicit none
      real(8), intent(in ) :: Xi
      integer, intent(in ) :: Nord
      integer, intent(out) :: NrdofH
      real(8), intent(out) :: ShapH(MAXP+1),GradH(MAXP+1)
c
      write(*,*) 'shape1H error: call shape1DH instead'
      stop
      !call shape1DH(Xi,Nord, NrdofH,ShapH,GradH)
c
      end subroutine shape1H
c
c
c----------------------------------------------------------------------
c
c   routine name       - shape1Q
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 14
c
c   purpose            - evaluate 1D L2 shape functions
c
c   arguments :
c     in:
c        Xi            - master element coordinates
c        Nord          - polynomial order of approximation (H1 sense)
c     out:
c        NrdofQ        - number of shape functions
c        ShapQ         - values of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape1Q(Xi,Nord, NrdofQ,ShapQ)
c
      use parameters , only : MAXP
c
      implicit none
      real(8), intent(in)  :: Xi
      integer, intent(in)  :: Nord
      integer, intent(out) :: NrdofQ
      real(8), intent(out) :: ShapQ(MAXP)
c
      write(*,*) 'shape1Q error: call shape1DQ instead'
      stop
      !call shape1DQ(Xi,Nord, NrdofQ,ShapQ)
c
      end subroutine shape1Q
c

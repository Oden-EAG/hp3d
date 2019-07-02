c----------------------------------------------------------------------
c----------------------------------------------------------------------
c In this file all the forking for SHAPE=0 or SHAPE=1 is performed for
c the case of the 1D sequence
c                           1D:    H1--->L2
c If SHAPE=0 the subroutines are located in shape_0/h1/shape1H.F
c                                           shape_0/l2/shape1Q.F  
c If SHAPE=1 the subroutines are located in shape_1/ContExactSequence.F
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
      real*8 , intent(in ) :: Xi
      integer, intent(in ) :: Nord
      integer, intent(out) :: NrdofH
      real*8 , intent(out) :: ShapH(MAXP+1),GradH(MAXP+1)
c
#if SHAPE_MODE==0
c
      call shape1H_0(Xi,Nord, NrdofH,ShapH,GradH)
c
#else
c
      call shape1DH(Xi,Nord, NrdofH,ShapH,GradH)
c  
#endif
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
      real*8 , intent(in)  :: Xi
      integer, intent(in)  :: Nord
      integer, intent(out) :: NrdofQ
      real*8 , intent(out) :: ShapQ(MAXP)
c
#if SHAPE_MODE==0
c
      call shape1Q_0(Xi,Nord, NrdofQ,ShapQ)
c
#else
c
      call shape1DQ(Xi,Nord, NrdofQ,ShapQ)
c
#endif
c
      end subroutine shape1Q
c

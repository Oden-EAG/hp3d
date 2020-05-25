c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c             IMPORTANT: This only applies for 2D problems!!!
c                 (when/if a 2D option is ever included)
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c In this file all the forking for SHAPE=0 or SHAPE=1 is performed for
c the traces of the 2D rotated sequence
c                   2D(Rot):    H1  -->  H(div)  -->  L2
c         Traces of 2D(Rot):  tr(H1)-->tr(H(div))
c The tr(H1) functions for this sequence are in Trace2DSequence.F
c If SHAPE=0 the subroutines are nowhere (rotated tr(H(curl) - shapEbe)
c If SHAPE=1 the subroutines are located in shape_1/TraceSequence.F
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
c   routine name       - shapVbe
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 14
c
c   purpose            - 1D edge trace of 2D H(div) edge functions
c                        (these are oriented 1D L2 functions -
c                        not 'bubbles')
c                        NOTE: only relevant in 2D problems!!
c
c   arguments :
c     in:
c       T              - local edge coordinate
c       Nord           - polynomial edge order (H1 sense)
c       Norient        - edge orientation
c     out:
c       NrdofV         - number of shape functions
c       ShapV          - values of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shapVbe(T,Nord,Norient, NrdofV,ShapV)
c
      use parameters , only : MAXP
c
      implicit none
      real(8), intent(in ) :: T
      integer, intent(in ) :: Nord
      integer, intent(in ) :: Norient
      integer, intent(out) :: NrdofV
      real(8), intent(out) :: ShapV(MAXP)
c
#if SHAPE_MODE==0
c
      write(*,*) 'shapVbe:This has NOT been implemented with SHAPE=0.'
      write(*,*) '        Please switch to SHAPE=1.'
      stop
c
#else
c
      call TraceEshapeV(T,Nord,Norient, NrdofV,ShapV)
c
#endif
c
      end subroutine shapVbe
c

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
c DEPRECATED INTERFACE:
c previously used to fork between old and new shape functions package.
c
c Use the new shape functions, located in:
c     shape_1/TraceSequence.F
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
      write(*,*) 'shapVbe error: call TraceEshapeV instead'
      stop
c     call TraceEshapeV(T,Nord,Norient, NrdofV,ShapV)
c
      end subroutine shapVbe
c

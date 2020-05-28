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
c   routine name       - shapHbe
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 14
c
c   purpose            - 1D edge trace of 2D H1 edge functions
c                        (these happen to be oriented 1D H1 'bubbles')
c
c   arguments :
c     in:
c       T              - local edge coordinate
c       Nord           - polynomial edge order (H1 sense)
c       Norient        - edge orientation
c     out:
c       NrdofH         - number of shape functions
c       ShapH          - values of shape functions
c       GradH          - local gradients of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shapHbe(T,Nord,Norient, NrdofH,ShapH,GradH)
c
      use parameters , only : MAXP
c
      implicit none
      real(8), intent(in ) :: T
      integer, intent(in ) :: Nord
      integer, intent(in ) :: Norient
      integer, intent(out) :: NrdofH
      real(8), intent(out) :: ShapH(MAXP-1),GradH(MAXP-1)
c
      write(*,*) 'shapHbe error: call TraceEshapeH instead'
      stop
c     call TraceEshapeH(T,Nord,Norient, NrdofH,ShapH,GradH)
c
      end subroutine shapHbe
c
c
c----------------------------------------------------------------------
c
c   routine name       - shapEbe
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 14
c
c   purpose            - 1D edge trace of 2D H(curl) edge functions
c                        (these are oriented 1D L2 functions -
c                        not 'bubbles')
c
c   arguments :
c     in:
c       T              - local edge coordinate
c       Nord           - polynomial edge order (H1 sense)
c       Norient        - edge orientation
c     out:
c       NrdofE         - number of shape functions
c       ShapE          - values of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shapEbe(T,Nord,Norient, NrdofE,ShapE)
c
      use parameters , only : MAXP
c
      implicit none
      real(8), intent(in ) :: T
      integer, intent(in ) :: Nord
      integer, intent(in ) :: Norient
      integer, intent(out) :: NrdofE
      real(8), intent(out) :: ShapE(MAXP)
c
      write(*,*) 'shapEbe error: call TraceEshapeE instead'
      stop
c     call TraceEshapeE(T,Nord,Norient, NrdofE,ShapE)
c
      end subroutine shapEbe
c

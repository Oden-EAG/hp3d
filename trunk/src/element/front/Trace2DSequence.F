c----------------------------------------------------------------------
c----------------------------------------------------------------------
c In this file all the forking for SHAPE=0 or SHAPE=1 is performed for
c the traces of the 2D sequence
c                      2D:    H1  -->  H(curl)  -->  L2
c            Traces of 2D:  tr(H1)-->tr(H(curl))
c If SHAPE=0 the subroutines are located in shape_0/h1/Hkernels.F
c                                           shape_0/hcurl/shapEbe.F
c If SHAPE=1 the subroutines are located in shape_1/TraceSequence.F
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
      real*8 , intent(in ) :: T
      integer, intent(in ) :: Nord
      integer, intent(in ) :: Norient
      integer, intent(out) :: NrdofH
      real*8 , intent(out) :: ShapH(MAXP-1),GradH(MAXP-1)
c
#if SHAPE_MODE==0
c      
      call shapHbe_0(T,Nord,Norient, NrdofH,ShapH,GradH)
c      
#else
c
      call TraceEshapeH(T,Nord,Norient, NrdofH,ShapH,GradH)
c      
#endif
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
      real*8 , intent(in ) :: T
      integer, intent(in ) :: Nord
      integer, intent(in ) :: Norient
      integer, intent(out) :: NrdofE
      real*8 , intent(out) :: ShapE(MAXP)
c     
#if SHAPE_MODE==0
c
      call shapEbe_0(T,Nord,Norient, NrdofE,ShapE)
c 
#else
c
      call TraceEshapeE(T,Nord,Norient, NrdofE,ShapE)
c         
#endif
c
      end subroutine shapEbe
c

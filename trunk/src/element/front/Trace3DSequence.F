c----------------------------------------------------------------------
c----------------------------------------------------------------------
c In this file all the forking for SHAPE=0 or SHAPE=1 is performed for
c the traces of the 3D sequence
c              3D:    H1  -->  H(curl)  -->  H(div)  -->  L2
c    Traces of 3D:  tr(H1)-->tr(H(curl))-->tr(H(div))
c If SHAPE=0 the subroutines are located in shape_0/h1/Hkernels.F
c                                           shape_0/hcurl/shapEb2.F
c If SHAPE=1 the subroutines are located in shape_1/TraceSequence.F
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
c   routine name       - shapHb2
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 14
c
c   purpose            - 2D face trace of 3D H1 shape functions
c                        (these happen to be oriented 2D H1 'bubbles')
c
c   arguments :
c     in:
c       Type           - face type (quad or tria)
c       T              - local face coordinates
c       Nord           - polynomial face order (H1 sense)
c       Norient        - face orientation
c     out:
c       NrdofH         - number of shape functions
c       ShapH          - values of shape functions
c       GradH          - values of derivatives of the shape functions
c                        wrt to master element coordinates      
c
c----------------------------------------------------------------------
c
      subroutine shapHb2(Type,T,Nord,Norient, NrdofH,ShapH,GradH)
c      
      use parameters , only : MAXmdlqH
c
      implicit none
      character(len=4), intent(in)  :: Type
      real(8), intent(in ) :: T(2)
      integer, intent(in ) :: Nord
      integer, intent(in ) :: Norient
      integer, intent(out) :: NrdofH
      real(8), intent(out) :: ShapH(MAXmdlqH),GradH(2,MAXmdlqH)
c
#if SHAPE_MODE==0
c
      call shapHb2_0(Type,T,Nord,Norient, NrdofH,ShapH,GradH)
c
#else
c
      call TraceFshapeH(Type,T,Nord,Norient, NrdofH,ShapH,GradH)
c  
#endif
c
      end subroutine shapHb2
c
c
c----------------------------------------------------------------------
c     
c   routine name       - shapEb2
c     
c----------------------------------------------------------------------
c     
c   latest revision    - Nov 14
c     
c   purpose            - 2D face trace of 3D H(curl) face functions
c                        (these are oriented 2D H(curl) 'bubbles')
c     
c   arguments :
c     in:
c       Type           - face type (quad or tria)
c       T              - local face coordinates
c       Nord           - polynomial face order (H1 sense)
c       Norient        - face orientation
c     out:
c       NrdofE         - number of shape functions
c       ShapE          - values of shape functions
c       CurlE          - curl of shape functions      
c     
c----------------------------------------------------------------------
c     
      subroutine shapEb2(Type,T,Nord,Norient, NrdofE,ShapE,CurlE)
c   
      use parameters , only : MAXmdlqE
c
      implicit none      
      character(len=4), intent(in)  :: Type
      real(8), intent(in ) :: T(2)
      integer, intent(in ) :: Nord
      integer, intent(in ) :: Norient
      integer, intent(out) :: NrdofE
      real(8), intent(out) :: ShapE(2,MAXmdlqE),CurlE(MAXmdlqE)
c     
#if SHAPE_MODE==0
c
      call shapEb2_0(Type,T,Nord,Norient, NrdofE,ShapE,CurlE)
c 
#else
c
      call TraceFshapeE(Type,T,Nord,Norient, NrdofE,ShapE,CurlE)
c         
#endif
c
      end subroutine shapEb2
c
c
c----------------------------------------------------------------------
c     
c   routine name       - shapVb2
c     
c----------------------------------------------------------------------
c     
c   latest revision    - Nov 14
c     
c   purpose            - 2D face trace of 3D H(div) face functions
c                        (these are oriented 2D L2 functions -
c                        not 'bubbles')
c     
c   arguments :
c     in:
c       Type           - face type (quad or tria)
c       T              - local face coordinates
c       Nord           - polynomial face order (H1 sense)
c       Norient        - face orientation
c     out:
c       NrdofV         - number of shape functions
c       ShapV          - values of shape functions   
c     
c----------------------------------------------------------------------
c     
      subroutine shapVb2(Type,T,Nord,Norient, NrdofV,ShapV)
c   
      use parameters , only : MAXmdlqQ
c
      implicit none      
      character(len=4), intent(in)  :: Type
      real(8), intent(in ) :: T(2)
      integer, intent(in ) :: Nord
      integer, intent(in ) :: Norient
      integer, intent(out) :: NrdofV
      real(8), intent(out) :: ShapV(MAXmdlqQ)
c     
#if SHAPE_MODE==0
c
      write(*,*) 'shapeVb2:This has NOT been implemented with SHAPE=0.'
      write(*,*) '         Please switch to SHAPE=1.'
      stop
c
#else
c
      call TraceFshapeV(Type,T,Nord,Norient, NrdofV,ShapV)
c         
#endif
c
      end subroutine shapVb2
c

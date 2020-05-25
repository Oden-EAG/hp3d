c----------------------------------------------------------------------
c----------------------------------------------------------------------
c In this file all the forking for SHAPE=0 or SHAPE=1 is performed for
c the case of the 3D sequence
c                 3D:    H1--->H(curl)--->H(div)--->L2
c If SHAPE=0 the subroutines are located in shape_0/h1/shape3H.F
c                                           shape_0/hcurl/shape3E.F 
c                                           shape_0/l2/shape3Q.F
c If SHAPE=1 the subroutines are located in shape_1/ContExactSequence.F
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
c   routine name       - shape3H
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 14
c
c   purpose            - routine evaluates H1 shape functions for 3D
c                        elements of various types
c
c   arguments :
c     in:
c       Type           - element type
c       Xi             - master element coordinates
c       Norder         - polynomial order for the nodes (H1 sense)
c       Nedge_orient   - edge orientations
c       Nface_orient   - face orientations
c     out:
c       NrdofH         - number of the element shape functions      
c       ShapH          - values of shape functions
c       GradH          - values of derivatives of the shape functions
c                        wrt to master element coordinates      
c
c----------------------------------------------------------------------
c
      subroutine shape3H(Type,Xi,Norder,Nedge_orient,Nface_orient, 
     .                   NrdofH,ShapH,GradH)
c
      use parameters , only : MAXbrickH
c 
      implicit none
      character(len=4), intent(in)  :: Type
      integer, intent(in)  :: Nedge_orient(12)
      integer, intent(in)  :: Nface_orient(6)
      integer, intent(in)  :: Norder(19)
      real(8), intent(in)  :: Xi(3)
      integer, intent(out) :: NrdofH
      real(8), intent(out) :: ShapH(MAXbrickH),GradH(3,MAXbrickH)
c
#if SHAPE_MODE==0
c      
      call shape3H_0(Type,Xi,Norder,Nedge_orient,Nface_orient, 
     .               NrdofH,ShapH,GradH)
c      
#else
c
      call shape3DH(Type,Xi,Norder,Nedge_orient,Nface_orient, 
     .               NrdofH,ShapH,GradH)
c      
#endif
c
      end subroutine shape3H
c
c
c----------------------------------------------------------------------
c
c   routine name       - shape3E
c
c----------------------------------------------------------------------
c
c   latest revision    - Oct 14
c
c   purpose            - routine evaluates H(curl) shape functions for 
c                        3D elements of various types
c
c   arguments :
c     in:
c       Type           - element type
c       Xi             - master element coordinates
c       Norder         - polynomial order for the nodes (H1 sense)
c       Nedge_orient   - edge orientations
c       Nface_orient   - face orientations
c     out:
c       NrdofE         - number of the element shape functions      
c       ShapE          - shape functions
c       CurlE          - curl of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape3E(Type,Xi,Norder,Nedge_orient,Nface_orient, 
     .                   NrdofE,ShapE,CurlE)    
c
      use parameters , only : MAXbrickE
c 
      implicit none
      character(len=4), intent(in)  :: Type
      integer, intent(in)  :: Nedge_orient(12)
      integer, intent(in)  :: Nface_orient(6)
      integer, intent(in)  :: Norder(19)
      real(8), intent(in)  :: Xi(3)
      integer, intent(out) :: NrdofE
      real(8), intent(out) :: ShapE(3,MAXbrickE),CurlE(3,MAXbrickE)
c     
#if SHAPE_MODE==0
c
      call shape3E_0(Type,Xi,Norder,Nedge_orient,Nface_orient, 
     .               NrdofE,ShapE,CurlE)
c 
#else
c
      call shape3DE(Type,Xi,Norder,Nedge_orient,Nface_orient, 
     .               NrdofE,ShapE,CurlE)
c         
#endif
c
      end subroutine shape3E
c
c
c----------------------------------------------------------------------
c
c   routine name       - shape3V
c
c----------------------------------------------------------------------
c
c   latest revision    - Oct 14
c
c   purpose            - routine evaluates H(div) shape functions for 
c                        3D elements of various types
c
c   arguments :
c     in:
c       Type           - element type
c       Xi             - master element coordinates
c       Norder         - polynomial order for the nodes (H1 sense)
c       Nface_orient   - face orientations
c     out:
c       NrdofV         - number of the element shape functions      
c       ShapV          - shape functions
c       DivV           - divergence of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape3V(Type,Xi,Norder,Nface_orient, 
     .                   NrdofV,ShapV,DivV)    
c
      use parameters , only : MAXbrickV
c 
      implicit none
      character(len=4), intent(in)  :: Type
      integer, intent(in)  :: Nface_orient(6)
      integer, intent(in)  :: Norder(19)
      real(8), intent(in)  :: Xi(3)
      integer, intent(out) :: NrdofV
      real(8), intent(out) :: ShapV(3,MAXbrickV),DivV(MAXbrickV)
c     
#if SHAPE_MODE==0
c
      write(*,*) 'shape3V:This has NOT been implemented with SHAPE=0.'
      write(*,*) '        Please switch to SHAPE=1.'
      stop
c
#else
c
      call shape3DV(Type,Xi,Norder,Nface_orient, 
     .               NrdofV,ShapV,DivV)
c         
#endif
c
      end subroutine shape3V
c
c
c----------------------------------------------------------------------
c
c   routine name       - shape3Q
c
c----------------------------------------------------------------------
c
c   latest revision    - Oct 2014
c
c   purpose            - routine evaluates L2 shape functions for 3D
c                        elements of various types
c
c   arguments :
c     in:
c       Type           - element type
c       Xi             - master element coordinates
c       Norder         - polynomial order of interior node (H1 sense)
c     out:
c       NrdofQ         - number of the element shape functions
c       ShapQ          - values of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape3Q(Type,Xi,Norder, NrdofQ,ShapQ)
c
      use parameters , only : MAXbrickQ
c
      implicit none
      character(len=4), intent(in)  :: Type
      real(8), intent(in)  :: Xi(3)
      integer, intent(in)  :: Norder(19)
      integer, intent(out) :: NrdofQ
      real(8), intent(out) :: ShapQ(MAXbrickQ)
c
#if SHAPE_MODE==0
c
      call shape3Q_0(Type,Xi,Norder, NrdofQ,ShapQ)
c
#else
c
      call shape3DQ(Type,Xi,Norder, NrdofQ,ShapQ)
c
#endif
c
      end subroutine shape3Q
c      

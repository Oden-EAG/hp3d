c----------------------------------------------------------------------
c----------------------------------------------------------------------
c DEPRECATED INTERFACE:
c previously used to fork between old and new shape functions package.
c
c Use the new shape functions, located in:
c     shape_1/ContExactSequence.F
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
      NrdofH = 0; ShapH = 0.d0; GradH = 0.d0
      write(*,*) 'shape3H error: call shape3DH instead'
      stop
c     call shape3DH(Type,Xi,Norder,Nedge_orient,Nface_orient,
c    .              NrdofH,ShapH,GradH)
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
      NrdofE = 0; ShapE = 0.d0; CurlE = 0.d0
      write(*,*) 'shape3E error: call shape3DE instead'
      stop
c     call shape3DE(Type,Xi,Norder,Nedge_orient,Nface_orient,
c    .              NrdofE,ShapE,CurlE)
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
      NrdofV = 0; ShapV = 0.d0; DivV = 0.d0
      write(*,*) 'shape3V error: call shape3DV instead'
      stop
c     call shape3DV(Type,Xi,Norder,Nface_orient,
c    .              NrdofV,ShapV,DivV)
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
      NrdofQ = 0; ShapQ = 0.d0
      write(*,*) 'shape3Q error: call shape3DQ instead'
      stop
c     call shape3DQ(Type,Xi,Norder, NrdofQ,ShapQ)
c
      end subroutine shape3Q
c

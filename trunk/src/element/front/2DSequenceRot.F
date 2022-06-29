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
c     shape_1/ContExactSequence.F
c     shape_1/Quadrilateral.F
c     shape_1/Triangle.F
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
c     routine name      - shape2V
c
c----------------------------------------------------------------------
c
c     latest revision:  - Nov 14
c
c     purpose:          - routine returns values of a 2D element
c                         H(div) shape functions and their divergence
c                         NOTE: only relevant in 2D problems!!
c
c     arguments:
c       in:
c         Type          - element type
c         Xi            - master element coordinates
c         Norder        - polynomial order for the nodes (H1 sense)
c         Nedge_orient  - edge orientations
c       out:
c         NrdofV        - number of dof
c         ShapV         - values of the shape functions at the point
c         DivV          - divergences of the shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape2V(Type,Xi,Norder,Nedge_orient,
     .                                              NrdofV,ShapV,DivV)
c
      use parameters , only : MAXquadV
c
      implicit none
      character(len=4), intent(in)  :: Type
      real(8), intent(in ) :: Xi(2)
      integer, intent(in ) :: Norder(5)
      integer, intent(in ) :: Nedge_orient(4)
      integer, intent(out) :: NrdofV
      real(8), intent(out) :: ShapV(2,MAXquadV),DivV(MAXquadV)
c
      NrdofV = 0; ShapV = 0.d0; DivV = 0.d0
      write(*,*) 'shape2V error: call shape2DV instead'
      stop
c     call shape2DV(Type,Xi,Norder,Nedge_orient, NrdofV,ShapV,DivV)
c
      end subroutine shape2V
c
c
c----------------------------------------------------------------------
c
c   routine name       - shapeVq
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 14
c
c   purpose            - evaluate 2D quad H(div) shape functions
c                        NOTE: only relevant in 2D problems!!
c
c     arguments :
c     in:
c       Xi             - master element coordinates
c       Norder         - polynomial order for the nodes (H1 sense)
c       Nedge_orient   - edge orientations
c     out:
c       NrdofV         - number of element dof
c       ShapV          - values of shape functions
c       DivV           - values of derivatives of the shape functions
c                        wrt to master element coordinates
c
c----------------------------------------------------------------------
c
      subroutine shapeVq(Xi,Norder,Nedge_orient, NrdofV,ShapV,DivV)
c
      use parameters , only : MAXquadV
c
      implicit none
      real(8), intent(in ) :: Xi(2)
      integer, intent(in ) :: Norder(5)
      integer, intent(in ) :: Nedge_orient(4)
      integer, intent(out) :: NrdofV
      real(8), intent(out) :: ShapV(2,MAXquadV),DivV(MAXquadV)
c
      NrdofV = 0; ShapV = 0.d0; DivV = 0.d0
      write(*,*) 'shapeVq error: call shape2DV instead'
      stop
c     call shape2V('quad',Xi,Norder,Nedge_orient, NrdofV,ShapV,DivV)
c
      end subroutine shapeVq
c
c
c----------------------------------------------------------------------
c
c   routine name       - shapeVt
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 14
c
c   purpose            - evaluate 2D triangle H(div) shape functions
c                        NOTE: only relevant in 2D problems!!
c
c     arguments :
c     in:
c       Xi             - master element coordinates
c       Norder         - polynomial order for the nodes (H1 sense)
c       Nedge_orient   - edge orientaions
c     out:
c       NrdofV         - number of element dof
c       ShapV          - values of shape functions
c       DivV           - values of derivatives of the shape functions
c                        wrt to master element coordinates
c
c----------------------------------------------------------------------
c
      subroutine shapeVt(Xi,Norder,Nedge_orient, NrdofV,ShapV,DivV)
c
      use parameters , only : MAXtriaV
c
      implicit none
      real(8), intent(in ) :: Xi(2)
      integer, intent(in ) :: Norder(4)
      integer, intent(in ) :: Nedge_orient(3)
      integer, intent(out) :: NrdofV
      real(8), intent(out) :: ShapV(2,MAXtriaV),DivV(MAXtriaV)
c
      NrdofV = 0; ShapV = 0.d0; DivV = 0.d0
      write(*,*) 'shapeVt error: call shape2DV instead'
      stop
c     call shape2V('tria',Xi,Norder,Nedge_orient, NrdofV,ShapV,DivV)
c
      end subroutine shapeVt
c

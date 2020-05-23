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
c the case of the 2D H(div), which is part of the 2D rotated sequence
c                     2D(Rot):    H1--->H(div)--->L2
c The H1 and L2 functions for this sequence are already in 2DSequence.F
c If SHAPE=0 the subroutines are nowhere (they are the rotated H(curl))
c If SHAPE=1 the subroutines are located in shape_1/ContExactSequence.F
c                                           shape_1/Quadrilateral.F
c                                           shape_1/Triangle.F
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
      real*8 , intent(in ) :: Xi(2)
      integer, intent(in ) :: Norder(5)
      integer, intent(in ) :: Nedge_orient(4)
      integer, intent(out) :: NrdofV
      real*8 , intent(out) :: ShapV(2,MAXquadV),DivV(MAXquadV)
c
#if SHAPE_MODE==0
c
      write(*,*) 'shape2V:This has NOT been implemented with SHAPE=0.'
      write(*,*) '        Please switch to SHAPE=1.'
      stop
c 
#else
c
      call shape2DV(Type,Xi,Norder,Nedge_orient, NrdofV,ShapV,DivV)
c
#endif
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
      real*8 , intent(in ) :: Xi(2)
      integer, intent(in ) :: Norder(5)
      integer, intent(in ) :: Nedge_orient(4)
      integer, intent(out) :: NrdofV
      real*8 , intent(out) :: ShapV(2,MAXquadV),DivV(MAXquadV)
c
      call shape2V('quad',Xi,Norder,Nedge_orient, NrdofV,ShapV,DivV)
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
      use parameters , only : MAXtriaV,MAXquadV
c
      implicit none      
      real*8 , intent(in ) :: Xi(2)
      integer, intent(in ) :: Norder(4)
      integer, intent(in ) :: Nedge_orient(3)
      integer, intent(out) :: NrdofV
      real*8 , intent(out) :: ShapV(2,MAXtriaV),DivV(MAXtriaV)
c
c  ...locals
      integer   :: norder_loc(5)
      integer   :: nedge_orient_loc(4)
      real*8    :: shapV_loc(2,MAXquadV),divV_loc(MAXquadV)
c
      norder_loc(1:4) = Norder(1:4);norder_loc(5)=1
      nedge_orient_loc(1:3) = Nedge_orient(1:3);nedge_orient_loc(4)=0
      call shape2V('tria',Xi,norder_loc,nedge_orient_loc, 
     .              NrdofV,shapV_loc,divV_loc)
      ShapV(1:2,1:NrdofV) = shapV_loc(1:2,1:NrdofV)
      DivV(1:NrdofV) = divV_loc(1:NrdofV)
c
      end subroutine shapeVt
c

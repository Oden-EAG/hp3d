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
c     routine name      - shape2H
c
c----------------------------------------------------------------------
c
c     latest revision:  - Nov 14
c
c     purpose:          - routine returns values of a 2D element
c                         H1 shape functions and their derivatives
c
c     arguments:
c       in:
c         Type          - element type
c         Xi            - master element coordinates
c         Norder        - polynomial order for the nodes (H1 sense)
c         Nedge_orient  - edge orientations
c
c       out:
c         NrdofH        - number of dof
c         ShapH         - values of the shape functions at the point
c         GradH         - gradients of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape2H(Type,Xi,Norder,Nedge_orient, 
     .                                              NrdofH,ShapH,GradH)
c
      use parameters , only : MAXquadH
c
      implicit none
      character(len=4), intent(in)  :: Type      
      real(8), intent(in ) :: Xi(2)
      integer, intent(in ) :: Norder(5)
      integer, intent(in ) :: Nedge_orient(4)
      integer, intent(out) :: NrdofH
      real(8), intent(out) :: ShapH(MAXquadH),GradH(2,MAXquadH)
c
      write(*,*) 'shape2H error: call shape2DH instead'
      stop
c     call shape2DH(Type,Xi,Norder,Nedge_orient, NrdofH,ShapH,GradH)
c
      end subroutine shape2H
c
c
c----------------------------------------------------------------------
c
c   routine name       - shapeHq
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 14
c
c   purpose            - evaluate quad H1 shape functions
c
c   arguments :
c     in:
c        Xi            - master element coordinates
c        Norder        - polynomial order for the nodes (H1 sense)
c        Nedge_orient  - edge orientations
c     out:
c        NrdofH        - number of dof
c        ShapH         - values of shape functions
c        GradH         - values of derivatives of the shape functions
c                        wrt to master element coordinates      
c
c----------------------------------------------------------------------
c
      subroutine shapeHq(Xi,Norder,Nedge_orient, NrdofH,ShapH,GradH)
c
      use parameters , only : MAXquadH
c
      implicit none      
      real(8), intent(in ) :: Xi(2)
      integer, intent(in ) :: Norder(5)
      integer, intent(in ) :: Nedge_orient(4)
      integer, intent(out) :: NrdofH
      real(8), intent(out) :: ShapH(MAXquadH),GradH(2,MAXquadH)
c
      write(*,*) 'shapeHq error: call shape2DH instead'
      stop
c     call shape2H('quad',Xi,Norder,Nedge_orient, NrdofH,ShapH,GradH)
c
      end subroutine shapeHq
c
c
c----------------------------------------------------------------------
c
c   routine name       - shapeHt
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 14
c
c   purpose            - evaluate triangle H1 shape functions
c
c   arguments :
c     in:
c        Xi            - master element coordinates
c        Norder        - polynomial order for the nodes (H1 sense)
c        Nedge_orient  - edge orientations
c     out:
c        NrdofH        - number of dof
c        ShapH         - values of shape functions
c        GradH         - values of derivatives of the shape functions
c                        wrt to master element coordinates       
c
c----------------------------------------------------------------------
c
      subroutine shapeHt(Xi,Norder,Nedge_orient, NrdofH,ShapH,GradH)
c
      use parameters , only : MAXtriaH
c
      implicit none      
      real(8), intent(in ) :: Xi(2)
      integer, intent(in ) :: Norder(4)
      integer, intent(in ) :: Nedge_orient(3)
      integer, intent(out) :: NrdofH
      real(8), intent(out) :: ShapH(MAXtriaH),GradH(2,MAXtriaH)
c
      write(*,*) 'shapeHt error: call shape2DH instead'
      stop
c     call shape2H('tria',Xi,Norder,Nedge_orient, NrdofH,ShapH,GradH)
c
c
      end subroutine shapeHt
c
c
c----------------------------------------------------------------------
c
c     routine name      - shape2E
c
c----------------------------------------------------------------------
c
c     latest revision:  - Nov 14
c
c     purpose:          - routine returns values of a 2D element
c                         H(curl) shape functions and their curl
c
c     arguments:
c       in:
c         Type          - element type
c         Xi            - master element coordinates
c         Norder        - polynomial order for the nodes (H1 sense)
c         Nedge_orient  - edge orientations
c       out:
c         NrdofE        - number of dof
c         ShapE         - values of the shape functions at the point
c         CurlE         - curl of the shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape2E(Type,Xi,Norder,Nedge_orient, 
     .                                              NrdofE,ShapE,CurlE)
c     
      use parameters , only : MAXquadE
c
      implicit none
      character(len=4), intent(in)  :: Type
      real(8), intent(in ) :: Xi(2)
      integer, intent(in ) :: Norder(5)
      integer, intent(in ) :: Nedge_orient(4)
      integer, intent(out) :: NrdofE
      real(8), intent(out) :: ShapE(2,MAXquadE),CurlE(MAXquadE)
c
      write(*,*) 'shape2E error: call shape2DE instead'
      stop
c     call shape2DE(Type,Xi,Norder,Nedge_orient, NrdofE,ShapE,CurlE)
c
      end subroutine shape2E
c
c
c----------------------------------------------------------------------
c
c   routine name       - shapeEq
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 14
c
c   purpose            - evaluate quad H(curl) shape functions
c
c   arguments :
c     in:
c       Xi             - master element coordinates
c       Norder         - polynomial order for the nodes (H1 sense)
c       Nedge_orient   - edge orientations
c     out:
c       NrdofE         - number of element dof
c       ShapE          - values of shape functions
c       CurlE          - curl of the shape functions      
c     
c----------------------------------------------------------------------
c
      subroutine shapeEq(Xi,Norder,Nedge_orient, NrdofE,ShapE,CurlE)
c     
      use parameters , only : MAXquadE
c
      implicit none      
      real(8), intent(in ) :: Xi(2)
      integer, intent(in ) :: Norder(5)
      integer, intent(in ) :: Nedge_orient(4)
      integer, intent(out) :: NrdofE
      real(8), intent(out) :: ShapE(2,MAXquadE),CurlE(MAXquadE)
c
      write(*,*) 'shapeEq error: call shape2DE instead'
      stop
c     call shape2E('quad',Xi,Norder,Nedge_orient, NrdofE,ShapE,CurlE)
c
      end subroutine shapeEq
c
c
c----------------------------------------------------------------------
c
c   routine name       - shapeEt
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 14
c
c   purpose            - evaluate triangle H(curl) shape functions
c
c   arguments :
c     in:
c       Xi             - master element coordinates
c       Norder         - polynomial order for the nodes (H1 sense)
c       Nedge_orient   - edge orientations
c     out:
c       NrdofE         - number of element dof
c       ShapE          - values of shape functions
c       CurlE          - curl of the shape functions      
c     
c----------------------------------------------------------------------
c     
      subroutine shapeEt(Xi,Norder,Nedge_orient, NrdofE,ShapE,CurlE)
c     
      use parameters , only : MAXtriaE
c
      implicit none
c
      real(8), intent(in ) :: Xi(2)
      integer, intent(in ) :: Norder(4)
      integer, intent(in ) :: Nedge_orient(3)
      integer, intent(out) :: NrdofE
      real(8), intent(out) :: ShapE(2,MAXtriaE),CurlE(MAXtriaE)
c
      write(*,*) 'shapeEt error: call shape2DE instead'
      stop
c     call shape2E('tria',Xi,Norder,Nedge_orient, NrdofE,ShapE,CurlE)
c
      end subroutine shapeEt
c
c
c----------------------------------------------------------------------
c
c     routine name      - shape2Q
c
c----------------------------------------------------------------------
c
c     latest revision:  - Nov 14
c
c     purpose:          - routine returns values of a 2D element
c                         L2 shape functions
c
c     arguments:
c       in:
c         Type          - element type
c         Xi            - master element coordinates
c         Norder        - polynomial order for the nodes (H1 sense)
c       out:
c         NrdofQ        - number of dof
c         ShapQ         - values of the shape functions at the point
c
c----------------------------------------------------------------------
c
      subroutine shape2Q(Type,Xi,Norder, NrdofQ,ShapQ)
c
      use parameters , only : MAXquadQ
c
      implicit none
      character(len=4), intent(in)  :: Type
      real(8), intent(in)  :: Xi(2)
      integer, intent(in)  :: Norder(5)
      integer, intent(out) :: NrdofQ
      real(8), intent(out) :: ShapQ(MAXquadQ)
c
      write(*,*) 'shape2Q error: call shape2DQ instead'
      stop
c     call shape2DQ(Type,Xi,Norder, NrdofQ,ShapQ)
c
      end subroutine shape2Q
c
c
c----------------------------------------------------------------------
c
c     routine name      - shapeQq
c
c----------------------------------------------------------------------
c
c     latest revision:  - Nov 14
c
c     purpose:          - routine returns values of quad L2
c                         shape functions
c
c     arguments:
c       in:
c         Xi            - master element coordinates
c         Norder        - polynomial order of face node (H1 sense)
c       out:
c         NrdofQ        - number of dof
c         ShapQ         - values of the shape functions at the point
c
c----------------------------------------------------------------------
c
      subroutine shapeQq(Xi,Norder, NrdofQ,ShapQ)
c
      use parameters , only : MAXquadQ
c
      implicit none
      real(8), intent(in)  :: Xi(2)
      integer, intent(in)  :: Norder
      integer, intent(out) :: NrdofQ
      real(8), intent(out) :: ShapQ(MAXquadQ)
c
      write(*,*) 'shapeQq error: call shape2DQ instead'
      stop
c     call shape2Q('quad',Xi,Norder, NrdofQ,ShapQ)
c
      end subroutine shapeQq
c
c
c----------------------------------------------------------------------
c
c     routine name      - shapeQt
c
c----------------------------------------------------------------------
c
c     latest revision:  - Nov 14
c
c     purpose:          - routine returns values of triangle L2
c                         shape functions
c
c     arguments:
c       in:
c         Xi            - master element coordinates
c         Norder        - polynomial order of face node (H1 sense)
c       out:
c         NrdofQ        - number of dof
c         ShapQ         - values of the shape functions at the point
c
c----------------------------------------------------------------------
c
      subroutine shapeQt(Xi,Norder, NrdofQ,ShapQ)
c
      use parameters , only : MAXtriaQ
c
      implicit none
      real(8), intent(in)  :: Xi(2)
      integer, intent(in)  :: Norder
      integer, intent(out) :: NrdofQ
      real(8), intent(out) :: ShapQ(MAXtriaQ)
c
      write(*,*) 'shapeQt error: call shape2DQ instead'
      stop
c     call shape2Q('tria',Xi,Norder, NrdofQ,ShapQ)
c
      end subroutine shapeQt
c

c----------------------------------------------------------------------
c----------------------------------------------------------------------
c In this file all the forking for SHAPE=0 or SHAPE=1 is performed for
c the case of the 2D sequence
c                      2D:    H1--->H(curl)--->L2
c If SHAPE=0 the subroutines are located in shape_0/h1/shape2H.F
c                                           shape_0/hcurl/shape2E.F 
c If SHAPE=1 the subroutines are located in shape_1/ContExactSequence.F
c                                           shape_1/Quadrilateral.F
c                                           shape_1/Triangle.F
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
      real*8 , intent(in ) :: Xi(2)
      integer, intent(in ) :: Norder(5)
      integer, intent(in ) :: Nedge_orient(4)
      integer, intent(out) :: NrdofH
      real*8 , intent(out) :: ShapH(MAXquadH),GradH(2,MAXquadH)
c
#if SHAPE_MODE==0
c
      write(*,*) 'shape2H:This has NOT been implemented with SHAPE=0.'
      write(*,*) '        Use shapeHq and shapeHt directly.'
      stop
c
#else
c
      call shape2DH(Type,Xi,Norder,Nedge_orient, NrdofH,ShapH,GradH)
c  
#endif
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
      real*8 , intent(in ) :: Xi(2)
      integer, intent(in ) :: Norder(5)
      integer, intent(in ) :: Nedge_orient(4)
      integer, intent(out) :: NrdofH
      real*8 , intent(out) :: ShapH(MAXquadH),GradH(2,MAXquadH)
c
#if SHAPE_MODE==0
c
      call shapeHq_0(Xi,Norder,Nedge_orient, NrdofH,ShapH,GradH)
c
#else
c
ccc      call shape2H('quad',Xi,Norder,Nedge_orient, NrdofH,ShapH,GradH)
      call shape2DH('quad',Xi,Norder,Nedge_orient, NrdofH,ShapH,GradH)
c  
#endif
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
      use parameters , only : MAXtriaH, MAXquadH
c
      implicit none      
      real*8 , intent(in ) :: Xi(2)
      integer, intent(in ) :: Norder(4)
      integer, intent(in ) :: Nedge_orient(3)
      integer, intent(out) :: NrdofH
      real*8 , intent(out) :: ShapH(MAXtriaH),GradH(2,MAXtriaH)
c
c  ...locals
      integer :: norder_loc(5), nedge_orient_loc(4)
      real*8  :: shapH_loc(MAXquadH),gradH_loc(2,MAXquadH)
c
#if SHAPE_MODE==0
c
      call shapeHt_0(Xi,Norder,Nedge_orient, NrdofH,ShapH,GradH)
c
#else
      norder_loc(1:4) = Norder(1:4); norder_loc(5) = 11
      nedge_orient_loc(1:3) = Nedge_orient(1:3);nedge_orient_loc(4)=0 
c
      call shape2H('tria',Xi,norder_loc,nedge_orient_loc, 
     .              NrdofH,shapH_loc,gradH_loc)
      ShapH(1:NrdofH) = shapH_loc(1:NrdofH)
      GradH(1:2,1:NrdofH) = gradH_loc(1:2,1:NrdofH)
c  
#endif
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
      real*8 , intent(in ) :: Xi(2)
      integer, intent(in ) :: Norder(5)
      integer, intent(in ) :: Nedge_orient(4)
      integer, intent(out) :: NrdofE
      real*8 , intent(out) :: ShapE(2,MAXquadE),CurlE(MAXquadE)
c
#if SHAPE_MODE==0
c
      write(*,*) 'shape2E:This has NOT been implemented with SHAPE=0.'
      write(*,*) '        Use shapeEq and shapeEt directly.'
      stop
c
#else
c
      call shape2DE(Type,Xi,Norder,Nedge_orient, NrdofE,ShapE,CurlE)
c  
#endif
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
      real*8 , intent(in ) :: Xi(2)
      integer, intent(in ) :: Norder(5)
      integer, intent(in ) :: Nedge_orient(4)
      integer, intent(out) :: NrdofE
      real*8 , intent(out) :: ShapE(2,MAXquadE),CurlE(MAXquadE)
c
#if SHAPE_MODE==0
c
      write(*,*) 'shapeEq:This has NOT been implemented with SHAPE=0.'
      write(*,*) '        Please switch to SHAPE=1.'
      stop
c 
#else
c
ccc      call shape2E('quad',Xi,Norder,Nedge_orient, NrdofE,ShapE,CurlE)
      call shape2DE('quad',Xi,Norder,Nedge_orient, NrdofE,ShapE,CurlE)
c
#endif
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
      use parameters , only : MAXtriaE,MaxquadE
c
      implicit none      
      real*8 , intent(in ) :: Xi(2)
      integer, intent(in ) :: Norder(4)
      integer, intent(in ) :: Nedge_orient(3)
      integer, intent(out) :: NrdofE
      real*8 , intent(out) :: ShapE(2,MAXtriaE),CurlE(MAXtriaE)
c
c  ...locals
      integer  :: norder_loc(5)
      integer  :: nedge_orient_loc(4)
      real*8   :: shapE_loc(2,MAXquadE),curlE_loc(MAXquadE)

c
#if SHAPE_MODE==0
c
      call shapeEt_0(Xi,Norder,Nedge_orient, NrdofE,ShapE,CurlE)
c
#else
c
      norder_loc(1:4) = Norder(1:4); norder_loc(5)=11
      nedge_orient_loc(1:3) = Nedge_orient(1:3);nedge_orient_loc(4)=0
ccc      call shape2E('tria',Xi,norder_loc,nedge_orient_loc, 
      call shape2DE('tria',Xi,norder_loc,nedge_orient_loc, 
     .              NrdofE,shapE_loc,curlE_loc)
      ShapE(1:2,1:NrdofE) = shapE_loc(1:2,1:NrdofE)
      CurlE(1:NrdofE) = curlE_loc(1:NrdofE)
c  
#endif
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
      real*8 , intent(in)  :: Xi(2)
      integer, intent(in)  :: Norder(5)
      integer, intent(out) :: NrdofQ
      real*8 , intent(out) :: ShapQ(MAXquadQ)
c
#if SHAPE_MODE==0
c
      write(*,*) 'shape2Q:This has NOT been implemented with SHAPE=0.'
      write(*,*) '        Please switch to SHAPE=1.'
      stop
c
#else
c
      call shape2DQ(Type,Xi,Norder, NrdofQ,ShapQ)
c
#endif
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
c         Nord          - polynomial order of face node (H1 sense)
c       out:
c         NrdofQ        - number of dof
c         ShapQ         - values of the shape functions at the point
c
c----------------------------------------------------------------------
c
      subroutine shapeQq(Xi,Nord, NrdofQ,ShapQ)
c
      use parameters , only : MAXquadQ
c
      implicit none
      real*8 , intent(in)  :: Xi(2)
      integer, intent(in)  :: Nord
      integer, intent(out) :: NrdofQ
      real*8 , intent(out) :: ShapQ(MAXquadQ)
c
c  ...locals
      integer :: norder(5)
c
      norder(1:4)=1
      norder(5) = Nord
      call shape2Q('quad',Xi,norder, NrdofQ,ShapQ)
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
c         Norder        - polynomial order  (H1 sense)
c       out:
c         NrdofQ        - number of dof
c         ShapQ         - values of the shape functions at the point
c
c----------------------------------------------------------------------
c
      subroutine shapeQt(Xi,Norder, NrdofQ,ShapQ)
c
      use parameters , only : MAXtriaQ,MAXquadQ
c
      implicit none
      real*8 , intent(in)  :: Xi(2)
      integer, intent(in)  :: Norder(4)
      integer, intent(out) :: NrdofQ
      real*8 , intent(out) :: ShapQ(MAXtriaQ)
c
c  ...locals
      integer :: norder_loc(5)
      real*8  :: shapQ_loc(MAXquadQ)
c
      norder_loc(1:4) = Norder(1:4); norder_loc(5)=11
      call shape2Q('tria',Xi,norder_loc, NrdofQ,shapQ_loc)
      ShapQ(1:NrdofQ) = shapQ_loc(1:NrdofQ)
c
      end subroutine shapeQt
c

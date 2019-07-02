c Routines:
c  - shape2DHBrokenTri
c  - shape2DEBrokenTri
c  - shape2DVBrokenTri
c  - shape2DQBrokenTri
c--------------------------------------------------------------------
c
c     routine name      - shape2DHBrokenTri
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 2D triangle
c                         BROKEN H1 shape functions
c
c     arguments:
c
c     in:
c          X            - master triangle coordinates from (0,1)^2
c          NordM        - polynomial order for middle node (H1 sense)
c          Nsize        - relevant sizes of local arrays
c
c     out:
c          NrdofH       - number of dof
c          ShapH        - values of the shape functions at the point
c          GradH        - gradients of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape2DHBrokenTri(Xi,NordM,Nsize, NrdofH,ShapH,GradH)
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofH
      integer :: noriE(3),norder(4)
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapH(Nsize(2))
      double precision, intent(out) :: GradH(1:2,Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      noriE(1:3)=0
      norder(1:4)=NordM
      call shape2DHTri(Xi,norder,noriE,Nsize, NrdofH,ShapH,GradH)
c
c
c  ...Option 2: Write more efficient routine for enriched functions
c
c
      end subroutine shape2DHBrokenTri
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape2DEBrokenTri
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 2D triangle
c                         BROKEN H(curl) shape functions
c
c     arguments:
c
c     in:
c          X            - master triangle coordinates from (0,1)^2
c          NordM        - polynomial order for middle node (H1 sense)
c          Nsize        - relevant sizes of local arrays
c
c     out:
c          NrdofE       - number of dof
c          ShapE        - values of the shape functions at the point
c          CurlE        - curl of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape2DEBrokenTri(Xi,NordM,Nsize, NrdofE,ShapE,CurlE)
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofE
      integer :: noriE(3),norder(4)
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapE(1:2,Nsize(2))
      double precision, intent(out) :: CurlE(Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      noriE(1:3)=0
      norder(1:4)=NordM
      call shape2DETri(Xi,norder,noriE,Nsize, NrdofE,ShapE,CurlE)
c
c
c  ...Option 2: Write more efficient routine for enriched functions
c
c
      end subroutine shape2DEBrokenTri
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape2DVBrokenTri
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 2D triangle
c                         BROKEN H(div) shape functions
c
c     arguments:
c
c     in:
c          X            - master triangle coordinates from (0,1)^2
c          NordM        - polynomial order for middle node (H1 sense)
c          Nsize        - relevant sizes of local arrays
c
c     out:
c          NrdofV       - number of dof
c          ShapV        - values of the shape functions at the point
c          DivV         - divergence of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape2DVBrokenTri(Xi,NordM,Nsize, NrdofV,ShapV,DivV)
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofV
      integer :: noriE(3),norder(4)
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapV(1:2,Nsize(2))
      double precision, intent(out) :: DivV(Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      noriE(1:3)=0
      norder(1:4)=NordM
      call shape2DVTri(Xi,norder,noriE,Nsize, NrdofV,ShapV,DivV)
c
c  ...Option 2: Write more efficient routine for enriched functions
c
c
      end subroutine shape2DVBrokenTri
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape2DQBrokenTri
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 2D triangle
c                         BROKEN L2 shape functions
c
c     arguments:
c
c     in:
c          X            - master triangle coordinates from (0,1)^2
c          NordM        - polynomial order for middle node (H1 sense)
c          Nsize        - relevant sizes of local arrays
c
c     out:
c          NrdofQ       - number of dof
c          ShapQ        - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape2DQBrokenTri(Xi,NordM,Nsize, NrdofQ,ShapQ)
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofQ
      integer :: norder(4)
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapQ(Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      norder(1:3)=1
      norder(4)=NordM
      call shape2DQTri(Xi,norder,Nsize, NrdofQ,ShapQ)
c
c  ...Option 2: Write a separate routine for enriched functions
c
c
      end subroutine shape2DQBrokenTri
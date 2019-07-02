c Routines:
c  - shape3DHBrokenTet
c  - shape3DEBrokenTet
c  - shape3DVBrokenTet
c  - shape3DQBrokenTet
c--------------------------------------------------------------------
c
c     routine name      - shape3DHBrokenTet
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 3D tetrahedron
c                         BROKEN H1 shape functions
c
c     arguments:
c
c     in:
c          X            - master tetrahedron coordinates from (0,1)^3
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
      subroutine shape3DHBrokenTet(Xi,NordM,Nsize, NrdofH,ShapH,GradH)
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofH
      integer :: noriE(6),noriF(4),norder(11)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapH(Nsize(2))
      double precision, intent(out) :: GradH(1:3,Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      noriE(1:6)=0
      noriF(1:4)=0
      norder(1:11)=NordM
      call shape3DHTet(Xi,norder,noriE,noriF,Nsize,
     .                                             NrdofH,ShapH,GradH)
c
c
c  ...Option 2: Write more efficient routine for enriched functions
c
c
      end subroutine shape3DHBrokenTet
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape3DEBrokenTet
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 3D tetrahedron
c                         BROKEN H(curl) shape functions
c
c     arguments:
c
c     in:
c          X            - master tetrahedron coordinates from (0,1)^3
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
      subroutine shape3DEBrokenTet(Xi,NordM,Nsize, NrdofE,ShapE,CurlE)
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofE
      integer :: noriE(6),noriF(4),norder(11)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapE(1:3,Nsize(2))
      double precision, intent(out) :: CurlE(1:3,Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      noriE(1:6)=0
      noriF(1:4)=0
      norder(1:11)=NordM
      call shape3DETet(Xi,norder,noriE,noriF,Nsize,
     .                                             NrdofE,ShapE,CurlE)
c
c
c  ...Option 2: Write more efficient routine for enriched functions
c
c
      end subroutine shape3DEBrokenTet
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape3DVBrokenTet
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 3D tetrahedron
c                         BROKEN H(div) shape functions
c
c     arguments:
c
c     in:
c          X            - master tetrahedron coordinates from (0,1)^3
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
      subroutine shape3DVBrokenTet(Xi,NordM,Nsize, NrdofV,ShapV,DivV)
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofV
      integer :: noriF(4),norder(11)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapV(1:3,Nsize(2))
      double precision, intent(out) :: DivV(Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      noriF(1:4)=0
      norder(1:6)=1
      norder(7:11)=NordM
      call shape3DVTet(Xi,norder,noriF,Nsize, NrdofV,ShapV,DivV)
c
c  ...Option 2: Write more efficient routine for enriched functions
c
c
      end subroutine shape3DVBrokenTet
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape3DQBrokenTet
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 3D tetrahedron
c                         BROKEN L2 shape functions
c
c     arguments:
c
c     in:
c          X            - master tetrahedron coordinates from (0,1)^3
c          NordM        - polynomial order for middle node (H1 sense)
c          Nsize        - relevant sizes of local arrays
c
c     out:
c          NrdofQ       - number of dof
c          ShapQ        - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape3DQBrokenTet(Xi,NordM,Nsize, NrdofQ,ShapQ)
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofQ
      integer :: norder(11)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapQ(Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      norder(1:10)=1
      norder(11)=NordM
      call shape3DQTet(Xi,norder,Nsize, NrdofQ,ShapQ)
c
c  ...Option 2: Write a separate routine for enriched functions
c
c
      end subroutine shape3DQBrokenTet
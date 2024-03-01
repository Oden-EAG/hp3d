! Routines:
!  - shape3DHBrokenTet
!  - shape3DEBrokenTet
!  - shape3DVBrokenTet
!  - shape3DQBrokenTet
!--------------------------------------------------------------------
!
!     routine name      - shape3DHBrokenTet
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!     purpose:          - routine returns values of 3D tetrahedron
!                         BROKEN H1 shape functions
!
!     arguments:
!
!     in:
!          X            - master tetrahedron coordinates from (0,1)^3
!          NordM        - polynomial order for middle node (H1 sense)
!          Nsize        - relevant sizes of local arrays
!
!     out:
!          NrdofH       - number of dof
!          ShapH        - values of the shape functions at the point
!          GradH        - gradients of the shape functions
!
!-----------------------------------------------------------------------
!
      subroutine shape3DHBrokenTet(Xi,NordM,Nsize, NrdofH,ShapH,GradH)
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofH
      integer :: noriE(6),noriF(4),norder(11)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapH(Nsize(2))
      double precision, intent(out) :: GradH(1:3,Nsize(2))
!
!  ...Option 1: Simply call the usual shape functions with enrichment
      noriE(1:6)=0
      noriF(1:4)=0
      norder(1:11)=NordM
      call shape3DHTet(Xi,norder,noriE,noriF,Nsize, &
                                                   NrdofH,ShapH,GradH)
!
!
!  ...Option 2: Write more efficient routine for enriched functions
!
!
      end subroutine shape3DHBrokenTet
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape3DEBrokenTet
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!     purpose:          - routine returns values of 3D tetrahedron
!                         BROKEN H(curl) shape functions
!
!     arguments:
!
!     in:
!          X            - master tetrahedron coordinates from (0,1)^3
!          NordM        - polynomial order for middle node (H1 sense)
!          Nsize        - relevant sizes of local arrays
!
!     out:
!          NrdofE       - number of dof
!          ShapE        - values of the shape functions at the point
!          CurlE        - curl of the shape functions
!
!-----------------------------------------------------------------------
!
      subroutine shape3DEBrokenTet(Xi,NordM,Nsize, NrdofE,ShapE,CurlE)
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofE
      integer :: noriE(6),noriF(4),norder(11)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapE(1:3,Nsize(2))
      double precision, intent(out) :: CurlE(1:3,Nsize(2))
!
!  ...Option 1: Simply call the usual shape functions with enrichment
      noriE(1:6)=0
      noriF(1:4)=0
      norder(1:11)=NordM
      call shape3DETet(Xi,norder,noriE,noriF,Nsize, &
                                                   NrdofE,ShapE,CurlE)
!
!
!  ...Option 2: Write more efficient routine for enriched functions
!
!
      end subroutine shape3DEBrokenTet
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape3DVBrokenTet
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!     purpose:          - routine returns values of 3D tetrahedron
!                         BROKEN H(div) shape functions
!
!     arguments:
!
!     in:
!          X            - master tetrahedron coordinates from (0,1)^3
!          NordM        - polynomial order for middle node (H1 sense)
!          Nsize        - relevant sizes of local arrays
!
!     out:
!          NrdofV       - number of dof
!          ShapV        - values of the shape functions at the point
!          DivV         - divergence of the shape functions
!
!-----------------------------------------------------------------------
!
      subroutine shape3DVBrokenTet(Xi,NordM,Nsize, NrdofV,ShapV,DivV)
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofV
      integer :: noriF(4),norder(11)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapV(1:3,Nsize(2))
      double precision, intent(out) :: DivV(Nsize(2))
!
!  ...Option 1: Simply call the usual shape functions with enrichment
      noriF(1:4)=0
      norder(1:6)=1
      norder(7:11)=NordM
      call shape3DVTet(Xi,norder,noriF,Nsize, NrdofV,ShapV,DivV)
!
!  ...Option 2: Write more efficient routine for enriched functions
!
!
      end subroutine shape3DVBrokenTet
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape3DQBrokenTet
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!     purpose:          - routine returns values of 3D tetrahedron
!                         BROKEN L2 shape functions
!
!     arguments:
!
!     in:
!          X            - master tetrahedron coordinates from (0,1)^3
!          NordM        - polynomial order for middle node (H1 sense)
!          Nsize        - relevant sizes of local arrays
!
!     out:
!          NrdofQ       - number of dof
!          ShapQ        - values of the shape functions at the point
!
!-----------------------------------------------------------------------
!
      subroutine shape3DQBrokenTet(Xi,NordM,Nsize, NrdofQ,ShapQ)
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofQ
      integer :: norder(11)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapQ(Nsize(2))
!
!  ...Option 1: Simply call the usual shape functions with enrichment
      norder(1:10)=1
      norder(11)=NordM
      call shape3DQTet(Xi,norder,Nsize, NrdofQ,ShapQ)
!
!  ...Option 2: Write a separate routine for enriched functions
!
!
      end subroutine shape3DQBrokenTet

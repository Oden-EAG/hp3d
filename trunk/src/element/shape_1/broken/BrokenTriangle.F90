! Routines:
!  - shape2DHBrokenTri
!  - shape2DEBrokenTri
!  - shape2DVBrokenTri
!  - shape2DQBrokenTri
!--------------------------------------------------------------------
!
!     routine name      - shape2DHBrokenTri
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!     purpose:          - routine returns values of 2D triangle
!                         BROKEN H1 shape functions
!
!     arguments:
!
!     in:
!          X            - master triangle coordinates from (0,1)^2
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
   subroutine shape2DHBrokenTri(Xi,NordM,Nsize, NrdofH,ShapH,GradH)
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofH
      integer :: noriE(3),norder(4)
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapH(Nsize(2))
      double precision, intent(out) :: GradH(1:2,Nsize(2))
!
!  ...Option 1: Simply call the usual shape functions with enrichment
      noriE(1:3)=0
      norder(1:4)=NordM
      call shape2DHTri(Xi,norder,noriE,Nsize, NrdofH,ShapH,GradH)
!
!
!  ...Option 2: Write more efficient routine for enriched functions
!
!
   end subroutine shape2DHBrokenTri
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape2DEBrokenTri
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!     purpose:          - routine returns values of 2D triangle
!                         BROKEN H(curl) shape functions
!
!     arguments:
!
!     in:
!          X            - master triangle coordinates from (0,1)^2
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
   subroutine shape2DEBrokenTri(Xi,NordM,Nsize, NrdofE,ShapE,CurlE)
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofE
      integer :: noriE(3),norder(4)
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapE(1:2,Nsize(2))
      double precision, intent(out) :: CurlE(Nsize(2))
!
!  ...Option 1: Simply call the usual shape functions with enrichment
      noriE(1:3)=0
      norder(1:4)=NordM
      call shape2DETri(Xi,norder,noriE,Nsize, NrdofE,ShapE,CurlE)
!
!
!  ...Option 2: Write more efficient routine for enriched functions
!
!
   end subroutine shape2DEBrokenTri
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape2DVBrokenTri
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!     purpose:          - routine returns values of 2D triangle
!                         BROKEN H(div) shape functions
!
!     arguments:
!
!     in:
!          X            - master triangle coordinates from (0,1)^2
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
   subroutine shape2DVBrokenTri(Xi,NordM,Nsize, NrdofV,ShapV,DivV)
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofV
      integer :: noriE(3),norder(4)
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapV(1:2,Nsize(2))
      double precision, intent(out) :: DivV(Nsize(2))
!
!  ...Option 1: Simply call the usual shape functions with enrichment
      noriE(1:3)=0
      norder(1:4)=NordM
      call shape2DVTri(Xi,norder,noriE,Nsize, NrdofV,ShapV,DivV)
!
!  ...Option 2: Write more efficient routine for enriched functions
!
!
   end subroutine shape2DVBrokenTri
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape2DQBrokenTri
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!     purpose:          - routine returns values of 2D triangle
!                         BROKEN L2 shape functions
!
!     arguments:
!
!     in:
!          X            - master triangle coordinates from (0,1)^2
!          NordM        - polynomial order for middle node (H1 sense)
!          Nsize        - relevant sizes of local arrays
!
!     out:
!          NrdofQ       - number of dof
!          ShapQ        - values of the shape functions at the point
!
!-----------------------------------------------------------------------
!
   subroutine shape2DQBrokenTri(Xi,NordM,Nsize, NrdofQ,ShapQ)
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofQ
      integer :: norder(4)
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapQ(Nsize(2))
!
!  ...Option 1: Simply call the usual shape functions with enrichment
      norder(1:3)=1
      norder(4)=NordM
      call shape2DQTri(Xi,norder,Nsize, NrdofQ,ShapQ)
!
!  ...Option 2: Write a separate routine for enriched functions
!
!
   end subroutine shape2DQBrokenTri

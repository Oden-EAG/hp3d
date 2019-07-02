c Routines:
c  - shape3DHBrokenPris
c  - shape3DEBrokenPris
c  - shape3DVBrokenPris
c  - shape3DQBrokenPris
c--------------------------------------------------------------------
c
c     routine name      - shape3DHBrokenPris
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 3D prism
c                         BROKEN H1 shape functions
c
c     arguments:
c
c     in:
c          X            - master prism coordinates from (0,1)^3
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
      subroutine shape3DHBrokenPris(Xi,NordM,Nsize, NrdofH,ShapH,GradH)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofH
      integer :: noriE(9),noriF(5),norder(15),nordB(2)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapH(Nsize(2))
      double precision, intent(out) :: GradH(1:3,Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      noriE(1:9)=0
      noriF(1:5)=0
      call decod(NordM,MODORDER,2, nordB)
      norder(1:6)=nordB(1)
      norder(7:9)=nordB(2)
      norder(10:11)=nordB(1)
      norder(12:14)=NordM
      norder(15)=NordM
      call shape3DHPris(Xi,norder,noriE,noriF,Nsize,
     .                                             NrdofH,ShapH,GradH)
c
c
c  ...Option 2: Write more efficient routine for enriched functions
c
c
      end subroutine shape3DHBrokenPris
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape3DEBrokenPris
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 3D prism
c                         BROKEN H(curl) shape functions
c
c     arguments:
c
c     in:
c          X            - master prism coordinates from (0,1)^3
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
      subroutine shape3DEBrokenPris(Xi,NordM,Nsize, NrdofE,ShapE,CurlE)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofE
      integer :: noriE(9),noriF(5),norder(15),nordB(2)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapE(1:3,Nsize(2))
      double precision, intent(out) :: CurlE(1:3,Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      noriE(1:9)=0
      noriF(1:5)=0
      call decod(NordM,MODORDER,2, nordB)
      norder(1:6)=nordB(1)
      norder(7:9)=nordB(2)
      norder(10:11)=nordB(1)
      norder(12:14)=NordM
      norder(15)=NordM
      call shape3DEPris(Xi,norder,noriE,noriF,Nsize,
     .                                             NrdofE,ShapE,CurlE)
c
c
c  ...Option 2: Write more efficient routine for enriched functions
c
c
      end subroutine shape3DEBrokenPris
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape3DVBrokenPris
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 3D prism
c                         BROKEN H(div) shape functions
c
c     arguments:
c
c     in:
c          X            - master prism coordinates from (0,1)^3
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
      subroutine shape3DVBrokenPris(Xi,NordM,Nsize, NrdofV,ShapV,DivV)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofV
      integer :: noriF(5),norder(15),nordB(2)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapV(1:3,Nsize(2))
      double precision, intent(out) :: DivV(Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      noriF(1:5)=0
      call decod(NordM,MODORDER,2, nordB)
      norder(1:9)=1
      norder(10:11)=nordB(1)
      norder(12:14)=NordM
      norder(15)=NordM
      call shape3DVPris(Xi,norder,noriF,Nsize, NrdofV,ShapV,DivV)
c
c  ...Option 2: Write more efficient routine for enriched functions
c
c
      end subroutine shape3DVBrokenPris
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape3DQBrokenPris
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 3D prism
c                         BROKEN L2 shape functions
c
c     arguments:
c
c     in:
c          X            - master prism coordinates from (0,1)^3
c          NordM        - polynomial order for middle node (H1 sense)
c          Nsize        - relevant sizes of local arrays
c
c     out:
c          NrdofQ       - number of dof
c          ShapQ        - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape3DQBrokenPris(Xi,NordM,Nsize, NrdofQ,ShapQ)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofQ
      integer :: norder(15),nordF
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapQ(Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      norder(1:11)=1
      call encod((/1,1/),MODORDER,2, nordF)
      norder(12:14)=nordF
      norder(15)=NordM
      call shape3DQPris(Xi,norder,Nsize, NrdofQ,ShapQ)
c
c  ...Option 2: Write a separate routine for enriched functions
c
c
      end subroutine shape3DQBrokenPris
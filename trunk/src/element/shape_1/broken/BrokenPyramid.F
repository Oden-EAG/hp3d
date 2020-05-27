c Routines:
c  - shape3DHBrokenPyra
c  - shape3DEBrokenPyra
c  - shape3DVBrokenPyra
c  - shape3DQBrokenPyra
c--------------------------------------------------------------------
c
c     routine name      - shape3DHBrokenPyra
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 3D pyramid
c                         BROKEN H1 shape functions
c
c     arguments:
c
c     in:
c          X            - master pyramid coordinates from (0,1)^3
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
      subroutine shape3DHBrokenPyra(Xi,NordM,Nsize, NrdofH,ShapH,GradH)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofH
      integer :: noriE(8),noriF(5),norder(14),nordF
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapH(Nsize(2))
      double precision, intent(out) :: GradH(1:3,Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      noriE(1:8)=0
      noriF(1:5)=0
      norder(1:8)=NordM
      call encod((/NordM,NordM/),MODORDER,2, nordF)
      norder(9)=nordF
      norder(10:13)=NordM
      norder(14)=NordM
      call shape3DHPyra(Xi,norder,noriE,noriF,Nsize,
     .                                             NrdofH,ShapH,GradH)
c
c
c  ...Option 2: Write more efficient routine for enriched functions
c
c
      end subroutine shape3DHBrokenPyra
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape3DEBrokenPyra
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 3D pyramid
c                         BROKEN H(curl) shape functions
c
c     arguments:
c
c     in:
c          X            - master pyramid coordinates from (0,1)^3
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
      subroutine shape3DEBrokenPyra(Xi,NordM,Nsize, NrdofE,ShapE,CurlE)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofE
      integer :: noriE(8),noriF(5),norder(14),nordF
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapE(1:3,Nsize(2))
      double precision, intent(out) :: CurlE(1:3,Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      noriE(1:8)=0
      noriF(1:5)=0
      norder(1:8)=NordM
      call encod((/NordM,NordM/),MODORDER,2, nordF)
      norder(9)=nordF
      norder(10:13)=NordM
      norder(14)=NordM
      call shape3DEPyra(Xi,norder,noriE,noriF,Nsize,
     .                                             NrdofE,ShapE,CurlE)
c
c
c  ...Option 2: Write more efficient routine for enriched functions
c
c
      end subroutine shape3DEBrokenPyra
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape3DVBrokenPyra
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 3D pyramid
c                         BROKEN H(div) shape functions
c
c     arguments:
c
c     in:
c          X            - master pyramid coordinates from (0,1)^3
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
      subroutine shape3DVBrokenPyra(Xi,NordM,Nsize, NrdofV,ShapV,DivV)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofV
      integer :: noriF(5),norder(14),nordF
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapV(1:3,Nsize(2))
      double precision, intent(out) :: DivV(Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      noriF(1:5)=0
      norder(1:8)=1
      call encod((/NordM,NordM/),MODORDER,2, nordF)
      norder(9)=nordF
      norder(10:13)=NordM
      norder(14)=NordM
      call shape3DVPyra(Xi,norder,noriF,Nsize, NrdofV,ShapV,DivV)
c
c  ...Option 2: Write more efficient routine for enriched functions
c
c
      end subroutine shape3DVBrokenPyra
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape3DQBrokenPyra
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 3D pyramid
c                         BROKEN L2 shape functions
c
c     arguments:
c
c     in:
c          X            - master pyramid coordinates from (0,1)^3
c          NordM        - polynomial order for middle node (H1 sense)
c          Nsize        - relevant sizes of local arrays
c
c     out:
c          NrdofQ       - number of dof
c          ShapQ        - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape3DQBrokenPyra(Xi,NordM,Nsize, NrdofQ,ShapQ)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofQ
      integer :: norder(14),nordF
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapQ(Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      norder(1:8)=1
      call encod((/1,1/),MODORDER,2, nordF)
      norder(9)=nordF
      norder(10:13)=1
      norder(14)=NordM
      call shape3DQPyra(Xi,norder,Nsize, NrdofQ,ShapQ)
c
c  ...Option 2: Write a separate routine for enriched functions
c
c
      end subroutine shape3DQBrokenPyra

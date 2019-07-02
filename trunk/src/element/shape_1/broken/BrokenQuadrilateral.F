c Routines:
c  - shape2DHBrokenQuad
c  - shape2DEBrokenQuad
c  - shape2DVBrokenQuad
c  - shape2DQBrokenQuad
c--------------------------------------------------------------------
c
c     routine name      - shape2DHBrokenQuad
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 2D quadrilateral
c                         BROKEN H1 shape functions
c
c     arguments:
c
c     in:
c          X            - master quadrilateral coordinates from (0,1)^2
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
      subroutine shape2DHBrokenQuad(Xi,NordM,Nsize, NrdofH,ShapH,GradH)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofH
      integer :: noriE(4),norder(5),nordB(2)
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapH(Nsize(2))
      double precision, intent(out) :: GradH(1:2,Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      noriE(1:4)=0
      call decod(NordM,MODORDER,2, nordB)
      norder(1:4)=(/nordB(1),nordB(2),nordB(1),nordB(2)/)
      norder(5)=NordM
      call shape2DHQuad(Xi,norder,noriE,Nsize, NrdofH,ShapH,GradH)
c
c
c  ...Option 2: Write more efficient routine for enriched functions
c
c
      end subroutine shape2DHBrokenQuad
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape2DEBrokenQuad
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 2D quadrilateral
c                         BROKEN H(curl) shape functions
c
c     arguments:
c
c     in:
c          X            - master quadrilateral coordinates from (0,1)^2
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
      subroutine shape2DEBrokenQuad(Xi,NordM,Nsize, NrdofE,ShapE,CurlE)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofE
      integer :: noriE(4),norder(5),nordB(2)
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapE(1:2,Nsize(2))
      double precision, intent(out) :: CurlE(Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      noriE(1:4)=0
      call decod(NordM,MODORDER,2, nordB)
      norder(1:4)=(/nordB(1),nordB(2),nordB(1),nordB(2)/)
      norder(5)=NordM
      call shape2DEQuad(Xi,norder,noriE,Nsize, NrdofE,ShapE,CurlE)
c
c
c  ...Option 2: Write more efficient routine for enriched functions
c
c
      end subroutine shape2DEBrokenQuad
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape2DVBrokenQuad
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 2D quadrilateral
c                         BROKEN H(div) shape functions
c
c     arguments:
c
c     in:
c          X            - master quadrilateral coordinates from (0,1)^2
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
      subroutine shape2DVBrokenQuad(Xi,NordM,Nsize, NrdofV,ShapV,DivV)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofV
      integer :: noriE(4),norder(5),nordB(2)
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapV(1:2,Nsize(2))
      double precision, intent(out) :: DivV(Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      noriE(1:4)=0
      call decod(NordM,MODORDER,2, nordB)
      norder(1:4)=(/nordB(1),nordB(2),nordB(1),nordB(2)/)
      norder(5)=NordM
      call shape2DVQuad(Xi,norder,noriE,Nsize, NrdofV,ShapV,DivV)
c
c  ...Option 2: Write more efficient routine for enriched functions
c
c
      end subroutine shape2DVBrokenQuad
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape2DQBrokenQuad
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 2D quadrilateral
c                         BROKEN L2 shape functions
c
c     arguments:
c
c     in:
c          X            - master quadrilateral coordinates from (0,1)^2
c          NordM        - polynomial order for middle node (H1 sense)
c          Nsize        - relevant sizes of local arrays
c
c     out:
c          NrdofQ       - number of dof
c          ShapQ        - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape2DQBrokenQuad(Xi,NordM,Nsize, NrdofQ,ShapQ)
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofQ
      integer :: norder(5)
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapQ(Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      norder(1:4)=1
      norder(5)=NordM
      call shape2DQQuad(Xi,norder,Nsize, NrdofQ,ShapQ)
c
c  ...Option 2: Write a separate routine for enriched functions
c
c
      end subroutine shape2DQBrokenQuad

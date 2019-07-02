c Routines:
c  - shape1DHBrokenSeg
c  - shape1DQBrokenSeg
c--------------------------------------------------------------------
c
c     routine name      - shape1DHBrokenSeg
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 1D segment
c                         BROKEN H1 shape functions
c
c     arguments:
c
c     in:
c          X            - master segment coordinates from (0,1)
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
      subroutine shape1DHBrokenSeg(Xi,NordM,Nsize, NrdofH,ShapH,GradH)
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofH
      double precision, intent(in)  :: Xi
      double precision, intent(out) :: ShapH(Nsize(2)),GradH(Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      call shape1DHSeg(Xi,NordM,Nsize, NrdofH,ShapH,GradH)
c
c
c  ...Option 2: Write more efficient routine for enriched functions
c
c
      end subroutine shape1DHBrokenSeg
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape1DQBrokenSeg
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 1D segment
c                         BROKEN L2 shape functions
c
c     arguments:
c
c     in:
c          X            - master segment coordinates from (0,1)
c          NordM        - polynomial order for middle node (H1 sense)
c          Nsize        - relevant sizes of local arrays
c
c     out:
c          NrdofQ       - number of dof
c          ShapQ        - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape1DQBrokenSeg(Xi,NordM,Nsize, NrdofQ,ShapQ)
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofQ
      double precision, intent(in)  :: Xi
      double precision, intent(out) :: ShapQ(Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      call shape1DQSeg(Xi,NordM,Nsize, NrdofQ,ShapQ)
c
c  ...Option 2: Write a separate routine for enriched functions
c
c
      end subroutine shape1DQBrokenSeg

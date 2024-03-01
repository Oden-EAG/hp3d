! Routines:
!  - shape1DHBrokenSeg
!  - shape1DQBrokenSeg
!--------------------------------------------------------------------
!
!     routine name      - shape1DHBrokenSeg
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!     purpose:          - routine returns values of 1D segment
!                         BROKEN H1 shape functions
!
!     arguments:
!
!     in:
!          X            - master segment coordinates from (0,1)
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
   subroutine shape1DHBrokenSeg(Xi,NordM,Nsize, NrdofH,ShapH,GradH)
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofH
      double precision, intent(in)  :: Xi
      double precision, intent(out) :: ShapH(Nsize(2)),GradH(Nsize(2))
!
!  ...Option 1: Simply call the usual shape functions with enrichment
      call shape1DHSeg(Xi,NordM,Nsize, NrdofH,ShapH,GradH)
!
!
!  ...Option 2: Write more efficient routine for enriched functions
!
!
   end subroutine shape1DHBrokenSeg
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape1DQBrokenSeg
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!     purpose:          - routine returns values of 1D segment
!                         BROKEN L2 shape functions
!
!     arguments:
!
!     in:
!          X            - master segment coordinates from (0,1)
!          NordM        - polynomial order for middle node (H1 sense)
!          Nsize        - relevant sizes of local arrays
!
!     out:
!          NrdofQ       - number of dof
!          ShapQ        - values of the shape functions at the point
!
!-----------------------------------------------------------------------
!
   subroutine shape1DQBrokenSeg(Xi,NordM,Nsize, NrdofQ,ShapQ)
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofQ
      double precision, intent(in)  :: Xi
      double precision, intent(out) :: ShapQ(Nsize(2))
!
!  ...Option 1: Simply call the usual shape functions with enrichment
      call shape1DQSeg(Xi,NordM,Nsize, NrdofQ,ShapQ)
!
!  ...Option 2: Write a separate routine for enriched functions
!
!
   end subroutine shape1DQBrokenSeg

! Routines:
!  - shape3DHBrokenPyra
!  - shape3DEBrokenPyra
!  - shape3DVBrokenPyra
!  - shape3DQBrokenPyra
!--------------------------------------------------------------------
!
!     routine name      - shape3DHBrokenPyra
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!> @brief         - routine returns values of 3D pyramid
!                         BROKEN H1 shape functions
!
!     arguments:
!
!     in:
!          X            - master pyramid coordinates from (0,1)^3
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
   subroutine shape3DHBrokenPyra(Xi,NordM,Nsize, NrdofH,ShapH,GradH)
!
      use parameters , only : MODORDER
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofH
      integer :: noriE(8),noriF(5),norder(14),nordF
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapH(Nsize(2))
      double precision, intent(out) :: GradH(1:3,Nsize(2))
!
!  ...Option 1: Simply call the usual shape functions with enrichment
      noriE(1:8)=0
      noriF(1:5)=0
      norder(1:8)=NordM
      call encod((/NordM,NordM/),MODORDER,2, nordF)
      norder(9)=nordF
      norder(10:13)=NordM
      norder(14)=NordM
      call shape3DHPyra(Xi,norder,noriE,noriF,Nsize, &
                                                   NrdofH,ShapH,GradH)
!
!
!  ...Option 2: Write more efficient routine for enriched functions
!
!
   end subroutine shape3DHBrokenPyra
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape3DEBrokenPyra
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!> @brief         - routine returns values of 3D pyramid
!                         BROKEN H(curl) shape functions
!
!     arguments:
!
!     in:
!          X            - master pyramid coordinates from (0,1)^3
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
   subroutine shape3DEBrokenPyra(Xi,NordM,Nsize, NrdofE,ShapE,CurlE)
!
      use parameters , only : MODORDER
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofE
      integer :: noriE(8),noriF(5),norder(14),nordF
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapE(1:3,Nsize(2))
      double precision, intent(out) :: CurlE(1:3,Nsize(2))
!
!  ...Option 1: Simply call the usual shape functions with enrichment
      noriE(1:8)=0
      noriF(1:5)=0
      norder(1:8)=NordM
      call encod((/NordM,NordM/),MODORDER,2, nordF)
      norder(9)=nordF
      norder(10:13)=NordM
      norder(14)=NordM
      call shape3DEPyra(Xi,norder,noriE,noriF,Nsize, &
                                                   NrdofE,ShapE,CurlE)
!
!
!  ...Option 2: Write more efficient routine for enriched functions
!
!
   end subroutine shape3DEBrokenPyra
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape3DVBrokenPyra
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!> @brief         - routine returns values of 3D pyramid
!                         BROKEN H(div) shape functions
!
!     arguments:
!
!     in:
!          X            - master pyramid coordinates from (0,1)^3
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
   subroutine shape3DVBrokenPyra(Xi,NordM,Nsize, NrdofV,ShapV,DivV)
!
      use parameters , only : MODORDER
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofV
      integer :: noriF(5),norder(14),nordF
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapV(1:3,Nsize(2))
      double precision, intent(out) :: DivV(Nsize(2))
!
!  ...Option 1: Simply call the usual shape functions with enrichment
      noriF(1:5)=0
      norder(1:8)=1
      call encod((/NordM,NordM/),MODORDER,2, nordF)
      norder(9)=nordF
      norder(10:13)=NordM
      norder(14)=NordM
      call shape3DVPyra(Xi,norder,noriF,Nsize, NrdofV,ShapV,DivV)
!
!  ...Option 2: Write more efficient routine for enriched functions
!
!
   end subroutine shape3DVBrokenPyra
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape3DQBrokenPyra
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!> @brief         - routine returns values of 3D pyramid
!                         BROKEN L2 shape functions
!
!     arguments:
!
!     in:
!          X            - master pyramid coordinates from (0,1)^3
!          NordM        - polynomial order for middle node (H1 sense)
!          Nsize        - relevant sizes of local arrays
!
!     out:
!          NrdofQ       - number of dof
!          ShapQ        - values of the shape functions at the point
!
!-----------------------------------------------------------------------
!
   subroutine shape3DQBrokenPyra(Xi,NordM,Nsize, NrdofQ,ShapQ)
!
      use parameters , only : MODORDER
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofQ
      integer :: norder(14),nordF
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapQ(Nsize(2))
!
!  ...Option 1: Simply call the usual shape functions with enrichment
      norder(1:8)=1
      call encod((/1,1/),MODORDER,2, nordF)
      norder(9)=nordF
      norder(10:13)=1
      norder(14)=NordM
      call shape3DQPyra(Xi,norder,Nsize, NrdofQ,ShapQ)
!
!  ...Option 2: Write a separate routine for enriched functions
!
!
   end subroutine shape3DQBrokenPyra

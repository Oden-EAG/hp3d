! Routines:
!  - shape2DHBrokenQuad
!  - shape2DEBrokenQuad
!  - shape2DVBrokenQuad
!  - shape2DQBrokenQuad
!--------------------------------------------------------------------
!
!     routine name      - shape2DHBrokenQuad
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!     purpose:          - routine returns values of 2D quadrilateral
!                         BROKEN H1 shape functions
!
!     arguments:
!
!     in:
!          X            - master quadrilateral coordinates from (0,1)^2
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
   subroutine shape2DHBrokenQuad(Xi,NordM,Nsize, NrdofH,ShapH,GradH)
!
      use parameters , only : MODORDER
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofH
      integer :: i,j,m
!      integer :: noriE(4),norder(5)
      integer :: nordB(2),ndofH(3)
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapH(Nsize(2))
      double precision, intent(out) :: GradH(1:2,Nsize(2))
      double precision :: shapH1(Nsize(1)+1),dshapH1(Nsize(1)+1)
      double precision :: shapH2(Nsize(1)+1),dshapH2(Nsize(1)+1)
!
!  ...Option 1: Simply call the usual shape functions with enrichment
!      noriE(1:4)=0
!      call decod(NordM,MODORDER,2, nordB)
!      norder(1:4)=(/nordB(1),nordB(2),nordB(1),nordB(2)/)
!      norder(5)=NordM
!      call shape2DHQuad(Xi,norder,noriE,Nsize, NrdofH,ShapH,GradH)
!
!
!  ...Option 2: Write more efficient routine for enriched functions
!
!  ...initiate counter for shape functions
      m=0
!
!  ...shape functions are tensor products of 1D shape functions
      call decod(NordM,MODORDER,2, nordB)
      call shape1HH(Xi(1),nordB(1), ndofH(1),shapH1,dshapH1)
      call shape1HH(Xi(2),nordB(2), ndofH(2),shapH2,dshapH2)
!
      do j=1,ndofH(2)
         do i=1,ndofH(1)
            m=m+1
            ShapH(m)   =  shapH1(i)* shapH2(j)
            GradH(1,m) = dshapH1(i)* shapH2(j)
            GradH(2,m) =  shapH1(i)*dshapH2(j)
         enddo
      enddo
!
!  ...give total degrees of freedom
      NrdofH = m
!
   end subroutine shape2DHBrokenQuad
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape2DEBrokenQuad
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!     purpose:          - routine returns values of 2D quadrilateral
!                         BROKEN H(curl) shape functions
!
!     arguments:
!
!     in:
!          X            - master quadrilateral coordinates from (0,1)^2
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
   subroutine shape2DEBrokenQuad(Xi,NordM,Nsize, NrdofE,ShapE,CurlE)
!
      use parameters , only : MODORDER
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofE
!      integer :: noriE(4),norder(5)
      integer :: i,j,m
      integer :: nordB(2),ndofH(3),ndofQ(3)
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapE(1:2,Nsize(2))
      double precision, intent(out) :: CurlE(Nsize(2))
      double precision :: shapH1(Nsize(1)+1),dshapH1(Nsize(1)+1)
      double precision :: shapH2(Nsize(1)+1),dshapH2(Nsize(1)+1)
      double precision :: shapQ1(Nsize(1))
      double precision :: shapQ2(Nsize(1))
!
!  ...Option 1: Simply call the usual shape functions with enrichment
!      noriE(1:4)=0
!      call decod(NordM,MODORDER,2, nordB)
!      norder(1:4)=(/nordB(1),nordB(2),nordB(1),nordB(2)/)
!      norder(5)=NordM
!      call shape2DEQuad(Xi,norder,noriE,Nsize, NrdofE,ShapE,CurlE)
!
!
!  ...Option 2: Write more efficient routine for enriched functions
!
!  ...initiate counter for shape functions
      m=0
!
!  ...shape functions are tensor products of 1D shape functions
      call decod(NordM,MODORDER,2, nordB)
      call shape1HH(Xi(1),nordB(1), ndofH(1),shapH1,dshapH1)
      call shape1HH(Xi(2),nordB(2), ndofH(2),shapH2,dshapH2)
      call shape1QQ(Xi(1),nordB(1), ndofQ(1),shapQ1)
      call shape1QQ(Xi(2),nordB(2), ndofQ(2),shapQ2)
!
!  ...shape functions with values along the x-axis
      do j=1,ndofH(2)
         do i=1,ndofQ(1)
            m=m+1
            ShapE(1,m) =  shapQ1(i)*shapH2(j)
            ShapE(2,m) =  0.d0
            CurlE(m) = -shapQ1(i)*dshapH2(j)
         enddo
      enddo
!
!  ...shape functions with values along the y-axis
      do j=1,ndofQ(2)
         do i=1,ndofH(1)
            m=m+1
            ShapE(1,m) =  0.d0
            ShapE(2,m) =  shapH1(i)*shapQ2(j)
            CurlE(m) = dshapH1(i)*shapQ2(j)
         enddo
      enddo
!
!  ...give total degrees of freedom
      NrdofE = m
!
!
   end subroutine shape2DEBrokenQuad
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape2DVBrokenQuad
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!     purpose:          - routine returns values of 2D quadrilateral
!                         BROKEN H(div) shape functions
!
!     arguments:
!
!     in:
!          X            - master quadrilateral coordinates from (0,1)^2
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
   subroutine shape2DVBrokenQuad(Xi,NordM,Nsize, NrdofV,ShapV,DivV)
!
      use parameters , only : MODORDER
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofV
      integer :: i,j,m
!      integer :: noriE(4),norder(5),
      integer :: nordB(2),ndofH(3),ndofQ(3)
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapV(1:2,Nsize(2))
      double precision, intent(out) :: DivV(Nsize(2))
      double precision :: shapH1(Nsize(1)+1),dshapH1(Nsize(1)+1)
      double precision :: shapH2(Nsize(1)+1),dshapH2(Nsize(1)+1)
      double precision :: shapQ1(Nsize(1))
      double precision :: shapQ2(Nsize(1))
!
!  ...Option 1: Simply call the usual shape functions with enrichment
!      noriE(1:4)=0
!      call decod(NordM,MODORDER,2, nordB)
!      norder(1:4)=(/nordB(1),nordB(2),nordB(1),nordB(2)/)
!      norder(5)=NordM
!      call shape2DVQuad(Xi,norder,noriE,Nsize, NrdofV,ShapV,DivV)
!
!  ...Option 2: Write more efficient routine for enriched functions
!
!  ...initiate counter for shape functions
      m=0
!
!  ...shape functions are tensor products of 1D shape functions
      call decod(NordM,MODORDER,2, nordB)
      call shape1HH(Xi(1),nordB(1), ndofH(1),shapH1,dshapH1)
      call shape1HH(Xi(2),nordB(2), ndofH(2),shapH2,dshapH2)
      call shape1QQ(Xi(1),nordB(1), ndofQ(1),shapQ1)
      call shape1QQ(Xi(2),nordB(2), ndofQ(2),shapQ2)
!
!  ...shape functions with values along the x-axis
      do j=1,ndofQ(2)
         do i=1,ndofH(1)
            m=m+1
            ShapV(1,m) =  shapH1(i)*shapQ2(j)
            ShapV(2,m) = 0.d0
            DivV(m)    = dshapH1(i)*shapQ2(j)
         enddo
      enddo
!
!  ...shape functions with values along the y-axis
      do j=1,ndofH(2)
         do i=1,ndofQ(1)
            m=m+1
            ShapV(1,m) = 0.d0
            ShapV(2,m) = shapQ1(i)* shapH2(j)
            DivV(m)    = shapQ1(i)*dshapH2(j)
         enddo
      enddo
!
!  ...give total degrees of freedom
      NrdofV = m
!
!
   end subroutine shape2DVBrokenQuad
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape2DQBrokenQuad
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!     purpose:          - routine returns values of 2D quadrilateral
!                         BROKEN L2 shape functions
!
!     arguments:
!
!     in:
!          X            - master quadrilateral coordinates from (0,1)^2
!          NordM        - polynomial order for middle node (H1 sense)
!          Nsize        - relevant sizes of local arrays
!
!     out:
!          NrdofQ       - number of dof
!          ShapQ        - values of the shape functions at the point
!
!-----------------------------------------------------------------------
!
   subroutine shape2DQBrokenQuad(Xi,NordM,Nsize, NrdofQ,ShapQ)
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofQ
      integer :: norder(5)
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapQ(Nsize(2))
!
!  ...Option 1: Simply call the usual shape functions with enrichment
      norder(1:4)=1
      norder(5)=NordM
      call shape2DQQuad(Xi,norder,Nsize, NrdofQ,ShapQ)
!
!  ...Option 2: Write a separate routine for enriched functions
!
!
   end subroutine shape2DQBrokenQuad

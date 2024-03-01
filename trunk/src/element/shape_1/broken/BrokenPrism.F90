! Routines:
!  - shape3DHBrokenPris
!  - shape3DEBrokenPris
!  - shape3DVBrokenPris
!  - shape3DQBrokenPris
!--------------------------------------------------------------------
!
!     routine name      - shape3DHBrokenPris
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!     purpose:          - routine returns values of 3D prism
!                         BROKEN H1 shape functions
!
!     arguments:
!
!     in:
!          X            - master prism coordinates from (0,1)^3
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
   subroutine shape3DHBrokenPris(Xi,NordM,Nsize, NrdofH,ShapH,GradH)
!
      use parameters , only : MODORDER
      use node_types
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofH
      integer :: i12,i3,m
!      integer :: noriE(9),noriF(5),norder(15)
      integer :: nordB(2),ndofH12,ndofH3
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapH(Nsize(2))
      double precision, intent(out) :: GradH(1:3,Nsize(2))
      double precision :: shapH12((Nsize(1)+2)*(Nsize(1)+1)/2), &
                         dshapH12(1:2,(Nsize(1)+2)*(Nsize(1)+1)/2)
      double precision :: shapH3(Nsize(1)+1),dshapH3(Nsize(1)+1)
!
!  ...Option 1: Simply call the usual shape functions with enrichment
!      noriE(1:9)=0
!      noriF(1:5)=0
!      call decod(NordM,MODORDER,2, nordB)
!      norder(1:6)=nordB(1)
!      norder(7:9)=nordB(2)
!      norder(10:11)=nordB(1)
!      norder(12:14)=NordM
!      norder(15)=NordM
!      call shape3DHPris(Xi,norder,noriE,noriF,Nsize, &
!                                                   NrdofH,ShapH,GradH)
!
!
!  ...Option 2: Write more efficient routine for enriched functions
!
!  ...initiate counter for shape functions
      m = 0
!
!  ...Shape functions are tensor product of H1 triangles and segments
      call decod(NordM,MODORDER,2, nordB)
      call shape2HH(MDLT,Xi(1:2),nordB(1), ndofH12,shapH12,dshapH12)
      call shape1HH(Xi(3),nordB(2), ndofH3,shapH3,dshapH3)
!
      do i3=1,ndofH3
         do i12=1,ndofH12
            m=m+1
            ShapH(m)   = shapH12(i12)*   shapH3(i3)
            GradH(1,m) = dshapH12(1,i12)*shapH3(i3)
            GradH(2,m) = dshapH12(2,i12)*shapH3(i3)
            GradH(3,m) = shapH12(i12)*  dshapH3(i3)
         enddo
      enddo
!
!  ...total degrees of freedom
      NrdofH = m
!
   end subroutine shape3DHBrokenPris
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape3DEBrokenPris
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!     purpose:          - routine returns values of 3D prism
!                         BROKEN H(curl) shape functions
!
!     arguments:
!
!     in:
!          X            - master prism coordinates from (0,1)^3
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
   subroutine shape3DEBrokenPris(Xi,NordM,Nsize, NrdofE,ShapE,CurlE)
!
      use parameters , only : MODORDER
      use node_types
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofE
      integer :: i12,i3,m
!      integer :: noriE(9),noriF(5),norder(15)
      integer :: nordB(2),ndofE12,ndofH12,ndofH3,ndofQ3
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapE(1:3,Nsize(2))
      double precision, intent(out) :: CurlE(1:3,Nsize(2))
!  ...Family 1 shape arrays
      double precision :: shapE12(1:2,Nsize(1)*(Nsize(1)+2)), &
                          curlE12(Nsize(1)*(Nsize(1)+2))
      double precision :: shapH3(Nsize(1)+1),dshapH3(Nsize(1)+1)
!  ...Family 2 shape arrays
      double precision :: shapH12((Nsize(1)+2)*(Nsize(1)+1)/2), &
                          dshapH12(1:2,(Nsize(1)+2)*(Nsize(1)+1)/2)
      double precision :: shapQ3(Nsize(1))
!
!  ...Option 1: Simply call the usual shape functions with enrichment
!      noriE(1:9)=0
!      noriF(1:5)=0
!      call decod(NordM,MODORDER,2, nordB)
!      norder(1:6)=nordB(1)
!      norder(7:9)=nordB(2)
!      norder(10:11)=nordB(1)
!      norder(12:14)=NordM
!      norder(15)=NordM
!      call shape3DEPris(Xi,norder,noriE,noriF,Nsize, &
!                                                   NrdofE,ShapE,CurlE)
!
!
!  ...Option 2: Write more efficient routine for enriched functions
!
!  ...initiate counter for shape functions
      m = 0
!
!  ...Decode order
      call decod(NordM,MODORDER,2, nordB)
!
!  ...Shape functions split into 2 families
!     Family 1: E12 x H3
      call shape2EE(MDLT,Xi(1:2),nordB(1), ndofE12,shapE12,curlE12)
      call shape1HH(Xi(3),nordB(2), ndofH3,shapH3,dshapH3)
!
!     Family 2: H12 x Q3
      call shape2HH(MDLT,Xi(1:2),nordB(1), ndofH12,shapH12,dshapH12)
      call shape1QQ(Xi(3),nordB(2), ndofQ3,shapQ3)
!
!  ...Evaluate Family 1
      do i3=1,ndofH3
         do i12=1,ndofE12
            m=m+1
!        ...value
            ShapE(1,m) = shapE12(1,i12)*shapH3(i3)
            ShapE(2,m) = shapE12(2,i12)*shapH3(i3)
            ShapE(3,m) = 0.d0
!        ...curl
            CurlE(1,m) = -shapE12(2,i12)*dshapH3(i3)
            CurlE(2,m) =  shapE12(1,i12)*dshapH3(i3)
            CurlE(3,m) =  curlE12(i12)*   shapH3(i3)
         enddo
      enddo
!
!  ...Evaluate Family 2
      do i3=1,ndofQ3
         do i12=1,ndofH12
            m=m+1
!        ...value
            ShapE(1,m) = 0.d0
            ShapE(2,m) = 0.d0
            ShapE(3,m) = shapH12(i12)*shapQ3(i3)
!        ...curl
            CurlE(1,m) =  dshapH12(2,i12)*shapQ3(i3)
            CurlE(2,m) = -dshapH12(1,i12)*shapQ3(i3)
            CurlE(3,m) =  0.d0
         enddo
      enddo
!
!  ...total degrees of freedom
      NrdofE = m
!
   end subroutine shape3DEBrokenPris
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape3DVBrokenPris
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!     purpose:          - routine returns values of 3D prism
!                         BROKEN H(div) shape functions
!
!     arguments:
!
!     in:
!          X            - master prism coordinates from (0,1)^3
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
   subroutine shape3DVBrokenPris(Xi,NordM,Nsize, NrdofV,ShapV,DivV)
!
      use parameters , only : MODORDER
      use node_types
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofV
      integer :: i12,i3,m
!      integer :: noriF(5),norder(15)
      integer :: nordB(2),ndofV12,ndofQ12,ndofH3,ndofQ3
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapV(1:3,Nsize(2))
      double precision, intent(out) :: DivV(Nsize(2))
!  ...Family 1 shape arrays
      double precision :: shapV12(1:2,Nsize(1)*(Nsize(1)+2)), &
                          divV12(Nsize(1)*(Nsize(1)+2))
      double precision :: shapQ3(Nsize(1))
!  ...Family 2 shape arrays
      double precision :: shapQ12(Nsize(1)*(Nsize(1)+1)/2)
      double precision :: shapH3(Nsize(1)+1),dshapH3(Nsize(1)+1)
!
!  ...Option 1: Simply call the usual shape functions with enrichment
!      noriF(1:5)=0
!      call decod(NordM,MODORDER,2, nordB)
!      norder(1:9)=1
!      norder(10:11)=nordB(1)
!      norder(12:14)=NordM
!      norder(15)=NordM
!      call shape3DVPris(Xi,norder,noriF,Nsize, NrdofV,ShapV,DivV)
!
!  ...Option 2: Write more efficient routine for enriched functions
!
!  ...initiate counter for shape functions
      m = 0
!
!  ...Decode order
      call decod(NordM,MODORDER,2, nordB)
!
!  ...Shape functions split into 2 families
!     Family 1: V12 x Q3
      call shape2VV(MDLT,Xi(1:2),nordB(1), ndofV12,shapV12,divV12)
      call shape1QQ(Xi(3),nordB(2), ndofQ3,shapQ3)
!
!     Family 2: Q12 x H3
      call shape2QQ(MDLT,Xi(1:2),nordB(1), ndofQ12,shapQ12)
      call shape1HH(Xi(3),nordB(2), ndofH3,shapH3,dshapH3)
!
!  ...Evaluate Family 1
      do i3=1,ndofQ3
         do i12=1,ndofV12
            m=m+1
!        ...value
            ShapV(1,m) = shapV12(1,i12)*shapQ3(i3)
            ShapV(2,m) = shapV12(2,i12)*shapQ3(i3)
            ShapV(3,m) = 0.d0
!        ...divergence
            DivV(m) = divV12(i12)*shapQ3(i3)
         enddo
      enddo
!
!  ...Evaluate Family 2
      do i3=1,ndofH3
         do i12=1,ndofQ12
            m=m+1
!        ...value
            ShapV(1,m) = 0.d0
            ShapV(2,m) = 0.d0
            ShapV(3,m) = shapQ12(i12)*shapH3(i3)
!        ...divergence
            DivV(m) = shapQ12(i12)*dshapH3(i3)
         enddo
      enddo
!
!  ...total degrees of freedom
      NrdofV = m
!
   end subroutine shape3DVBrokenPris
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape3DQBrokenPris
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!     purpose:          - routine returns values of 3D prism
!                         BROKEN L2 shape functions
!
!     arguments:
!
!     in:
!          X            - master prism coordinates from (0,1)^3
!          NordM        - polynomial order for middle node (H1 sense)
!          Nsize        - relevant sizes of local arrays
!
!     out:
!          NrdofQ       - number of dof
!          ShapQ        - values of the shape functions at the point
!
!-----------------------------------------------------------------------
!
   subroutine shape3DQBrokenPris(Xi,NordM,Nsize, NrdofQ,ShapQ)
!
      use parameters , only : MODORDER
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofQ
      integer :: norder(15),nordF
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapQ(Nsize(2))
!
!  ...Option 1: Simply call the usual shape functions with enrichment
      norder(1:11)=1
      call encod((/1,1/),MODORDER,2, nordF)
      norder(12:14)=nordF
      norder(15)=NordM
      call shape3DQPris(Xi,norder,Nsize, NrdofQ,ShapQ)
!
!  ...Option 2: Write a separate routine for enriched functions
!
!
   end subroutine shape3DQBrokenPris

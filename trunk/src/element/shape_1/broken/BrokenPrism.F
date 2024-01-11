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
      use node_types
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofH
      integer :: i12,i3,m
c      integer :: noriE(9),noriF(5),norder(15)
      integer :: nordB(2),ndofH12,ndofH3
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapH(Nsize(2))
      double precision, intent(out) :: GradH(1:3,Nsize(2))
      double precision :: shapH12((Nsize(1)+2)*(Nsize(1)+1)/2),
     .                   dshapH12(1:2,(Nsize(1)+2)*(Nsize(1)+1)/2)
      double precision :: shapH3(Nsize(1)+1),dshapH3(Nsize(1)+1)
c
c  ...Option 1: Simply call the usual shape functions with enrichment
c      noriE(1:9)=0
c      noriF(1:5)=0
c      call decod(NordM,MODORDER,2, nordB)
c      norder(1:6)=nordB(1)
c      norder(7:9)=nordB(2)
c      norder(10:11)=nordB(1)
c      norder(12:14)=NordM
c      norder(15)=NordM
c      call shape3DHPris(Xi,norder,noriE,noriF,Nsize,
c     .                                             NrdofH,ShapH,GradH)
c
c
c  ...Option 2: Write more efficient routine for enriched functions
c
c  ...initiate counter for shape functions
      m = 0
c
c  ...Shape functions are tensor product of H1 triangles and segments
      call decod(NordM,MODORDER,2, nordB)
      call shape2HH(MDLT,Xi(1:2),nordB(1), ndofH12,shapH12,dshapH12)
      call shape1HH(Xi(3),nordB(2), ndofH3,shapH3,dshapH3)
c
      do i3=1,ndofH3
         do i12=1,ndofH12
            m=m+1
            ShapH(m)   = shapH12(i12)*   shapH3(i3)
            GradH(1,m) = dshapH12(1,i12)*shapH3(i3)
            GradH(2,m) = dshapH12(2,i12)*shapH3(i3)
            GradH(3,m) = shapH12(i12)*  dshapH3(i3)
         enddo
      enddo
c
c  ...total degrees of freedom
      NrdofH = m
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
      use node_types
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofE
      integer :: i12,i3,m
c      integer :: noriE(9),noriF(5),norder(15)
      integer :: nordB(2),ndofE12,ndofH12,ndofH3,ndofQ3
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapE(1:3,Nsize(2))
      double precision, intent(out) :: CurlE(1:3,Nsize(2))
c  ...Family 1 shape arrays
      double precision :: shapE12(1:2,Nsize(1)*(Nsize(1)+2)),
     .                    curlE12(Nsize(1)*(Nsize(1)+2))
      double precision :: shapH3(Nsize(1)+1),dshapH3(Nsize(1)+1)
c  ...Family 2 shape arrays
      double precision :: shapH12((Nsize(1)+2)*(Nsize(1)+1)/2),
     .                    dshapH12(1:2,(Nsize(1)+2)*(Nsize(1)+1)/2)
      double precision :: shapQ3(Nsize(1))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
c      noriE(1:9)=0
c      noriF(1:5)=0
c      call decod(NordM,MODORDER,2, nordB)
c      norder(1:6)=nordB(1)
c      norder(7:9)=nordB(2)
c      norder(10:11)=nordB(1)
c      norder(12:14)=NordM
c      norder(15)=NordM
c      call shape3DEPris(Xi,norder,noriE,noriF,Nsize,
c     .                                             NrdofE,ShapE,CurlE)
c
c
c  ...Option 2: Write more efficient routine for enriched functions
c
c  ...initiate counter for shape functions
      m = 0
c
c  ...Decode order
      call decod(NordM,MODORDER,2, nordB)
c
c  ...Shape functions split into 2 families
c     Family 1: E12 x H3
      call shape2EE(MDLT,Xi(1:2),nordB(1), ndofE12,shapE12,curlE12)
      call shape1HH(Xi(3),nordB(2), ndofH3,shapH3,dshapH3)
c
c     Family 2: H12 x Q3
      call shape2HH(MDLT,Xi(1:2),nordB(1), ndofH12,shapH12,dshapH12)
      call shape1QQ(Xi(3),nordB(2), ndofQ3,shapQ3)
c
c  ...Evaluate Family 1
      do i3=1,ndofH3
         do i12=1,ndofE12
            m=m+1
c        ...value
            ShapE(1,m) = shapE12(1,i12)*shapH3(i3)
            ShapE(2,m) = shapE12(2,i12)*shapH3(i3)
            ShapE(3,m) = 0.d0
c        ...curl
            CurlE(1,m) = -shapE12(2,i12)*dshapH3(i3)
            CurlE(2,m) =  shapE12(1,i12)*dshapH3(i3)
            CurlE(3,m) =  curlE12(i12)*   shapH3(i3)
         enddo
      enddo
c
c  ...Evaluate Family 2
      do i3=1,ndofQ3
         do i12=1,ndofH12
            m=m+1
c        ...value
            ShapE(1,m) = 0.d0
            ShapE(2,m) = 0.d0
            ShapE(3,m) = shapH12(i12)*shapQ3(i3)
c        ...curl
            CurlE(1,m) =  dshapH12(2,i12)*shapQ3(i3)
            CurlE(2,m) = -dshapH12(1,i12)*shapQ3(i3)
            CurlE(3,m) =  0.d0
         enddo
      enddo
c
c  ...total degrees of freedom
      NrdofE = m
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
      use node_types
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofV
      integer :: i12,i3,m
c      integer :: noriF(5),norder(15)
      integer :: nordB(2),ndofV12,ndofQ12,ndofH3,ndofQ3
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapV(1:3,Nsize(2))
      double precision, intent(out) :: DivV(Nsize(2))
c  ...Family 1 shape arrays
      double precision :: shapV12(1:2,Nsize(1)*(Nsize(1)+2)),
     .                    divV12(Nsize(1)*(Nsize(1)+2))
      double precision :: shapQ3(Nsize(1))
c  ...Family 2 shape arrays
      double precision :: shapQ12(Nsize(1)*(Nsize(1)+1)/2)
      double precision :: shapH3(Nsize(1)+1),dshapH3(Nsize(1)+1)
c
c  ...Option 1: Simply call the usual shape functions with enrichment
c      noriF(1:5)=0
c      call decod(NordM,MODORDER,2, nordB)
c      norder(1:9)=1
c      norder(10:11)=nordB(1)
c      norder(12:14)=NordM
c      norder(15)=NordM
c      call shape3DVPris(Xi,norder,noriF,Nsize, NrdofV,ShapV,DivV)
c
c  ...Option 2: Write more efficient routine for enriched functions
c
c  ...initiate counter for shape functions
      m = 0
c
c  ...Decode order
      call decod(NordM,MODORDER,2, nordB)
c
c  ...Shape functions split into 2 families
c     Family 1: V12 x Q3
      call shape2VV(MDLT,Xi(1:2),nordB(1), ndofV12,shapV12,divV12)
      call shape1QQ(Xi(3),nordB(2), ndofQ3,shapQ3)
c
c     Family 2: Q12 x H3
      call shape2QQ(MDLT,Xi(1:2),nordB(1), ndofQ12,shapQ12)
      call shape1HH(Xi(3),nordB(2), ndofH3,shapH3,dshapH3)
c
c  ...Evaluate Family 1
      do i3=1,ndofQ3
         do i12=1,ndofV12
            m=m+1
c        ...value
            ShapV(1,m) = shapV12(1,i12)*shapQ3(i3)
            ShapV(2,m) = shapV12(2,i12)*shapQ3(i3)
            ShapV(3,m) = 0.d0
c        ...divergence
            DivV(m) = divV12(i12)*shapQ3(i3)
         enddo
      enddo
c
c  ...Evaluate Family 2
      do i3=1,ndofH3
         do i12=1,ndofQ12
            m=m+1
c        ...value
            ShapV(1,m) = 0.d0
            ShapV(2,m) = 0.d0
            ShapV(3,m) = shapQ12(i12)*shapH3(i3)
c        ...divergence
            DivV(m) = shapQ12(i12)*dshapH3(i3)
         enddo
      enddo
c
c  ...total degrees of freedom
      NrdofV = m
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

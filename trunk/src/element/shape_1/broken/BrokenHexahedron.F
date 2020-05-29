c Routines:
c  - shape3DHBrokenHexa
c  - shape3DEBrokenHexa
c  - shape3DVBrokenHexa
c  - shape3DQBrokenHexa
c--------------------------------------------------------------------
c
c     routine name      - shape3DHBrokenHexa
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 3D hexahedron
c                         BROKEN H1 shape functions
c
c     arguments:
c
c     in:
c          X            - master hexahedron coordinates from (0,1)^3
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
      subroutine shape3DHBrokenHexa(Xi,NordM,Nsize, NrdofH,ShapH,GradH)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofH
      !integer :: noriE(12),noriF(6),norder(19),nordF(3)
      integer :: i,j,k,m
      integer :: nordB(3),ndofH(3),iprint
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapH(Nsize(2))
      double precision, intent(out) :: GradH(1:3,Nsize(2))
      double precision :: shapH1(Nsize(1)+1),dshapH1(Nsize(1)+1)
      double precision :: shapH2(Nsize(1)+1),dshapH2(Nsize(1)+1)
      double precision :: shapH3(Nsize(1)+1),dshapH3(Nsize(1)+1)
c
c  ...Option 1: Simply call the usual shape functions with enrichment
c      noriE(1:12)=0
c      noriF(1:6)=0
c      call decod(NordM,MODORDER,2, (/nordF(1),nordB(3)/))
c      call decod(nordF(1),MODORDER,2, nordB(1:2))
c      call encod((/nordB(1),nordB(3)/),MODORDER,2, nordF(2))
c      call encod(nordB(2:3),MODORDER,2, nordF(3))
c      norder(1:4)=(/nordB(1),nordB(2),nordB(1),nordB(2)/)
c      norder(5:8)=(/nordB(1),nordB(2),nordB(1),nordB(2)/)
c      norder(9:12)=nordB(3)
c      norder(13:14)=nordF(1)
c      norder(15:18)=(/nordF(2),nordF(3),nordF(2),nordF(3)/)
c      norder(19)=NordM
c      call shape3DHHexa(Xi,norder,noriE,noriF,Nsize,
c     .                                             NrdofH,ShapH,GradH)
c
c
c  ...Option 2: Write more efficient routine for enriched functions
c
c  ...debugging flag
      iprint=0
c  ...initiate counter for shape functions
      m=0
c
c  ...shape functions are tensor products of 1D shape functions
      call decod(NordM,MODORDER,3, nordB)
      call shape1HH(Xi(1),nordB(1), ndofH(1),shapH1,dshapH1)
      call shape1HH(Xi(2),nordB(2), ndofH(2),shapH2,dshapH2)
      call shape1HH(Xi(3),nordB(3), ndofH(3),shapH3,dshapH3)
c
      do k=1,ndofH(3)
        do j=1,ndofH(2)
          do i=1,ndofH(1)
            m=m+1
            ShapH(m)   =  shapH1(i)* shapH2(j)* shapH3(k)
            GradH(1,m) = dshapH1(i)* shapH2(j)* shapH3(k)
            GradH(2,m) =  shapH1(i)*dshapH2(j)* shapH3(k)
            GradH(3,m) =  shapH1(i)* shapH2(j)*dshapH3(k)
          enddo
        enddo
      enddo
c
c  ...give total degrees of freedom
      NrdofH = m
c
c  ...print this when debugging
      if (iprint.ge.1) then
        write(*,7001) Xi(1:3)
 7001   format('shape3DHBrokenHexa: Xi = ',3f8.3)
        do m=1,NrdofH
          write(*,7002) m,ShapH(m),GradH(1:3,m)
 7002     format('k = ',i3,' ShapH, GradH = ',e12.5,3x,3e12.5)
        enddo
        call pause
      endif
c
c
      end subroutine shape3DHBrokenHexa
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape3DEBrokenHexa
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 3D hexahedron
c                         BROKEN H(curl) shape functions
c
c     arguments:
c
c     in:
c          X            - master hexahedron coordinates from (0,1)^3
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
      subroutine shape3DEBrokenHexa(Xi,NordM,Nsize, NrdofE,ShapE,CurlE)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofE
      !integer :: noriE(12),noriF(6),norder(19),nordF(3)
      integer :: i,j,k,m
      integer :: nordB(3),ndofH(3),ndofQ(3),iprint
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapE(1:3,Nsize(2))
      double precision, intent(out) :: CurlE(1:3,Nsize(2))
      double precision :: shapH1(Nsize(1)+1),dshapH1(Nsize(1)+1)
      double precision :: shapH2(Nsize(1)+1),dshapH2(Nsize(1)+1)
      double precision :: shapH3(Nsize(1)+1),dshapH3(Nsize(1)+1)
      double precision :: shapQ1(Nsize(1))
      double precision :: shapQ2(Nsize(1))
      double precision :: shapQ3(Nsize(1))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
c      noriE(1:12)=0
c      noriF(1:6)=0
c      call decod(NordM,MODORDER,2, (/nordF(1),nordB(3)/))
c      call decod(nordF(1),MODORDER,2, nordB(1:2))
c      call encod((/nordB(1),nordB(3)/),MODORDER,2, nordF(2))
c      call encod(nordB(2:3),MODORDER,2, nordF(3))
c      norder(1:4)=(/nordB(1),nordB(2),nordB(1),nordB(2)/)
c      norder(5:8)=(/nordB(1),nordB(2),nordB(1),nordB(2)/)
c      norder(9:12)=nordB(3)
c      norder(13:14)=nordF(1)
c      norder(15:18)=(/nordF(2),nordF(3),nordF(2),nordF(3)/)
c      norder(19)=NordM
c      call shape3DEHexa(Xi,norder,noriE,noriF,Nsize,
c     .                                             NrdofE,ShapE,CurlE)
c
c
c  ...Option 2: Write more efficient routine for enriched functions
c
c  ...debugging flag
      iprint=0
c  ...initiate counter for shape functions
      m=0
c
c  ...shape functions are tensor products of 1D shape functions
      call decod(NordM,MODORDER,3, nordB)
      call shape1HH(Xi(1),nordB(1), ndofH(1),shapH1,dshapH1)
      call shape1HH(Xi(2),nordB(2), ndofH(2),shapH2,dshapH2)
      call shape1HH(Xi(3),nordB(3), ndofH(3),shapH3,dshapH3)
      call shape1QQ(Xi(1),nordB(1), ndofQ(1),shapQ1)
      call shape1QQ(Xi(2),nordB(2), ndofQ(2),shapQ2)
      call shape1QQ(Xi(3),nordB(3), ndofQ(3),shapQ3)
c
c  ...shape functions with values along the x-axis
      do k=1,ndofH(3)
        do j=1,ndofH(2)
          do i=1,ndofQ(1)
            m=m+1
            ShapE(1,m) =  shapQ1(i)*shapH2(j)* shapH3(k)
            ShapE(2,m) =  0.d0
            ShapE(3,m) =  0.d0
            CurlE(1,m) =  0.d0
            CurlE(2,m) =  shapQ1(i)* shapH2(j)*dshapH3(k)
            CurlE(3,m) = -shapQ1(i)*dshapH2(j)* shapH3(k)
          enddo
        enddo
      enddo
c
c  ...shape functions with values along the y-axis
      do k=1,ndofH(3)
        do j=1,ndofQ(2)
          do i=1,ndofH(1)
            m=m+1
            ShapE(1,m) =  0.d0
            ShapE(2,m) =  shapH1(i)*shapQ2(j)* shapH3(k)
            ShapE(3,m) =  0.d0
            CurlE(1,m) = -shapH1(i)*shapQ2(j)*dshapH3(k)
            CurlE(2,m) =  0.d0
            CurlE(3,m) = dshapH1(i)*shapQ2(j)* shapH3(k)
          enddo
        enddo
      enddo
c
c  ...shape functions with values along the z-axis
      do k=1,ndofQ(3)
        do j=1,ndofH(2)
          do i=1,ndofH(1)
            m=m+1
            ShapE(1,m) = 0.d0
            ShapE(2,m) = 0.d0
            ShapE(3,m) =  shapH1(i)* shapH2(j)* shapQ3(k)
            CurlE(1,m) =  shapH1(i)*dshapH2(j)* shapQ3(k)
            CurlE(2,m) =-dshapH1(i)* shapH2(j)* shapQ3(k)
            CurlE(3,m) = 0.d0
          enddo
        enddo
      enddo
c
c  ...give total degrees of freedom
      NrdofE = m
c
c  ...print this when debugging
      if (iprint.ge.1) then
        write(*,7001) Xi(1:3)
 7001   format('shape3DEBrokenHexa: Xi = ',3f8.3)
        do m=1,NrdofE
          write(*,7002) k,ShapE(1:3,m),CurlE(1:3,m)
 7002     format('k = ',i3,' ShapE, CurlE = ',3e12.5,3x,3e12.5)
        enddo
        call pause
      endif
c
c
      end subroutine shape3DEBrokenHexa
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape3DVBrokenHexa
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 3D hexahedron
c                         BROKEN H(div) shape functions
c
c     arguments:
c
c     in:
c          X            - master hexahedron coordinates from (0,1)^3
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
      subroutine shape3DVBrokenHexa(Xi,NordM,Nsize, NrdofV,ShapV,DivV)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofV
      !integer :: noriF(6),norder(19),nordF(3)
      integer :: i,j,k,m
      integer :: nordB(3),ndofH(3),ndofQ(3),iprint
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapV(1:3,Nsize(2))
      double precision, intent(out) :: DivV(Nsize(2))
      double precision :: shapH1(Nsize(1)+1),dshapH1(Nsize(1)+1)
      double precision :: shapH2(Nsize(1)+1),dshapH2(Nsize(1)+1)
      double precision :: shapH3(Nsize(1)+1),dshapH3(Nsize(1)+1)
      double precision :: shapQ1(Nsize(1))
      double precision :: shapQ2(Nsize(1))
      double precision :: shapQ3(Nsize(1))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
c      noriF(1:6)=0
c      call decod(NordM,MODORDER,2, (/nordF(1),nordB(3)/))
c      call decod(nordF(1),MODORDER,2, nordB(1:2))
c      call encod((/nordB(1),nordB(3)/),MODORDER,2, nordF(2))
c      call encod(nordB(2:3),MODORDER,2, nordF(3))
c      norder(1:12)=1
c      norder(13:14)=nordF(1)
c      norder(15:18)=(/nordF(2),nordF(3),nordF(2),nordF(3)/)
c      norder(19)=NordM
c      call shape3DVHexa(Xi,norder,noriF,Nsize, NrdofV,ShapV,DivV)
c
c  ...Option 2: Write more efficient routine for enriched functions
c
c  ...debugging flag
      iprint=0
c  ...initiate counter for shape functions
      m=0
c
c  ...shape functions are tensor products of 1D shape functions
      call decod(NordM,MODORDER,3, nordB)
      call shape1HH(Xi(1),nordB(1), ndofH(1),shapH1,dshapH1)
      call shape1HH(Xi(2),nordB(2), ndofH(2),shapH2,dshapH2)
      call shape1HH(Xi(3),nordB(3), ndofH(3),shapH3,dshapH3)
      call shape1QQ(Xi(1),nordB(1), ndofQ(1),shapQ1)
      call shape1QQ(Xi(2),nordB(2), ndofQ(2),shapQ2)
      call shape1QQ(Xi(3),nordB(3), ndofQ(3),shapQ3)
c
c  ...shape functions with values along the x-axis
      do k=1,ndofQ(3)
        do j=1,ndofQ(2)
          do i=1,ndofH(1)
            m=m+1
            ShapV(1,m) =  shapH1(i)* shapQ2(j)* shapQ3(k)
            ShapV(2,m) = 0.d0
            ShapV(3,m) = 0.d0
            DivV(m)    = dshapH1(i)* shapQ2(j)* shapQ3(k)
          enddo
        enddo
      enddo
c
c  ...shape functions with values along the y-axis
      do k=1,ndofQ(3)
        do j=1,ndofH(2)
          do i=1,ndofQ(1)
            m=m+1
            ShapV(1,m) = 0.d0
            ShapV(2,m) =  shapQ1(i)* shapH2(j)* shapQ3(k)
            ShapV(3,m) = 0.d0
            DivV(m)    =  shapQ1(i)*dshapH2(j)* shapQ3(k)
          enddo
        enddo
      enddo
c
c  ...shape functions with values along the z-axis
      do k=1,ndofH(3)
        do j=1,ndofQ(2)
          do i=1,ndofQ(1)
            m=m+1
            ShapV(1,m) = 0.d0
            ShapV(2,m) = 0.d0
            ShapV(3,m) =  shapQ1(i)* shapQ2(j)* shapH3(k)
            DivV(m)    =  shapQ1(i)* shapQ2(j)*dshapH3(k)
          enddo
        enddo
      enddo
c
c  ...give total degrees of freedom
      NrdofV = m
c
c  ...print this when debugging
      if (iprint.ge.1) then
        write(*,7001) Xi(1:3)
 7001   format('shape3DVBrokenHexa: Xi = ',3f8.3)
        do m=1,NrdofV
          write(*,7002) m,ShapV(1:3,m),DivV(m)
 7002     format('k = ',i3,' ShapV, DivV = ',3e12.5,3x,e12.5)
        enddo
        call pause
      endif
c
c
      end subroutine shape3DVBrokenHexa
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape3DQBrokenHexa
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns values of 3D hexahedron
c                         BROKEN L2 shape functions
c
c     arguments:
c
c     in:
c          X            - master hexahedron coordinates from (0,1)^3
c          NordM        - polynomial order for middle node (H1 sense)
c          Nsize        - relevant sizes of local arrays
c
c     out:
c          NrdofQ       - number of dof
c          ShapQ        - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape3DQBrokenHexa(Xi,NordM,Nsize, NrdofQ,ShapQ)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofQ
      integer :: norder(19),nordF
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapQ(Nsize(2))
c
c  ...Option 1: Simply call the usual shape functions with enrichment
      norder(1:12)=1
      call encod((/1,1/),MODORDER,2, nordF)
      norder(13:18)=nordF
      norder(19)=NordM
      call shape3DQHexa(Xi,norder,Nsize, NrdofQ,ShapQ)
c
c  ...Option 2: Write a separate routine for enriched functions
c
c
      end subroutine shape3DQBrokenHexa


! Routines:
!  - shape3DHBrokenHexa
!  - shape3DEBrokenHexa
!  - shape3DVBrokenHexa
!  - shape3DQBrokenHexa
!--------------------------------------------------------------------
!
!     routine name      - shape3DHBrokenHexa
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!     purpose:          - routine returns values of 3D hexahedron
!                         BROKEN H1 shape functions
!
!     arguments:
!
!     in:
!          X            - master hexahedron coordinates from (0,1)^3
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
   subroutine shape3DHBrokenHexa(Xi,NordM,Nsize, NrdofH,ShapH,GradH)
!
      use parameters , only : MODORDER
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofH
!      integer :: noriE(12),noriF(6),norder(19),nordF(3)
      integer :: i,j,k,m
      integer :: nordB(3),ndofH(3)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapH(Nsize(2))
      double precision, intent(out) :: GradH(1:3,Nsize(2))
      double precision :: shapH1(Nsize(1)+1),dshapH1(Nsize(1)+1)
      double precision :: shapH2(Nsize(1)+1),dshapH2(Nsize(1)+1)
      double precision :: shapH3(Nsize(1)+1),dshapH3(Nsize(1)+1)
!
#if HP3D_DEBUG
!  ...debugging flag
      integer :: iprint
      iprint=0
#endif
!
!  ...Option 1: Simply call the usual shape functions with enrichment
!      noriE(1:12)=0
!      noriF(1:6)=0
!      call decod(NordM,MODORDER,2, (/nordF(1),nordB(3)/))
!      call decod(nordF(1),MODORDER,2, nordB(1:2))
!      call encod((/nordB(1),nordB(3)/),MODORDER,2, nordF(2))
!      call encod(nordB(2:3),MODORDER,2, nordF(3))
!      norder(1:4)=(/nordB(1),nordB(2),nordB(1),nordB(2)/)
!      norder(5:8)=(/nordB(1),nordB(2),nordB(1),nordB(2)/)
!      norder(9:12)=nordB(3)
!      norder(13:14)=nordF(1)
!      norder(15:18)=(/nordF(2),nordF(3),nordF(2),nordF(3)/)
!      norder(19)=NordM
!      call shape3DHHexa(Xi,norder,noriE,noriF,Nsize, &
!                                                   NrdofH,ShapH,GradH)
!
!
!  ...Option 2: Write more efficient routine for enriched functions
!
!  ...initiate counter for shape functions
      m=0
!
!  ...shape functions are tensor products of 1D shape functions
      call decod(NordM,MODORDER,3, nordB)
      call shape1HH(Xi(1),nordB(1), ndofH(1),shapH1,dshapH1)
      call shape1HH(Xi(2),nordB(2), ndofH(2),shapH2,dshapH2)
      call shape1HH(Xi(3),nordB(3), ndofH(3),shapH3,dshapH3)
!
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
!
!  ...give total degrees of freedom
      NrdofH = m
!
#if HP3D_DEBUG
!  ...print this when debugging
      if (iprint.ge.1) then
        write(*,7001) Xi(1:3)
 7001   format('shape3DHBrokenHexa: Xi = ',3f8.3)
        do m=1,NrdofH
          write(*,7002) m,ShapH(m),GradH(1:3,m)
 7002     format('k = ',i3,' ShapH, GradH = ',e12.5,3x,3e12.5)
        enddo
        call pause
      endif
#endif
!
!
   end subroutine shape3DHBrokenHexa
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape3DEBrokenHexa
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!     purpose:          - routine returns values of 3D hexahedron
!                         BROKEN H(curl) shape functions
!
!     arguments:
!
!     in:
!          X            - master hexahedron coordinates from (0,1)^3
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
   subroutine shape3DEBrokenHexa(Xi,NordM,Nsize, NrdofE,ShapE,CurlE)
!
      use parameters , only : MODORDER
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofE
!      integer :: noriE(12),noriF(6),norder(19),nordF(3)
      integer :: i,j,k,m
      integer :: nordB(3),ndofH(3),ndofQ(3)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapE(1:3,Nsize(2))
      double precision, intent(out) :: CurlE(1:3,Nsize(2))
      double precision :: shapH1(Nsize(1)+1),dshapH1(Nsize(1)+1)
      double precision :: shapH2(Nsize(1)+1),dshapH2(Nsize(1)+1)
      double precision :: shapH3(Nsize(1)+1),dshapH3(Nsize(1)+1)
      double precision :: shapQ1(Nsize(1))
      double precision :: shapQ2(Nsize(1))
      double precision :: shapQ3(Nsize(1))
!
#if HP3D_DEBUG
!  ...debugging flag
      integer :: iprint
      iprint=0
#endif
!
!  ...Option 1: Simply call the usual shape functions with enrichment
!      noriE(1:12)=0
!      noriF(1:6)=0
!      call decod(NordM,MODORDER,2, (/nordF(1),nordB(3)/))
!      call decod(nordF(1),MODORDER,2, nordB(1:2))
!      call encod((/nordB(1),nordB(3)/),MODORDER,2, nordF(2))
!      call encod(nordB(2:3),MODORDER,2, nordF(3))
!      norder(1:4)=(/nordB(1),nordB(2),nordB(1),nordB(2)/)
!      norder(5:8)=(/nordB(1),nordB(2),nordB(1),nordB(2)/)
!      norder(9:12)=nordB(3)
!      norder(13:14)=nordF(1)
!      norder(15:18)=(/nordF(2),nordF(3),nordF(2),nordF(3)/)
!      norder(19)=NordM
!      call shape3DEHexa(Xi,norder,noriE,noriF,Nsize, &
!                                                   NrdofE,ShapE,CurlE)
!
!
!  ...Option 2: Write more efficient routine for enriched functions
!
!  ...initiate counter for shape functions
      m=0
!
!  ...shape functions are tensor products of 1D shape functions
      call decod(NordM,MODORDER,3, nordB)
      call shape1HH(Xi(1),nordB(1), ndofH(1),shapH1,dshapH1)
      call shape1HH(Xi(2),nordB(2), ndofH(2),shapH2,dshapH2)
      call shape1HH(Xi(3),nordB(3), ndofH(3),shapH3,dshapH3)
      call shape1QQ(Xi(1),nordB(1), ndofQ(1),shapQ1)
      call shape1QQ(Xi(2),nordB(2), ndofQ(2),shapQ2)
      call shape1QQ(Xi(3),nordB(3), ndofQ(3),shapQ3)
!
!  ...shape functions with values along the x-axis
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
!
!  ...shape functions with values along the y-axis
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
!
!  ...shape functions with values along the z-axis
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
!
!  ...give total degrees of freedom
      NrdofE = m
!
#if HP3D_DEBUG
!  ...print this when debugging
      if (iprint.ge.1) then
        write(*,7001) Xi(1:3)
 7001   format('shape3DEBrokenHexa: Xi = ',3f8.3)
        do m=1,NrdofE
          write(*,7002) k,ShapE(1:3,m),CurlE(1:3,m)
 7002     format('k = ',i3,' ShapE, CurlE = ',3e12.5,3x,3e12.5)
        enddo
        call pause
      endif
#endif
!
!
   end subroutine shape3DEBrokenHexa
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape3DVBrokenHexa
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!     purpose:          - routine returns values of 3D hexahedron
!                         BROKEN H(div) shape functions
!
!     arguments:
!
!     in:
!          X            - master hexahedron coordinates from (0,1)^3
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
   subroutine shape3DVBrokenHexa(Xi,NordM,Nsize, NrdofV,ShapV,DivV)
!
      use parameters , only : MODORDER
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofV
!      integer :: noriF(6),norder(19),nordF(3)
      integer :: i,j,k,m
      integer :: nordB(3),ndofH(3),ndofQ(3)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapV(1:3,Nsize(2))
      double precision, intent(out) :: DivV(Nsize(2))
      double precision :: shapH1(Nsize(1)+1),dshapH1(Nsize(1)+1)
      double precision :: shapH2(Nsize(1)+1),dshapH2(Nsize(1)+1)
      double precision :: shapH3(Nsize(1)+1),dshapH3(Nsize(1)+1)
      double precision :: shapQ1(Nsize(1))
      double precision :: shapQ2(Nsize(1))
      double precision :: shapQ3(Nsize(1))
!
#if HP3D_DEBUG
!  ...debugging flag
      integer :: iprint
      iprint=0
#endif
!
!  ...Option 1: Simply call the usual shape functions with enrichment
!      noriF(1:6)=0
!      call decod(NordM,MODORDER,2, (/nordF(1),nordB(3)/))
!      call decod(nordF(1),MODORDER,2, nordB(1:2))
!      call encod((/nordB(1),nordB(3)/),MODORDER,2, nordF(2))
!      call encod(nordB(2:3),MODORDER,2, nordF(3))
!      norder(1:12)=1
!      norder(13:14)=nordF(1)
!      norder(15:18)=(/nordF(2),nordF(3),nordF(2),nordF(3)/)
!      norder(19)=NordM
!      call shape3DVHexa(Xi,norder,noriF,Nsize, NrdofV,ShapV,DivV)
!
!
!  ...Option 2: Write more efficient routine for enriched functions
!
!  ...initiate counter for shape functions
      m=0
!
!  ...shape functions are tensor products of 1D shape functions
      call decod(NordM,MODORDER,3, nordB)
      call shape1HH(Xi(1),nordB(1), ndofH(1),shapH1,dshapH1)
      call shape1HH(Xi(2),nordB(2), ndofH(2),shapH2,dshapH2)
      call shape1HH(Xi(3),nordB(3), ndofH(3),shapH3,dshapH3)
      call shape1QQ(Xi(1),nordB(1), ndofQ(1),shapQ1)
      call shape1QQ(Xi(2),nordB(2), ndofQ(2),shapQ2)
      call shape1QQ(Xi(3),nordB(3), ndofQ(3),shapQ3)
!
!  ...shape functions with values along the x-axis
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
!
!  ...shape functions with values along the y-axis
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
!
!  ...shape functions with values along the z-axis
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
!
!  ...give total degrees of freedom
      NrdofV = m
!
#if HP3D_DEBUG
!  ...print this when debugging
      if (iprint.ge.1) then
        write(*,7001) Xi(1:3)
 7001   format('shape3DVBrokenHexa: Xi = ',3f8.3)
        do m=1,NrdofV
          write(*,7002) m,ShapV(1:3,m),DivV(m)
 7002     format('k = ',i3,' ShapV, DivV = ',3e12.5,3x,e12.5)
        enddo
        call pause
      endif
#endif
!
!
   end subroutine shape3DVBrokenHexa
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape3DQBrokenHexa
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 17
!
!     purpose:          - routine returns values of 3D hexahedron
!                         BROKEN L2 shape functions
!
!     arguments:
!
!     in:
!          X            - master hexahedron coordinates from (0,1)^3
!          NordM        - polynomial order for middle node (H1 sense)
!          Nsize        - relevant sizes of local arrays
!
!     out:
!          NrdofQ       - number of dof
!          ShapQ        - values of the shape functions at the point
!
!-----------------------------------------------------------------------
!
   subroutine shape3DQBrokenHexa(Xi,NordM,Nsize, NrdofQ,ShapQ)
!
      use parameters , only : MODORDER
!
      implicit none
      integer, intent(in)  :: NordM,Nsize(2)
      integer, intent(out) :: NrdofQ
      integer :: norder(19),nordF
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapQ(Nsize(2))
!
!  ...Option 1: Simply call the usual shape functions with enrichment
      norder(1:12)=1
      call encod((/1,1/),MODORDER,2, nordF)
      norder(13:18)=nordF
      norder(19)=NordM
      call shape3DQHexa(Xi,norder,Nsize, NrdofQ,ShapQ)
!
!  ...Option 2: Write a separate routine for enriched functions
!
!
   end subroutine shape3DQBrokenHexa


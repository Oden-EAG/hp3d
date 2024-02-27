c Routines:
c  - shape2DHTri
c  - shape2DETri
c  - shape2DVTri
c  - shape2DQTri
c
c----------------------------------------------------------------------
c
c     routine name      - shape2DHTri
c
c----------------------------------------------------------------------
c
c     latest revision:  - Nov 14, Apr 17, Jul 21
c
c     purpose:          - evaluate triangle H1 shape functions and
c                         their gradient
c
c     arguments:
c
c     in:
c       X               - master element coordinates
c       Nord            - polynomial order for the nodes (H1 sense)
c       NoriE           - edge orientations
c       Nsize           - relevant sizes of local arrays
c
c     out:
c       NrdofH          - number of dof
c       ShapH           - values of the shape functions
c       GradH           - gradients of the shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape2DHTri(X,Nord,NoriE,Nsize, NrdofH,ShapH,GradH)
c
      implicit none
      integer, intent(in)  :: Nord(4),NoriE(3),Nsize(2)
      integer, intent(out) :: NrdofH
      integer :: i,j,nij,m,v,e,N,nordE,ndofE,nordF,ndofF
      integer :: minI,minJ,minIJ,maxI,maxJ,maxIJ
      logical :: IdecE,IdecF
      double precision, intent(in)  :: X(2)
      double precision, intent(out) :: ShapH(Nsize(2))
      double precision, intent(out) :: GradH(2,Nsize(2))
      double precision :: Nu(0:2),DNu(2,0:2)
      double precision :: NubV(3),DNubV(2,3)
      double precision :: NupE(0:1,3),DNupE(2,0:1,3)
      double precision :: GNupE(0:1),GDNupE(2,0:1)
      double precision :: phiE(2:Nsize(1)),DphiE(2,2:Nsize(1))
      double precision :: phiTri(2:Nsize(1)-1,1:Nsize(1)-2)
      double precision :: DphiTri(2,2:Nsize(1)-1,1:Nsize(1)-2)
c
#if HP3D_DEBUG
c  ...debugging flag
      integer :: iprint
      iprint=0
#endif
c
c  ...spatial dimensions
      N=2
c
c  ...initiate counter for shape functions
      m=0
c
c  ...local parameters
      minI = 2; minJ = 1
      minIJ = minI+minJ
c
c  ...Define affine coordinates and gradients
      call AffineTriangle(X, Nu,DNu)
c
c  ...VERTEX SHAPE FUNCTIONS
      call BlendTriV(Nu,DNu, NubV,DNubV)
      do v=1,3
        m=m+1
c
        ShapH(m)     = NubV(v)
        GradH(1:N,m) = DNubV(1:N,v)
      enddo
c
c  ...EDGE SHAPE FUNCTIONS
      call ProjectTriE(Nu,DNu, NupE,DNupE,IdecE)
c  ...loop over edges
      do e=1,3
c    ...local parameters
        nordE = Nord(e)
        ndofE = nordE-1
        if (ndofE.gt.0) then
c      ...local parameters (again)
          maxI = nordE
c      ...orient
          call OrientE(NupE(0:1,e),DNupE(1:N,0:1,e),NoriE(e),N,
     .                                                    GNupE,GDNupE)
c      ...construct the shape functions
          call AncPhiE(GNupE,GDNupE,nordE,IdecE,N,
     .                            phiE(minI:maxI),DphiE(1:N,minI:maxI))
          do i=minI,maxI
            m=m+1
c
            ShapH(m)     = phiE(i)
            GradH(1:N,m) = DphiE(1:N,i)
          enddo
        endif
      enddo
c
c  ...FACE BUBBLE FUNCTIONS
c  ...local parameters
      nordF = Nord(4)
      ndofF = (nordF-1)*(nordF-2)/2
      IdecF = .TRUE.
      if (ndofF.gt.0) then
c    ...local parameters (again)
        maxIJ = nordF
        maxI = maxIJ-minJ
        maxJ = maxIJ-minI
c    ...construct the shape functions
        call AncPhiTri(Nu,DNu,nordF,IdecF,N,
     .                                     phiTri(minI:maxI,minJ:maxJ),
     .                                DphiTri(1:N,minI:maxI,minJ:maxJ))
        do nij=minIJ,maxIJ
          do i=minI,nij-minJ
            j=nij-i
            m=m+1
c
            ShapH(m)     = phiTri(i,j)
            GradH(1:N,m) = DphiTri(1:N,i,j)
          enddo
        enddo
      endif
c
c  ...give total degrees of freedom
      NrdofH = m
c
#if HP3D_DEBUG
c  ...print this when debugging
      if (iprint.eq.1) then
        write(*,7001) X(1:2),Nord(1:4),NoriE(1:3)
 7001   format('shape2DHTri: Xi = ',2f8.3,/,
     .         'Norder  = ',3i2,2x,i2,/,
     .         'Norient = ',3i2)
        write(*,7002)
 7002   format('VERTEX SHAPE FUNCTIONS = ')
        do v=1,3
          m=v
          write(*,7003) m,ShapH(m),GradH(1:2,m)
 7003     format('k = ',i3,' ShapH, GradH = ',e12.5,3x,2e12.5)
        enddo
        do e=1,3
          ndofE = Nord(e)-1
          if (ndofE.gt.0) then
            write(*,7004) e
 7004       format('EDGE SHAPE FUNCTIONS = ',i2)
            do j=1,ndofE
              m=m+1
              write(*,7003) m,ShapH(m),GradH(1:2,m)
            enddo
          endif
        enddo
        if (ndofF.gt.0) then
          write(*,7005)
 7005     format('FACE BUBBLES = ')
          do j=1,ndofF
            m=m+1
            write(*,7003) m,ShapH(m),GradH(1:2,m)
          enddo
        endif
        call pause
      endif
#endif
c
      end subroutine shape2DHTri
c
c
c----------------------------------------------------------------------
c
c     routine name      - shape2DETri
c
c----------------------------------------------------------------------
c
c     latest revision:  - Nov 14, Apr 17, Jul 21
c
c     purpose:          - evaluate triangle H(curl) shape functions and
c                         their curls
c
c     arguments:
c
c     in:
c          X            - master triangle coordinates from (0,1)^2
c          Nord         - polynomial order for the nodes (H1 sense)
c          NoriE        - edge orientations
c          Nsize        - relevant sizes of local arrays
c
c     out:
c          NrdofE       - number of dof
c          ShapE        - values of the shape functions
c          CurlE        - curls of the shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape2DETri(X,Nord,NoriE,Nsize, NrdofE,ShapE,CurlE)
c
      implicit none
      integer, intent(in)  :: Nord(4),NoriE(3),Nsize(2)
      integer, intent(out) :: NrdofE
      integer :: i,j,nij,m,e,N,nordE,ndofE,nordF,ndofF
      integer :: minI,minJ,minIJ,maxI,maxJ,maxIJ,abc(3),fam,famctr
      logical :: IdecE,IdecF
      double precision, intent(in)  :: X(2)
      double precision, intent(out) :: ShapE(2,Nsize(2))
      double precision, intent(out) :: CurlE(Nsize(2))
      double precision :: Nu(0:2),DNu(2,0:2)
      double precision :: NupE(0:1,3),DNupE(2,0:1,3)
      double precision :: GNupE(0:1),GDNupE(2,0:1)
      double precision :: EE(2,0:Nsize(1)-1),CurlEE(0:Nsize(1)-1)
      double precision :: ETri(2,0:Nsize(1)-2,1:Nsize(1)-1)
      double precision :: CurlETri(0:Nsize(1)-2,1:Nsize(1)-1)
c
#if HP3D_DEBUG
c  ...debugging flag
      integer :: iprint
      iprint=0
#endif
c
c  ...spatial dimensions
      N=2
c
c  ...initiate counter for shape functions
      m=0
c
c  ...local parameters
      minI = 0; minJ = 1
      minIJ = minI+minJ
c
c  ...Define affine coordinates and gradients
      call AffineTriangle(X, Nu,DNu)
c
c  ...EDGE SHAPE FUNCTIONS
      call ProjectTriE(Nu,DNu, NupE,DNupE,IdecE)
c  ...loop over edges
      do e=1,3
c    ...local parameters
        nordE = Nord(e)
        ndofE = nordE
        if (ndofE.gt.0) then
c      ...local parameters (again)
          maxI = nordE-1
c      ...orient first
          call OrientE(NupE(0:1,e),DNupE(1:N,0:1,e),NoriE(e),N,
     .                                                    GNupE,GDNupE)
c      ...construct the shape functions
          call AncEE(GNupE,GDNupE,nordE,IdecE,N,
     .                             EE(1:N,minI:maxI),CurlEE(minI:maxI))
          do i=minI,maxI
            m=m+1
c
            ShapE(1:N,m) = EE(1:N,i)
            CurlE(m)     = CurlEE(i)
          enddo
        endif
      enddo
c
c  ...FACE BUBBLE FUNCTIONS
c  ...local parameters
      nordF = Nord(4)
      ndofF = nordF*(nordF-1)/2
      IdecF = .TRUE.
      if (ndofF.gt.0) then
c    ...local parameters (again)
        maxIJ = nordF-1
        maxI = maxIJ-minJ
        maxJ = maxIJ-minI
c    ...loop over families
        famctr=m
        do fam=0,1
          m=famctr+fam-1
          abc = cshift((/0,1,2/),fam)
c      ...construct the shape functions
          call AncETri(Nu(abc),DNu(1:N,abc),nordF,IdecF,N,
     .                                   ETri(1:N,minI:maxI,minJ:maxJ),
     .                                   CurlETri(minI:maxI,minJ:maxJ))
          do nij=minIJ,maxIJ
            do i=minI,nij-minJ
              j=nij-i
              m=m+2
c
              ShapE(1:N,m) = ETri(1:N,i,j)
              CurlE(m)     = CurlETri(i,j)
            enddo
          enddo
        enddo
      endif
c
c  ...give total degrees of freedom
      NrdofE = m
c
#if HP3D_DEBUG
c  ...print this when debugging
      if (iprint.eq.1) then
        write(*,7001) X(1:2),Nord(1:4),NoriE(1:3)
 7001   format('shape2DETri: Xi = ',2f8.3,/,
     .         'Norder  = ',3i2,2x,i2,/,
     .         'Norient = ',3i2)
        m=0
        do e=1,3
          ndofE = Nord(e)
          if (ndofE.gt.0) then
            write(*,7002) e
 7002       format('SHAPE FUNCTIONS FOR EDGE = ',i2)
            do j=1,ndofE
              m=m+1
              write(*,7003) m,ShapE(1:2,m),CurlE(m)
 7003         format('k = ',i3,' ShapE, CurlE = ',2e12.5,3x,e12.5)
            enddo
          endif
        enddo
        if (ndofF.gt.0) then
          write(*,7004)
 7004     format('FACE BUBBLES = ')
          famctr=m
          do fam=0,1
            m=famctr+fam-1
            write(*,7005) fam
 7005       format('family = ',i2)
            do j=1,ndofF
              m=m+2
              write(*,7003) m,ShapE(1:2,m),CurlE(m)
            enddo
          enddo
        endif
        call pause
      endif
#endif
c
      end subroutine shape2DETri
c
c
c----------------------------------------------------------------------
c
c     routine name      - shape2DVTri
c
c----------------------------------------------------------------------
c
c     latest revision:  - Nov 14, Apr 17
c
c     purpose:          - evaluate triangle H(div) shape functions and
c                         their divergences
c
c     arguments :
c
c     in:
c        Xi             - master element coordinates
c        Nord           - polynomial order for the nodes (H1 sense)
c        NoriE          - edge orientations
c        Nsize          - relevant sizes of local arrays
c
c     out:
c        NrdofV         - number of dof
c        ShapV          - values of shape functions
c        DivV           - divergences of the shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape2DVTri(X,Nord,NoriE,Nsize, NrdofV,ShapV,DivV)
c
      implicit none
      integer, intent(in)  :: Nord(4),NoriE(3),Nsize(2)
      integer, intent(out) :: NrdofV
c
      double precision, intent(in)  :: X(2)
      double precision, intent(out) :: ShapV(2,Nsize(2))
      double precision, intent(out) :: DivV(Nsize(2))
c
      double precision :: shapE(2,Nsize(2))
      integer          :: m
c
#if HP3D_DEBUG
      integer :: j,e,ndofE,nordF,ndofF,famctr,fam
c  ...debugging flag
      integer :: iprint
      iprint=0
#endif
c
c  ...compute H(curl) shape functions
c  ...remember that NrdofE = NrdofV, div(V) = curl(E)
      call shape2DETri(X,Nord,NoriE,Nsize, NrdofV,shapE,DivV)
c
c  ...'rotate' shape functions
      do m=1,NrdofV
        ShapV(1,m) = shapE(2,m)
        ShapV(2,m) = -shapE(1,m)
      end do
c
#if HP3D_DEBUG
c  ...print this when debugging
      if (iprint.eq.1) then
        write(*,7001) X(1:2),Nord(1:4),NoriE(1:3)
 7001   format('shape2DVTri: Xi = ',2f8.3,/,
     .         'Norder  = ',3i2,2x,i2,/,
     .         'Norient = ',3i2)
        m=0
        do e=1,3
          ndofE = Nord(e)
          if (ndofE.gt.0) then
            write(*,7002) e
 7002       format('SHAPE FUNCTIONS FOR EDGE = ',i2)
            do j=1,ndofE
              m=m+1
              write(*,7003) m,ShapV(1:2,m),DivV(m)
 7003         format('k = ',i3,' ShapV, DivV = ',2e12.5,3x,e12.5)
            enddo
          endif
        enddo
        nordF = Nord(4)
        ndofF = nordF*(nordF-1)/2
        if (ndofF.gt.0) then
          write(*,7004)
 7004     format('FACE BUBBLES = ')
          famctr=m
          do fam=0,1
            m=famctr+fam-1
            write(*,7005) fam
 7005       format('family = ',i2)
            do j=1,ndofF
              m=m+2
              write(*,7003) m,ShapV(1:2,m),DivV(m)
            enddo
          enddo
        endif
        call pause
      endif
#endif
c
      end subroutine shape2DVTri
c
c
c----------------------------------------------------------------------
c
c     routine name      - shape2DQTri
c
c----------------------------------------------------------------------
c
c     latest revision:  - Nov 14, Apr 17
c
c     purpose:          - evaluate triangle L2 shape functions
c
c     arguments :
c
c     in:
c        Xi             - master element coordinates
c        Nord           - polynomial order of face node (H1 sense)
c        Nsize          - relevant sizes of local arrays
c
c     out:
c        NrdofQ         - number of dof
c        ShapQ          - values of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape2DQTri(X,Nord,Nsize, NrdofQ,ShapQ)
c
      implicit none
      integer, intent(in)  :: Nord(4),Nsize(2)
      integer, intent(out) :: NrdofQ
      integer :: i,j,nij,m,N,ndofF
      integer :: minalpha,minI,minJ,minIJ,maxI,maxJ,maxIJ
      double precision, intent(in)  :: X(2)
      double precision, intent(out) :: ShapQ(Nsize(2))
      double precision :: Nu(0:2),DNu(2,0:2)
      double precision :: homP(0:Nsize(1)-1)
      double precision :: homPal(0:Nsize(1)-1,0:Nsize(1)-1)
c
#if HP3D_DEBUG
c  ...debugging flag
      integer :: iprint
      iprint=0
#endif
c
c  ...spatial dimensions
      N=2
c
c  ...initiate counter for shape functions
      m=0
c
c  ...Define affine coordinates and gradients
      call AffineTriangle(X, Nu,DNu)
c
c  ...order and dof
      ndofF = (Nord(4)+1)*Nord(4)/2
      if (ndofF.gt.0) then
c
c    ...local parameters
        minI = 0; minJ = 0
        minIJ = minI+minJ
        maxIJ = Nord(4)-1
        maxI = maxIJ-minJ
        maxJ = maxIJ-minI
        minalpha = 2*minI+1
c    ...construct shape functions with homogenized Legendre and Jacobi
c    ...polynomials: homP and homPal respectively
        call HomLegendre(Nu(0:1),maxI, homP(minI:maxI))
        call HomJacobi((/Nu(0)+Nu(1),Nu(2)/),maxJ,minalpha,
     .                                     homPal(minI:maxI,minJ:maxJ))
c    ...construct the shape functions
        do nij=minIJ,maxIJ
          do i=minI,nij-minJ
            j=nij-i
            m=m+1
c
            ShapQ(m) = homP(i)*homPal(i,j)
          enddo
        enddo
      endif
c
c  ...give total degrees of freedom
      NrdofQ = m
c
#if HP3D_DEBUG
c  ...print this when debugging
      if (iprint.eq.1) then
        write(*,7002) X(1:2),Nord(4)
 7002   format('shape2DQTri: Xi = ',2f8.3,/,
     .         'Norder  = ',i2)
        if (ndofF.gt.0) then
          write(*,7003)
 7003     format('FACE FUNCTIONS = ')
          do m=1,ndofF
            write(*,7004) m,ShapQ(m)
 7004       format('k = ',i3,' ShapQ = ',e12.5)
          enddo
        endif
        call pause
      endif
#endif
c
      end subroutine shape2DQTri
c

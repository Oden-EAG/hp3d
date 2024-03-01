! Routines:
!  - shape2DHTri
!  - shape2DETri
!  - shape2DVTri
!  - shape2DQTri
!
!----------------------------------------------------------------------
!
!     routine name      - shape2DHTri
!
!----------------------------------------------------------------------
!
!     latest revision:  - Nov 14, Apr 17, Jul 21
!
!     purpose:          - evaluate triangle H1 shape functions and
!                         their gradient
!
!     arguments:
!
!     in:
!       X               - master element coordinates
!       Nord            - polynomial order for the nodes (H1 sense)
!       NoriE           - edge orientations
!       Nsize           - relevant sizes of local arrays
!
!     out:
!       NrdofH          - number of dof
!       ShapH           - values of the shape functions
!       GradH           - gradients of the shape functions
!
!----------------------------------------------------------------------
!
   subroutine shape2DHTri(X,Nord,NoriE,Nsize, NrdofH,ShapH,GradH)
!
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
!
#if HP3D_DEBUG
!  ...debugging flag
      integer :: iprint
      iprint=0
#endif
!
!  ...spatial dimensions
      N=2
!
!  ...initiate counter for shape functions
      m=0
!
!  ...local parameters
      minI = 2; minJ = 1
      minIJ = minI+minJ
!
!  ...Define affine coordinates and gradients
      call AffineTriangle(X, Nu,DNu)
!
!  ...VERTEX SHAPE FUNCTIONS
      call BlendTriV(Nu,DNu, NubV,DNubV)
      do v=1,3
        m=m+1
!
        ShapH(m)     = NubV(v)
        GradH(1:N,m) = DNubV(1:N,v)
      enddo
!
!  ...EDGE SHAPE FUNCTIONS
      call ProjectTriE(Nu,DNu, NupE,DNupE,IdecE)
!  ...loop over edges
      do e=1,3
!    ...local parameters
        nordE = Nord(e)
        ndofE = nordE-1
        if (ndofE.gt.0) then
!      ...local parameters (again)
          maxI = nordE
!      ...orient
          call OrientE(NupE(0:1,e),DNupE(1:N,0:1,e),NoriE(e),N, &
                                                          GNupE,GDNupE)
!      ...construct the shape functions
          call AncPhiE(GNupE,GDNupE,nordE,IdecE,N, &
                                  phiE(minI:maxI),DphiE(1:N,minI:maxI))
          do i=minI,maxI
            m=m+1
!
            ShapH(m)     = phiE(i)
            GradH(1:N,m) = DphiE(1:N,i)
          enddo
        endif
      enddo
!
!  ...FACE BUBBLE FUNCTIONS
!  ...local parameters
      nordF = Nord(4)
      ndofF = (nordF-1)*(nordF-2)/2
      IdecF = .TRUE.
      if (ndofF.gt.0) then
!    ...local parameters (again)
        maxIJ = nordF
        maxI = maxIJ-minJ
        maxJ = maxIJ-minI
!    ...construct the shape functions
        call AncPhiTri(Nu,DNu,nordF,IdecF,N, &
                                           phiTri(minI:maxI,minJ:maxJ), &
                                      DphiTri(1:N,minI:maxI,minJ:maxJ))
        do nij=minIJ,maxIJ
          do i=minI,nij-minJ
            j=nij-i
            m=m+1
!
            ShapH(m)     = phiTri(i,j)
            GradH(1:N,m) = DphiTri(1:N,i,j)
          enddo
        enddo
      endif
!
!  ...give total degrees of freedom
      NrdofH = m
!
#if HP3D_DEBUG
!  ...print this when debugging
      if (iprint.eq.1) then
        write(*,7001) X(1:2),Nord(1:4),NoriE(1:3)
 7001   format('shape2DHTri: Xi = ',2f8.3,/, &
               'Norder  = ',3i2,2x,i2,/, &
               'Norient = ',3i2)
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
!
   end subroutine shape2DHTri
!
!
!----------------------------------------------------------------------
!
!     routine name      - shape2DETri
!
!----------------------------------------------------------------------
!
!     latest revision:  - Nov 14, Apr 17, Jul 21
!
!     purpose:          - evaluate triangle H(curl) shape functions and
!                         their curls
!
!     arguments:
!
!     in:
!          X            - master triangle coordinates from (0,1)^2
!          Nord         - polynomial order for the nodes (H1 sense)
!          NoriE        - edge orientations
!          Nsize        - relevant sizes of local arrays
!
!     out:
!          NrdofE       - number of dof
!          ShapE        - values of the shape functions
!          CurlE        - curls of the shape functions
!
!----------------------------------------------------------------------
!
   subroutine shape2DETri(X,Nord,NoriE,Nsize, NrdofE,ShapE,CurlE)
!
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
!
#if HP3D_DEBUG
!  ...debugging flag
      integer :: iprint
      iprint=0
#endif
!
!  ...spatial dimensions
      N=2
!
!  ...initiate counter for shape functions
      m=0
!
!  ...local parameters
      minI = 0; minJ = 1
      minIJ = minI+minJ
!
!  ...Define affine coordinates and gradients
      call AffineTriangle(X, Nu,DNu)
!
!  ...EDGE SHAPE FUNCTIONS
      call ProjectTriE(Nu,DNu, NupE,DNupE,IdecE)
!  ...loop over edges
      do e=1,3
!    ...local parameters
        nordE = Nord(e)
        ndofE = nordE
        if (ndofE.gt.0) then
!      ...local parameters (again)
          maxI = nordE-1
!      ...orient first
          call OrientE(NupE(0:1,e),DNupE(1:N,0:1,e),NoriE(e),N, &
                                                          GNupE,GDNupE)
!      ...construct the shape functions
          call AncEE(GNupE,GDNupE,nordE,IdecE,N, &
                                   EE(1:N,minI:maxI),CurlEE(minI:maxI))
          do i=minI,maxI
            m=m+1
!
            ShapE(1:N,m) = EE(1:N,i)
            CurlE(m)     = CurlEE(i)
          enddo
        endif
      enddo
!
!  ...FACE BUBBLE FUNCTIONS
!  ...local parameters
      nordF = Nord(4)
      ndofF = nordF*(nordF-1)/2
      IdecF = .TRUE.
      if (ndofF.gt.0) then
!    ...local parameters (again)
        maxIJ = nordF-1
        maxI = maxIJ-minJ
        maxJ = maxIJ-minI
!    ...loop over families
        famctr=m
        do fam=0,1
          m=famctr+fam-1
          abc = cshift((/0,1,2/),fam)
!      ...construct the shape functions
          call AncETri(Nu(abc),DNu(1:N,abc),nordF,IdecF,N, &
                                         ETri(1:N,minI:maxI,minJ:maxJ), &
                                         CurlETri(minI:maxI,minJ:maxJ))
          do nij=minIJ,maxIJ
            do i=minI,nij-minJ
              j=nij-i
              m=m+2
!
              ShapE(1:N,m) = ETri(1:N,i,j)
              CurlE(m)     = CurlETri(i,j)
            enddo
          enddo
        enddo
      endif
!
!  ...give total degrees of freedom
      NrdofE = m
!
#if HP3D_DEBUG
!  ...print this when debugging
      if (iprint.eq.1) then
        write(*,7001) X(1:2),Nord(1:4),NoriE(1:3)
 7001   format('shape2DETri: Xi = ',2f8.3,/, &
               'Norder  = ',3i2,2x,i2,/, &
               'Norient = ',3i2)
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
!
   end subroutine shape2DETri
!
!
!----------------------------------------------------------------------
!
!     routine name      - shape2DVTri
!
!----------------------------------------------------------------------
!
!     latest revision:  - Nov 14, Apr 17
!
!     purpose:          - evaluate triangle H(div) shape functions and
!                         their divergences
!
!     arguments :
!
!     in:
!        Xi             - master element coordinates
!        Nord           - polynomial order for the nodes (H1 sense)
!        NoriE          - edge orientations
!        Nsize          - relevant sizes of local arrays
!
!     out:
!        NrdofV         - number of dof
!        ShapV          - values of shape functions
!        DivV           - divergences of the shape functions
!
!----------------------------------------------------------------------
!
   subroutine shape2DVTri(X,Nord,NoriE,Nsize, NrdofV,ShapV,DivV)
!
      implicit none
      integer, intent(in)  :: Nord(4),NoriE(3),Nsize(2)
      integer, intent(out) :: NrdofV
!
      double precision, intent(in)  :: X(2)
      double precision, intent(out) :: ShapV(2,Nsize(2))
      double precision, intent(out) :: DivV(Nsize(2))
!
      double precision :: shapE(2,Nsize(2))
      integer          :: m
!
#if HP3D_DEBUG
      integer :: j,e,ndofE,nordF,ndofF,famctr,fam
!  ...debugging flag
      integer :: iprint
      iprint=0
#endif
!
!  ...compute H(curl) shape functions
!  ...remember that NrdofE = NrdofV, div(V) = curl(E)
      call shape2DETri(X,Nord,NoriE,Nsize, NrdofV,shapE,DivV)
!
!  ...'rotate' shape functions
      do m=1,NrdofV
        ShapV(1,m) = shapE(2,m)
        ShapV(2,m) = -shapE(1,m)
      end do
!
#if HP3D_DEBUG
!  ...print this when debugging
      if (iprint.eq.1) then
        write(*,7001) X(1:2),Nord(1:4),NoriE(1:3)
 7001   format('shape2DVTri: Xi = ',2f8.3,/, &
               'Norder  = ',3i2,2x,i2,/, &
               'Norient = ',3i2)
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
!
   end subroutine shape2DVTri
!
!
!----------------------------------------------------------------------
!
!     routine name      - shape2DQTri
!
!----------------------------------------------------------------------
!
!     latest revision:  - Nov 14, Apr 17
!
!     purpose:          - evaluate triangle L2 shape functions
!
!     arguments :
!
!     in:
!        Xi             - master element coordinates
!        Nord           - polynomial order of face node (H1 sense)
!        Nsize          - relevant sizes of local arrays
!
!     out:
!        NrdofQ         - number of dof
!        ShapQ          - values of shape functions
!
!----------------------------------------------------------------------
!
   subroutine shape2DQTri(X,Nord,Nsize, NrdofQ,ShapQ)
!
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
!
#if HP3D_DEBUG
!  ...debugging flag
      integer :: iprint
      iprint=0
#endif
!
!  ...spatial dimensions
      N=2
!
!  ...initiate counter for shape functions
      m=0
!
!  ...Define affine coordinates and gradients
      call AffineTriangle(X, Nu,DNu)
!
!  ...order and dof
      ndofF = (Nord(4)+1)*Nord(4)/2
      if (ndofF.gt.0) then
!
!    ...local parameters
        minI = 0; minJ = 0
        minIJ = minI+minJ
        maxIJ = Nord(4)-1
        maxI = maxIJ-minJ
        maxJ = maxIJ-minI
        minalpha = 2*minI+1
!    ...construct shape functions with homogenized Legendre and Jacobi
!    ...polynomials: homP and homPal respectively
        call HomLegendre(Nu(0:1),maxI, homP(minI:maxI))
        call HomJacobi((/Nu(0)+Nu(1),Nu(2)/),maxJ,minalpha, &
                                           homPal(minI:maxI,minJ:maxJ))
!    ...construct the shape functions
        do nij=minIJ,maxIJ
          do i=minI,nij-minJ
            j=nij-i
            m=m+1
!
            ShapQ(m) = homP(i)*homPal(i,j)
          enddo
        enddo
      endif
!
!  ...give total degrees of freedom
      NrdofQ = m
!
#if HP3D_DEBUG
!  ...print this when debugging
      if (iprint.eq.1) then
        write(*,7002) X(1:2),Nord(4)
 7002   format('shape2DQTri: Xi = ',2f8.3,/, &
               'Norder  = ',i2)
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
!
   end subroutine shape2DQTri
!

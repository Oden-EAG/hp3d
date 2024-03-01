! Routines:
!  - TraceEshapeH
!  - TraceEshapeE
!  - TraceEshapeV
!  - TraceFshapeH
!  - TraceFshapeE
!  - TraceFshapeV
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!                    2D:   H1  --->  Hcurl  --->  L2
!             Trace(2D): tr(H1)--->tr(Hcurl)
!
!                    2D:   H1  --->  Hdiv   --->  L2  (rotated)
!             Trace(2D): tr(H1)--->tr(Hdiv)           (rotated)
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!                             Trace(2D H1)
!----------------------------------------------------------------------
!
!     routine name      - TraceEshapeH
!
!----------------------------------------------------------------------
!
!     latest revision:  - Nov 14
!
!     purpose:          - 1D edge trace of 2D H1 edge functions
!                         (these happen to be oriented 1D H1 'bubbles')
!
!     arguments:
!
!     in:
!          T            - local edge coordinate
!          Nord         - polynomial edge order (H1 sense)
!          Nori         - edge orientation
!
!     out:
!          NrdofH       - number of trace shape functions
!          ShapH        - values of trace shape functions
!          GradH        - local gradients of trace shape functions
!
!----------------------------------------------------------------------
!
   subroutine TraceEshapeH(T,Nord,Nori, NrdofH,ShapH,GradH)
!
      use parameters , only : MAXP
!
      implicit none
      integer, intent(in)  :: Nord,Nori
      integer, intent(out) :: NrdofH
      integer :: N,m,ndofE,minI,maxI,i
      logical :: IdecE
      double precision, intent(in)  :: T
      double precision, intent(out) :: ShapH(MAXP-1),GradH(MAXP-1)
      double precision :: Mu(0:1),DMu(0:1),GMu(0:1),GDMu(0:1)
      double precision :: phiE(2:Nord),DphiE(2:Nord)
!
#if HP3D_DEBUG
!  ...debugging flag
      integer :: iprint
      iprint=0
#endif
!
!  ...spatial dimensions
      N=1
!
!  ...initiate counter for shape functions
      m=0
!
!  ...Define affine coordinates and gradients
      call AffineSegment(T, Mu,DMu)
!
!  ...TRACE OF 2D EDGE FUNCTIONS
      call checkpolyorder(Nord)
      ndofE = Nord-1
      if (ndofE.gt.0) then
!    ...local parameters
        minI  = 2
        maxI  = Nord
        IdecE = .true.
!    ...orient -- this is important for traces!!!
        call OrientE(Mu,DMu,Nori,N, GMu,GDMu)
!    ...construct the shape functions
        call AncPhiE(GMu,GDMu,Nord,IdecE,N, phiE,DphiE)
        do i=minI,maxI
          m=m+1
!
          ShapH(m) = phiE(i)
          GradH(m) = DphiE(i)
        enddo
      endif
!
!  ...give total degrees of freedom
      NrdofH = m
!
#if HP3D_DEBUG
!  ...print this when debugging
      if (iprint.eq.1) then
        write(*,7001) T,Nord,Nori
 7001   format('TraceEshapeH: T = ',f8.3,/, &
               'Norder  = ',i2,/, &
               'Norient = ',i2)
!
        if (ndofE.gt.0) then
          write(*,*) 'TRACE OF 2D H1 EDGE FUNCTIONS = '
          do m=1,ndofE
            write(*,7002) m,ShapH(m),GradH(m)
 7002       format('k = ',i3,' ShapH, GradH = ',e12.5,3x,e12.5)
          enddo
        endif
        call pause
      endif
#endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Check polynomial order
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine checkpolyorder(Norder)
!
        integer, intent(in) :: Norder
!
        if ((Norder.lt.1).or.(Norder.gt.MAXP)) then
          write(*,7003) Norder
          write(*,7004) MAXP
 7003     format('TraceEshapeH: Polynomial order = ',i3,' is either ')
 7004     format('              less than 1 or more than MAXP = ',i3)
          stop 1
        endif
!
      end subroutine checkpolyorder
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   end subroutine TraceEshapeH
!
!----------------------------------------------------------------------
!                           Trace(2D Hcurl)
!----------------------------------------------------------------------
!
!     routine name      - TraceEshapeE
!
!----------------------------------------------------------------------
!
!     latest revision:  - Nov 14
!
!     purpose:          - 1D edge trace of 2D H(curl) edge functions
!                         (these are oriented 1D L2 functions -
!                         not 'bubbles')
!
!     arguments:
!
!     in:
!          T            - local edge coordinate
!          Nord         - polynomial edge order (H1 sense)
!          Nori         - edge orientation
!
!     out:
!          NrdofE       - number of trace shape functions
!          ShapE        - values of trace shape functions
!
!-----------------------------------------------------------------------
!
   subroutine TraceEshapeE(T,Nord,Nori, NrdofE,ShapE)
!
      use parameters , only : MAXP
!
      implicit none
      integer, intent(in ) :: Nord,Nori
      integer, intent(out) :: NrdofE
      integer :: N,m,ndofE,minI,maxI,i
      double precision, intent(in ) :: T
      double precision, intent(out) :: ShapE(MAXP)
      double precision :: Mu(0:1),DMu(0:1),GMu(0:1),GDMu(0:1)
      double precision :: homP(0:Nord-1)
      double precision :: jac
!
#if HP3D_DEBUG
!  ...debugging flag
      integer :: iprint
      iprint=0
#endif
!
!  ...spatial dimensions
      N=1
!
!  ...initiate counter for shape functions
      m=0
!
!  ...Define affine coordinates and gradients
      call AffineSegment(T, Mu,DMu)
!
!  ...TRACE OF 2D EDGE FUNCTIONS
      call checkpolyorder(Nord)
      ndofE = Nord
      if (ndofE.gt.0) then
!    ...local parameters (again)
        minI  = 0
        maxI  = Nord-1
!    ...orient -- this is important for traces!!!
        call OrientE(Mu,DMu,Nori,N, GMu,GDMu)
!    ...construct the shape functions
        call HomLegendre(GMu,maxI, homP)
        do i=minI,maxI
          m=m+1
!
          jac = GDMu(1)
          ShapE(m) = homP(i)*jac
        enddo
      endif
!
!  ...give total degrees of freedom
      NrdofE = m
!
#if HP3D_DEBUG
!  ...print this when debugging
      if (iprint.eq.1) then
        write(*,7001) T,Nord,Nori
 7001   format('TraceEshapeE: T = ',f8.3,/, &
               'Norder  = ',i2,/, &
               'Norient = ',i2)
!
        if (ndofE.gt.0) then
          write(*,*) 'TRACE OF 2D H(curl) EDGE FUNCTIONS = '
          do m=1,ndofE
            write(*,7002) m,ShapE(m)
 7002       format('k = ',i3,'ShapE = ',e12.5)
          enddo
        endif
        call pause
      endif
#endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Check polynomial order
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine checkpolyorder(Norder)
!
        integer, intent(in) :: Norder
!
        if ((Norder.lt.1).or.(Norder.gt.MAXP)) then
          write(*,7003) Norder
          write(*,7004) MAXP
 7003     format('TraceEshapeE: Polynomial order = ',i3,' is either ')
 7004     format('              less than 1 or more than MAXP = ',i3)
          stop 1
        endif
!
      end subroutine checkpolyorder
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   end subroutine TraceEshapeE
!
!----------------------------------------------------------------------
!                           Trace(2D Hdiv)
!----------------------------------------------------------------------
!
!     routine name      - TraceEshapeV
!
!----------------------------------------------------------------------
!
!     latest revision:  - Nov 14
!
!     purpose:          - 1D edge trace of 2D H(div) edge functions
!                         (these are oriented 1D L2 functions -
!                         not 'bubbles')
!
!     arguments:
!
!     in:
!          T            - local edge coordinate
!          Nord         - polynomial edge order (H1 sense)
!          Nori         - edge orientation
!
!     out:
!          NrdofV       - number of trace shape functions
!          ShapV        - values of trace shape functions
!
!-----------------------------------------------------------------------
!
   subroutine TraceEshapeV(T,Nord,Nori, NrdofV,ShapV)
!
      use parameters , only : MAXP
!
      implicit none
      integer, intent(in ) :: Nord,Nori
      integer, intent(out) :: NrdofV
!
      double precision, intent(in ) :: T
      double precision, intent(out) :: ShapV(MAXP)
!
#if HP3D_DEBUG
      integer :: m,ndofE
!  ...debugging flag
      integer :: iprint
      iprint=0
#endif
!
!  ...the traces are the same as the H(curl) traces
      call TraceEshapeE(T,Nord,Nori, NrdofV,ShapV)
!
#if HP3D_DEBUG
!  ...print this when debugging
      if (iprint.eq.1) then
        write(*,7001) T,Nord,Nori
 7001   format('TraceEshapeV: T = ',f8.3,/, &
               'Norder  = ',i2,/, &
               'Norient = ',i2)
!
        ndofE = Nord
        if (ndofE.gt.0) then
          write(*,*) 'TRACE OF 2D H(div) EDGE FUNCTIONS = '
          do m=1,ndofE
            write(*,7002) m,ShapV(m)
 7002       format('k = ',i3,'ShapV = ',e12.5)
          enddo
        endif
        call pause
      endif
#endif
!
   end subroutine TraceEshapeV
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!              3D:   H1  --->  Hcurl  --->  Hdiv  --->  L2
!       Trace(3D): tr(H1)--->tr(Hcurl)--->tr(Hdiv)
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!                             Trace(3D H1)
!----------------------------------------------------------------------
!
!     routine name      - TraceFshapeH
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - 2D face trace of 3D H1 face functions
!                         (these happen to be oriented 2D H1 'bubbles')
!
!     arguments:
!
!     in:
!          Ftype        - face type (quad or triangle)
!          T            - local face coordinate
!          Nord         - polynomial face order (H1 sense)
!          Nori         - face orientation
!
!     out:
!          NrdofH       - number of trace shape functions
!          ShapH        - values of trace shape functions
!          GradH        - local gradients of trace shape functions
!
!----------------------------------------------------------------------
!
   subroutine TraceFshapeH(Ftype,T,Nord,Nori, NrdofH,ShapH,GradH)
!
      use parameters , only : MODORDER,MAXP,MAXmdlqH
      use node_types
!
      implicit none
      integer, intent(in)  :: Ftype,Nord,Nori
      integer, intent(out) :: NrdofH
      double precision, intent(in)  :: T(2)
      double precision, intent(out) :: ShapH(MAXmdlqH)
      double precision, intent(out) :: GradH(2,MAXmdlqH)
!
      select case(Ftype)
      case(TRIA,MDLT)
        call traceTriFshapeH
      case(QUAD,MDLQ,RECT)
        call traceQuadFshapeH
      case default
        write(*,*)'TraceFshapeH: Type = ', S_Type(Ftype)
        stop 1
      end select
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Triangle traces
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine traceTriFshapeH
!
        integer :: i,j,nij,m,N,ndofF
        integer :: minI,minJ,minIJ,maxI,maxJ,maxIJ
        logical :: IdecF
        double precision :: Nu(0:2),DNu(2,0:2)
        double precision :: GNu(0:2),GDNu(2,0:2)
        double precision :: phiTri(2:Nord-1,1:Nord-2)
        double precision :: DphiTri(2,2:Nord-1,1:Nord-2)
!
#if HP3D_DEBUG
!    ...debugging flag
        integer :: iprint
        iprint=0
#endif
!
!    ...spatial dimensions
        N=2
!
!    ...initiate counter for shape functions
        m=0
!
!    ...Define affine coordinates and gradients
        call AffineTriangle(T, Nu,DNu)
!
!    ...TRACE OF 3D FACE FUNCTIONS
        call checkpolyorder(Nord)
        ndofF = (Nord-1)*(Nord-2)/2
        IdecF = .true.
        if (ndofF.gt.0) then
!      ...local parameters
          minI  = 2
          minJ  = 1
          minIJ = minI+minJ
          maxIJ = Nord
          maxI  = maxIJ-minJ
          maxJ  = maxIJ-minI
!      ...orient -- this is important for traces!!!
          call OrientTri(Nu,DNu,Nori,N, GNu,GDNu)
!      ...construct the shape functions
          call AncPhiTri(GNu,GDNu,Nord,IdecF,N, phiTri,DphiTri)
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
!    ...give total degrees of freedom
        NrdofH = m
!
#if HP3D_DEBUG
!    ...print this when debugging
        if (iprint.eq.1) then
          write(*,7001) T(1:2),Nord,Nori
 7001     format('traceTriFshapeH: T = ',2f8.3,/, &
                 'Norder  = ',i2,/, &
                 'Norient = ',i2)
          if (ndofF.gt.0) then
            write(*,*) 'TRACE OF 3D H1 TRIANGLE FACE FUNCTIONS = '
            do m=1,ndofF
              write(*,7003) m,ShapH(m),GradH(1:2,m)
 7003         format('k = ',i3,' ShapH, GradH = ',e12.5,3x,2e12.5)
            enddo
          endif
          call pause
        endif
#endif
!
      end subroutine traceTriFshapeH
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Quadrilateral traces
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine traceQuadFshapeH
!
        integer :: N,m,i,j,nordF(2),ndofF,minI,minJ,maxI,maxJ
        logical :: IdecF(2),GIdecF(2)
        double precision :: Mu(0:1,2),DMu(2,0:1,2)
        double precision :: GMu(0:1,2),GDMu(2,0:1,2)
        double precision :: phiQuad(2:MAXP,2:MAXP)
        double precision :: DphiQuad(2,2:MAXP,2:MAXP)
!
#if HP3D_DEBUG
!    ...debugging flag
        integer :: iprint
        iprint=0
#endif
!
!    ...spatial dimensions
        N=2
!
!    ...initiate counter for shape functions
        m=0
!
!    ...Define affine coordinates and gradients
        call AffineQuadrilateral(T, Mu,DMu)
!
!    ...TRACE OF 3D FACE FUNCTIONS
        IdecF(1:2) = .true.
        call decod(Nord,MODORDER,2, nordF)
        call checkpolyorder(nordF(1))
        call checkpolyorder(nordF(2))
        ndofF = (nordF(1)-1)*(nordF(2)-1)
        if (ndofF.gt.0) then
!      ...local parameters
          minI = 2
          minJ = 2
          maxI = nordF(1)
          maxJ = nordF(2)
!      ...orient
          call OrientQuad(Mu,DMu,Nori,IdecF,N, GMu,GDMu,GIdecF)
!      ...construct the shape functions
          call AncPhiQuad(GMu,GDMu,nordF,GIdecF,N, &
                                          phiQuad(minI:maxI,minJ:maxJ), &
                                     DphiQuad(1:N,minI:maxI,minJ:maxJ))
          do j=minJ,maxJ
            do i=minI,maxI
              m=m+1
!
              ShapH(m)     = phiQuad(i,j)
              GradH(1:N,m) = DphiQuad(1:N,i,j)
            enddo
          enddo
        endif
!
!    ...give total degrees of freedom
        NrdofH = m
!
#if HP3D_DEBUG
!    ...print this when debugging
        if (iprint.eq.1) then
          write(*,7001) T(1:2),Nord,Nori
 7001     format('traceQuadFshapeH: T = ',2f8.3,/, &
                 'Norder  = ',i2,/, &
                 'Norient = ',i2)
          if (ndofF.gt.0) then
            write(*,*) 'TRACE OF 3D H1 QUAD FACE FUNCTIONS = '
            do m=1,ndofF
              write(*,7003) m,ShapH(m),GradH(1:2,m)
 7003         format('k = ',i3,' ShapH, GradH = ',e12.5,3x,2e12.5)
            enddo
          endif
          call pause
        endif
#endif
!
      end subroutine traceQuadFshapeH
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Check polynomial order
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine checkpolyorder(Norder)
!
        integer, intent(in) :: Norder
!
        if ((Norder.lt.1).or.(Norder.gt.MAXP)) then
          write(*,7003) Norder
          write(*,7004) MAXP
 7003     format('TraceFshapeH: Polynomial order = ',i3,' is either ')
 7004     format('              less than 1 or more than MAXP = ',i3)
          stop 1
        endif
!
      end subroutine checkpolyorder
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   end subroutine TraceFshapeH
!
!----------------------------------------------------------------------
!                          Trace(3D H(curl))
!----------------------------------------------------------------------
!
!     routine name      - TraceFshapeE
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - 2D face trace of 3D H(curl) face functions
!                         (these happen to be oriented 2D H(curl)
!                         'bubbles')
!
!     arguments:
!
!     in:
!          Ftype        - face type (quad or triangle)
!          T            - local face coordinate
!          Nord         - polynomial face order (H1 sense)
!          Nori         - face orientation
!
!     out:
!          NrdofE       - number of trace shape functions
!          ShapE        - values of trace shape functions
!          CurlE        - local curls of trace shape functions
!
!----------------------------------------------------------------------
!
   subroutine TraceFshapeE(Ftype,T,Nord,Nori, NrdofE,ShapE,CurlE)
!
      use parameters , only : MODORDER,MAXP,MAXmdlqE
      use node_types
!
      implicit none
      integer, intent(in)  :: Ftype,Nord,Nori
      integer, intent(out) :: NrdofE
      double precision, intent(in)  :: T(2)
      double precision, intent(out) :: ShapE(2,MAXmdlqE)
      double precision, intent(out) :: CurlE(MAXmdlqE)
!
      select case(Ftype)
      case(TRIA,MDLT)
        call traceTriFshapeE
      case(QUAD,MDLQ,RECT)
        call traceQuadFshapeE
      case default
        write(*,*)'TraceFshapeE: Type = ', S_Type(Ftype)
        stop 1
      end select
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Triangle traces
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine traceTriFshapeE
!
        integer :: i,j,nij,m,N,ndofF
        integer :: minI,minJ,minIJ,maxI,maxJ,maxIJ,abc(3),fam,famctr
        logical :: IdecF
        double precision :: Nu(0:2),DNu(2,0:2)
        double precision :: GNu(0:2),GDNu(2,0:2)
        double precision :: ETri(2,0:Nord-2,1:Nord-1)
        double precision :: CurlETri(0:Nord-2,1:Nord-1)
!
#if HP3D_DEBUG
!    ...debugging flag
        integer :: iprint
        iprint=0
#endif
!
!    ...spatial dimensions
        N=2
!
!    ...initiate counter for shape functions
        m=0
!
!    ...Define affine coordinates and gradients
        call AffineTriangle(T, Nu,DNu)
!
!    ...TRACE OF 3D FACE FUNCTIONS
        call checkpolyorder(Nord)
        ndofF = Nord*(Nord-1)/2
        IdecF = .true.
        if (ndofF.gt.0) then
!    ...local parameters
          minI  = 0
          minJ  = 1
          minIJ = minI+minJ
          maxIJ = Nord-1
          maxI  = maxIJ-minJ
          maxJ  = maxIJ-minI
!      ...orient -- this is important for traces!!!
          call OrientTri(Nu,DNu,Nori,N, GNu,GDNu)
!      ...loop over families
          famctr=m
          do fam=0,1
            m=famctr+fam-1
            abc = cshift((/0,1,2/),fam)
!        ...construct the shape functions
            call AncETri(GNu(abc),GDNu(1:N,abc),Nord,IdecF,N, &
                                                         ETri,CurlETri)
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
!    ...give total degrees of freedom
        NrdofE = m
!
#if HP3D_DEBUG
!    ...print this when debugging
        if (iprint.eq.1) then
          write(*,7001) T(1:2),Nord,Nori
 7001     format('traceTriFshapeE: T = ',2f8.3,/, &
                 'Norder  = ',i2,/, &
                 'Norient = ',i2)
          if (ndofF.gt.0) then
            write(*,*) 'TRACE OF 3D H(curl) TRIANGLE FACE FUNCTIONS = '
            famctr=0
            do fam=0,1
              m=famctr+fam-1
              write(*,7005) fam
 7005         format('family = ',i2)
              do j=1,ndofF
                m=m+2
                write(*,7003) m,ShapE(1:2,m),CurlE(m)
 7003           format('k = ',i3,' ShapE, CurlE = ',2e12.5,3x,e12.5)
              enddo
            enddo
          endif
          call pause
        endif
#endif
!
      end subroutine traceTriFshapeE
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Quadrilateral traces
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine traceQuadFshapeE
!
        integer :: N,m,i,j,ij(2),ig,jg,a,b,ab(2),fam
        integer :: nordF(2),ndofF(0:1),minF(2),maxF(2)
        logical :: IdecF(2),GIdecF(2)
        double precision :: Mu(0:1,2),DMu(2,0:1,2)
        double precision :: GMu(0:1,2),GDMu(2,0:1,2)
        double precision :: EQuad(2,0:MAXP-1,2:MAXP)
        double precision :: curlEQuad(0:MAXP-1,2:MAXP)
!
#if HP3D_DEBUG
!    ...debugging flag
        integer :: iprint
        iprint=0
#endif
!
!    ...spatial dimensions
        N=2
!
!    ...initiate counter for shape functions
        m=0
!
!    ...Define affine coordinates and gradients
        call AffineQuadrilateral(T, Mu,DMu)
!
!    ...TRACE OF 3D FACE FUNCTIONS
        IdecF(1:2) = .true.
        call decod(Nord,MODORDER,2, nordF)
        call checkpolyorder(nordF(1))
        call checkpolyorder(nordF(2))
!    ...orient -- this is important for traces!!!
        call OrientQuad(Mu,DMu,Nori,IdecF,N, GMu,GDMu,GIdecF)
!    ...loop over families
        do fam=0,1
          ab = cshift((/1,2/),fam)
          a = ab(1); b = ab(2)
          ndofF(fam) = nordF(a)*(nordF(b)-1)
          if (ndofF(fam).gt.0) then
!        ...local parameters
            minF(1) = 0
            minF(2) = 2
            maxF(1) = nordF(a)-1
            maxF(2) = nordF(b)
!        ...construct the shape functions
            call AncEQuad(GMu(0:1,ab),GDMu(1:N,0:1,ab), &
                                             nordF(ab),GIdecF(ab),N, &
                            EQuad(1:N,minF(1):maxF(1),minF(2):maxF(2)), &
                            curlEQuad(minF(1):maxF(1),minF(2):maxF(2)))
!        ...in the code the outer loop always is
!           numbered wrt the second global face axis
            minF = cshift(minF,-fam); maxF = cshift(maxF,-fam)
            do jg=minF(2),maxF(2)
              do ig=minF(1),maxF(1)
                ij = cshift((/ig,jg/),fam)
                i = ij(1); j = ij(2)
                m=m+1
!
                ShapE(1:N,m) = EQuad(1:N,i,j)
                CurlE(m)     = curlEQuad(i,j)
              enddo
            enddo
          endif
        enddo
!
!    ...give total degrees of freedom
        NrdofE = m
!
#if HP3D_DEBUG
!    ...print this when debugging
        if (iprint.eq.1) then
          write(*,7001) T(1:2),Nord,Nori
 7001     format('traceQuadFshapeE: T = ',2f8.3,/, &
                 'Norder  = ',i2,/, &
                 'Norient = ',i2)
          m=0
          if ((ndofF(0)+ndofF(1)).gt.0) then
            write(*,*) 'TRACE OF 3D H(curl) QUAD FACE FUNCTIONS = '
            do fam=0,1
              if (ndofF(fam).gt.0) then
                write(*,7005) fam
 7005           format('family = ',i2)
                do j=1,ndofF(fam)
                  m=m+1
                  write(*,7003) m,ShapE(1:2,m),CurlE(m)
 7003             format('k = ',i3,' ShapE, CurlE = ',2e12.5,3x,e12.5)
                enddo
              endif
            enddo
          endif
          call pause
        endif
#endif
!
      end subroutine traceQuadFshapeE
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Check polynomial order
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine checkpolyorder(Norder)
!
        integer, intent(in) :: Norder
!
        if ((Norder.lt.1).or.(Norder.gt.MAXP)) then
          write(*,7003) Norder
          write(*,7004) MAXP
 7003     format('TraceFshapeE: Polynomial order = ',i3,' is either ')
 7004     format('              less than 1 or more than MAXP = ',i3)
          stop 1
        endif
!
      end subroutine checkpolyorder
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   end subroutine TraceFshapeE
!
!----------------------------------------------------------------------
!                          Trace(3D H(div))
!----------------------------------------------------------------------
!
!     routine name      - TraceFshapeV
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - 2D face trace of 3D H(div) face functions
!                         these are oriented 2D L2 face funtions
!
!     arguments:
!
!     in:
!          Ftype        - face type (quad or triangle)
!          T            - local face coordinate
!          Nord         - polynomial face order (H1 sense)
!          Nori         - face orientation
!
!     out:
!          NrdofV       - number of trace shape functions
!          ShapV        - values of trace shape functions
!
!----------------------------------------------------------------------
!
   subroutine TraceFshapeV(Ftype,T,Nord,Nori, NrdofV,ShapV)
!
      use parameters , only : MODORDER,MAXP,MAXmdlqV
      use node_types
!
      implicit none
      integer, intent(in)  :: Ftype,Nord,Nori
      integer, intent(out) :: NrdofV
      double precision, intent(in)  :: T(2)
      double precision, intent(out) :: ShapV(MAXmdlqV)
!
      select case(Ftype)
      case(TRIA,MDLT)
        call traceTriFshapeV
      case(QUAD,MDLQ,RECT)
        call traceQuadFshapeV
      case default
        write(*,*)'TraceFshapeV: Type = ', S_Type(Ftype)
        stop 1
      end select
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Triangle traces
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine traceTriFshapeV
!
        integer :: i,j,nij,m,N,ndofF
        integer :: minalpha,minI,minJ,minIJ,maxI,maxJ,maxIJ
        double precision :: Nu(0:2),DNu(2,0:2)
        double precision :: GNu(0:2),GDNu(2,0:2)
        double precision :: homP(0:Nord-1),homPal(0:Nord-1,0:Nord-1)
        double precision :: jac(1)
!
#if HP3D_DEBUG
!    ...debugging flag
        integer :: iprint
        iprint=0
#endif
!
!    ...spatial dimensions
        N=2
!
!    ...initiate counter for shape functions
        m=0
!
!    ...Define affine coordinates and gradients
        call AffineTriangle(T, Nu,DNu)
!
!    ...TRACE OF 3D FACE FUNCTIONS
        call checkpolyorder(Nord)
        ndofF = (Nord+1)*Nord/2
        if (ndofF.gt.0) then
!      ...local parameters
          minI  = 0
          minJ  = 0
          minIJ = minI+minJ
          maxIJ = Nord-1
          maxI  = maxIJ-minJ
          maxJ  = maxIJ-minI
          minalpha = 2*minI+1
!      ...orient
          call OrientTri(Nu,DNu,Nori,N, GNu,GDNu)
!      ...construct the shape functions
!      ...get homogenized Legendre polynomials, homP
          call HomLegendre(GNu(0:1),maxI, homP)
!      ...get homogenized Jacobi polynomials, homPal
          call HomJacobi((/GNu(0)+GNu(1),GNu(2)/),maxJ,minalpha, &
                                                               homPal)
          do nij=minIJ,maxIJ
            do i=minI,nij-minJ
              j=nij-i
              m=m+1
!
              call cross(N,GDNu(1:N,1),GDNu(1:N,2), jac)
              ShapV(m) = homP(i)*homPal(i,j)*jac(1)
            enddo
          enddo
        endif
!
!    ...give total degrees of freedom
        NrdofV = m
!
#if HP3D_DEBUG
!    ...print this when debugging
        if (iprint.eq.1) then
          write(*,7001) T(1:2),Nord,Nori
 7001     format('traceTriFshapeV: T = ',2f8.3,/, &
                 'Norder  = ',i2,/, &
                 'Norient = ',i2)
          if (ndofF.gt.0) then
            write(*,*) 'TRACE OF 3D H(div) TRIANGLE FACE FUNCTIONS = '
            do m=1,ndofF
              write(*,7004) m,ShapV(m)
 7004         format('k = ',i3,' ShapV = ',e12.5)
            enddo
          endif
          call pause
        endif
#endif
!
      end subroutine traceTriFshapeV
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Quadrilateral traces
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine traceQuadFshapeV
!
        integer :: i,j,m,N,nordF(2),ndofF
        integer :: minI,minJ,maxI,maxJ
        logical :: IdecF(2),GIdecF(2)
        double precision :: Mu(0:1,2),DMu(2,0:1,2)
        double precision :: GMu(0:1,2),GDMu(2,0:1,2)
        double precision :: homP(0:MAXP-1,2)
        double precision :: jac(1)
!
#if HP3D_DEBUG
!    ...debugging flag
        integer :: iprint
        iprint=0
#endif
!
!    ...spatial dimensions
        N=2
!
!    ...initiate counter for shape functions
        m=0
!
!    ...Define affine coordinates and gradients
        call AffineQuadrilateral(T, Mu,DMu)
!
!    ...TRACE OF 3D FACE FUNCTIONS
        IdecF(1:2) = .true.
        call decod(Nord,MODORDER,2, nordF)
        call checkpolyorder(nordF(1))
        call checkpolyorder(nordF(2))
        ndofF = nordF(1)*nordF(2)
        if (ndofF.gt.0) then
!      ...local parameters (again)
          minI  = 0
          minJ  = 0
          maxI = nordF(1)-1
          maxJ = nordF(2)-1
!      ...orient
          call OrientQuad(Mu,DMu,Nori,IdecF,N, GMu,GDMu,GIdecF)
!      ...construct the shape functions
          call HomLegendre(GMu(0:1,1),maxI, homP(minI:maxI,1))
          call HomLegendre(GMu(0:1,2),maxJ, homP(minJ:maxJ,2))
          do j=minJ,maxJ
            do i=minI,maxI
              m=m+1
!
              call cross(N,GDMu(1:N,1,1),GDMu(1:N,1,2), jac)
              ShapV(m) = homP(i,1)*homP(j,2)*jac(1)
            enddo
          enddo
        endif
!
!    ...give total degrees of freedom
        NrdofV = m
!
#if HP3D_DEBUG
!    ...print this when debugging
        if (iprint.eq.1) then
          write(*,7001) T(1:2),Nord,Nori
 7001     format('traceQuadFshapeV: T = ',2f8.3,/, &
                 'Norder  = ',i2,/, &
                 'Norient = ',i2)
          if (ndofF.gt.0) then
            write(*,*) 'TRACE OF 3D H(div) QUAD FACE FUNCTIONS = '
            do m=1,ndofF
              write(*,7004) m,ShapV(m)
 7004         format('k = ',i3,' ShapV = ',e12.5)
            enddo
          endif
          call pause
        endif
#endif
!
      end subroutine traceQuadFshapeV
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Check polynomial order
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine checkpolyorder(Norder)
!
        integer, intent(in) :: Norder
!
        if ((Norder.lt.1).or.(Norder.gt.MAXP)) then
          write(*,7003) Norder
          write(*,7004) MAXP
 7003     format('TraceFshapeV: Polynomial order = ',i3,' is either ')
 7004     format('              less than 1 or more than MAXP = ',i3)
          stop 1
        endif
!
      end subroutine checkpolyorder
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   end subroutine TraceFshapeV
!

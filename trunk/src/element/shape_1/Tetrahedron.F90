! Routines:
!  - shape3DHTet
!  - shape3DETet
!  - shape3DVTet
!  - shape3DQTet
!--------------------------------------------------------------------
!
!     routine name      - shape3DHTet
!
!--------------------------------------------------------------------
!
!     latest revision:  - Oct 14, Apr 17, Jul 21
!
!     purpose:          - routine returns values of 3D tetrahedron element
!                         H1 shape functions and their derivatives
!
!     arguments:
!
!     in:
!          X            - master tetrahedron coordinates from (0,1)^3
!          Nord         - polynomial order for the nodes (H1 sense)
!          NoriE        - edge orientation
!          NoriF        - face orientation
!          Nsize        - relevant sizes of local arrays
!
!     out:
!          NrdofH       - number of dof
!          ShapH        - values of the shape functions at the point
!          GradH        - gradients of the shape functions
!
!-----------------------------------------------------------------------
!
   subroutine shape3DHTet(X,Nord,NoriE,NoriF,Nsize, &
                                                    NrdofH,ShapH,GradH)
!
      implicit none
      integer, intent(in)  :: Nord(11),NoriE(6),NoriF(4),Nsize(2)
      integer, intent(out) :: NrdofH
      integer :: i,j,k,nij,nijk,m,v,e,f,N,ndofE,nordE,nordF,ndofF
      integer :: minI,minJ,minK,minIJ,minIJK,maxI,maxJ,maxK,maxIJ,maxIJK
      integer :: nordB,ndofB,minbeta
      logical :: IdecE,IdecF,IdecB(2)
      double precision, intent(in)  :: X(3)
      double precision, intent(out) :: ShapH(Nsize(2))
      double precision, intent(out) :: GradH(1:3,Nsize(2))
      double precision :: Lam(0:3),DLam(1:3,0:3)
      double precision :: LambV(4),DLambV(1:3,4)
      double precision :: LampE(0:1,6),DLampE(3,0:1,6)
      double precision :: GLampE(0:1),GDLampE(3,0:1)
      double precision :: LampF(0:2,4),DLampF(3,0:2,4)
      double precision :: GLampF(0:2),GDLampF(3,0:2)
      double precision :: phiE(2:Nsize(1)),DphiE(1:3,2:Nsize(1))
      double precision :: phiTri(2:Nsize(1)-1,1:Nsize(1)-2)
      double precision :: DphiTri(3,2:Nsize(1)-1,1:Nsize(1)-2)
      double precision :: homLbet(3:Nsize(1)-1,1:Nsize(1)-3)
      double precision :: DhomLbet(3,3:Nsize(1)-1,1:Nsize(1)-3)
!
#if HP3D_DEBUG
!  ...debugging flag
      integer :: iprint
      iprint=0
#endif
!
!  ...spatial dimensions
      N=3
!
!  ...local parameters
      minI = 2; minJ = 1; minK = 1
      minIJ = minI+minJ;
      minIJK = minIJ+minK
!
!  ...initiate counter for shape functions
      m=0
!
!  ...Define affine coordinates and gradients
      call AffineTetrahedron(X, Lam,DLam)
!
!  ...VERTEX SHAPE FUNCTIONS
      call BlendTetV(Lam,DLam, LambV,DLambV)
      do v=1,4
       m=m+1
       ShapH(m) = LambV(v)
       GradH(1:N,m) = DLambV(1:N,v)
      enddo
!
!  ...EDGE SHAPE FUNCTIONS
      call ProjectTetE(Lam,DLam, LampE,DLampE,IdecE)
!  ...loop over edges
      do e=1,6
!    ...local parameters
        nordE = Nord(e)
        ndofE = nordE-1
        if (ndofE.gt.0) then
!      ...local parameters (again)
          maxI = nordE
!      ...orient first
          call OrientE(LampE(0:1,e),DLampE(1:N,0:1,e),NoriE(e),N, &
                                                         GLampE,GDLampE)
!      ...construct the shape functions
          call AncPhiE(GLampE,GDLampE,nordE,IdecE,N, &
                                   phiE(minI:maxI),DphiE(1:N,minI:maxI))
          do i=minI,maxI
            m=m+1
            ShapH(m) = phiE(i)
            GradH(1:N,m) = DphiE(1:N,i)
          enddo
        endif
      enddo
!
!  ...FACE SHAPE FUNCTIONS
      call ProjectTetF(Lam,DLam, LampF,DLampF,IdecF)
!  ...loop over faces
      do f=1,4
!    ...local parameters
        nordF = Nord(6+f)
        ndofF = (nordF-1)*(nordF-2)/2
        if (ndofF.gt.0) then
!      ...local parameters (again)
          maxIJ = nordF
          maxI = maxIJ-minJ
          maxJ = maxIJ-minI
!      ...orient first
          call OrientTri(LampF(0:2,f),DLampF(1:N,0:2,f),NoriF(f),N, &
                                                        GLampF,GDLampF)
!      ...construct the shape functions
          call AncPhiTri(GLampF,GDLampF,nordF,IdecF,N, &
                                           phiTri(minI:maxI,minJ:maxJ), &
                                      DphiTri(1:N,minI:maxI,minJ:maxJ))
            do nij=minIJ,maxIJ
              do i=minI,nij-minJ
                j=nij-i
                m=m+1
!
                ShapH(m) = phiTri(i,j)
                GradH(1:N,m) = DphiTri(1:N,i,j)
              enddo
            enddo
        endif
      enddo
!
!  ...BUBBLE FUNCTIONS
!  ...local parameters
      nordB = Nord(11)
      ndofB = (nordB-1)*(nordB-2)*(nordB-3)/6
!  ...if necessary, create bubbles
      if (ndofB.gt.0) then
!    ...local parameters (again)
        IdecB(1) = IdecF; IdecB(2) = .true.
        minbeta = 2*minIJ
        maxIJK = nordB
        maxIJ = maxIJK-minK
        maxI = maxIJ-minJ
        maxJ = maxIJ-minI
        maxK = maxIJK-minIJ
!    ...call phiTri and HomIJacobi - no need to orient
        call AncPhiTri(Lam(0:2),DLam(1:N,0:2),nordB-minK,IdecB(1),N, &
                                           phiTri(minI:maxI,minJ:maxJ), &
                                      DphiTri(1:N,minI:maxI,minJ:maxJ))
        call HomIJacobi((/1-Lam(3),Lam(3)/), &
                 (/-DLam(1:N,3),DLam(1:N,3)/),maxK,minbeta,IdecB(2),N, &
                                        homLbet(minIJ:maxIJ,minK:maxK), &
                                   DhomLbet(1:N,minIJ:maxIJ,minK:maxK))
        do nijk=minIJK,maxIJK
          do nij=minIJ,nijk-minK
            do i=minI,nij-minJ
                j=nij-i
                k=nijk-nij
                m=m+1
!
                ShapH(m) = phiTri(i,j)*homLbet(nij,k)
                GradH(1:N,m) = homLbet(nij,k)*DphiTri(1:N,i,j) &
                             + phiTri(i,j)*DhomLbet(1:N,nij,k)
              enddo
            enddo
          enddo
      endif
!
!  ...give total degrees of freedom
      NrdofH = m
!
#if HP3D_DEBUG
!  ...print this when debugging
      if (iprint.ge.1) then
        write(*,7001) X(1:3),Nord(1:11), &
                      NoriE(1:6),NoriF(1:4)
 7001   format('shape3DHTet: Xi = ',3f8.3,/, &
               'Norder = ',6i2,3x,4i2,3x,i2,/, &
               'orient = ',6i2,3x,4i2)
        write(*,7010)
 7010   format('VERTEX SHAPE FUNCTIONS = ')
        do v=1,4
          m=v
          write(*,7002) m,ShapH(m),GradH(1:3,m)
 7002     format('k = ',i3,' ShapH, GradH = ',e12.5,3x,3e12.5)
        enddo
        do e=1,6
          ndofE = Nord(e)-1
          if (ndofE.gt.0) then
            write(*,7011) e
 7011       format('SHAPE FUNCTIONS FOR EDGE = ',i2)
            do j=1,ndofE
              m=m+1
              write(*,7002) m,ShapH(m),GradH(1:3,m)
            enddo
          endif
        enddo
        do f=1,4
          nordF = Nord(6+f)
          ndofF = (nordF-1)*(nordF-2)/2
          if (ndofF.gt.0) then
            write(*,7012) f
 7012       format('SHAPE FUNCTIONS FOR FACE = ',i2)
            do j=1,ndofF
              m=m+1
              write(*,7002) m,ShapH(m),GradH(1:3,m)
            enddo
          endif
        enddo
        nordB = Nord(11)
        ndofB = (nordB-1)*(nordB-2)*(nordB-3)/6
        if (ndofB.gt.0) then
          write(*,7013)
 7013     format('BUBBLES = ')
          do j=1,ndofB
            m=m+1
            write(*,7002) m,ShapH(m),GradH(1:3,m)
          enddo
        endif
        call pause
      endif
#endif
!
!
   end subroutine shape3DHTet
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape3DETet
!
!--------------------------------------------------------------------
!
!     latest revision:  - Oct 14, Apr 17, Jul 21
!
!     purpose:          - routine returns values of 3D tetrahedron element
!                         H(curl) shape functions and their derivatives
!
!     arguments:
!
!     in:
!          X            - master tetrahedron coordinates from (0,1)^3
!          Nord         - polynomial order for the nodes (H1 sense)
!          NoriE        - edge orientation
!          NoriF        - face orientation
!          Nsize        - relevant sizes of local arrays
!
!     out:
!          NrdofE       - number of dof
!          ShapE        - values of the shape functions at the point
!          CurlE        - cur lof the shape functions
!
!-----------------------------------------------------------------------
!
   subroutine shape3DETet(X,Nord,NoriE,NoriF,Nsize, &
                                                    NrdofE,ShapE,CurlE)
!
      implicit none
      integer, intent(in)  :: Nord(11),NoriE(6),NoriF(4),Nsize(2)
      integer, intent(out) :: NrdofE
      integer :: i,j,k,nij,nijk,m,e,f,N,nordE,ndofE,nordF,ndofF
      integer :: minI,minJ,minK,minIJ,minIJK,maxI,maxJ,maxK,maxIJ,maxIJK
      integer :: nordB,ndofB,minbeta,famctr,fam,abc(3),abcd(4),d
      logical :: IdecE,IdecF,IdecB(2)
      double precision, intent(in)  :: X(3)
      double precision, intent(out) :: ShapE(3,Nsize(2))
      double precision, intent(out) :: CurlE(3,Nsize(2))
      double precision :: Lam(0:3),DLam(3,0:3)
      double precision :: LampE(0:1,6),DLampE(3,0:1,6)
      double precision :: GLampE(0:1),GDLampE(3,0:1)
      double precision :: LampF(0:2,4),DLampF(3,0:2,4)
      double precision :: GLampF(0:2),GDLampF(3,0:2)
      double precision :: EE(3,0:Nsize(1)-1),CurlEE(3,0:Nsize(1)-1)
      double precision :: ETri(3,0:Nsize(1)-2,1:Nsize(1)-1)
      double precision :: CurlETri(3,0:Nsize(1)-2,1:Nsize(1)-1)
      double precision :: homLbet(1:Nsize(1)-1,1:Nsize(1)-1)
      double precision :: DhomLbet(1:3,1:Nsize(1)-1,1:Nsize(1)-1)
      double precision :: DhomLbetxETri(3)
!
#if HP3D_DEBUG
!  ...debugging flag
      integer :: iprint
      iprint=0
#endif
!
!  ...spatial dimension
      N=3
!
!  ...local parameters
      minI = 0; minJ = 1; minK = 1
      minIJ = minI+minJ
      minIJK = minIJ+minK
!
!  ...initiate counter for shape functions
      m=0
!
!  ...Define affine coordinates and gradients
      call AffineTetrahedron(X, Lam,DLam)
!
!  ...EDGE SHAPE FUNCTIONS
      call ProjectTetE(Lam,DLam, LampE,DLampE,IdecE)
!  ...loop over edges
      do e=1,6
!    ...local parameters
        nordE = Nord(e)
        ndofE = nordE
        if (ndofE.gt.0) then
!      ...local parameters (again)
          maxI = nordE-1
!      ...orient
          call OrientE(LampE(0:1,e),DLampE(1:N,0:1,e),NoriE(e),N, &
                                                         GLampE,GDLampE)
!      ...construct the shape functions
          call AncEE(GLampE,GDLampE,nordE,IdecE,N, EE(1:N,minI:maxI), &
                                                  CurlEE(1:N,minI:maxI))
          do i=minI,maxI
            m=m+1
!
            ShapE(1:N,m) = EE(1:N,i)
            CurlE(1:N,m) = CurlEE(1:N,i)
          enddo
        endif
      enddo
!
!  ...FACE SHAPE FUNCTIONS
      call ProjectTetF(Lam,DLam, LampF,DLampF,IdecF)
!
!  ...loop over faces
      do f=1,4
!    ...local parameters
        nordF = Nord(6+f)
        ndofF = nordF*(nordF-1)/2
        if (ndofF.gt.0) then
!      ...local parameters (again)
          maxIJ = nordF-1
          maxI = maxIJ-minJ
          maxJ = maxIJ-minI
!      ...orient
          call OrientTri(LampF(0:2,f),DLampF(1:N,0:2,f),NoriF(f),N, &
                                                         GLampF,GDLampF)
!      ...loop over families
          famctr=m
          do fam=0,1
            m=famctr+fam-1
            abc = cshift((/0,1,2/),fam)
!        ...construct the shape functions
            call AncETri(GLampF(abc),GDLampF(1:N,abc),nordF,IdecF,N, &
                                          ETri(1:N,minI:maxI,minJ:maxJ), &
                                      CurlETri(1:N,minI:maxI,minJ:maxJ))
            do nij=minIJ,maxIJ
              do i=minI,nij-minJ
                j=nij-i
                m=m+2
!
                ShapE(1:N,m) = ETri(1:N,i,j)
                CurlE(1:N,m) = CurlETri(1:N,i,j)
              enddo
            enddo
          enddo
        endif
      enddo
!
!  ...BUBBLE FUNCTIONS
!  ...local parameters
      nordB = Nord(11)
      ndofB = nordB*(nordB-1)*(nordB-2)/6
!  ...if necessary, create bubbles
      if (ndofB.gt.0) then
!    ...local parameters (again)
        IdecB(1) = IdecF; IdecB(2) = .true.
        minbeta = 2*minIJ
        maxIJK = nordB-1
        maxIJ = maxIJK-minK
        maxI = maxIJ-minJ
        maxJ = maxIJ-minI
        maxK = maxIJK-minIJ
!    ...loop over families
        famctr=m
        do fam=0,2
          m=famctr+fam-2
          abcd = cshift((/0,1,2,3/),fam)
          abc = abcd(1:3)
          d = abcd(4)
!      ...now construct the shape functions (no need to orient)
          call AncETri(Lam(abc),DLam(1:N,abc),NordB-minK,IdecB(1),N, &
                                          ETri(1:N,minI:maxI,minJ:maxJ), &
                                      CurlETri(1:N,minI:maxI,minJ:maxJ))
          call HomIJacobi((/1-Lam(d),Lam(d)/), &
                   (/-DLam(1:N,d),DLam(1:N,d)/),maxK,minbeta,IdecB(2),N, &
                                         homLbet(minIJ:maxIJ,minK:maxK), &
                                    DhomLbet(1:N,minIJ:maxIJ,minK:maxK))

          do nijk=minIJK,maxIJK
            do nij=minIJ,nijk-minK
              do i=minI,nij-minJ
                j=nij-i
                k=nijk-nij
                m=m+3
!
                ShapE(1:N,m) = ETri(1:N,i,j)*homLbet(nij,k)
!
                call cross(N,DhomLbet(1:N,nij,k),ETri(1:N,i,j), &
                                                         DhomLbetxETri)
!
                CurlE(1:N,m) = homLbet(nij,k)*CurlETri(1:N,i,j) &
                             + DhomLbetxETri
              enddo
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
      if (iprint.ge.1) then
        write(*,7001) X(1:3),Nord(1:11), &
                      NoriE(1:6),NoriF(1:4)
 7001   format('shape3DETet: Xi = ',3f8.3,/, &
               'Norder = ',6i2,3x,4i2,3x,i2,/, &
               'orient = ',6i2,3x,4i2)
        m=0
        do e=1,6
          ndofE = Nord(e)
          if (ndofE.gt.0) then
            write(*,7011) e
 7011       format('SHAPE FUNCTIONS FOR EDGE = ',i2)
            do j=1,ndofE
              m=m+1
              write(*,7002) m,ShapE(1:N,m),CurlE(1:N,m)
 7002         format('k = ',i3,' ShapE, CurlE = ',3e12.5,3x,3e12.5)
            enddo
          endif
        enddo
        do f=1,4
          nordF = Nord(6+f)
          ndofF = nordF*(nordF-1)/2
          if (ndofF.gt.0) then
            write(*,7012) f
 7012       format('SHAPE FUNCTIONS FOR FACE = ',i2)
            famctr=m
            do fam=0,1
              m=famctr+fam-1
              write(*,7003) fam
 7003         format('family = ',i2)
              do j=1,ndofF
                m=m+2
                write(*,7002) m,ShapE(1:N,m),CurlE(1:N,m)
              enddo
            enddo
          endif
        enddo
        nordB = Nord(11)
        ndofB = nordB*(nordB-1)*(nordB-2)/6
        if (ndofB.gt.0) then
          write(*,7013)
 7013     format('BUBBLES = ')
          famctr=m
          do fam=0,2
            m=famctr-fam-2
            write(*,7003) fam
            do j=1,ndofB
              m=m+3
              write(*,7002) m,ShapE(1:N,m),CurlE(1:N,m)
            enddo
          enddo
        endif
        call pause
      endif
#endif
!
!
   end subroutine shape3DETet
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape3DVTet
!
!--------------------------------------------------------------------
!
!     latest revision:  - Oct 14, Apr 17, Jul 21
!
!     purpose:          - routine returns values of 3D tetrahedron element
!                         H(div) shape functions and their divergences
!
!     arguments:
!
!     in:
!          X            - master tetrahedron coordinates from (0,1)^3
!          Nord         - polynomial order for the nodes (H1 sense)
!          NoriF        - face orientation
!          Nsize        - relevant sizes of local arrays
!
!     out:
!          NrdofV       - number of dof
!          ShapV        - values of the shape functions at the point
!          DivV         - divergence of the shape functions
!
!-----------------------------------------------------------------------
!
   subroutine shape3DVTet(X,Nord,NoriF,Nsize, NrdofV,ShapV,DivV)
!
      implicit none
      integer, intent(in)  :: Nord(11),NoriF(4),Nsize(2)
      integer, intent(out) :: NrdofV
      integer :: i,j,k,nij,nijk,m,f,N,nordF,ndofF
      integer :: minI,minJ,minK,minIJ,minIJK,maxI,maxJ,maxK,maxIJ,maxIJK
      integer :: nordB,ndofB,minbeta,famctr,fam,abc(3),abcd(4),d
      logical :: IdecF,IdecB(2)
      double precision, intent(in)  :: X(3)
      double precision, intent(out) :: ShapV(3,Nsize(2))
      double precision, intent(out) :: DivV(Nsize(2))
      double precision :: Lam(0:3),DLam(3,0:3)
      double precision :: LampF(0:2,4),DLampF(3,0:2,4)
      double precision :: GLampF(0:2),GDLampF(3,0:2)
      double precision :: VTri(3,0:Nsize(1)-1,0:Nsize(1)-1)
      double precision :: DivVTri(0:Nsize(1)-1,0:Nsize(1)-1)
      double precision :: homLbet(0:Nsize(1)-2,1:Nsize(1)-1)
      double precision :: DhomLbet(3,0:Nsize(1)-2,1:Nsize(1)-1)
      double precision :: DhomLbetVTri
!
#if HP3D_DEBUG
!  ...debugging flag
      integer :: iprint
      iprint=0
#endif
!
!  ...spatial dimensions
      N=3
!
!  ...local parameters
      minI = 0; minJ = 0; minK = 1
      minIJ = minI+minJ;
      minIJK = minIJ+minK
!  ...initiate counter for shape functions
      m=0
!
!  ...Define affine coordinates and gradients
      call AffineTetrahedron(X, Lam,DLam)
!
!  ...FACE SHAPE FUNCTIONS
      call ProjectTetF(Lam,DLam, LampF,DLampF,IdecF)
      do f=1,4
!    ...local parameters
        nordF = Nord(6+f)
        ndofF = (nordF+1)*nordF/2
        if (ndofF.gt.0) then
!      ...local parameters (again)
          maxIJ = nordF-1
          maxI = maxIJ-minJ
          maxJ = maxIJ-minI
!      ...orient
          call OrientTri(LampF(0:2,f),DLampF(1:N,0:2,f),NoriF(f),N, &
                                                         GLampF,GDLampF)
!      ...construct the shape functions
          call AncVTri(GLampF,GDLampF,nordF,IdecF,N, &
                                          VTri(1:N,minI:maxI,minJ:maxJ), &
                                           DivVTri(minI:maxI,minJ:maxJ))
          do nij=minIJ,maxIJ
            do i=minI,nij-minJ
              j=nij-i
              m=m+1
!
              ShapV(1:N,m) = VTri(1:N,i,j)
              DivV(m) = DivVTri(i,j)
            enddo
          enddo
        endif
      enddo
!
!  ...BUBBLE FUNCTIONS
!  ...local parameters
      nordB = Nord(11)
      ndofB = (nordB+1)*nordB*(nordB-1)/6
!  ...if necessary, create bubbles
      if (ndofB.gt.0) then
!    ...local parameters (again)
        IdecB(1) = IdecF; IdecB(2) = .true.
        minbeta = 2*(minIJ+1)
        maxIJK = nordB-1
        maxIJ = maxIJK-minK
        maxI = maxIJ-minJ
        maxJ = maxIJ-minI
        maxK = maxIJK-minIJ
!    ...loop over families
        famctr=m
        do fam=0,2
          m=famctr+fam-2
          abcd = cshift((/0,1,2,3/),fam)
          abc = abcd(1:3)
          d = abcd(4)
!      ...construct the shape functions (no need to orient)
          call AncVTri(Lam(abc),DLam(1:N,abc),nordB-minK,IdecB(1),N, &
                                          VTri(1:N,minI:maxI,minJ:maxJ), &
                                           DivVTri(minI:maxI,minJ:maxJ))
          call HomIJacobi((/1-Lam(d),Lam(d)/), &
                   (/-DLam(1:N,d),DLam(1:N,d)/),maxK,minbeta,IdecB(2),N, &
                                         homLbet(minIJ:maxIJ,minK:maxK), &
                                    DhomLbet(1:N,minIJ:maxIJ,minK:maxK))
          do nijk=minIJK,maxIJK
            do nij=minIJ,nijk-minK
              do i=minI,nij-minJ
                j=nij-i
                k=nijk-nij
                m=m+3
!
                ShapV(1:N,m) = VTri(1:N,i,j)*homLbet(nij,k)
!
                call dot_product(DhomLbet(1:N,nij,k),VTri(1:N,i,j), &
                                                           DhomLbetVTri)
!
                DivV(m) = homLbet(nij,k)*DivVTri(i,j)+DhomLbetVTri
              enddo
            enddo
          enddo
        enddo
      endif
!
!  ...give total degrees of freedom
      NrdofV = m
!
#if HP3D_DEBUG
!  ...print this when debugging
      if (iprint.ge.1) then
        write(*,7001) X(1:3),Nord(7:11),NoriF(1:4)
 7001   format('shape3DVTet: Xi = ',3f8.3,/, &
               'Norder = ',3(4i2,2x),2i3,2x,4i3,3x,i4,/, &
               'orient = ',3(4i2,2x),2i3,2x,4i3)
        m=0
        do f=1,4
          nordF = Nord(6+f)
          ndofF = (nordF+1)*nordF/2
          if (ndofF.gt.0) then
            write(*,7012) f
 7012       format('SHAPE FUNCTIONS FOR FACE = ',i2)
            do j=1,ndofF
              m=m+1
              write(*,7002) m,ShapV(1:N,m),DivV(m)
 7002         format('k = ',i3,' ShapV, DivV= ',3e12.5,3x,e12.5)
            enddo
          endif
        enddo
        nordB = Nord(11)
        ndofB = (nordB+1)*nordB*(nordB-1)/6
        if (ndofB.gt.0) then
          write(*,7013)
 7013     format('BUBBLES = ')
          famctr=m
          do fam=0,2
            m=famctr+fam-2
            write(*,7003) fam
 7003       format('family ',i2)
            do j=1,ndofB
              m=m+3
              write(*,7002) m,ShapV(1:N,m),DivV(m)
            enddo
          enddo
        endif
        call pause
      endif
#endif
!
!
   end subroutine shape3DVTet
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape3DQTet
!
!--------------------------------------------------------------------
!
!     latest revision:  - Oct 14, Apr 17
!
!     purpose:          - routine returns values of 3D tetrahedron
!                         element L2 shape functions
!
!     arguments:
!
!     in:
!          X            - master tetrahedron coordinates from (0,1)^3
!          Nord         - polynomial order for the nodes (H1 sense)
!          Nsize        - relevant sizes of local arrays
!
!     out:
!          NrdofQ       - number of dof
!          ShapQ        - values of the shape functions at the point
!
!-----------------------------------------------------------------------
!
   subroutine shape3DQTet(X,Nord,Nsize, NrdofQ,ShapQ)
!
      implicit none
      integer, intent(in)  :: Nord(11),Nsize(2)
      integer, intent(out) :: NrdofQ
      integer :: i,j,k,nij,nijk,m,N,nordB,ndofB,minalpha,minbeta
      integer :: minI,minJ,minK,minIJ,minIJK,maxI,maxJ,maxK,maxIJ,maxIJK
      double precision, intent(in)  :: X(3)
      double precision, intent(out) :: ShapQ(Nsize(2))
      double precision :: Lam(0:3),DLam(1:3,0:3)
      double precision :: homP(0:Nsize(1)-1)
      double precision :: homPal(0:Nsize(1)-1,0:Nsize(1)-1)
      double precision :: homPbet(0:Nsize(1)-1,0:Nsize(1)-1)
!
#if HP3D_DEBUG
!  ...debugging flag
      integer :: iprint
      iprint=0
#endif
!
!  ...spatial dimensions
      N=3
!  ...initiate counter for shape functions
      m=0
!
!  ...Define affine coordinates and gradients
      call AffineTetrahedron(X, Lam,DLam)
!
!  ...local parameters
      nordB = Nord(11)
      ndofB = (nordB+2)*(nordB+1)*nordB/6
      minI = 0; minJ = 0; minK = 0
      minIJ = minI+minJ
      minIJK = minIJ+minK
      minalpha = 2*minI+1
      minbeta = 2*(minIJ+1)
      maxIJK = NordB-1
      maxIJ = maxIJK-minK
      maxI = maxIJ-minJ
      maxJ = maxIJ-minI
      maxK = maxIJK-minIJ
!
!  ...get homogenized Legendre polynomials, homP
      call HomLegendre(Lam(0:1),maxI, homP(minI:maxI))
!
!  ...get homogenized Jacobi polynomials, homPal
      call HomJacobi((/Lam(0)+Lam(1),Lam(2)/),maxIJ,minalpha, &
                                            homPal(minI:maxI,minJ:maxJ))
!  ...get homogenized Jacobi polynomials, homPbet
      call HomJacobi((/1-Lam(3),Lam(3)/),maxK,minbeta, &
                                         homPbet(minIJ:maxIJ,minK:maxK))
!
!  ...construct shape functions
      do nijk=minIJK,maxIJK
        do nij=minIJ,nijk-minK
          do i=minI,nij-minJ
            j=nij-i
            k=nijk-nij
            m=m+1
!
            ShapQ(m) = homP(i)*homPal(i,j)*homPbet(nij,k)
          enddo
        enddo
      enddo
!
!  ...give total degrees of freedom
      NrdofQ = m
!
#if HP3D_DEBUG
!  ...print this when debugging
      if (iprint.ge.1) then
        write(*,7001) X(1:3),Nord(11)
 7001   format('shape3DQTet: Xi = ',3f8.3,/, &
               'Norder = ',i2)
        nordB = Nord(11)
        ndofB = (nordB+2)*(nordB+1)*nordB/6
        if (ndofB.gt.0) then
          write(*,7013)
 7013     format('BUBBLES = ')
          m=0
          do j=1,ndofB
            m=m+1
            write(*,7002) m,ShapQ(m)
7002     format('k = ',i3,' ShapQ = ',e12.5)
          enddo
        endif
        call pause
      endif
#endif
!
   end subroutine shape3DQTet

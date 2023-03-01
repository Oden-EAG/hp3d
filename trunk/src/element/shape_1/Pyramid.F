c Routines:
c  - shape3DHPyra
c  - shape3DEPyra
c  - shape3DVPyra
c  - shape3DQPyra
c--------------------------------------------------------------------
c
c     routine name      - shape3DHPyra
c
c--------------------------------------------------------------------
c
c     latest revision:  - Jan 15, Apr 17, Jul 21
c
c     purpose:          - routine returns values of 3D pyramid element
c                         H1 shape functions and their derivatives
c
c     arguments:
c
c     in:
c          Xi           - master pyramid coordinates from (0,1)^3
c          Nord         - polynomial order for the nodes (H1 sense)
c          NoriE        - edge orientation
c          NoriF        - face orientation
c          Nsize        - relevant sizes of local arrays
c
c     out:
c          NrdofH       - number of dof
c          ShapH        - values of the shape functions at the point
c          GradH        - gradients of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3DHPyra(Xi,Nord,NoriE,NoriF,Nsize,
     .                                              NrdofH,ShapH,GradH)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: Nord(14),NoriE(8),NoriF(5),Nsize(2)
      integer, intent(out) :: NrdofH
      integer :: N,m,v,e,f,i,j,k,nij,nordME,ndofME,nordTE,ndofTE
      integer :: nordQF(2),ndofQF,nordTF,ndofTF,nordB,ndofB
      integer :: minI,minJ,minK,minIJ,maxI,maxJ,maxK,maxIJ
      logical :: IdecME,IdecTE,IdecQF(2),GIdecQF(2),IdecTF,IdecB(3)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapH(Nsize(2))
      double precision, intent(out) :: GradH(3,Nsize(2))
      double precision :: Lam(5),DLam(3,5)
      double precision :: Mu(0:1,2),DMu(3,0:1,2)
      double precision :: Nu(0:2,2),DNu(3,0:2,2)
      double precision :: MuZ(0:1),DMuZ(3,0:1)
      double precision :: LambV(5),DLambV(3,5)
      double precision :: MubME(4),DMubME(3,4)
      double precision :: NupME(0:1,4),DNupME(3,0:1,4)
      double precision :: GNupME(0:1),GDNupME(3,0:1)
      double precision :: LampTE(0:1,4),DLampTE(3,0:1,4)
      double precision :: GLampTE(0:1),GDLampTE(3,0:1)
      double precision :: MupQF(0:1,2),DMupQF(3,0:1,2)
      double precision :: GMupQF(0:1,2),GDMupQF(3,0:1,2)
      double precision :: MubTF(4),DMubTF(3,4)
      double precision :: NupTF(0:2,4),DNupTF(3,0:2,4)
      double precision :: GNupTF(0:2),GDNupTF(3,0:2)
      double precision :: phiE(2:Nsize(1)),DphiE(3,2:Nsize(1))
      double precision :: phiQuad(2:Nsize(1),2:Nsize(1))
      double precision :: DphiQuad(3,2:Nsize(1),2:Nsize(1))
      double precision :: phiTri(2:Nsize(1)-1,1:Nsize(1)-2)
      double precision :: DphiTri(3,2:Nsize(1)-1,1:Nsize(1)-2)
c
#if DEBUG_MODE
c  ...debugging flag
      integer :: iprint=0
#endif
c
c  ...spatial dimensions
      N=3
c
c  ...initiate counter for shape functions
      m=0
c
c  ...Define affine coordinates and gradients
      call AffinePyramid(Xi, Lam,DLam,Mu,DMu,Nu,DNu,MuZ,DMuZ)
c
c  ...VERTEX SHAPE FUNCTIONS
      call BlendPyraV(Lam,DLam, LambV,DLambV)
      do v=1,5
       m=m+1
       ShapH(m)     = LambV(v)
       GradH(1:N,m) = DLambV(1:N,v)
      enddo
c
c  ...MIXED EDGE SHAPE FUNCTIONS
      call BlendProjectPyraME(Mu,DMu,Nu,DNu,
     .                                 MubME,DMubME,NupME,DNupME,IdecME)
c  ...loop over edges
      do e=1,4
c    ...local parameters
        nordME = Nord(e)
        ndofME = nordME-1
        if (ndofME.gt.0) then
c      ...local parameters (again)
          minI = 2
          maxI = nordME
c      ...orient first
          call OrientE(NupME(0:1,e),DNupME(1:N,0:1,e),NoriE(e),N,
     .                                                   GNupME,GDNupME)
c      ...construct the shape functions
          call AncPhiE(GNupME,GDNupME,nordME,IdecME,N,
     .                             phiE(minI:maxI),DphiE(1:N,minI:maxI))
          do i=minI,maxI
            m=m+1
            ShapH(m)     = MubME(e)*phiE(i)
            GradH(1:N,m) = DMubME(1:N,e)*phiE(i)
     .                   + MubME(e)*DphiE(1:N,i)
          enddo
        endif
      enddo
c
c  ...TRIANGLE EDGE SHAPE FUNCTIONS
      call ProjectPyraTE(Lam,DLam, LampTE,DLampTE,IdecTE)
c  ...loop over edges
      do e=1,4
c    ...local parameters
        nordTE = Nord(e+4)
        ndofTE = nordTE-1
        if (ndofTE.gt.0) then
c      ...local parameters (again)
          minI = 2
          maxI = nordTE
c      ...orient first
          call OrientE(LampTE(0:1,e),DLampTE(1:N,0:1,e),NoriE(e+4),N,
     .                                                 GLampTE,GDLampTE)
c      ...construct the shape functions
          call AncPhiE(GLampTE,GDLampTE,nordTE,IdecTE,N,
     .                             phiE(minI:maxI),DphiE(1:N,minI:maxI))
          do i=minI,maxI
            m=m+1
            ShapH(m)     = phiE(i)
            GradH(1:N,m) = DphiE(1:N,i)
          enddo
        endif
      enddo
c
c  ...QUADRILATERAL FACE SHAPE FUNCTION
      call ProjectPyraQF(Mu,DMu, MupQF,DMupQF,IdecQF)
c  ...local parameters
      call decod(Nord(9),MODORDER,2, nordQF)
      ndofQF = (nordQF(1)-1)*(nordQF(2)-1)
      if (ndofQF.gt.0) then
c    ...local parameters (again)
        minI = 2
        minJ = 2
        maxI = nordQF(1)
        maxJ = nordQF(2)
c    ...orient first
        call OrientQuad(MupQF,DMupQF,NoriF(1),IdecQF,N,
     .                                           GMupQF,GDMupQF,GIdecQF)
c    ...construct the shape functions
        call AncPhiQuad(GMupQF,GDMupQF,nordQF,GIdecQF,N,
     .                                     phiQuad(minI:maxI,minJ:maxJ),
     .                                DphiQuad(1:N,minI:maxI,minJ:maxJ))
          do j=minJ,maxJ
            do i=minI,maxI
              m=m+1
c
              ShapH(m)     = MuZ(0)*phiQuad(i,j)
              GradH(1:N,m) = DMuZ(1:N,0)*phiQuad(i,j)
     .                     + MuZ(0)*DphiQuad(1:N,i,j)
            enddo
          enddo
        endif
c
c  ...TRIANGULAR FACE SHAPE FUNCTIONS
      call BlendProjectPyraTF(Mu,DMu,Nu,DNu,
     .                                 MubTF,DMubTF,NupTF,DNupTF,IdecTF)
c  ...loop over faces
      do f=1,4
c    ...local parameters
        nordTF = Nord(9+f)
        ndofTF = (nordTF-1)*(nordTF-2)/2
        if (ndofTF.gt.0) then
c      ...local parameters (again)
          minI = 2
          minJ = 1
          minIJ = minI+minJ
          maxIJ = nordTF
          maxI = maxIJ-minJ
          maxJ = maxIJ-minI
c      ...orient first
          call OrientTri(NupTF(0:2,f),DNupTF(1:N,0:2,f),NoriF(f+1),N,
     .                                                   GNupTF,GDNupTF)
c      ...construct the shape functions
          call AncPhiTri(GNupTF,GDNupTF,nordTF,IdecTF,N,
     .                                      phiTri(minI:maxI,minJ:maxJ),
     .                                 DphiTri(1:N,minI:maxI,minJ:maxJ))
            do nij=minIJ,maxIJ
              do i=minI,nij-minJ
                j=nij-i
                m=m+1
c
                ShapH(m)     = MubTF(f)*phiTri(i,j)
                GradH(1:N,m) = DMubTF(1:N,f)*phiTri(i,j)
     .                       + MubTF(f)*DphiTri(1:N,i,j)
              enddo
            enddo
        endif
      enddo
c
c  ...BUBBLE FUNCTIONS
c  ...local parameters
      nordB = Nord(14)
      ndofB = (nordB-1)**3
c  ...if necessary, create bubbles
      if (ndofB.gt.0) then
c    ...local parameters (again)
        IdecB(1:2) = IdecQF; IdecB(3) = .TRUE.
        minI = 2
        minJ = 2
        minK = 2
        maxI = nordB
        maxJ = nordB
        maxK = nordB
c    ...call phiQuad and phiE - no need to orient
        call AncPhiQuad(Mu,DMu,(/nordB,nordB/),IdecB(1:2),N,
     .                                     phiQuad(minI:maxI,minJ:maxJ),
     .                                DphiQuad(1:N,minI:maxI,minJ:maxJ))
        call AncPhiE(MuZ,DMuZ,nordB,IdecB(3),N,
     .                             phiE(minK:maxK),DphiE(1:N,minK:maxK))
        do k=minK,maxK
          do j=minJ,maxJ
            do i=minI,maxI
              m=m+1
c
              ShapH(m) = phiQuad(i,j)*phiE(k)
              GradH(1:N,m) = phiQuad(i,j)*DphiE(1:N,k)
     .                     + DphiQuad(1:N,i,j)*phiE(k)
            enddo
          enddo
        enddo
      endif
c
c  ...give total degrees of freedom
      NrdofH = m
c
#if DEBUG_MODE
c  ...print this when debugging
      if (iprint.ge.1) then
        write(*,7001) Xi(1:3),Nord(1:14),NoriE(1:8),NoriF(1:5)
 7001   format('shape3DHPyra: Xi = ',3f8.3,/,
     .         'Norder = ',8i2,3x,5i2,3x,i2,/,
     .         'orient = ',8i2,3x,5i2)
        write(*,7010)
 7010   format('VERTEX SHAPE FUNCTIONS = ')
        do v=1,5
          m=v
          write(*,7002) m,ShapH(m),GradH(1:3,m)
 7002     format('k = ',i3,' ShapH, GradH = ',e12.5,3x,3e12.5)
        enddo
        do e=1,4
          ndofME = Nord(e)-1
          if (ndofME.gt.0) then
            write(*,7011) e
 7011       format('SHAPE FUNCTIONS FOR MIXED EDGE = ',i2)
            do j=1,ndofME
              m=m+1
              write(*,7002) m,ShapH(m),GradH(1:3,m)
            enddo
          endif
        enddo
        do e=1,4
          ndofTE = Nord(e)-1
          if (ndofTE.gt.0) then
            write(*,7012) e
 7012       format('SHAPE FUNCTIONS FOR TRIANGLE EDGE = ',i2)
            do j=1,ndofTE
              m=m+1
              write(*,7002) m,ShapH(m),GradH(1:3,m)
            enddo
          endif
        enddo
        call decod(Nord(9),MODORDER,2, nordQF)
        ndofQF = (nordQF(1)-1)*(nordQF(2)-1)
        if (ndofQF.gt.0) then
          write(*,*) 'SHAPE FUNCTIONS FOR QUADRILATERAL FACE'
          do j=1,ndofQF
            m=m+1
            write(*,7002) m,ShapH(m),GradH(1:3,m)
          enddo
        endif
        do f=1,4
          nordTF = Nord(9+f)
          ndofTF = (nordTF-1)*(nordTF-2)/2
          if (ndofTF.gt.0) then
            write(*,7013) f
 7013       format('SHAPE FUNCTIONS FOR TRIANGLE FACE = ',i2)
            do j=1,ndofTF
              m=m+1
              write(*,7002) m,ShapH(m),GradH(1:3,m)
            enddo
          endif
        enddo
        nordB = Nord(14)
        ndofB = (nordB-1)**3
        if (ndofB.gt.0) then
          write(*,7015)
 7015     format('BUBBLES = ')
          do j=1,ndofB
            m=m+1
            write(*,7002) m,ShapH(m),GradH(1:3,m)
          enddo
        endif
        call pause
      endif
#endif
c
c
      end subroutine shape3DHPyra
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape3DEPyra
c
c--------------------------------------------------------------------
c
c     latest revision:  - Jan 15, Apr 17, Jul 21
c
c     purpose:          - routine returns values of 3D pyramid element
c                         H(curl) shape functions and their curls
c
c     arguments:
c
c     in:
c          Xi           - master pyramid coordinates from (0,1)^3
c          Nord         - polynomial order for the nodes (H1 sense)
c          NoriE        - edge orientation
c          NoriF        - face orientation
c          Nsize        - relevant sizes of local arrays
c
c     out:
c          NrdofE       - number of dof
c          ShapE        - values of the shape functions at the point
c          CurlE        - curl of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3DEPyra(Xi,Nord,NoriE,NoriF,Nsize,
     .                                              NrdofE,ShapE,CurlE)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: Nord(14),NoriE(8),NoriF(5),Nsize(2)
      integer, intent(out) :: NrdofE
      integer :: N,m,e,f,i,j,k,nij,nordME,ndofME,nordTE,ndofTE
      integer :: nordQF(2),ndofQF,nordTF,ndofTF,nordB,ndofB,ab(2),a,b
      integer :: minI,minJ,minK,minIJ,maxI,maxJ,maxK,maxIJ,fam,famctr
      integer :: ij(2),ig,jg,minF(2),maxF(2),abc(3),ijmax
      logical :: IdecME,IdecTE,IdecQF(2),GIdecQF(2),IdecTF,IdecB(3)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapE(3,Nsize(2))
      double precision, intent(out) :: CurlE(3,Nsize(2))
      double precision :: Lam(5),DLam(3,5)
      double precision :: Mu(0:1,2),DMu(3,0:1,2)
      double precision :: Nu(0:2,2),DNu(3,0:2,2)
      double precision :: MuZ(0:1),DMuZ(3,0:1)
      double precision :: MubME(4),DMubME(3,4)
      double precision :: NupME(0:1,4),DNupME(3,0:1,4)
      double precision :: GNupME(0:1),GDNupME(3,0:1)
      double precision :: LampTE(0:1,4),DLampTE(3,0:1,4)
      double precision :: GLampTE(0:1),GDLampTE(3,0:1)
      double precision :: MupQF(0:1,2),DMupQF(3,0:1,2)
      double precision :: GMupQF(0:1,2),GDMupQF(3,0:1,2)
      double precision :: MubTF(4),DMubTF(3,4)
      double precision :: NupTF(0:2,4),DNupTF(3,0:2,4)
      double precision :: GNupTF(0:2),GDNupTF(3,0:2)
      double precision :: EE(3,0:Nsize(1)-1),CurlEE(3,0:Nsize(1)-1)
      double precision :: EQuad(3,0:Nsize(1)-1,2:Nsize(1))
      double precision :: CurlEQuad(3,0:Nsize(1)-1,2:Nsize(1))
      double precision :: ETri(3,0:Nsize(1)-2,1:Nsize(1)-1)
      double precision :: CurlETri(3,0:Nsize(1)-2,1:Nsize(1)-1)
      double precision :: phiE(2:Nsize(1)),DphiE(3,2:Nsize(1))
      double precision :: phiQuad(2:Nsize(1),2:Nsize(1))
      double precision :: DphiQuad(3,2:Nsize(1),2:Nsize(1))
      double precision :: DMubMExEE(3),DMuZxEQuad(3),DMubTFxETri(3)
      double precision :: tmp,vectmp(3),DphiQuadxDMuZ(3)
c
#if DEBUG_MODE
c  ...debugging flag
      integer :: iprint=0
#endif
c
c  ...spatial dimensions
      N=3
c
c  ...initiate counter for shape functions
      m=0
c
c  ...Define affine coordinates and gradients
      call AffinePyramid(Xi, Lam,DLam,Mu,DMu,Nu,DNu,MuZ,DMuZ)
c
c  ...MIXED EDGE SHAPE FUNCTIONS
      call BlendProjectPyraME(Mu,DMu,Nu,DNu,
     .                                 MubME,DMubME,NupME,DNupME,IdecME)
c  ...loop over edges
      do e=1,4
c    ...local parameters
        nordME = Nord(e)
        ndofME = nordME
        if (ndofME.gt.0) then
c      ...local parameters (again)
          minI = 0
          maxI = nordME-1
c      ...orient first
          call OrientE(NupME(0:1,e),DNupME(1:N,0:1,e),NoriE(e),N,
     .                                                   GNupME,GDNupME)
c      ...construct the shape functions
          call AncEE(GNupME,GDNupME,nordME,IdecME,N,
     .                          EE(1:N,minI:maxI),CurlEE(1:N,minI:maxI))
          do i=minI,maxI
            m=m+1
c
            ShapE(1:N,m) = MubME(e)*EE(1:N,i)
            call cross(3,DMubME(1:N,e),EE(1:N,i), DMubMExEE)
            CurlE(1:N,m) = MubME(e)*CurlEE(1:N,i)
     .                   + DMubMExEE
          enddo
        endif
      enddo
c
c  ...TRIANGLE EDGE SHAPE FUNCTIONS
      call ProjectPyraTE(Lam,DLam, LampTE,DLampTE,IdecTE)
c  ...loop over edges
      do e=1,4
c    ...local parameters
        nordTE = Nord(e+4)
        ndofTE = nordTE
        if (ndofTE.gt.0) then
c      ...local parameters (again)
          minI = 0
          maxI = nordTE-1
c      ...orient first
          call OrientE(LampTE(0:1,e),DLampTE(1:N,0:1,e),NoriE(e+4),N,
     .                                                 GLampTE,GDLampTE)
c      ...construct the shape functions
          call AncEE(GLampTE,GDLampTE,nordTE,IdecTE,N,
     .                          EE(1:N,minI:maxI),CurlEE(1:N,minI:maxI))
          do i=minI,maxI
            m=m+1
c
            ShapE(1:N,m) = EE(1:N,i)
            CurlE(1:N,m) = CurlEE(1:N,i)
          enddo
        endif
      enddo
c
c  ...QUADRILATERAL FACE SHAPE FUNCTIONS
      call ProjectPyraQF(Mu,DMu, MupQF,DMupQF,IdecQF)
c  ...local parameters
      call decod(Nord(9),MODORDER,2, nordQF)
c  ...orient
      call OrientQuad(MupQF,DMupQF,NoriF(1),IdecQF,N,
     .                                           GMupQF,GDMupQF,GIdecQF)
c  ...loop over families
      do fam=0,1
        ab = cshift((/1,2/),fam);
        a = ab(1); b = ab(2)
        ndofQF = nordQF(a)*(nordQF(b)-1)
        if (ndofQF.gt.0) then
c      ...local parameters (again)
          minF(1) = 0
          minF(2) = 2
          maxF(1) = nordQF(a)-1
          maxF(2) = nordQF(b)
c      ...construct the shape functions
          call AncEQuad(GMupQF(0:1,ab),GDMupQF(1:N,0:1,ab),nordQF(ab),
     .                                                  GIdecQF(ab),N,
     .                       EQuad(1:N,minF(1):maxF(1),minF(2):maxF(2)),
     .                   CurlEQuad(1:N,minF(1):maxF(1),minF(2):maxF(2)))
c      ...in the code the outer loop always is
c      ...numbered wrt the second global face axis
          minF = cshift(minF,-fam); maxF = cshift(maxF,-fam)
          do jg=minF(2),maxF(2)
            do ig=minF(1),maxF(1)
              ij = cshift((/ig,jg/),fam)
              i = ij(1); j = ij(2)
              m=m+1
c
              ShapE(1:N,m) = EQuad(1:N,i,j)*MuZ(0)**2
              call cross(N,DMuZ(1:N,0),EQuad(1:N,i,j), DMuZxEQuad)
              CurlE(1:N,m) = CurlEQuad(1:N,i,j)*MuZ(0)**2
     .                     + 2*DMuZxEQuad*MuZ(0)
            enddo
          enddo
        endif
      enddo
c
c  ...TRIANGULAR FACE SHAPE FUNCTIONS
      call BlendProjectPyraTF(Mu,DMu,Nu,DNu,
     .                                 MubTF,DMubTF,NupTF,DNupTF,IdecTF)
c  ...loop over faces
      do f=1,4
c    ...local parameters
        nordTF = Nord(9+f)
        ndofTF = nordTF*(nordTF-1)/2
        if (ndofTF.gt.0) then
c      ...local parameters (again)
          minI  = 0
          minJ  = 1
          minIJ = minI+minJ
          maxIJ = nordTF-1
          maxI  = maxIJ-minJ
          maxJ  = maxIJ-minI
c      ...orient first
          call OrientTri(NupTF(0:2,f),DNupTF(1:N,0:2,f),NoriF(f+1),N,
     .                                                   GNupTF,GDNupTF)
c      ...loop over families
          famctr=m
          do fam=0,1
            m=famctr+fam-1
            abc = cshift((/0,1,2/),fam)
c        ...construct the shape functions
            call AncETri(GNupTF(abc),GDNupTF(1:N,abc),nordTF,IdecTF,N,
     .                                    ETri(1:N,minI:maxI,minJ:maxJ),
     .                                CurlETri(1:N,minI:maxI,minJ:maxJ))
            do nij=minIJ,maxIJ
              do i=minI,nij-minJ
                j=nij-i
                m=m+2
c
                ShapE(1:N,m) = ETri(1:N,i,j)*MubTF(f)
                call cross(N,DMubTF(1:N,f),ETri(1:N,i,j), DMubTFxETri)
                CurlE(1:N,m) = MubTF(f)*CurlETri(1:N,i,j)
     .                       + DMubTFxETri
              enddo
            enddo
          enddo
        endif
      enddo
c
c  ...BUBBLE FUNCTIONS
c  ...local parameters
      nordB = Nord(14)
c  ...FAMILY 1 (gradients of H1 bubbles)
      ndofB = (nordB-1)**3
c  ...if necessary, create bubbles
      if (ndofB.gt.0) then
c    ...local parameters (again)
        IdecB(1:2) = .TRUE.; IdecB(3) = .TRUE.
        minI = 2
        minJ = 2
        minK = 2
        maxI = nordB
        maxJ = nordB
        maxK = nordB
c    ...construct shape functions, no need to orient
        call AncPhiQuad(Mu,DMu,(/nordB,nordB/),IdecB(1:2),N,
     .                                     phiQuad(minI:maxI,minJ:maxJ),
     .                                DphiQuad(1:N,minI:maxI,minJ:maxJ))
        call AncPhiE(MuZ,DMuZ,nordB,IdecB(3),N,
     .                             phiE(minK:maxK),DphiE(1:N,minK:maxK))
        do k=minK,maxK
          do j=minJ,maxJ
            do i=minI,maxI
              m=m+1
c
              ShapE(1:N,m) = phiQuad(i,j)*DphiE(1:N,k)
     .                     + DphiQuad(1:N,i,j)*phiE(k)
              CurlE(1:N,m) = 0.d0
            enddo
          enddo
        enddo
      endif
c  ...FAMILY 2 AND 3 (induced from quad face functions)
      ndofB = nordB*(nordB-1)**2
c  ...if necessary, create bubbles
      if (ndofB.gt.0) then
c    ...local parameters
        IdecB(1:2) = .TRUE.; IdecB(3) = .TRUE.
        minF(1) = 0
        minF(2) = 2
        minK    = 2
        maxF(1) = nordB-1
        maxF(2) = nordB
        maxK    = nordB
c    ...loop over families
        do fam=0,1
          ab = cshift((/1,2/),fam);
          a = ab(1); b = ab(2)
c    ...construct shape functions, no need to orient
          call AncEQuad(Mu(0:1,ab),DMu(1:N,0:1,ab),(/nordB,nordB/),
     .                                                    IdecB(ab),N,
     .                       EQuad(1:N,minF(1):maxF(1),minF(2):maxF(2)),
     .                   CurlEQuad(1:N,minF(1):maxF(1),minF(2):maxF(2)))
          call AncPhiE(MuZ,DMuZ,nordB,IdecB(3),N,
     .                             phiE(minK:maxK),DphiE(1:N,minK:maxK))
c      ...in the code the outer loop always is
c      ...numbered wrt the second global face axis
          minF = cshift(minF,-fam); maxF = cshift(maxF,-fam)
          do k=minK,maxK
            do jg=minF(2),maxF(2)
              do ig=minF(1),maxF(1)
                ij = cshift((/ig,jg/),fam)
                i = ij(1); j = ij(2)
                m=m+1
c
                ShapE(1:N,m) = MuZ(0)*EQuad(1:N,i,j)*phiE(k)
                call cross(3,DMuZ(1:N,0)*phiE(k) + MuZ(0)*DphiE(1:N,k),
     .                                           EQuad(1:N,i,j), vectmp)
                CurlE(1:N,m) = MuZ(0)*CurlEQuad(1:N,i,j)*phiE(k)
     .                       + vectmp
              enddo
            enddo
          enddo
        enddo
      endif
c  ...FAMILY 4
      ndofB = (nordB-1)**2
c  ...if necessary, create bubbles
      if (ndofB.gt.0) then
c    ...local parameters (again)
        IdecB(1:2) = .TRUE.
        minI = 2
        minJ = 2
        maxI = nordB
        maxJ = nordB
c    ...construct shape functions, no need to orient
        call AncPhiQuad(Mu,DMu,(/nordB,nordB/),IdecB(1:2),N,
     .                                     phiQuad(minI:maxI,minJ:maxJ),
     .                                DphiQuad(1:N,minI:maxI,minJ:maxJ))
        do j=minJ,maxJ
          do i=minI,maxI
            m=m+1
            ijmax = max(i,j)
c
            tmp = phiQuad(i,j)*MuZ(0)**(ijmax-1)
            ShapE(1:N,m) = tmp*DMuZ(1:N,1)
            call cross(N,DphiQuad(1:N,i,j),DMuZ(1:N,1), DphiQuadxDMuZ)
            CurlE(1:N,m) = DphiQuadxDMuZ*MuZ(0)**(ijmax-1)
          enddo
        enddo
      endif

c
c  ...give total degrees of freedom
      NrdofE = m
c
#if DEBUG_MODE
c  ...print this when debugging
      if (iprint.ge.1) then
        write(*,7001)Xi(1:3),Nord(1:14),NoriE(1:8),NoriF(1:5),NrdofE
 7001   format('shape3DEPyra: Xi = ',3f8.3,/,
     .         'Norder = ',4i3,1x,4i3,3x,i3,1x,4i3,3x,i3,/,
     .         'orient = ',4i3,1x,4i3,3x,i3,1x,4i3,/,
     .         'NrdofE = ',2i3)
        m=0
        do e=1,4
          nordME = Nord(e)
          ndofME = nordME
          if (ndofME.gt.0) then
            write(*,7011) e
 7011       format('SHAPE FUNCTIONS FOR MIXED EDGE = ',i2)
            do j=1,ndofME
              m=m+1
              write(*,7002) m,ShapE(1:N,m),CurlE(1:N,m)
 7002         format('k = ',i3,' ShapE, CurlE = ',3e12.5,3x,3e12.5)
            enddo
          endif
        enddo
        do e=1,4
          nordTE = Nord(4+e)
          ndofTE = nordTE
          if (ndofTE.gt.0) then
            write(*,7012) e
 7012       format('SHAPE FUNCTIONS FOR TRIANGULAR EDGE = ',i2)
            do j=1,ndofTE
              m=m+1
              write(*,7003) m,ShapE(1:N,m),CurlE(1:N,m)
 7003         format('k = ',i3,' ShapE, CurlE = ',3e12.5,3x,3e12.5)
            enddo
          endif
        enddo
        call decod(Nord(9),MODORDER,2, nordQF)
        ndofQF = nordQF(1)*(nordQF(2)-1)
        if (ndofQF.gt.0) then
          write(*,*) 'SHAPE FUNCTIONS FOR QUADRILATERAL FACE'
          do fam=0,1
            write(*,7004) fam
 7004       format('Family ',i2)
            do j=1,ndofQF
              m=m+1
              write(*,7002) m,ShapE(1:N,m),CurlE(1:N,m)
            enddo
          enddo
        endif
        do f=1,4
          nordTF = Nord(9+f)
          ndofTF = nordTF*(nordTF-1)/2
          if (ndofTF.gt.0) then
            write(*,7014) f
 7014       format('SHAPE FUNCTIONS FOR TRIANGLAR FACE = ',i2)
            famctr=m
            do fam=0,1
              m=famctr+fam-1
              write(*,7004) fam
              do j=1,ndofTF
                m=m+2
                write(*,7002) m,ShapE(1:N,m),CurlE(1:N,m)
              enddo
            enddo
          endif
        enddo
        nordB = Nord(14)
        ndofB = (nordB-1)**3
        if (ndofB.gt.0) then
          write(*,*) 'SHAPE FUNCTIONS FOR BUBBLES'
          write(*,*) 'Family 1'
          do j=1,ndofB
            m=m+1
            write(*,7002) m,ShapE(1:N,m),CurlE(1:N,m)
          enddo
        endif
        ndofB = nordB*(nordB-1)**2
        if (ndofB.gt.0) then
          write(*,*) 'SHAPE FUNCTIONS FOR QUAD-TYPE BUBBLES'
          do fam=0,1
            write(*,7004) fam+2
            do j=1,ndofB
              m=m+1
              write(*,7002) m,ShapE(1:N,m),CurlE(1:N,m)
            enddo
          enddo
        endif
        ndofB = (nordB-1)**2
        if (ndofB.gt.0) then
          write(*,*) 'Family 4'
          do j=1,ndofB
            m=m+1
            write(*,7002) m,ShapE(1:N,m),CurlE(1:N,m)
          enddo
        endif
        call pause
      endif
#endif
c
c
      end subroutine shape3DEPyra
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape3DVPyra
c
c--------------------------------------------------------------------
c
c     latest revision:  - Jan 15, Apr 17, Jul 21
c
c     purpose:          - routine returns values of 3D pyramid element
c                         H(div) shape functions and their divergences
c
c     arguments:
c
c     in:
c          Xi           - master pyramid coordinates from (0,1)^3
c          Nord         - polynomial order for the nodes (H1 sense)
c          NoriF        - face orientation
c          Nsize        - relevant sizes of local arrays
c
c     out:
c          NrdofV       - number of dof
c          ShapV        - values of the shape functions at the point
c          DivV         - divergence of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3DVPyra(Xi,Nord,NoriF,Nsize, NrdofV,ShapV,DivV)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: Nord(14),NoriF(5),Nsize(2)
      integer, intent(out) :: NrdofV
      integer :: N,m,f,i,j,k,nij
      integer :: nordQF(2),ndofQF,nordTF,ndofTF,nordB,ndofB,ab(2),a,b
      integer :: minI,minJ,minK,minIJ,maxI,maxJ,maxK,maxIJ,fam
      integer :: ij(2),ig,jg,minF(2),maxF(2),ijmax
      logical :: IdecQF(2),GIdecQF(2),IdecTF,IdecLamTF
      logical :: IdecB(3)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapV(3,Nsize(2))
      double precision, intent(out) :: DivV(Nsize(2))
      double precision :: Lam(5),DLam(3,5)
      double precision :: Mu(0:1,2),DMu(3,0:1,2)
      double precision :: Nu(0:2,2),DNu(3,0:2,2)
      double precision :: MuZ(0:1),DMuZ(3,0:1)
      double precision :: MupQF(0:1,2),DMupQF(3,0:1,2)
      double precision :: GMupQF(0:1,2),GDMupQF(3,0:1,2)
      double precision :: MubTF(4),DMubTF(3,4)
      double precision :: NupTF(0:2,4),DNupTF(3,0:2,4)
      double precision :: LampTF(0:2,4),DLampTF(3,0:2,4)
      double precision :: GNupTF(0:2),GDNupTF(3,0:2)
      double precision :: GLampTF(0:2),GDLampTF(3,0:2)
      double precision :: VQuad(3,0:Nsize(1)-1,0:Nsize(1)-1)
      double precision :: DivVQuad(0:Nsize(1)-1,0:Nsize(1)-1)
      double precision :: VTri(3,0:Nsize(1)-1,0:Nsize(1)-1)
      double precision :: DivVTri(0:Nsize(1)-1,0:Nsize(1)-1)
      double precision :: VTri2(3,0:Nsize(1)-1,0:Nsize(1)-1)
      double precision :: DivVTri2(0:Nsize(1)-1,0:Nsize(1)-1)
      double precision :: EQuad(3,0:Nsize(1)-1,2:Nsize(1))
      double precision :: CurlEQuad(3,0:Nsize(1)-1,2:Nsize(1))
      double precision :: phiE(2:Nsize(1)),DphiE(3,2:Nsize(1))
      double precision :: phiE2(2:Nsize(1)),DphiE2(3,2:Nsize(1))
      double precision :: phiQuad(2:Nsize(1),2:Nsize(1))
      double precision :: DphiQuad(3,2:Nsize(1),2:Nsize(1))
      double precision :: tmp,tmp1,tmp2,tmp3
      double precision :: vectmp(3),vectmp1(3),vectmp2(3)
      double precision :: DMuZVQuad,DphiQuadxDMuZ(3)
c
#if DEBUG_MODE
c  ...debugging flag
      integer :: iprint=0
#endif
c
c  ...spatial dimensions
      N=3
c
c  ...initiate counter for shape functions
      m=0
c
c  ...Define affine coordinates and gradients
      call AffinePyramid(Xi, Lam,DLam,Mu,DMu,Nu,DNu,MuZ,DMuZ)
c
c  ...QUADRILATERAL FACE SHAPE FUNCTIONS
      call ProjectPyraQF(Mu,DMu, MupQF,DMupQF,IdecQF)
c  ...local parameters
      call decod(Nord(9),MODORDER,2, nordQF)
      ndofQF = nordQF(1)*nordQF(2)
      if (ndofQF.gt.0) then
        minI = 0
        minJ = 0
        maxI = NordQF(1)-1
        maxJ = NordQF(2)-1
c    ...orient
        call OrientQuad(MupQF,DMupQF,NoriF(1),IdecQF,N,
     .                                           GMupQF,GDMupQF,GIdecQF)
c    ...construct the shape functions
        call AncVQuad(GMupQF,GDMupQF,nordQF,GIdecQF,N,
     .                                   VQuad(1:N,minI:maxI,minJ:maxJ),
     .                                    DivVQuad(minI:maxI,minJ:maxJ))
        do j=minJ,maxJ
          do i=minI,maxI
            m=m+1
c
            ShapV(1:N,m) = MuZ(0)**3*VQuad(1:N,i,j)
            call dot_product(DMuZ(1:N,0),VQuad(1:N,i,j), DMuZVQuad)
            DivV(m)      = MuZ(0)**3*DivVQuad(i,j)
     .                   + 3*MuZ(0)**2*DMuZVQuad
          enddo
        enddo
      endif
c
c  ...TRIANGULAR FACE SHAPE FUNCTIONS
      call BlendProjectPyraTF(Mu,DMu,Nu,DNu,
     .                                 MubTF,DMubTF,NupTF,DNupTF,IdecTF)
      call ProjectPyraLamTF(Lam,DLam, LampTF,DLampTF,IdecLamTF)
c  ...loop over faces
      do f=1,4
c    ...local parameters
        nordTF = Nord(9+f)
        ndofTF = (nordTF+1)*nordTF/2
        if (ndofTF.gt.0) then
c      ...local parameters (again)
          minI  = 0
          minJ  = 0
          minIJ = minI+minJ
          maxIJ = nordTF-1
          maxI  = maxIJ-minJ
          maxJ  = maxIJ-minI
c      ...orient first
          call OrientTri(NupTF(0:2,f),DNupTF(1:N,0:2,f),NoriF(f+1),N,
     .                                                   GNupTF,GDNupTF)
          call OrientTri(LampTF(0:2,f),DLampTF(1:N,0:2,f),NoriF(f+1),N,
     .                                                 GLampTF,GDLampTF)
c        ...construct the shape functions
          call AncVTri(GNupTF,GDNupTF,nordTF,IdecTF,N,
     .                                    VTri(1:N,minI:maxI,minJ:maxJ),
     .                                     DivVTri(minI:maxI,minJ:maxJ))
          call AncVTri(GLampTF,GDLampTF,nordTF,IdecLamTF,N,
     .                                   VTri2(1:N,minI:maxI,minJ:maxJ),
     .                                    DivVTri2(minI:maxI,minJ:maxJ))
          do nij=minIJ,maxIJ
            do i=minI,nij-minJ
              j=nij-i
              m=m+1
c
              vectmp1 = MubTF(f)*VTri(1:N,i,j)
              vectmp2 = VTri2(1:N,i,j)/MubTF(f)
              ShapV(1:N,m) = (vectmp1+vectmp2)/2
c
              tmp1 = MubTF(f)*DivVTri(i,j)
              tmp2 = DivVTri2(i,j)/MubTF(f)
              vectmp = VTri(1:N,i,j)-VTri2(1:N,i,j)/(MubTF(f)**2)
              call dot_product(DMubTF(1:N,f),vectmp, tmp3)
              DivV(m)      = (tmp1+tmp2+tmp3)/2
            enddo
          enddo
        endif
      enddo
c
c  ...BUBBLE FUNCTIONS
      nordB = Nord(14)
c  ...FAMILY 1 AND 2 (curl of families 2 and 3 from H(curl))
      ndofB = nordB*(nordB-1)**2
c  ...if necessary, create bubbles
      if (ndofB.gt.0) then
c    ...local parameters
        IdecB(1:2) = .TRUE.; IdecB(3) = .TRUE.
        minF(1) = 0
        minF(2) = 2
        minK    = 2
        maxF(1) = nordB-1
        maxF(2) = nordB
        maxK    = nordB
c    ...loop over families
        do fam=0,1
          ab = cshift((/1,2/),fam);
          a = ab(1); b = ab(2)
c    ...construct shape functions, no need to orient
          call AncEQuad(Mu(0:1,ab),DMu(1:N,0:1,ab),(/nordB,nordB/),
     .                                                    IdecB(ab),N,
     .                       EQuad(1:N,minF(1):maxF(1),minF(2):maxF(2)),
     .                   CurlEQuad(1:N,minF(1):maxF(1),minF(2):maxF(2)))
          call AncPhiE(MuZ,DMuZ,nordB,IdecB(3),N,
     .                             phiE(minK:maxK),DphiE(1:N,minK:maxK))
c      ...in the code the outer loop always is
c      ...numbered wrt the second global face axis
          minF = cshift(minF,-fam); maxF = cshift(maxF,-fam)
          do k=minK,maxK
            do jg=minF(2),maxF(2)
              do ig=minF(1),maxF(1)
                ij = cshift((/ig,jg/),fam)
                i = ij(1); j = ij(2)
                m=m+1
c
                call cross(3,DMuZ(1:N,0)*phiE(k) + MuZ(0)*DphiE(1:N,k),
     .                                           EQuad(1:N,i,j), vectmp)
                ShapV(1:N,m) = MuZ(0)*CurlEQuad(1:N,i,j)*phiE(k)
     .                       + vectmp
                DivV(m)      = 0
              enddo
            enddo
          enddo
        enddo
      endif
c  ...FAMILY 3 (curl of family 4 from H(curl))
      ndofB = (nordB-1)**2
c  ...if necessary, create bubbles
      if (ndofB.gt.0) then
c    ...local parameters (again)
        IdecB(1:2) = .TRUE.; IdecB(3) = .TRUE.
        minI = 2
        minJ = 2
        maxI = nordB
        maxJ = nordB
c    ...construct shape functions, no need to orient
        call AncPhiQuad(Mu,DMu,(/nordB,nordB/),IdecB(1:2),N,
     .                                     phiQuad(minI:maxI,minJ:maxJ),
     .                                DphiQuad(1:N,minI:maxI,minJ:maxJ))
        do j=minJ,maxJ
          do i=minI,maxI
            m=m+1
            ijmax = max(i,j)
c
            call cross(N,DphiQuad(1:N,i,j),DMuZ(1:N,1), DphiQuadxDMuZ)
            ShapV(1:N,m) = DphiQuadxDMuZ*MuZ(0)**(ijmax-1)
            DivV(m)      = 0
          enddo
        enddo
      endif
c  ...FAMILY 4 (induced from quad face functions)
      ndofB = nordB**2*(nordB-1)
c  ...if necessary, create bubbles
      if (ndofB.gt.0) then
c    ...local parameters (again)
        IdecB(1:2) = .TRUE.
        minI = 0
        minJ = 0
        minK = 2
        maxI = nordB-1
        maxJ = nordB-1
        maxK = nordB
c    ...construct shape functions, no need to orient
        call AncVQuad(Mu,DMu,(/nordB,nordB/),IdecB(1:2),N,
     .                                   VQuad(1:N,minI:maxI,minJ:maxJ),
     .                                    DivVQuad(minI:maxI,minJ:maxJ))
        call AncPhiE(MuZ,DMuZ,nordB,IdecB(3),N,
     .                             phiE(minK:maxK),DphiE(1:N,minK:maxK))
        do k=minK,maxK
          do j=minJ,maxJ
            do i=minI,maxI
              m=m+1
c
              ShapV(1:N,m) = MuZ(0)**2*VQuad(1:N,i,j)*phiE(k)
              call dot_product(2*MuZ(0)*phiE(k)*DMuZ(1:N,0)
     .                     + MuZ(0)**2*DphiE(1:N,k),VQuad(1:N,i,j), tmp)
              DivV(m)      = MuZ(0)**2*DivVQuad(i,j)*phiE(k)
     .                     + tmp
            enddo
          enddo
        enddo
      endif
c  ...FAMILY 5
      ndofB = (nordB-1)**2
c  ...if necessary, create bubbles
      if (ndofB.gt.0) then
c    ...local parameters (again)
        IdecB(1:2) = .TRUE.
        minI = 2
        minJ = 2
        maxI = nordB
        maxJ = nordB
c    ...construct shape functions, no need to orient
        call AncPhiE(Mu(0:1,1),DMu(1:N,0:1,1),nordB,IdecB(1),N,
     .                             phiE(minI:maxI),DphiE(1:N,minI:maxI))
        call AncPhiE(Mu(0:1,2),DMu(1:N,0:1,2),nordB,IdecB(2),N,
     .                           phiE2(minJ:maxJ),DphiE2(1:N,minJ:maxJ))
        do j=minJ,maxJ
          do i=minI,maxI
            m=m+1
c
            call cross(N,DphiE(1:N,i),DphiE2(1:N,j), vectmp1)
            call cross(N,DMuZ(1:N,0),phiE(i)*DphiE2(1:N,j)
     .                                 - DphiE(1:N,i)*phiE2(j), vectmp2)
            vectmp = MuZ(0)*(MuZ(0)*vectmp1+vectmp2)
            ShapV(1:N,m) = MuZ(1)**(i-1)*vectmp
            call dot_product(DMuZ(1:N,1),vectmp, tmp)
            DivV(m)      = (i-1)*MuZ(1)**(i-2)*tmp
          enddo
        enddo
      endif
c  ...FAMILY 6 AND 7
      ndofB = nordB-1
c  ...if necessary, create bubbles
      if (ndofB.gt.0) then
c    ...local parameters (again)
        IdecB(1:2) = .TRUE.
        minI = 2
        maxI = nordB
        do fam=0,1
          ab = cshift((/1,2/),fam)
          a = ab(1); b = ab(2)
c      ...construct shape functions, no need to orient
          call AncPhiE(Mu(0:1,a),DMu(1:N,0:1,a),nordB,IdecB(a),N,
     .                             phiE(minI:maxI),DphiE(1:N,minI:maxI))
          do i=minI,maxI
            m=m+1
c
            call cross(N,DphiE(1:N,i)*MuZ(0)**2
     .              + 2*MuZ(0)*phiE(i)*DMuZ(1:N,0),DMu(1:N,b,1), vectmp)
            ShapV(1:N,m) = MuZ(1)**(i-1)*vectmp
            call dot_product(DMuZ(1:N,1),vectmp, tmp)
            DivV(m)      = (i-1)*MuZ(1)**(i-2)*tmp
          enddo
        enddo
      endif
c
c  ...total degrees of freedom
      NrdofV = m
c
#if DEBUG_MODE
c  ...print this when debugging
      if (iprint.ge.1) then
        write(*,7001)Xi(1:3),Nord(1:14),NoriF(1:5),NrdofV
 7001   format('shape3DVPyra: Xi = ',3f8.3,/,
     .         'Norder = ',4i3,1x,4i3,3x,i3,1x,4i3,3x,i3,/,
     .         'orient = ',i3,1x,4i3,/,
     .         'NrdofV = ',3i3)
        m=0
        call decod(Nord(9),MODORDER,2, nordQF)
        ndofQF = nordQF(1)*nordQF(2)
        if (ndofQF.gt.0) then
          write(*,*) 'SHAPE FUNCTIONS FOR QUADRILATERAL FACE'
          do j=1,ndofQF
            m=m+1
            write(*,7002) m,ShapV(1:N,m),DivV(m)
 7002       format('k = ',i3,' ShapV, DivV = ',3e12.5,3x,e12.5)
          enddo
        endif
        do f=1,4
          nordTF = Nord(9+f)
          ndofTF = (nordTF+1)*nordTF/2
          if (ndofTF.gt.0) then
            write(*,7014) f
 7014       format('SHAPE FUNCTIONS FOR TRIANGLAR FACE = ',i2)
            do j=1,ndofTF
              m=m+1
              write(*,7002) m,ShapV(1:N,m),DivV(m)
            enddo
          endif
        enddo
        nordB = Nord(14)
        ndofB = nordB*(nordB-1)**2
        if (ndofB.gt.0) then
          write(*,*) 'BUBBLES'
          do fam=0,1
            write(*,7004) fam+1
 7004       format('Family ',i2)
            do j=1,ndofB
              m=m+1
              write(*,7002) m,ShapV(1:N,m),DivV(m)
            enddo
          enddo
        endif
        ndofB = (nordB-1)**2
        if (ndofB.gt.0) then
          write(*,*) 'Family 3'
          do j=1,ndofB
            m=m+1
            write(*,7002) m,ShapV(1:N,m),DivV(m)
          enddo
        endif
        ndofB = nordB**2*(nordB-1)
        if (ndofB.gt.0) then
          write(*,*) 'Family 4'
          do j=1,ndofB
            m=m+1
            write(*,7002) m,ShapV(1:N,m),DivV(m)
          enddo
        endif
        ndofB = (nordB-1)**2
        if (ndofB.gt.0) then
          write(*,*) 'Family 5'
          do j=1,ndofB
            m=m+1
            write(*,7002) m,ShapV(1:N,m),DivV(m)
          enddo
        endif
        ndofB = nordB-1
        if (ndofB.gt.0) then
          do fam=0,1
            write(*,7004) fam+6
            do j=1,ndofB
              m=m+1
              write(*,7002) m,ShapV(1:N,m),DivV(m)
            enddo
          enddo
        endif
        call pause
      endif
#endif
c
c
      end subroutine shape3DVPyra


c--------------------------------------------------------------------
c
c     routine name      - shape3DQPyra
c
c--------------------------------------------------------------------
c
c     latest revision:  - Jan 15, Apr 17
c
c     purpose:          - routine returns values of 3D pyramid
c                         element L2 shape functions
c
c     arguments:
c
c     in:
c          Xi           - master pyramid coordinates from (0,1)^3
c          Nord         - polynomial order for the nodes (H1 sense)
c          Nsize        - relevant sizes of local arrays
c
c     out:
c          NrdofQ       - number of dof
c          ShapQ        - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape3DQPyra(Xi,Nord,Nsize, NrdofQ,ShapQ)
c
      implicit none
      integer, intent(in)  :: Nord(14),Nsize(2)
      integer, intent(out) :: NrdofQ
      integer :: N,m,i,j,k,nordB,ndofB,minI,minJ,minK,maxI,maxJ
      integer :: maxK
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapQ(Nsize(2))
      double precision :: Lam(5),DLam(3,5)
      double precision :: Mu(0:1,2),DMu(3,0:1,2)
      double precision :: Nu(0:2,2),DNu(3,0:2,2)
      double precision :: MuZ(0:1),DMuZ(3,0:1)
      double precision :: homP(0:Nsize(1)-1,3)
c
#if DEBUG_MODE
c  ...debugging flag
      integer :: iprint=0
#endif
c
c  ...spatial dimensions
      N=3
c
c  ...initiate counter for shape functions
      m=0
c
c  ...Define affine coordinates and gradients
      call AffinePyramid(Xi, Lam,DLam,Mu,DMu,Nu,DNu,MuZ,DMuZ)
c
c  ...local parameters
      nordB = Nord(14)
      ndofB = nordB**3
      minI = 0
      minJ = 0
      minK = 0
      maxI = Nord(14)-1
      maxJ = Nord(14)-1
      maxK = Nord(14)-1
c
c  ...get homogenized Legendre polynomials, homP
      call HomLegendre(Mu(0:1,1),maxI, homP(minI:maxI,1))
      call HomLegendre(Mu(0:1,2),maxJ, homP(minJ:maxJ,2))
      call HomLegendre(MuZ,maxK, homP(minK:maxK,3))
c
c  ...construct shape functions
      do k=minK,maxK
        do j=minJ,maxJ
          do i=minI,maxI
            m=m+1
c
            ShapQ(m) = homP(i,1)*homP(j,2)*homP(k,3)
          enddo
        enddo
      enddo
c
c  ...give total degrees of freedom
      NrdofQ = m
c
#if DEBUG_MODE
c  ...print this when debugging
      if (iprint.ge.1) then
        write(*,7001) Xi(1:3),Nord(14)
 7001   format('shape3DQPyra: Xi = ',3f8.3,/,
     .         'Norder = ',i2)
        nordB = Nord(14)
        ndofB = nordB**3
        if (ndofB.gt.0) then
          write(*,7010)
 7010     format('BUBBLES = ')
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
c
      end subroutine shape3DQPyra

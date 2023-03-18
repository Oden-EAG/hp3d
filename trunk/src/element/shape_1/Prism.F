c Routines:
c  - shape3DHPris
c  - shape3DEPris
c  - shape3DVPris
c  - shape3DQPris
c--------------------------------------------------------------------
c
c     routine name      - shape3DHPris
c
c--------------------------------------------------------------------
c
c     latest revision:  - Oct 14, Apr 17, Jul 21
c
c     purpose:          - routine returns values of 3D triangular prism
c                         element H1 shape functions and their derivatives
c
c     arguments:
c
c     in:
c          X            - master prism coordinates from (0,1)^3
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
      subroutine shape3DHPris(X,Nord,NoriE,NoriF,Nsize,
     .                                              NrdofH,ShapH,GradH)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: Nord(15),NoriE(9),NoriF(5),Nsize(2)
      integer, intent(out) :: NrdofH
      integer :: N,m,v,e,f,i,j,k,nij,nordME,ndofME,nordQE,ndofQE
      integer :: nordTF,ndofTF,nordQF(2),ndofQF,nordB(2),ndofB
      integer :: minI,minJ,minK,minIJ,maxI,maxJ,maxK,maxIJ
      logical :: IdecME,IdecQE,IdecTF,IdecQF(2,3),GIdecQF(2),IdecB(2)
      double precision, intent(in)  :: X(3)
      double precision, intent(out) :: ShapH(Nsize(2))
      double precision, intent(out) :: GradH(1:3,Nsize(2))
      double precision :: Mu(0:1),DMu(3,0:1)
      double precision :: Nu(0:2),DNu(3,0:2)
      double precision :: MubV(6),DMubV(3,6),NubV(6),DNubV(3,6)
      double precision :: MubME(6),DMubME(3,6)
      double precision :: NupME(0:1,6),DNupME(3,0:1,6)
      double precision :: GNupME(0:1),GDNupME(3,0:1)
      double precision :: NubQE(3),DNubQE(3,3)
      double precision :: MupQE(0:1,3),DMupQE(3,0:1,3)
      double precision :: GMupQE(0:1),GDMupQE(3,0:1)
      double precision :: MubTF(2),DMubTF(3,2)
      double precision :: NupTF(0:2,2),DNupTF(3,0:2,2)
      double precision :: GNupTF(0:2),GDNupTF(3,0:2)
      double precision :: STpQF(0:1,2,3),DSTpQF(3,0:1,2,3)
      double precision :: GSTpQF(0:1,2),GDSTpQF(3,0:1,2)
      double precision :: phiE(2:Nsize(1)),DphiE(3,2:Nsize(1))
      double precision :: phiTri(2:Nsize(1)-1,1:Nsize(1)-2)
      double precision :: DphiTri(3,2:Nsize(1)-1,1:Nsize(1)-2)
      double precision :: phiQuad(2:Nsize(1),2:Nsize(1))
      double precision :: DphiQuad(3,2:Nsize(1),2:Nsize(1))
c
#if DEBUG_MODE
c  ...debugging flag
      integer :: iprint
      iprint=0
#endif
c
c  ...spatial dimensions
      N=3
c
c  ...initiate counter for shape functions
      m=0
c
c  ...Define affine coordinates and gradients
      call AffinePrism(X, Mu,DMu,Nu,DNu)
c
c  ...VERTEX SHAPE FUNCTIONS
      call BlendPrisV(Mu,DMu,Nu,DNu, MubV,DMubV,NubV,DNubV)
      do v=1,6
        m=m+1
        ShapH(m) = NubV(v)*MubV(v)
        GradH(1:N,m) = DNubV(1:N,v)*MubV(v)
     .               + NubV(v)*DMubV(1:N,v)
      enddo
c
c  ...EDGE SHAPE FUNCTIONS
c  ...mixed edges
      call BlendProjectPrisME(Mu,DMu,Nu,DNu,
     .                                 MubME,DMubME,NupME,DNupME,IdecME)
c  ...loop over edges
      do e=1,6
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
            ShapH(m) = phiE(i)*MubME(e)
            GradH(1:N,m) = DphiE(1:N,i)*MubME(e)
     .                   + phiE(i)*DMubME(1:N,e)
          enddo
        endif
      enddo
c  ...quadrilateral edges
      call BlendProjectPrisQE(Mu,DMu,Nu,DNu,
     .                                 NubQE,DNubQE,MupQE,DMupQE,IdecQE)
c  ...loop over edges
      do e=1,3
c    ...local parameters
        nordQE = Nord(6+e)
        ndofQE = nordQE-1
        if (ndofQE.gt.0) then
c      ...local parameters (again)
          minI = 2
          maxI = nordQE
c      ...orient first
          call OrientE(MupQE(0:1,e),DMupQE(1:N,0:1,e),NoriE(6+e),N,
     .                                                   GMupQE,GDMupQE)
c      ...construct the shape functions
          call AncPhiE(GMupQE,GDMupQE,nordQE,IdecQE,N,
     .                             phiE(minI:maxI),DphiE(1:N,minI:maxI))
          do i=minI,maxI
            m=m+1
            ShapH(m) = phiE(i)*NubQE(e)
            GradH(1:N,m) = DphiE(1:N,i)*NubQE(e)
     .                   + phiE(i)*DNubQE(1:N,e)
          enddo
        endif
      enddo
c
c  ...FACE SHAPE FUNCTIONS
c  ...triangle faces
      call BlendProjectPrisTF(Mu,DMu,Nu,DNu,
     .                                 MubTF,DMubTF,NupTF,DNupTF,IdecTF)
c  ...loop over faces
      do f=1,2
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
c      ...orient
          call OrientTri(NupTF(0:2,f),DNupTF(1:N,0:2,f),NoriF(f),N,
     .                                                   GNupTF,GDNupTF)
c      ...construct the shape functions
          call AncPhiTri(GNupTF,GDNupTF,NordTF,IdecTF,N,
     .                                      phiTri(minI:maxI,minJ:maxJ),
     .                                 DphiTri(1:N,minI:maxI,minJ:maxJ))
          do nij=minIJ,maxIJ
            do i=minI,nij-minJ
              j=nij-i
              m=m+1
c
              ShapH(m) = phiTri(i,j)*MubTF(f)
              GradH(1:N,m) = DphiTri(1:N,i,j)*MubTF(f)
     .                     + phiTri(i,j)*DMubTF(1:N,f)
            enddo
          enddo
        endif
      enddo
c  ...quadrilateral faces
      call ProjectPrisQF(Mu,DMu,Nu,DNu, STpQF,DSTpQF,IdecQF)
c  ...loop over faces
      do f=1,3
c    ...local parameters
        call decod(Nord(11+f),MODORDER,2, nordQF)
        ndofQF = (nordQF(1)-1)*(nordQF(2)-1)
        if (ndofQF.gt.0) then
c      ...local parameters (again)
          minI = 2
          minJ = 2
          maxI = nordQF(1)
          maxJ = nordQF(2)
c      ...orient
          call OrientQuad(STpQF(0:1,1:2,f),DSTpQF(1:N,0:1,1:2,f),
     .               NoriF(f+2),IdecQF(1:2,f),N, GSTpQF,GDSTpQF,GIdecQF)
c      ...construct the shape functions
          call AncPhiQuad(GSTpQF,GDSTpQF,nordQF,GIdecQF,N,
     .                                     phiQuad(minI:maxI,minJ:maxJ),
     .                                DphiQuad(1:N,minI:maxI,minJ:maxJ))
          do j=minJ,maxJ
            do i=minI,maxI
              m=m+1
              ShapH(m) = phiQuad(i,j)
              GradH(1:N,m) = DphiQuad(1:N,i,j)
            enddo
          enddo
        endif
      enddo
c
c  ...BUBBLE FUNCTIONS
c  ...local parameters
      call decod(Nord(15),MODORDER,2, nordB)
      ndofB = (nordB(1)-1)*(nordB(1)-2)*(nordB(2)-1)/2
c  ...if necessary, create bubbles
      if (ndofB.gt.0) then
c    ...local parameters (again)
        IdecB(1) = IdecTF
        IdecB(2) = IdecQE
        minI = 2
        minJ = 1
        minK = 2
        minIJ = minI+minJ
        maxIJ = nordB(1)
        maxI = maxIJ-minJ
        maxJ = maxIJ-minI
        maxK = nordB(2)
c    ...call phiTri and phiE - no need to orient
        call AncPhiTri(Nu,DNu,nordB(1),IdecB(1),N,
     .                                      phiTri(minI:maxI,minJ:maxJ),
     .                                 DphiTri(1:N,minI:maxI,minJ:maxJ))
        call AncPhiE(Mu,DMu,nordB(2),IdecB(2),N,
     .                             phiE(minK:maxK),DphiE(1:N,minK:maxK))
        do k=minK,maxK
          do nij=minIJ,maxIJ
            do i=minI,nij-minJ
              j=nij-i
              m=m+1
c
              ShapH(m) = phiTri(i,j)*phiE(k)
              GradH(1:N,m) = DphiTri(1:N,i,j)*phiE(k)
     .                     + phiTri(i,j)*DphiE(1:N,k)
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
        write(*,7001) X(1:3),Nord(1:15),NoriE(1:9),NoriF(1:5),NrdofH
 7001   format('shape3DHPris: Xi = ',3f8.3,/,
     .         'Norder = ',3i3,1x,3i3,2x,3i3,3x,2i3,2x,3i3,1x,i3,/,
     .         'orient = ',3i3,1x,3i3,2x,3i3,3x,2i3,2x,3i3,/,
     .         'NrdofH = ',i3)
        write(*,7010)
 7010   format('VERTEX SHAPE FUNCTIONS = ')
        do v=1,6
          m=v
          write(*,7002) m,ShapH(m),GradH(1:3,m)
 7002     format('k = ',i3,' ShapH, GradH = ',e12.5,3x,3e12.5)
        enddo
        do e=1,6
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
        do e=1,3
          ndofQE = Nord(6+e)-1
          if (ndofQE.gt.0) then
            write(*,7012) e
 7012       format('SHAPE FUNCTIONS FOR QUAD EDGE = ',i2)
            do j=1,ndofQE
              m=m+1
              write(*,7002) m,ShapH(m),GradH(1:3,m)
            enddo
          endif
        enddo
        do f=1,2
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
        do f=1,3
          call decod(Nord(11+f),MODORDER,2, nordQF)
          ndofQF = (nordQF(1)-1)*(nordQF(2)-1)
          if (ndofQF.gt.0) then
            write(*,7014) f
 7014       format('SHAPE FUNCTIONS FOR QUAD FACE = ',i2)
            do j=1,ndofQF
              m=m+1
              write(*,7002) m,ShapH(m),GradH(1:3,m)
            enddo
          endif
        enddo
        call decod(Nord(15),MODORDER,2, nordB)
        ndofB = (nordB(1)-1)*(nordB(1)-2)*(nordB(2)-1)/2
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
      end subroutine shape3DHPris

c--------------------------------------------------------------------
c
c     routine name      - shape3EPris
c
c--------------------------------------------------------------------
c
c     latest revision:  - Oct 14, Apr 17, Jul 21
c
c     purpose:          - routine returns values of 3D triangular prism
c                         element H(curl) shape functions and their
c                         derivatives
c
c     arguments:
c
c     in:
c          X            - master prism coordinates from (0,1)^3
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

      subroutine shape3DEPris(X,Nord,NoriE,NoriF,Nsize,
     .                                              NrdofE,ShapE,CurlE)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: Nord(15),NoriE(9),NoriF(5),Nsize(2)
      integer, intent(out) :: NrdofE
      integer :: N,m,e,f,i,j,k,nij,nordME,ndofME,nordQE,ndofQE
      integer :: nordTF,ndofTF,nordQF(2),ndofQF,nordB(2),ndofB
      integer :: minI,minJ,minK,minIJ,maxI,maxJ,maxK,maxIJ
      integer :: famctr,fam,a,b,ab(2),abc(3),ij(2),ig,jg,minF(2),maxF(2)
      logical :: IdecME,IdecQE,IdecTF,IdecQF(2,3),GIdecQF(2),IdecB(2)
      double precision, intent(in)  :: X(3)
      double precision, intent(out) :: ShapE(3,Nsize(2))
      double precision, intent(out) :: CurlE(3,Nsize(2))
      double precision :: Mu(0:1),DMu(3,0:1)
      double precision :: Nu(0:2),DNu(3,0:2)
      double precision :: MubME(6),DMubME(3,6)
      double precision :: NupME(0:1,6),DNupME(3,0:1,6)
      double precision :: GNupME(0:1),GDNupME(3,0:1)
      double precision :: NubQE(3),DNubQE(3,3)
      double precision :: MupQE(0:1,3),DMupQE(3,0:1,3)
      double precision :: GMupQE(0:1),GDMupQE(3,0:1)
      double precision :: MubTF(2),DMubTF(3,2)
      double precision :: NupTF(0:2,2),DNupTF(3,0:2,2)
      double precision :: GNupTF(0:2),GDNupTF(3,0:2)
      double precision :: STpQF(0:1,2,3),DSTpQF(3,0:1,2,3)
      double precision :: GSTpQF(0:1,2),GDSTpQF(3,0:1,2)
      double precision :: EE(3,0:Nsize(1)-1),CurlEE(3,0:Nsize(1)-1)
      double precision :: ETri(3,0:Nsize(1)-2,1:Nsize(1)-1)
      double precision :: CurlETri(3,0:Nsize(1)-2,1:Nsize(1)-1)
      double precision :: EQuad(3,0:Nsize(1)-1,2:Nsize(1))
      double precision :: CurlEQuad(3,0:Nsize(1)-1,2:Nsize(1))
      double precision :: PhiE(2:Nsize(1)),DPhiE(3,2:Nsize(1))
      double precision :: PhiTri(2:Nsize(1)-1,1:Nsize(1)-2)
      double precision :: DPhiTri(3,2:Nsize(1)-1,1:Nsize(1)-2)
      double precision :: DMubMExEE(3),DMubTFxETri(3),DPhiExETri(3)
c
#if DEBUG_MODE
c  ...debugging flag
      integer :: iprint
      iprint=0
#endif
c
c  ...spatial dimensions
      N=3
c
c  ...initiate counter for shape functions
      m=0
c
c  ...Define affine coordinates and gradients
      call AffinePrism(X, Mu,DMu,Nu,DNu)
c
c  ...EDGE SHAPE FUNCTIONS
c  ...MIXED EDGES
      call BlendProjectPrisME(Mu,DMu,Nu,DNu,
     .                                 MubME,DMubME,NupME,DNupME,IdecME)
c  ...loop over edges
      do e=1,6
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
c  ...QUADRILATERAL EDGES
      call BlendProjectPrisQE(Mu,DMu,Nu,DNu,
     .                                 NubQE,DNubQE,MupQE,DMupQE,IdecQE)
c  ...loop over edges
      do e=1,3
c    ...local parameters
        nordQE = Nord(6+e)
        ndofQE = nordQE
        if (ndofQE.gt.0) then
c      ...local parameters (again)
          minI = 0
          maxI = nordQE-1
c      ...orient first
          call OrientE(MupQE(0:1,e),DMupQE(1:N,0:1,e),NoriE(6+e),N,
     .                                                   GMupQE,GDMupQE)
c      ...construct the shape functions (CurlEE should be returned as all 0)
          call AncEE(GMupQE,GDMupQE,nordQE,IdecQE,N,
     .                          EE(1:N,minI:maxI),CurlEE(1:N,minI:maxI))
          do i=minI,maxI
            m=m+1
c
            ShapE(1:N,m) = NubQE(e)*EE(1:3,i)
            call cross(3,DNubQE(1:N,e),EE(1:N,i), CurlE(1:N,m))
          enddo
        endif
      enddo
c
c  ...FACE SHAPE FUNCTIONS
c  ...triangle faces
      call BlendProjectPrisTF(Mu,DMu,Nu,DNu,
     .                                 MubTF,DMubTF,NupTF,DNupTF,IdecTF)
c  ...loop over faces
      do f=1,2
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
c      ...orient
          call OrientTri(NupTF(0:2,f),DNupTF(1:N,0:2,f),NoriF(f),N,
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
c
                call cross(3,DMubTF(1:N,f),ETri(1:N,i,j), DMubTFxETri)
                CurlE(1:N,m) = MubTF(f)*CurlETri(1:N,i,j)
     .                       + DMubTFxETri
              enddo
            enddo
          enddo
        endif
      enddo
c  ...quadrilateral faces
      call ProjectPrisQF(Mu,DMu,Nu,DNu, STpQF,DSTpQF,IdecQF)
c  ...loop over faces
      do f=1,3
c    ...local parameters
        call decod(Nord(11+f),MODORDER,2, nordQF)
c    ...orient
        call OrientQuad(STpQF(0:1,1:2,f),DSTpQF(1:N,0:1,1:2,f),
     .               NoriF(2+f),IdecQF(1:2,f),N, GSTpQF,GDSTpQF,GIdecQF)
c    ...loop over families
        do fam=0,1
          ab = cshift((/1,2/),fam);
          a = ab(1); b = ab(2)
          ndofQF = nordQF(a)*(nordQF(b)-1)
          if (ndofQF.gt.0) then
c        ...local parameters (again)
            minF(1) = 0
            minF(2) = 2
            maxF(1) = nordQF(a)-1
            maxF(2) = nordQF(b)
c        ...construct the shape functions
            call AncEQuad(GSTpQF(0:1,ab),GDSTpQF(1:N,0:1,ab),
     .                    nordQF(ab),GIdecQF(ab),N,
     .                       EQuad(1:N,minF(1):maxF(1),minF(2):maxF(2)),
     .                   CurlEQuad(1:N,minF(1):maxF(1),minF(2):maxF(2)))
c        ...in the code the outer loop always is
c        ...numbered wrt the second global face axis
            minF = cshift(minF,-fam); maxF = cshift(maxF,-fam)
            do jg=minF(2),maxF(2)
              do ig=minF(1),maxF(1)
                ij = cshift((/ig,jg/),fam)
                i = ij(1); j = ij(2)
                m=m+1
c
                ShapE(1:N,m) = EQuad(1:N,i,j)
                CurlE(1:N,m) = CurlEQuad(1:N,i,j)
              enddo
            enddo
          endif
        enddo
      enddo
c
c  ...BUBBLE FUNCTIONS
c  ...Families 1 and 2 (Triangle type)
c  ...local parameters
      call decod(Nord(15),MODORDER,2, nordB)
      ndofB = nordB(1)*(nordB(1)-1)*(nordB(2)-1)/2
c  ...if necessary, create bubbles
      if (ndofB.gt.0) then
c    ...local parameters (again)
        IdecB(1) = IdecTF
        IdecB(2) = IdecQE
        minI  = 0
        minJ  = 1
        minK  = 2
        minIJ = minI+minJ
        maxIJ = nordB(1)-1
        maxI  = maxIJ-minJ
        maxJ  = maxIJ-minI
        maxK  = nordB(2)
c    ...loop over families
        famctr=m
        do fam=0,1
          m=famctr+fam-1
          abc = cshift((/0,1,2/),fam)
c      ...now construct the shape functions (no need to orient)
          call AncETri(Nu(abc),DNu(1:N,abc),nordB(1),IdecB(1),N,
     .                                    ETri(1:N,minI:maxI,minJ:maxJ),
     .                                CurlETri(1:N,minI:maxI,minJ:maxJ))
          call AncPhiE(Mu,DMu,nordB(2),IdecB(2),N,
     .                             PhiE(minK:maxK),DPhiE(1:N,minK:maxK))
          do k=minK,maxK
            do nij=minIJ,maxIJ
              do i=minI,nij-minJ
                j=nij-i
                m=m+2
c
                ShapE(1:N,m) = ETri(1:N,i,j)*PhiE(k)
c
                call cross(N,DPhiE(1:N,k),ETri(1:N,i,j), DPhiExETri)
                CurlE(1:N,m) = PhiE(k)*CurlETri(1:N,i,j)
     .                       + DPhiExETri
              enddo
            enddo
          enddo
        enddo
      endif
c  ...Family 3 (Quadrilateral type)
c  ...local parameters
      ndofB = (nordB(1)-1)*(nordB(1)-2)*nordB(2)/2
c  ...if necessary, create bubbles
      if (ndofB.gt.0) then
c    ...local parameters (again)
        IdecB(1) = IdecTF
        IdecB(2) = IdecQE
        minI  = 2
        minJ  = 1
        minK  = 0
        minIJ = minI+minJ
        maxIJ = nordB(1)
        maxI  = maxIJ-minJ
        maxJ  = maxIJ-minI
        maxK  = nordB(2)-1
c    ...now construct the shape functions (no need to orient)
        call AncphiTri(Nu,DNu,nordB(1),IdecB(1),N,
     .                                      PhiTri(minI:maxI,minJ:maxJ),
     .                                 DPhiTri(1:N,minI:maxI,minJ:maxJ))
        call AncEE(Mu,DMu,nordB(2),IdecB(2),N,
     .                          EE(1:N,minK:maxK),CurlEE(1:N,minK:maxK))
        do k=minK,maxK
          do nij=minIJ,maxIJ
            do i=minI,nij-minJ
              j=nij-i
              m=m+1
c
              ShapE(1:N,m) = PhiTri(i,j)*EE(1:N,k)
              call cross(3,DPhiTri(1:N,i,j),EE(1:N,k), CurlE(1:N,m))
            enddo
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
        write(*,7001) X(1:3),Nord(1:15),NoriE(1:9),NoriF(1:5),NrdofE
 7001   format('shape3DEPris: Xi = ',3f8.3,/,
     .         'Norder = ',3i3,1x,3i3,2x,3i3,3x,2i3,2x,3i3,1x,i3,/,
     .         'orient = ',3i3,1x,3i3,2x,3i3,3x,2i3,2x,3i3,/,
     .         'NrdofE = ',i3)
        m=0
        do e=1,6
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
        do e=1,3
          nordQE = Nord(6+e)
          ndofQE = nordQE
          if (ndofQE.gt.0) then
            write(*,7012) e
 7012       format('SHAPE FUNCTIONS FOR QUAD EDGE = ',i2)
            do j=1,ndofQE
              m=m+1
              write(*,7003) m,ShapE(1:N,m),CurlE(1:N,m)
 7003         format('k = ',i3,' ShapE, CurlE = ',3e12.5,3x,3e12.5)
            enddo
          endif
        enddo
        do f=1,2
          nordTF = Nord(9+f)
          ndofTF = nordTF*(nordTF-1)/2
          if (ndofTF.gt.0) then
            write(*,7013) f
 7013       format('SHAPE FUNCTIONS FOR TRIANGLE FACE = ',i2)
            famctr=m
            do fam=0,1
              m=famctr+fam-1
              write(*,7004) fam
 7004         format('family = ',i2)
              do j=1,ndofTF
                m=m+2
                write(*,7002) m,ShapE(1:N,m),CurlE(1:N,m)
              enddo
            enddo
          endif
        enddo
        do f=1,3
          call decod(Nord(11+f),MODORDER,2, nordQF)
          ndofQF = nordQF(a)*(nordQF(b)-1)
          if (ndofQF.gt.0) then
            write(*,7014) f
 7014       format('SHAPE FUNCTIONS FOR QUAD FACE = ',i2)
            do fam=0,1
              write(*,7004) fam
              do j=1,ndofQF
                m=m+1
                write(*,7002) m,ShapE(1:N,m),CurlE(1:N,m)
              enddo
            enddo
          endif
        enddo
        call decod(Nord(15),MODORDER,2, nordB)
        ndofB = nordB(1)*(nordB(1)-1)*(nordB(2)-1)/2
        if (ndofB.gt.0) then
          write(*,*) 'SHAPE FUNCTIONS FOR TRI-TYPE BUBBLES'
          famctr=m
          do fam=0,1
            m=famctr+fam-1
            write(*,7004) fam
            do j=1,ndofB
              m=m+2
              write(*,7002) m,ShapE(1:N,m),CurlE(1:N,m)
            enddo
          enddo
        endif
        ndofB = (nordB(1)-1)*(nordB(1)-2)*nordB(2)/2
        if (ndofB.gt.0) then
          write(*,*) 'SHAPE FUNCTIONS FOR QUAD-TYPE BUBBLES'
          write(*,7004) 2
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
      end subroutine shape3DEPris
c
c
c-----------------------------------------------------------------------
c
c     routine name      - shape3DVPris
c
c-----------------------------------------------------------------------
c
c     latest revision:  - Nov 14, Apr 17, Jul 21
c
c     purpose:          - routine returns values of 3D Prism element
c                         H(div) shape functions and their divergences
c
c     arguments:
c
c     in:
c          X            - master prism coordinates from (0,1)^3
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
      subroutine shape3DVPris(X,Nord,NoriF,Nsize, NrdofV,ShapV,DivV)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: Nord(15),NoriF(5),Nsize(2)
      integer, intent(out) :: NrdofV
      integer :: N,m,f,i,j,k,nij
      integer :: nordTF,ndofTF,nordQF(2),ndofQF,nordB(2),ndofB
      integer :: minI,minJ,minK,minIJ,maxI,maxJ,maxK,maxIJ
      integer :: famctr,fam,abc(3)
      logical :: IdecTF,IdecQF(2,3),GIdecQF(2),IdecB(2)
      double precision, intent(in)  :: X(3)
      double precision, intent(out) :: ShapV(3,Nsize(2))
      double precision, intent(out) :: DivV(Nsize(2))
      double precision :: Mu(0:1),DMu(3,0:1)
      double precision :: Nu(0:2),DNu(3,0:2)
      double precision :: MubTF(2),DMubTF(3,2)
      double precision :: NupTF(0:2,2),DNupTF(3,0:2,2)
      double precision :: GNupTF(0:2),GDNupTF(3,0:2)
      double precision :: STpQF(0:1,2,3),DSTpQF(3,0:1,2,3)
      double precision :: GSTpQF(0:1,2),GDSTpQF(3,0:1,2)
      double precision :: VTri(3,0:Nsize(1)-1,0:Nsize(1)-1)
      double precision :: DivVTri(0:Nsize(1)-1,0:Nsize(1)-1)
      double precision :: VQuad(3,0:Nsize(1)-1,0:Nsize(1)-1)
      double precision :: DivVQuad(0:Nsize(1)-1,0:Nsize(1)-1)
      double precision :: ETri(3,0:Nsize(1)-2,1:Nsize(1)-1)
      double precision :: CurlETri(3,0:Nsize(1)-2,1:Nsize(1)-1)
      double precision :: EE(3,0:Nsize(1)-1),CurlEE(3,0:Nsize(1)-1)
      double precision :: PhiE(2:Nsize(1)),DPhiE(3,2:Nsize(1))
c
#if DEBUG_MODE
c  ...debugging flag
      integer :: iprint
      iprint=0
#endif
c
c  ...spatial dimensions
      N=3
c
c  ...initiate counter for shape functions
      m=0
c
c  ...Define affine coordinates and gradients
      call AffinePrism(X, Mu,DMu,Nu,DNu)
c
c  ...FACE SHAPE FUNCTIONS
c  ...triangle faces
      call BlendProjectPrisTF(Mu,DMu,Nu,DNu,
     .                                 MubTF,DMubTF,NupTF,DNupTF,IdecTF)
      do f=1,2
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
c      ...orient
          call OrientTri(NupTF(0:2,f),DNupTF(1:N,0:2,f),NoriF(f),N,
     .                                                   GNupTF,GDNupTF)
c      ...construct the shape functions (DivVTri should be 0)
          call AncVTri(GNupTF,GDNupTF,nordTF,IdecTF,N,
     .                                    VTri(1:N,minI:maxI,minJ:maxJ),
     .                                     DivVTri(minI:maxI,minJ:maxJ))
          do nij=minIJ,maxIJ
            do i=minI,nij-minJ
              j=nij-i
              m=m+1
c
              ShapV(1:N,m) = MubTF(f)*VTri(1:N,i,j)
              call dot_product(DMubTF(1:N,f),VTri(1:N,i,j), DivV(m))
            enddo
          enddo
        endif
      enddo
c  ...quadrilateral faces
      call ProjectPrisQF(Mu,DMu,Nu,DNu, STpQF,DSTpQF,IdecQF)
c  ...loop over faces
      do f=1,3
c    ...local parameters
        call decod(Nord(11+f),MODORDER,2, nordQF)
        ndofQF = nordQF(1)*nordQF(2)
        if (ndofQF.gt.0) then
c      ...local parameters (again)
          minI = 0
          minJ = 0
          maxI = nordQF(1)-1
          maxJ = nordQF(2)-1
c      ...orient
          call OrientQuad(STpQF(0:1,1:2,f),DSTpQF(1:N,0:1,1:2,f),
     .               NoriF(2+f),IdecQF(1:2,f),N, GSTpQF,GDSTpQF,GIdecQF)
c      ...construct the shape functions
          call AncVQuad(GSTpQF,GDSTpQF,nordQF,GIdecQF,N,
     .                                   VQuad(1:N,minI:maxI,minJ:maxJ),
     .                                    DivVQuad(minI:maxI,minJ:maxJ))
          do j=minJ,maxJ
            do i=minI,maxI
              m=m+1
c
              ShapV(1:N,m) = VQuad(1:N,i,j)
              DivV(m) = DivVQuad(i,j)
            enddo
          enddo
        endif
      enddo
c
c  ...BUBBLE FUNCTIONS
c  ...Families 1 and 2 (Triangle type)
c  ...local parameters
      call decod(Nord(15),MODORDER,2, nordB)
      ndofB = nordB(1)*(nordB(1)-1)*nordB(2)/2
c  ...if necessary, create bubbles
      if (ndofB.gt.0) then
c    ...local parameters (again)
        IdecB(1) = IdecTF
        IdecB(2) = .TRUE.
        minI  = 0
        minJ  = 1
        minK  = 0
        minIJ = minI+minJ
        maxIJ = nordB(1)-1
        maxI  = maxIJ-minJ
        maxJ  = maxIJ-minI
        maxK  = nordB(2)-1
c    ...loop over families
        famctr=m
        do fam=0,1
          m=famctr+fam-1
          abc = cshift((/0,1,2/),fam)
c      ...now construct the shape functions (no need to orient)
          call AncETri(Nu(abc),DNu(1:N,abc),nordB(1)-minK,IdecB(1),N,
     .                                    ETri(1:N,minI:maxI,minJ:maxJ),
     .                                CurlETri(1:N,minI:maxI,minJ:maxJ))
          call AncEE(Mu,DMu,nordB(2),IdecB(2),N,
     .                          EE(1:N,minK:maxK),CurlEE(1:N,minK:maxK))
          do k=minK,maxK
            do nij=minIJ,maxIJ
              do i=minI,nij-minJ
                j=nij-i
                m=m+2
c
                call cross(3,ETri(1:N,i,j),EE(1:N,k), ShapV(1:N,m))
                call dot_product(EE(1:N,k),CurlETri(1:N,i,j), DivV(m))
              enddo
            enddo
          enddo
        enddo
      endif
c  ...Family 3 (Quadrilateral type)
c  ...local parameters
      ndofB = (nordB(1)+1)*nordB(1)*(nordB(2)-1)/2
c  ...if necessary, create bubbles
      if (ndofB.gt.0) then
c    ...local parameters (again)
        IdecB(1) = IdecTF
        IdecB(2) = .TRUE.
        minI  = 0
        minJ  = 0
        minK  = 2
        minIJ = minI+minJ
        maxIJ = nordB(1)-1
        maxI  = maxIJ-minJ
        maxJ  = maxIJ-minI
        maxK  = nordB(2)
c      ...construct the shape functions (DivVTri should be 0)
        call AncVTri(Nu,DNu,nordB(1),IdecB(1),N,
     .                                    VTri(1:N,minI:maxI,minJ:maxJ),
     .                                     DivVTri(minI:maxI,minJ:maxJ))
        call AncPhiE(Mu,DMu,nordB(2),IdecB(2),N,
     .                             PhiE(minK:maxK),DPhiE(1:N,minK:maxK))
        do k=minK,maxK
          do nij=minIJ,maxIJ
            do i=minI,nij-minJ
              j=nij-i
              m=m+1
c
              ShapV(1:N,m) = PhiE(k)*VTri(1:N,i,j)
              call dot_product(DPhiE(1:N,k),VTri(1:N,i,j),  DivV(m))
            enddo
          enddo
        enddo
      endif
c
c  ...give total degrees of freedom
      NrdofV = m
c
#if DEBUG_MODE
c  ...print this when debugging
      if (iprint.ge.1) then
        write(*,7001) X(1:3),Nord(1:15),NoriF(1:5),NrdofV
 7001   format('shape3DVPris: Xi = ',3f8.3,/,
     .         'Norder = ',3i3,1x,3i3,2x,3i3,3x,2i3,2x,3i3,1x,i3,/,
     .         'orient = ',2i3,2x,3i3,/,
     .         'NrdofV = ',i3)
        m=0
        do f=1,2
          nordTF = Nord(9+f)
          ndofTF = (nordTF+1)*nordTF/2
          if (ndofTF.gt.0) then
            write(*,7012) f
 7012       format('SHAPE FUNCTIONS FOR TRIANGLE FACE = ',i2)
            do j=1,ndofTF
              m=m+1
              write(*,7002) m,ShapV(1:N,m),DivV(m)
 7002         format('k = ',i3,' ShapV, DivV = ',3e12.5,3x,e12.5)
            enddo
          endif
        enddo
        do f=1,3
          call decod(Nord(11+f),MODORDER,2, nordQF)
          ndofQF = nordQF(1)*nordQF(2)
          if (ndofQF.gt.0) then
            write(*,7013) f
 7013       format('SHAPE FUNCTIONS FOR QUAD FACE = ',i2)
            do j=1,ndofQF
              m=m+1
              write(*,7002) m,ShapV(1:N,m),DivV(m)
            enddo
          endif
        enddo
        call decod(Nord(15),MODORDER,2, nordB)
        ndofB = nordB(1)*(nordB(1)-1)*nordB(2)/2
        if (ndofB.gt.0) then
          write(*,*) 'BUBBLES : '
          write(*,*) 'SHAPE FUNCTIONS FOR TRIANGLE-TYPE BUBBLES'
          famctr=m
          do fam=0,1
            m=famctr+fam-1
            write(*,7003) fam
 7003       format('family ',i2)
            do j=1,ndofB
              m=m+2
              write(*,7002) m,ShapV(1:N,m),DivV(m)
            enddo
          enddo
        endif
        ndofB = (nordB(1)+1)*nordB(1)*(nordB(2)-1)/2
        if (ndofB.gt.0) then
          write(*,*) 'SHAPE FUNCTIONS FOR QUAD-TYPE BUBBLES'
          write(*,7003) 2
          do j=1,ndofB
            m=m+1
            write(*,7002) m,ShapV(1:N,m),DivV(m)
          enddo
        endif
        call pause
      endif
#endif
c
c
      end subroutine shape3DVPris
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape3DQPris
c
c--------------------------------------------------------------------
c
c     latest revision:  - Nov 14, Apr 17
c
c     purpose:          - routine returns values of 3D Prism
c                         element L2 shape functions
c
c     arguments:
c
c     in:
c          X            - master prism coordinates from (0,1)^3
c          Nord         - polynomial order for the nodes (H1 sense)
c          Nsize        - relevant sizes of local arrays
c
c     out:
c          NrdofQ       - number of dof
c          ShapQ        - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape3DQPris(X,Nord,Nsize, NrdofQ,ShapQ)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: Nord(15),Nsize(2)
      integer, intent(out) :: NrdofQ
      integer :: i,j,k,nij,m,N,nordB(2),ndofB,minalpha
      integer :: minI,minJ,minK,minIJ,maxI,maxJ,maxK,maxIJ
      double precision, intent(in)  :: X(3)
      double precision, intent(out) :: ShapQ(Nsize(2))
      double precision :: Mu(0:1),DMu(3,0:1),Nu(0:2),DNu(3,0:2)
      double precision :: homP(0:Nsize(1)-1)
      double precision :: homPal(0:Nsize(1)-1,0:Nsize(1)-1)
      double precision :: homPz(0:Nsize(1)-1)
c
#if DEBUG_MODE
c  ...debugging flag
      integer :: iprint
      iprint=0
#endif
c
c  ...spatial dimensions
      N=3
c  ...initiate counter for shape functions
      m=0
c
c  ...Define affine coordinates and gradients
      call AffinePrism(X, Mu,DMu,Nu,DNu)
c  ...local parameters
      call decod(Nord(15),MODORDER,2, nordB)
      ndofB = (nordB(1)+1)*nordB(1)*nordB(2)/2
      minI  = 0
      minJ  = 0
      minK  = 0
      minIJ = minI+minJ
      maxIJ = nordB(1)-1
      maxI  = maxIJ-minJ
      maxJ  = maxIJ-minI
      maxK  = nordB(2)-1
      minalpha = 2*minI+1
c
c  ...get homogenized Legendre polynomials, homP
      call HomLegendre(Nu(0:1),maxI, homP(minI:maxI))
c
c  ...get homogenized Jacobi polynomials, homPal
      call HomJacobi((/Nu(0)+Nu(1),Nu(2)/),maxIJ,minalpha,
     .                                      homPal(minI:maxI,minJ:maxJ))
c  ...get homogenized Legendre polynomials in z-drection, homPz
      call HomLegendre(Mu(0:1),maxK, homPz(minK:maxK))
c
c  ...construct shape functions
      do k=minK,maxK
        do nij=minIJ,maxIJ
          do i=minI,nij-minJ
            j=nij-i
            m=m+1
c
            ShapQ(m) = homP(i)*homPal(i,j)*homPz(k)
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
        write(*,7001) X(1:3),Nord(15)
 7001   format('shape3DQPris: Xi = ',3f8.3,/,
     .         'Norder = ',i2)
        call decod(Nord(15),MODORDER,2, nordB)
        ndofB = (nordB(1)+1)*nordB(1)*nordB(2)/2
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
c
      end subroutine shape3DQPris

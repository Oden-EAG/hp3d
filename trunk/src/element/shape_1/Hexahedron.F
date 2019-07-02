c Routines:
c  - shape3DHHexa
c  - shape3DEHexa
c  - shape3DVHexa
c  - shape3DQHexa
c--------------------------------------------------------------------
c
c     routine name      - shape3DHHexa
c
c--------------------------------------------------------------------
c
c     latest revision:  - Oct 14, Apr 17
c
c     purpose:          - routine returns values of 3D hexahedron element
c                         H1 shape functions and their derivatives
c
c     arguments:
c
c     in:
c          X            - master hexahedron coordinates from (0,1)^3
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
      subroutine shape3DHHexa(Xi,Nord,NoriE,NoriF,Nsize, 
     .                                              NrdofH,ShapH,GradH)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: Nord(19),NoriE(12),NoriF(6),Nsize(2)
      integer, intent(out) :: NrdofH
      integer :: i,j,k,m,v,e,f,N,ndofE,nordF(2),ndofF
      integer :: nordB(3),ndofB,iprint
      logical :: IdecE,IdecF(2),GIdecF(2),IdecB(3)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapH(Nsize(2))
      double precision, intent(out) :: GradH(1:3,Nsize(2))
      double precision :: Mu(1:3,0:1),DMu(1:3,1:3,0:1)
      double precision :: MubV(1:8,1:3),DMubV(1:8,1:3,1:3)
      double precision :: MubE(1:12,1:2),DMubE(1:12,1:3,1:2)
      double precision :: MupE(1:12,0:1),DMupE(1:12,1:3,0:1)
      double precision :: GMupE(0:1),GDMupE(1:3,0:1)
      double precision :: MubF(1:6),DMubF(1:6,1:3)
      double precision :: MupF(1:6,1:2,0:1),DMupF(1:6,1:3,1:2,0:1)
      double precision :: GMupF(1:2,0:1),GDMupF(1:3,1:2,0:1)
      double precision :: phiE(2:Nsize(1)),DphiE(1:3,2:Nsize(1))
      double precision :: phiQuad(2:Nsize(1),2:Nsize(1))
      double precision :: DphiQuad(1:3,2:Nsize(1),2:Nsize(1))
c
c  ...debugging flag
      iprint=0
c  ...spatial dimensions
      N=3
c  ...initiate counter for shape functions
      m=0
c
c  ...Define affine coordinates
      call AffineHexahedron(Xi, Mu,DMu)
c
c  ...First the vertices
c  ...call the blending functions
      call BlendHexaV(Mu,DMu, MubV,DMubV)
      do v=1,8
        m=m+1
        ShapH(m) = MubV(v,1)*MubV(v,2)*MubV(v,3)
        GradH(1:N,m) = MubV(v,1)*MubV(v,2)*DMubV(v,1:N,3)
     .               + MubV(v,1)*DMubV(v,1:N,2)*MubV(v,3)
     .               + DMubV(v,1:N,1)*MubV(v,2)*MubV(v,3)
      enddo
c
c  ...Second the edges
c  ...call the blending and projections
      call BlendProjectHexaE(Mu,DMu, MubE,DMubE,MupE,DMupE,IdecE)
      do e=1,12
        ndofE=Nord(e)-1
        if (ndofE.gt.0) then        
c      ...orient first
          call OrientE(MupE(e,0:1),DMupE(e,1:N,0:1),NoriE(e),N, 
     .                                                 GMupE,GDMupE)
c      ...construct the shape functions
          call AncPhiE(GMupE,GDMupE,Nord(e),IdecE,N, 
     .                         phiE(2:Nord(e)),DphiE(1:N,2:Nord(e)))
          do i=2,Nord(e)
            m=m+1
            ShapH(m) = MubE(e,1)*MubE(e,2)*phiE(i)
            GradH(1:N,m) = MubE(e,1)*MubE(e,2)*DphiE(1:N,i)
     .                   + MubE(e,1)*DMubE(e,1:N,2)*phiE(i)
     .                   + DMubE(e,1:N,1)*MubE(e,2)*phiE(i)
          enddo
        endif
      enddo
c
c  ...Third the faces
c  ...call the blending and projections
      call BlendProjectHexaF(Mu,DMu, MubF,DMubF,MupF,DMupF,IdecF)
      do f=1,6
        call decod(Nord(12+f),MODORDER,2, nordF)
        ndofF = (nordF(1)-1)*(nordF(2)-1)
        if (ndofF.gt.0) then        
c      ...orient first
          call OrientQuad(MupF(f,1:2,0:1),DMupF(f,1:N,1:2,0:1),
     .                         NoriF(f),IdecF,N, GMupF,GDMupF,GIdecF)
c      ...orders already take into account the orientations, so
c      ...no need for swapping nordF
c      ...now construct the shape functions
          call AncPhiQuad(GMupF,GDMupF,nordF,GIdecF,N, 
     .                          phiQuad(2:nordF(1),2:nordF(2)),
     .                           DphiQuad(1:N,2:nordF(1),2:nordF(2)))
          do j=2,nordF(2)
            do i=2,nordF(1)
              m=m+1
              ShapH(m) = MubF(f)*phiQuad(i,j)
              GradH(1:N,m) = MubF(f)*DphiQuad(1:N,i,j)
     .                     + DMubF(f,1:N)*phiQuad(i,j)
            enddo
          enddo
        endif
      enddo
c
c  ...Finally the bubbles
c  ...find order
      call decod(Nord(19),MODORDER,3, nordB)
      ndofB = (nordB(1)-1)*(nordB(2)-1)*(nordB(3)-1)
      IdecB(1) = .TRUE.; IdecB(2) = .TRUE.; IdecB(3) = .TRUE.
c  ...if necessary, create bubbles
      if (ndofB.gt.0) then
c    ...call phiQuad and phiE - no need to orient
        call AncPhiQuad(Mu(1:2,0:1),DMu(1:N,1:2,0:1),nordB(1:2),
     .                 IdecB(1:2),N, phiQuad(2:nordB(1),2:nordB(2)),
     .                            DphiQuad(1:N,2:nordB(1),2:nordB(2)))
        call AncPhiE(Mu(3,0:1),DMu(1:N,3,0:1),nordB(3),IdecB(3),N, 
     .                        phiE(2:nordB(3)),DphiE(1:N,2:nordB(3)))
        do k=2,nordB(3)
          do j=2,nordB(2)
            do i=2,nordB(1)
              m=m+1
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
c  ...print this when debugging
      if (iprint.ge.1) then
        write(*,7001) Xi(1:3),Nord(1:19),NoriE(1:12),NoriF(1:6),NrdofH
 7001   format('shape3DHHexa: Xi = ',3f8.3,/,
     .         'Norder = ',3(4i2,2x),2i3,2x,4i3,3x,i4,/,
     .         'orient = ',3(4i2,2x),2i3,2x,4i3,/,'NrdofH = ',i3)
        write(*,7010)
 7010   format('VERTEX SHAPE FUNCTIONS = ')
        do v=1,8
          m=v
          write(*,7002) m,ShapH(m),GradH(1:3,m)
 7002     format('k = ',i3,' ShapH, GradH = ',e12.5,3x,3e12.5)
        enddo
        do e=1,12
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
        do f=1,6
          call decod(Nord(12+f),MODORDER,2, nordF)
          ndofF = (nordF(1)-1)*(nordF(2)-1)
          if (ndofF.gt.0) then
            write(*,7012) f
 7012       format('SHAPE FUNCTIONS FOR FACE = ',i2)
            do j=1,ndofF
              m=m+1
              write(*,7002) m,ShapH(m),GradH(1:3,m)
            enddo
          endif
        enddo
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
c
      end subroutine shape3DHHexa
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape3DEHexa
c
c--------------------------------------------------------------------
c
c     latest revision:  - Oct 14, Apr 17
c
c     purpose:          - routine returns values of 3D hexahedron element
c                         H(curl) shape functions and their derivatives
c
c     arguments:
c
c     in:
c          X            - master hexahedron coordinates from (0,1)^3
c          Nord         - polynomial order for the nodes (H1 sense)
c          NoriE        - edge orientation
c          NoriF        - face orientation
c          Nsize        - relevant sizes of local arrays
c
c     out:
c          NrdofE       - number of dof
c          ShapE        - values of the shape functions at the point
c          CurlE        - cur lof the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3DEHexa(Xi,Nord,NoriE,NoriF,Nsize,
     .                                              NrdofE,ShapE,CurlE)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: Nord(19),NoriE(12),NoriF(6),Nsize(2)
      integer, intent(out) :: NrdofE
      integer :: i,j,k,ig,jg,kg,m,e,f,fam,a,b,c,ab(2),abc(3),N,iprint
      integer :: ndofE,nordF(2),ndofF(0:1),minF(2),maxF(2),ij(2)
      integer :: nordB(3),ndofB(0:2),minB(3),maxB(3),ijk(3)
      logical :: IdecE,IdecF(2),GIdecF(2),IdecB(3)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapE(1:3,Nsize(2))
      double precision, intent(out) :: CurlE(1:3,Nsize(2))
      double precision :: Mu(1:3,0:1),DMu(1:3,1:3,0:1)
      double precision :: MubE(1:12,1:2),DMubE(1:12,1:3,1:2)
      double precision :: MupE(1:12,0:1),DMupE(1:12,1:3,0:1)
      double precision :: GMupE(0:1),GDMupE(1:3,0:1)
      double precision :: MubF(1:6),DMubF(1:6,1:3)
      double precision :: MupF(1:6,1:2,0:1),DMupF(1:6,1:3,1:2,0:1)
      double precision :: GMupF(1:2,0:1),GDMupF(1:3,1:2,0:1)
      double precision :: EE(1:3,0:Nsize(1)-1),curlEE(1:3,0:Nsize(1)-1)
      double precision :: EQuad(1:3,0:Nsize(1)-1,2:Nsize(1))
      double precision :: curlEQuad(1:3,0:Nsize(1)-1,2:Nsize(1))
      double precision :: phiE(2:Nsize(1)),DphiE(1:3,2:Nsize(1))
      double precision :: DTemp(1:3),CTemp(1:3)
c
c  ...debugging flag
      iprint=0
c  ...spatial dimensions
      N=3
c  ...initiate counter for shape functions
      m=0
c
c  ...Define affine coordinates
      call AffineHexahedron(Xi, Mu,DMu)
c
c  ...First the edges
c  ...call the blending and projections
      call BlendProjectHexaE(Mu,DMu, MubE,DMubE,MupE,DMupE,IdecE)
      do e=1,12
        ndofE=Nord(e)
        if (ndofE.gt.0) then        
c      ...orient first
          call OrientE(MupE(e,0:1),DMupE(e,1:N,0:1),NoriE(e),N, 
     .                                                 GMupE,GDMupE)
c      ...construct the shape functions
          call AncEE(GMupE,GDMupE,Nord(e),IdecE,N, 
     .                    EE(1:N,0:Nord(e)-1),curlEE(1:N,0:Nord(e)-1))
          do i=0,Nord(e)-1
            m=m+1
            DTemp = MubE(e,1)*DMubE(e,1:N,2)+DMubE(e,1:N,1)*MubE(e,2)
            call cross(N,DTemp,EE(1:N,i), CTemp)
            ShapE(1:N,m) = MubE(e,1)*MubE(e,2)*EE(1:N,i)
            CurlE(1:N,m) = CTemp
          enddo
        endif
      enddo
c
c  ...Second the faces
c  ...call the blending and projections
      call BlendProjectHexaF(Mu,DMu, MubF,DMubF,MupF,DMupF,IdecF)
      do f=1,6
c    ...find order
c    ...these already account for orientations
        call decod(Nord(12+f),MODORDER,2, nordF)
c    ...orient the variables first (except the order)
        call OrientQuad(MupF(f,1:2,0:1),DMupF(f,1:N,1:2,0:1),
     .                         NoriF(f),IdecF,N, GMupF,GDMupF,GIdecF)
c    ...loop over the two families
        do fam=0,1
c      ...get the (global) face axis indexing for the family (a,b)
c      ...fam=0->(1,2), fam=1->(2,1)
          ab = cshift((/1,2/),fam); 
          a = ab(1); b = ab(2)
c      ...degrees of freedom (dof) for this family
          ndofF(fam) = nordF(a)*(nordF(b)-1)
c      ...now construct the shape functions if necessary
          if (ndofF(fam).gt.0) then
            call AncEQuad(GMupF(ab,0:1),GDMupF(1:N,ab,0:1),
     .                       nordF(ab),GIdecF(ab),N, 
     .                        EQuad(1:N,0:nordF(a)-1,2:nordF(b)),
     .                         curlEQuad(1:N,0:nordF(a)-1,2:nordF(b)))
c        ...the following manipulations are necessary due to
c        ...some conventions in the code: the outer loop always is
c        ...numbered wrt the second global face axis
            minF(1) = 0; minF(2) = 2
            maxF(1) = nordF(a)-1; maxF(2) = nordF(b)
            minF = cshift(minF,-fam); maxF = cshift(maxF,-fam)
            do jg=minF(2),maxF(2)
              do ig=minF(1),maxF(1)
                ij = cshift((/ig,jg/),fam); 
                i = ij(1); j = ij(2)
                m=m+1
                call cross(N,DMubF(f,1:N),EQuad(1:N,i,j), CTemp)
                ShapE(1:N,m) = MubF(f)*EQuad(1:N,i,j)
                CurlE(1:N,m) = MubF(f)*curlEQuad(1:N,i,j)+CTemp
              enddo
            enddo
          endif
        enddo
      enddo
c
c  ...Finally the bubbles
c  ...find order
      call decod(Nord(19),MODORDER,3, nordB)
      IdecB(1) = .TRUE.; IdecB(2) = .TRUE.; IdecB(3) = .TRUE.
c  ...loop over the three families
      do fam=0,2
c    ...get the interior axis indexing for the family (a,b,c)
c    ...fam=0->(1,2,3), fam=1->(2,3,1), fam=2->(3,1,2)
        abc = cshift((/1,2,3/),fam); 
        a = abc(1); b = abc(2); c = abc(3); ab(1) = a; ab(2) = b
c    ...degrees of freedom (dof) for this family
        ndofB(fam) = nordB(a)*(nordB(b)-1)*(nordB(c)-1)
c    ...create the bubbles for this family if necessary
        if (ndofB(fam).gt.0) then
c      ...call EQuad and phiE with appropriate indexing
          call AncEQuad(Mu(ab,0:1),DMu(1:N,ab,0:1),
     .                                nordB(ab),IdecB(ab),N,
     .                     EQuad(1:N,0:nordB(a)-1,2:nordB(b)),
     .                   curlEQuad(1:N,0:nordB(a)-1,2:nordB(b)))
          call AncPhiE(Mu(c,0:1),DMu(1:N,c,0:1),nordB(c),IdecB(c),N, 
     .                        phiE(2:nordB(c)),DphiE(1:N,2:nordB(c)))
c      ...the following manipulations are necessary due to
c      ...some conventions in the code: the outer loop always wrt the 
c      ...third axis, the inner loop wrt to the first axis.
          minB(1) = 0; minB(2) = 2; minB(3) = 2
          maxB(1) = nordB(a)-1; maxB(2) = nordB(b); maxB(3) = nordB(c)
          minB = cshift(minB,-fam); maxB = cshift(maxB,-fam)
          do kg=minB(3),maxB(3)
            do jg=minB(2),maxB(2)
              do ig=minB(1),maxB(1) 
                ijk = cshift((/ig,jg,kg/),fam); 
                i = ijk(1); j = ijk(2); k = ijk(3)
                m=m+1
                call cross(N,DphiE(1:N,k),EQuad(1:N,i,j), CTemp)
                ShapE(1:N,m) = EQuad(1:N,i,j)*phiE(k)
                CurlE(1:N,m) = phiE(k)*curlEQuad(1:N,i,j)+CTemp
              enddo
            enddo
          enddo
        endif
      enddo
c
c  ...give total degrees of freedom
      NrdofE = m
c
c  ...print this when debugging
      if (iprint.ge.1) then
        write(*,7001) Xi(1:3),Nord(1:19),
     .                NoriE(1:12),NoriF(1:6),NrdofE
 7001   format('shape3DEHexa: Xi = ',3f8.3,/,
     .         'Norder = ',3(4i2,2x),2i3,2x,4i3,3x,i4,/,
     .         'orient = ',3(4i2,2x),2i3,2x,4i3,/,'NrdofE = ',i3)
        m=0
        do e=1,12
          ndofE = Nord(e)
          if (ndofE.gt.0) then
            write(*,7011) e
 7011       format('SHAPE FUNCTIONS FOR EDGE = ',i2)
            do j=1,ndofE
              m=m+1
              write(*,7003) m,ShapE(1:N,m),CurlE(1:N,m)
 7003         format('k = ',i3,' ShapE, CurlE = ',3e12.5,3x,3e12.5)
            enddo
          endif
        enddo
        do f=1,6
          call decod(Nord(12+f),MODORDER,2, nordF)
          ndofF(0) = nordF(1)*(nordF(2)-1)
          ndofF(1) = (nordF(1)-1)*nordF(2)
          if ((ndofF(0)+ndofF(1)).gt.0) then
            write(*,7012) f
 7012       format('SHAPE FUNCTIONS FOR FACE = ',i2)
            do j=1,ndofF(0)+ndofF(1)
              m=m+1
              write(*,7003) m,ShapE(1:N,m),CurlE(1:N,m)
            enddo
          endif
        enddo
        if ((ndofB(0)+ndofB(1)+ndofB(2)).gt.0) then
          write(*,7013)
 7013     format('BUBBLES = ')
          do j=1,ndofB(0)+ndofB(1)+ndofB(2)
            m=m+1
            write(*,7003) m,ShapE(1:N,m),CurlE(1:N,m)
          enddo
        endif
        call pause
      endif
c
      end subroutine shape3DEHexa
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape3DVHexa
c
c--------------------------------------------------------------------
c
c     latest revision:  - Oct 14, Apr 17
c
c     purpose:          - routine returns values of 3D hexahedron element
c                         H(div) shape functions and their divergences
c
c     arguments:
c
c     in:
c          X            - master hexahedron coordinates from (0,1)^3
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
      subroutine shape3DVHexa(Xi,Nord,NoriF,Nsize, NrdofV,ShapV,DivV)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: Nord(19),NoriF(6),Nsize(2)
      integer, intent(out) :: NrdofV
      integer :: i,j,k,ig,jg,kg,m,f,fam,a,b,c,ab(2),abc(3),N,iprint
      integer :: nordF(2),ndofF
      integer :: nordB(3),ndofB(0:2),minB(3),maxB(3),ijk(3)
      logical :: IdecF(2),GIdecF(2),IdecB(3)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapV(1:3,Nsize(2))
      double precision, intent(out) :: DivV(Nsize(2))
      double precision :: Mu(1:3,0:1),DMu(1:3,1:3,0:1)
      double precision :: MubF(1:6),DMubF(1:6,1:3)
      double precision :: MupF(1:6,1:2,0:1),DMupF(1:6,1:3,1:2,0:1)
      double precision :: GMupF(1:2,0:1),GDMupF(1:3,1:2,0:1)
      double precision :: VQuad(1:3,0:Nsize(1)-1,0:Nsize(1)-1)
      double precision :: divVQuad(0:Nsize(1)-1,0:Nsize(1)-1)
      double precision :: phiE(2:Nsize(1)),DphiE(1:3,2:Nsize(1))
c
c  ...debugging flag
      iprint=0
c  ...spatial dimensions
      N=3
c  ...initiate counter for shape functions
      m=0
c
c  ...Define affine coordinates
      call AffineHexahedron(Xi, Mu,DMu)
c
c  ...First the faces
c  ...call the blending and projections
      call BlendProjectHexaF(Mu,DMu, MubF,DMubF,MupF,DMupF,IdecF)
      do f=1,6
c    ...find order
c    ...these already account for orientations
        call decod(Nord(12+f),MODORDER,2, nordF)
c    ...orient the variables first (except the order)
        call OrientQuad(MupF(f,1:2,0:1),DMupF(f,1:N,1:2,0:1),
     .                         NoriF(f),IdecF,N, GMupF,GDMupF,GIdecF)
        ndofF = nordF(1)*nordF(2)
c    ...now construct the shape functions if necessary
        if (ndofF.gt.0) then
          call AncVQuad(GMupF,GDMupF,nordF,GIdecF,N, 
     .                   VQuad(1:N,0:nordF(1)-1,0:nordF(2)-1),
     .                       divVQuad(0:nordF(1)-1,0:nordF(2)-1))
          do j=0,nordF(2)-1
            do i=0,nordF(1)-1
              m=m+1
              ShapV(1:N,m) = MubF(f)*VQuad(1:N,i,j)
              call dot_product(DMubF(f,1:N),VQuad(1:N,i,j), DivV(m))
            enddo
          enddo
        endif
      enddo
c
c  ...Finally the bubbles
c  ...find order
      call decod(Nord(19),MODORDER,3, nordB)
      IdecB(1) = .TRUE.; IdecB(2) = .TRUE.; IdecB(3) = .TRUE.
c  ...loop over the three families
      do fam=0,2
c    ...get the interior axis indexing for the family (a,b,c)
c    ...fam=0->(1,2,3), fam=1->(2,3,1), fam=2->(3,1,2)
        abc = cshift((/1,2,3/),fam); 
        a = abc(1); b = abc(2); c = abc(3); ab(1) = a; ab(2) = b
c    ...degrees of freedom (dof) for this family
        ndofB(fam) = nordB(a)*nordB(b)*(nordB(c)-1)
c    ...create the bubbles for this family if necessary
        if (ndofB(fam).gt.0) then
c      ...call VQuad and phiE with appropriate indexing
          call AncVQuad(Mu(ab,0:1),DMu(1:N,ab,0:1),
     .                                nordB(ab),IdecB(ab),N,
     .                   VQuad(1:N,0:nordB(a)-1,0:nordB(b)-1),
     .                         divVQuad(0:nordB(a)-1,0:nordB(b)-1))
          call AncPhiE(Mu(c,0:1),DMu(1:N,c,0:1),nordB(c),IdecB(c),N, 
     .                        phiE(2:nordB(c)),DphiE(1:N,2:nordB(c)))
c      ...the following manipulations are necessary due to
c      ...some conventions in the code: the outer loop always wrt the 
c      ...third axis, the inner loop wrt to the first axis.
          minB(1) = 0; minB(2) = 0; minB(3) = 2
          maxB(1) = nordB(a)-1; maxB(2) = nordB(b)-1; maxB(3) = nordB(c)
          minB = cshift(minB,-fam); maxB = cshift(maxB,-fam)
          do kg=minB(3),maxB(3)
            do jg=minB(2),maxB(2)
              do ig=minB(1),maxB(1)
                ijk = cshift((/ig,jg,kg/),fam); 
                i = ijk(1); j = ijk(2); k = ijk(3)
                m=m+1
                ShapV(1:N,m) = phiE(k)*VQuad(1:N,i,j)
                call dot_product(DphiE(1:N,k),VQuad(1:N,i,j), DivV(m))
              enddo
            enddo
          enddo
        endif
      enddo
c
c  ...give total degrees of freedom
      NrdofV = m
c
c  ...print this when debugging
      if (iprint.ge.1) then
        write(*,7001) Xi(1:3),Nord(13:19),NoriF(1:6),NrdofV
 7001   format('shape3DVHexa: Xi = ',3f8.3,/,
     .         'Norder = ',2i3,2x,4i3,3x,i4,/,
     .         'orient = ',2i3,2x,4i3,/,'NrdofV = ',i3)
        m=0
        do f=1,6
          call decod(Nord(12+f),MODORDER,2, nordF)
          ndofF = nordF(1)*nordF(2)
          if (ndofF.gt.0) then
            write(*,7012) f
 7012       format('SHAPE FUNCTIONS FOR FACE = ',i2)
            do j=1,ndofF
              m=m+1
              write(*,7003) m,ShapV(1:N,m),DivV(m)
 7003         format('k = ',i3,' ShapV, DivV = ',3e12.5,3x,e12.5)
            enddo
          endif
        enddo
        if ((ndofB(0)+ndofB(1)+ndofB(2)).gt.0) then
          write(*,7013)
 7013     format('BUBBLES = ')
          do j=1,ndofB(0)+ndofB(1)+ndofB(2)
            m=m+1
            write(*,7003) m,ShapV(1:N,m),DivV(m)
          enddo
        endif
        call pause
      endif
c
      end subroutine shape3DVHexa
c
c
c--------------------------------------------------------------------
c
c     routine name      - shape3DQHexa
c
c--------------------------------------------------------------------
c
c     latest revision:  - Oct 14, Apr 17
c
c     purpose:          - routine returns values of 3D hexahedron
c                         element L2 shape functions
c
c     arguments:
c
c     in:
c          X            - master hexahedron coordinates from (0,1)^3
c          Nord         - polynomial order for the nodes (H1 sense)
c          Nsize        - relevant sizes of local arrays
c
c     out:
c          NrdofQ       - number of dof
c          ShapQ        - values of the shape functions at the point
c
c-----------------------------------------------------------------------
c
      subroutine shape3DQHexa(Xi,Nord,Nsize, NrdofQ,ShapQ)
c
      use parameters , only : MODORDER
c
      implicit none
      integer, intent(in)  :: Nord(19),Nsize(2)
      integer, intent(out) :: NrdofQ
      integer :: i,j,k,m,N,nordB(3),ndofB,iprint
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapQ(Nsize(2))
      double precision :: Mu(1:3,0:1),DMu(1:3,1:3,0:1)
      double precision :: homP(1:3,0:Nsize(1)-1)
c
c  ...debugging flag
      iprint=0
c  ...spatial dimensions
      N=3
c  ...initiate counter for shape functions
      m=0
c
c  ...Define affine coordinates
      call AffineHexahedron(Xi, Mu,DMu)
c
c  ...There are only bubbles
c  ...find order
      call decod(Nord(19),MODORDER,3, nordB)
      ndofB = nordB(1)*nordB(2)*nordB(3)
c  ...if necessary, create bubbles - always necessary if p,q,r>=1
      if (ndofB.gt.0) then
c    ...call Legendre polynomials - no need to orient
        call HomLegendre(Mu(1,0:1),nordB(1)-1, homP(1,0:nordB(1)-1))
        call HomLegendre(Mu(2,0:1),nordB(2)-1, homP(2,0:nordB(2)-1))
        call HomLegendre(Mu(3,0:1),nordB(3)-1, homP(3,0:nordB(3)-1))
        do k=0,nordB(3)-1
          do j=0,nordB(2)-1
            do i=0,nordB(1)-1
              m=m+1
              ShapQ(m) = homP(1,i)*homP(2,j)*homP(3,k)
            enddo
          enddo
        enddo
      endif
c
c  ...give total degrees of freedom
      NrdofQ = m
c
c  ...print this when debugging
      if (iprint.ge.1) then
        write(*,7001) Xi(1:3),Nord(19),NrdofQ
 7001   format('shap3Q_bric: Xi = ',3f8.3,' Nord = ',i3,/,
     .               'NrdofQ = ',i3)
        do m=1,NrdofQ
          write(*,7002) m,ShapQ(m)
 7002     format('k = ',i3,' ShapQ, = ',e12.5)
        enddo
        call pause
      endif
c
      end subroutine shape3DQHexa

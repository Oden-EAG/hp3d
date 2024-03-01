! Routines:
!  - shape3DHHexa
!  - shape3DEHexa
!  - shape3DVHexa
!  - shape3DQHexa
!--------------------------------------------------------------------
!
!     routine name      - shape3DHHexa
!
!--------------------------------------------------------------------
!
!     latest revision:  - Oct 14, Apr 17, Jul 21
!
!> @brief         - routine returns values of 3D hexahedron element
!                         H1 shape functions and their derivatives
!
!     arguments:
!
!     in:
!          X            - master hexahedron coordinates from (0,1)^3
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
   subroutine shape3DHHexa(Xi,Nord,NoriE,NoriF,Nsize, &
                                                    NrdofH,ShapH,GradH)
!
      use parameters , only : MODORDER
!
      implicit none
      integer, intent(in)  :: Nord(19),NoriE(12),NoriF(6),Nsize(2)
      integer, intent(out) :: NrdofH
      integer :: i,j,k,m,v,e,f,N,ndofE,nordF(2),ndofF
      integer :: nordB(3),ndofB
      logical :: IdecE,IdecF(2),GIdecF(2),IdecB(3)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapH(Nsize(2))
      double precision, intent(out) :: GradH(1:3,Nsize(2))
      double precision :: Mu(0:1,1:3),DMu(1:3,0:1,1:3)
      double precision :: MubV(1:3,1:8),DMubV(1:3,1:3,1:8)
      double precision :: MubE(1:2,1:12),DMubE(1:3,1:2,1:12)
      double precision :: MupE(0:1,1:12),DMupE(1:3,0:1,1:12)
      double precision :: GMupE(0:1),GDMupE(1:3,0:1)
      double precision :: MubF(1:6),DMubF(1:3,1:6)
      double precision :: MupF(0:1,1:2,1:6),DMupF(1:3,0:1,1:2,1:6)
      double precision :: GMupF(1:2,0:1),GDMupF(1:3,1:2,0:1)
      double precision :: phiE(2:Nsize(1)),DphiE(1:3,2:Nsize(1))
      double precision :: phiQuad(2:Nsize(1),2:Nsize(1))
      double precision :: DphiQuad(1:3,2:Nsize(1),2:Nsize(1))
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
!  ...Define affine coordinates
      call AffineHexahedron(Xi, Mu,DMu)
!
!  ...First the vertices
!  ...call the blending functions
      call BlendHexaV(Mu,DMu, MubV,DMubV)
      do v=1,8
        m=m+1
        ShapH(m) = MubV(1,v)*MubV(2,v)*MubV(3,v)
        GradH(1:N,m) = MubV(1,v)*MubV(2,v)*DMubV(1:N,3,v) &
                     + MubV(1,v)*DMubV(1:N,2,v)*MubV(3,v) &
                     + DMubV(1:N,1,v)*MubV(2,v)*MubV(3,v)
      enddo
!
!  ...Second the edges
!  ...call the blending and projections
      call BlendProjectHexaE(Mu,DMu, MubE,DMubE,MupE,DMupE,IdecE)
      do e=1,12
        ndofE=Nord(e)-1
        if (ndofE.gt.0) then
!      ...orient first
          call OrientE(MupE(0:1,e),DMupE(1:N,0:1,e),NoriE(e),N, &
                                                       GMupE,GDMupE)
!      ...construct the shape functions
          call AncPhiE(GMupE,GDMupE,Nord(e),IdecE,N, &
                               phiE(2:Nord(e)),DphiE(1:N,2:Nord(e)))
          do i=2,Nord(e)
            m=m+1
            ShapH(m) = MubE(1,e)*MubE(2,e)*phiE(i)
            GradH(1:N,m) = MubE(1,e)*MubE(2,e)*DphiE(1:N,i) &
                         + MubE(1,e)*DMubE(1:N,2,e)*phiE(i) &
                         + DMubE(1:N,1,e)*MubE(2,e)*phiE(i)
          enddo
        endif
      enddo
!
!  ...Third the faces
!  ...call the blending and projections
      call BlendProjectHexaF(Mu,DMu, MubF,DMubF,MupF,DMupF,IdecF)
      do f=1,6
        call decod(Nord(12+f),MODORDER,2, nordF)
        ndofF = (nordF(1)-1)*(nordF(2)-1)
        if (ndofF.gt.0) then
!      ...orient first
          call OrientQuad(MupF(0:1,1:2,f),DMupF(1:N,0:1,1:2,f), &
                               NoriF(f),IdecF,N, GMupF,GDMupF,GIdecF)
!      ...orders already take into account the orientations, so
!      ...no need for swapping nordF
!      ...now construct the shape functions
          call AncPhiQuad(GMupF,GDMupF,nordF,GIdecF,N, &
                                phiQuad(2:nordF(1),2:nordF(2)), &
                                 DphiQuad(1:N,2:nordF(1),2:nordF(2)))
          do j=2,nordF(2)
            do i=2,nordF(1)
              m=m+1
              ShapH(m) = MubF(f)*phiQuad(i,j)
              GradH(1:N,m) = MubF(f)*DphiQuad(1:N,i,j) &
                           + DMubF(1:N,f)*phiQuad(i,j)
            enddo
          enddo
        endif
      enddo
!
!  ...Finally the bubbles
!  ...find order
      call decod(Nord(19),MODORDER,3, nordB)
      ndofB = (nordB(1)-1)*(nordB(2)-1)*(nordB(3)-1)
      IdecB(1) = .true.; IdecB(2) = .true.; IdecB(3) = .true.
!  ...if necessary, create bubbles
      if (ndofB.gt.0) then
!    ...call phiQuad and phiE - no need to orient
        call AncPhiQuad(Mu(0:1,1:2),DMu(1:N,0:1,1:2),nordB(1:2), &
                       IdecB(1:2),N, phiQuad(2:nordB(1),2:nordB(2)), &
                                  DphiQuad(1:N,2:nordB(1),2:nordB(2)))
        call AncPhiE(Mu(0:1,3),DMu(1:N,0:1,3),nordB(3),IdecB(3),N, &
                              phiE(2:nordB(3)),DphiE(1:N,2:nordB(3)))
        do k=2,nordB(3)
          do j=2,nordB(2)
            do i=2,nordB(1)
              m=m+1
              ShapH(m) = phiQuad(i,j)*phiE(k)
              GradH(1:N,m) = phiQuad(i,j)*DphiE(1:N,k) &
                           + DphiQuad(1:N,i,j)*phiE(k)
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
        write(*,7001) Xi(1:3),Nord(1:19),NoriE(1:12),NoriF(1:6),NrdofH
 7001   format('shape3DHHexa: Xi = ',3f8.3,/, &
               'Norder = ',3(4i2,2x),2i3,2x,4i3,3x,i4,/, &
               'orient = ',3(4i2,2x),2i3,2x,4i3,/,'NrdofH = ',i3)
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
#endif
!
   end subroutine shape3DHHexa
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape3DEHexa
!
!--------------------------------------------------------------------
!
!     latest revision:  - Oct 14, Apr 17, Jul 21
!
!> @brief         - routine returns values of 3D hexahedron element
!                         H(curl) shape functions and their derivatives
!
!     arguments:
!
!     in:
!          X            - master hexahedron coordinates from (0,1)^3
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
   subroutine shape3DEHexa(Xi,Nord,NoriE,NoriF,Nsize, &
                                                    NrdofE,ShapE,CurlE)
!
      use parameters , only : MODORDER
!
      implicit none
      integer, intent(in)  :: Nord(19),NoriE(12),NoriF(6),Nsize(2)
      integer, intent(out) :: NrdofE
      integer :: i,j,k,ig,jg,kg,m,e,f,fam,a,b,c,ab(2),abc(3),N
      integer :: ndofE,nordF(2),ndofF(0:1),minF(2),maxF(2),ij(2)
      integer :: nordB(3),ndofB(0:2),minB(3),maxB(3),ijk(3)
      logical :: IdecE,IdecF(2),GIdecF(2),IdecB(3)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapE(1:3,Nsize(2))
      double precision, intent(out) :: CurlE(1:3,Nsize(2))
      double precision :: Mu(0:1,1:3),DMu(1:3,0:1,1:3)
      double precision :: MubE(1:2,1:12),DMubE(1:3,1:2,1:12)
      double precision :: MupE(0:1,1:12),DMupE(1:3,0:1,1:12)
      double precision :: GMupE(0:1),GDMupE(1:3,0:1)
      double precision :: MubF(1:6),DMubF(1:3,1:6)
      double precision :: MupF(0:1,1:2,1:6),DMupF(1:3,0:1,1:2,1:6)
      double precision :: GMupF(0:1,1:2),GDMupF(1:3,0:1,1:2)
      double precision :: EE(1:3,0:Nsize(1)-1),curlEE(1:3,0:Nsize(1)-1)
      double precision :: EQuad(1:3,0:Nsize(1)-1,2:Nsize(1))
      double precision :: curlEQuad(1:3,0:Nsize(1)-1,2:Nsize(1))
      double precision :: phiE(2:Nsize(1)),DphiE(1:3,2:Nsize(1))
      double precision :: DTemp(1:3),CTemp(1:3)
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
!  ...Define affine coordinates
      call AffineHexahedron(Xi, Mu,DMu)
!
!  ...First the edges
!  ...call the blending and projections
      call BlendProjectHexaE(Mu,DMu, MubE,DMubE,MupE,DMupE,IdecE)
      do e=1,12
        ndofE=Nord(e)
        if (ndofE.gt.0) then
!      ...orient first
          call OrientE(MupE(0:1,e),DMupE(1:N,0:1,e),NoriE(e),N, &
                                                       GMupE,GDMupE)
!      ...construct the shape functions
          call AncEE(GMupE,GDMupE,Nord(e),IdecE,N, &
                          EE(1:N,0:Nord(e)-1),curlEE(1:N,0:Nord(e)-1))
          do i=0,Nord(e)-1
            m=m+1
            DTemp = MubE(1,e)*DMubE(1:N,2,e)+DMubE(1:N,1,e)*MubE(2,e)
            call cross(N,DTemp,EE(1:N,i), CTemp)
            ShapE(1:N,m) = MubE(1,e)*MubE(2,e)*EE(1:N,i)
            CurlE(1:N,m) = CTemp
          enddo
        endif
      enddo
!
!  ...Second the faces
!  ...call the blending and projections
      call BlendProjectHexaF(Mu,DMu, MubF,DMubF,MupF,DMupF,IdecF)
      do f=1,6
!    ...find order
!    ...these already account for orientations
        call decod(Nord(12+f),MODORDER,2, nordF)
!    ...orient the variables first (except the order)
        call OrientQuad(MupF(0:1,1:2,f),DMupF(1:N,0:1,1:2,f), &
                               NoriF(f),IdecF,N, GMupF,GDMupF,GIdecF)
!    ...loop over the two families
        do fam=0,1
!      ...get the (global) face axis indexing for the family (a,b)
!      ...fam=0->(1,2), fam=1->(2,1)
          ab = cshift((/1,2/),fam);
          a = ab(1); b = ab(2)
!      ...degrees of freedom (dof) for this family
          ndofF(fam) = nordF(a)*(nordF(b)-1)
!      ...now construct the shape functions if necessary
          if (ndofF(fam).gt.0) then
            call AncEQuad(GMupF(0:1,ab),GDMupF(1:N,0:1,ab), &
                             nordF(ab),GIdecF(ab),N, &
                              EQuad(1:N,0:nordF(a)-1,2:nordF(b)), &
                               curlEQuad(1:N,0:nordF(a)-1,2:nordF(b)))
!        ...the following manipulations are necessary due to
!        ...some conventions in the code: the outer loop always is
!        ...numbered wrt the second global face axis
            minF(1) = 0; minF(2) = 2
            maxF(1) = nordF(a)-1; maxF(2) = nordF(b)
            minF = cshift(minF,-fam); maxF = cshift(maxF,-fam)
            do jg=minF(2),maxF(2)
              do ig=minF(1),maxF(1)
                ij = cshift((/ig,jg/),fam);
                i = ij(1); j = ij(2)
                m=m+1
                call cross(N,DMubF(1:N,f),EQuad(1:N,i,j), CTemp)
                ShapE(1:N,m) = MubF(f)*EQuad(1:N,i,j)
                CurlE(1:N,m) = MubF(f)*curlEQuad(1:N,i,j)+CTemp
              enddo
            enddo
          endif
        enddo
      enddo
!
!  ...Finally the bubbles
!  ...find order
      call decod(Nord(19),MODORDER,3, nordB)
      IdecB(1) = .true.; IdecB(2) = .true.; IdecB(3) = .true.
!  ...loop over the three families
      do fam=0,2
!    ...get the interior axis indexing for the family (a,b,c)
!    ...fam=0->(1,2,3), fam=1->(2,3,1), fam=2->(3,1,2)
        abc = cshift((/1,2,3/),fam);
        a = abc(1); b = abc(2); c = abc(3); ab(1) = a; ab(2) = b
!    ...degrees of freedom (dof) for this family
        ndofB(fam) = nordB(a)*(nordB(b)-1)*(nordB(c)-1)
!    ...create the bubbles for this family if necessary
        if (ndofB(fam).gt.0) then
!      ...call EQuad and phiE with appropriate indexing
          call AncEQuad(Mu(0:1,ab),DMu(1:N,0:1,ab), &
                                      nordB(ab),IdecB(ab),N, &
                           EQuad(1:N,0:nordB(a)-1,2:nordB(b)), &
                         curlEQuad(1:N,0:nordB(a)-1,2:nordB(b)))
          call AncPhiE(Mu(0:1,c),DMu(1:N,0:1,c),nordB(c),IdecB(c),N, &
                              phiE(2:nordB(c)),DphiE(1:N,2:nordB(c)))
!      ...the following manipulations are necessary due to
!      ...some conventions in the code: the outer loop always wrt the
!      ...third axis, the inner loop wrt to the first axis.
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
!
!  ...give total degrees of freedom
      NrdofE = m
!
#if HP3D_DEBUG
!  ...print this when debugging
      if (iprint.ge.1) then
        write(*,7001) Xi(1:3),Nord(1:19), &
                      NoriE(1:12),NoriF(1:6),NrdofE
 7001   format('shape3DEHexa: Xi = ',3f8.3,/, &
               'Norder = ',3(4i2,2x),2i3,2x,4i3,3x,i4,/, &
               'orient = ',3(4i2,2x),2i3,2x,4i3,/,'NrdofE = ',i3)
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
#endif
!
   end subroutine shape3DEHexa
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape3DVHexa
!
!--------------------------------------------------------------------
!
!     latest revision:  - Oct 14, Apr 17, Jul 21
!
!> @brief         - routine returns values of 3D hexahedron element
!                         H(div) shape functions and their divergences
!
!     arguments:
!
!     in:
!          X            - master hexahedron coordinates from (0,1)^3
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
   subroutine shape3DVHexa(Xi,Nord,NoriF,Nsize, NrdofV,ShapV,DivV)
!
      use parameters , only : MODORDER
!
      implicit none
      integer, intent(in)  :: Nord(19),NoriF(6),Nsize(2)
      integer, intent(out) :: NrdofV
      integer :: i,j,k,ig,jg,kg,m,f,fam,a,b,c,ab(2),abc(3),N
      integer :: nordF(2),ndofF
      integer :: nordB(3),ndofB(0:2),minB(3),maxB(3),ijk(3)
      logical :: IdecF(2),GIdecF(2),IdecB(3)
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapV(1:3,Nsize(2))
      double precision, intent(out) :: DivV(Nsize(2))
      double precision :: Mu(0:1,1:3),DMu(1:3,0:1,1:3)
      double precision :: MubF(1:6),DMubF(1:3,1:6)
      double precision :: MupF(0:1,1:2,1:6),DMupF(1:3,0:1,1:2,1:6)
      double precision :: GMupF(1:2,0:1),GDMupF(1:3,1:2,0:1)
      double precision :: VQuad(1:3,0:Nsize(1)-1,0:Nsize(1)-1)
      double precision :: divVQuad(0:Nsize(1)-1,0:Nsize(1)-1)
      double precision :: phiE(2:Nsize(1)),DphiE(1:3,2:Nsize(1))
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
!  ...Define affine coordinates
      call AffineHexahedron(Xi, Mu,DMu)
!
!  ...First the faces
!  ...call the blending and projections
      call BlendProjectHexaF(Mu,DMu, MubF,DMubF,MupF,DMupF,IdecF)
      do f=1,6
!    ...find order
!    ...these already account for orientations
        call decod(Nord(12+f),MODORDER,2, nordF)
!    ...orient the variables first (except the order)
        call OrientQuad(MupF(0:1,1:2,f),DMupF(1:N,0:1,1:2,f), &
                               NoriF(f),IdecF,N, GMupF,GDMupF,GIdecF)
        ndofF = nordF(1)*nordF(2)
!    ...now construct the shape functions if necessary
        if (ndofF.gt.0) then
          call AncVQuad(GMupF,GDMupF,nordF,GIdecF,N, &
                         VQuad(1:N,0:nordF(1)-1,0:nordF(2)-1), &
                             divVQuad(0:nordF(1)-1,0:nordF(2)-1))
          do j=0,nordF(2)-1
            do i=0,nordF(1)-1
              m=m+1
              ShapV(1:N,m) = MubF(f)*VQuad(1:N,i,j)
              call dot_product(DMubF(1:N,f),VQuad(1:N,i,j), DivV(m))
            enddo
          enddo
        endif
      enddo
!
!  ...Finally the bubbles
!  ...find order
      call decod(Nord(19),MODORDER,3, nordB)
      IdecB(1) = .true.; IdecB(2) = .true.; IdecB(3) = .true.
!  ...loop over the three families
      do fam=0,2
!    ...get the interior axis indexing for the family (a,b,c)
!    ...fam=0->(1,2,3), fam=1->(2,3,1), fam=2->(3,1,2)
        abc = cshift((/1,2,3/),fam);
        a = abc(1); b = abc(2); c = abc(3); ab(1) = a; ab(2) = b
!    ...degrees of freedom (dof) for this family
        ndofB(fam) = nordB(a)*nordB(b)*(nordB(c)-1)
!    ...create the bubbles for this family if necessary
        if (ndofB(fam).gt.0) then
!      ...call VQuad and phiE with appropriate indexing
          call AncVQuad(Mu(0:1,ab),DMu(1:N,0:1,ab), &
                                      nordB(ab),IdecB(ab),N, &
                         VQuad(1:N,0:nordB(a)-1,0:nordB(b)-1), &
                               divVQuad(0:nordB(a)-1,0:nordB(b)-1))
          call AncPhiE(Mu(0:1,c),DMu(1:N,0:1,c),nordB(c),IdecB(c),N, &
                              phiE(2:nordB(c)),DphiE(1:N,2:nordB(c)))
!      ...the following manipulations are necessary due to
!      ...some conventions in the code: the outer loop always wrt the
!      ...third axis, the inner loop wrt to the first axis.
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
!
!  ...give total degrees of freedom
      NrdofV = m
!
#if HP3D_DEBUG
!  ...print this when debugging
      if (iprint.ge.1) then
        write(*,7001) Xi(1:3),Nord(13:19),NoriF(1:6),NrdofV
 7001   format('shape3DVHexa: Xi = ',3f8.3,/, &
               'Norder = ',2i3,2x,4i3,3x,i4,/, &
               'orient = ',2i3,2x,4i3,/,'NrdofV = ',i3)
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
#endif
!
   end subroutine shape3DVHexa
!
!
!--------------------------------------------------------------------
!
!     routine name      - shape3DQHexa
!
!--------------------------------------------------------------------
!
!     latest revision:  - Oct 14, Apr 17
!
!> @brief         - routine returns values of 3D hexahedron
!                         element L2 shape functions
!
!     arguments:
!
!     in:
!          X            - master hexahedron coordinates from (0,1)^3
!          Nord         - polynomial order for the nodes (H1 sense)
!          Nsize        - relevant sizes of local arrays
!
!     out:
!          NrdofQ       - number of dof
!          ShapQ        - values of the shape functions at the point
!
!-----------------------------------------------------------------------
!
   subroutine shape3DQHexa(Xi,Nord,Nsize, NrdofQ,ShapQ)
!
      use parameters , only : MODORDER
!
      implicit none
      integer, intent(in)  :: Nord(19),Nsize(2)
      integer, intent(out) :: NrdofQ
      integer :: i,j,k,m,N,nordB(3),ndofB
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapQ(Nsize(2))
      double precision :: Mu(0:1,1:3),DMu(1:3,0:1,1:3)
      double precision :: homP(0:Nsize(1)-1,1:3)
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
!  ...Define affine coordinates
      call AffineHexahedron(Xi, Mu,DMu)
!
!  ...There are only bubbles
!  ...find order
      call decod(Nord(19),MODORDER,3, nordB)
      ndofB = nordB(1)*nordB(2)*nordB(3)
!  ...if necessary, create bubbles - always necessary if p,q,r>=1
      if (ndofB.gt.0) then
!    ...call Legendre polynomials - no need to orient
        call HomLegendre(Mu(0:1,1),nordB(1)-1, homP(0:nordB(1)-1,1))
        call HomLegendre(Mu(0:1,2),nordB(2)-1, homP(0:nordB(2)-1,2))
        call HomLegendre(Mu(0:1,3),nordB(3)-1, homP(0:nordB(3)-1,3))
        do k=0,nordB(3)-1
          do j=0,nordB(2)-1
            do i=0,nordB(1)-1
              m=m+1
              ShapQ(m) = homP(i,1)*homP(j,2)*homP(k,3)
            enddo
          enddo
        enddo
      endif
!
!  ...give total degrees of freedom
      NrdofQ = m
!
#if HP3D_DEBUG
!  ...print this when debugging
      if (iprint.ge.1) then
        write(*,7001) Xi(1:3),Nord(19),NrdofQ
 7001   format('shap3Q_bric: Xi = ',3f8.3,' Nord = ',i3,/, &
                     'NrdofQ = ',i3)
        do m=1,NrdofQ
          write(*,7002) m,ShapQ(m)
 7002     format('k = ',i3,' ShapQ, = ',e12.5)
        enddo
        call pause
      endif
#endif
!
   end subroutine shape3DQHexa

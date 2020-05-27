c----------------------------------------------------------------------
c> Purpose : H^1 shape functions for 3D elements (version developed by
c            P. Gatto, A.D. 2008)
c
c  @param[in ] Type         - element type
c  @param[in ] Xi           - master element coordinates
c  @param[in ] Norder       - order of approximation for the nodes
c  @param[in ] Nedge_orient - edge orientations
c  @param[in ] Nface_orient - face orientations
c  @param[out] Nrdof        - number of the element shape functions
c  @param[out] ShapH        - values of shape functions
c  @param[out] DshapH       - values of derivatives of shape functions
c                             wrt to master element coordinates
c
c  rev@Sep 14
c----------------------------------------------------------------------
c
      subroutine shape3H_0(Type,Xi,Norder,Nedge_orient,Nface_orient,
     .                     Nrdof,ShapH,DshapH)
c
      use parameters
#include "syscom.blk"
c
      character(len=4)                 , intent(in)  :: Type
      integer , dimension(12)          , intent(in)  :: Nedge_orient
      integer , dimension(6)           , intent(in)  :: Nface_orient
      integer , dimension(19)          , intent(in)  :: Norder
      real*8  , dimension(3)           , intent(in)  :: Xi
      integer ,                          intent(out) :: Nrdof
      real*8  , dimension(  MAXbrickH) , intent(out) :: ShapH
      real*8  , dimension(3,MAXbrickH) , intent(out) :: DshapH
c
c----------------------------------------------------------------------
c
      iprint=0
c
      if (iprint.eq.1) then
        write(*,7001)Type,Xi(1:3)
 7001   format(' shape3H: Type,Xi = ',a5,2x,3(e12.5,2x))
      endif
c
      select case(Type)
      case('pris','mdlp')
        call shapeHp(Xi,Norder,Nedge_orient,Nface_orient,
     .               Nrdof,ShapH,DshapH)
      case('bric','mdlb')
        call shapeHb(Xi,Norder,Nedge_orient,Nface_orient,
     .               Nrdof,ShapH,DshapH)
      case('tetr','mdln')
        call shapeHn(Xi,Norder,Nedge_orient,Nface_orient,
     .               Nrdof,ShapH,DshapH)
      case('pyra','mdld')
        call shapeHd(Xi,Norder,Nedge_orient,Nface_orient,
     .               Nrdof,ShapH,DshapH)
c
      case default
        write(*,*) 'shape3H_0: error, Type = ', Type
        stop
      endselect
c
c
      endsubroutine shape3H_0
c
c
c
c----------------------------------------------------------------------
c
c   routine name       - shapeHp
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 08
c
c   purpose            - evaluate H1 shape functions for a prismatic
c                        element
c
c   arguments :
c     in:
c         Nedge_orient - edge orientations
c         Nface_orient - face orientations
c         Xi           - master element coordinates
c         Norder       - order of approximation for the nodes
c     out:
c         Nrdof        - number of dof
c         ShapH        - values of shape functions
c         DshapH       - values of derivatives of the shape functions
c                        wrt to master element coordinates
c
c----------------------------------------------------------------------
c
      subroutine shapeHp(Xi,Norder, Nedge_orient,Nface_orient,
     .                   Nrdof,ShapH,DshapH)
c
      use parameters
      use element_data
#include "syscom.blk"
c
      dimension Nedge_orient(9), Nface_orient(5),Norder(15)
      dimension Xi(3),ShapH(MAXprismH),DshapH(3,MAXprismH)
c
c  ...derivatives of edge coordinate
      dimension dtedxi(3)
c
c  ...face coordinates
      dimension tf(2),dtfdxi(2,3)
c
c  ...edge bubble shape functions
      dimension shapHb(MAXP-1),dshapHb(MAXP-1)
c
c  ...edge kernel shape functions
      dimension shapHk(MAXP-1),dshapHk(MAXP-1)
c
c  ...triangular mid-face node bubble shape functions
      dimension shapHt(MAXmdltH),dshapHt(2,MAXmdltH)
c
c  ...affine coordinates for triangular faces
      dimension vshapt(3),dvshapt(2,3)
c
c  ...1D vertex shape functions in the xi_3 direction
      dimension vshap(2),dvshap(2)
c-----------------------------------------------------------------------
c
      iprint=0
c
c  ...affine coordinates for the triangular faces
      vshapt(1) = 1.d0 - (Xi(1) + Xi(2))
      dvshapt(1,1) = -1.d0; dvshapt(2,1) = -1.d0
      vshapt(2) = Xi(1)
      dvshapt(1,2) = 1.d0; dvshapt(2,2) = 0.d0
      vshapt(3) = Xi(2)
      dvshapt(1,3) = 0.d0; dvshapt(2,3) = 1.d0
c
c  ...1D shape functions in the xi_3 direction
      vshap(1) = 1.d0 - Xi(3); dvshap(1) = -1.d0
      vshap(2) = Xi(3); dvshap(2) = 1.d0
c
c
      k=0
c
c  ...vertex shape functions
      do j=1,2
      do i=1,3
        k=k+1
        ShapH(k) = vshapt(i)*vshap(j)
        DshapH(1:2,k) = dvshapt(1:2,i)*vshap(j)
        DshapH(3,k) = vshapt(i)*dvshap(j)
      enddo
      enddo
c
c  ...mid-edge nodes shape functions
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c  ...loop through horizontal edges
      ie=0
      do j=1,2
      do i=1,3
        ie=ie+1
c
c  .....get the edge vertices specifying the local edge orientation
        iv1=TRIAN_EDGE_TO_VERT(1,i) ; iv2=TRIAN_EDGE_TO_VERT(2,i)
c
c  .....project Xi(1:2) onto the edge
        call proj_t2e(iv1,iv2,vshapt,dvshapt, te,dtedxi(1:2))
c
c  .....evaluate edge kernel shape functions
        call shapHke(te,Norder(ie),Nedge_orient(ie),
     .               nrdofHe,shapHk,dshapHk)
c
c  .....evaluate edge shape functions
        do l=1,nrdofHe
          k=k+1
          ShapH(k) = vshapt(iv1)*vshapt(iv2)*vshap(j)*shapHk(l)
          DshapH(1:2,k)
     .    = dvshapt(1:2,iv1)*vshapt(iv2)*vshap(j)*shapHk(l)
     .    + vshapt(iv1)*dvshapt(1:2,iv2)*vshap(j)*shapHk(l)
     .    + vshapt(iv1)*vshapt(iv2)*vshap(j)*dshapHk(l)*dtedxi(1:2)
          DshapH(3,k) = vshapt(iv1)*vshapt(iv2)*dvshap(j)*shapHk(l)
        enddo
      enddo
      enddo
c
c  ...loop through vertical edges
      do i=1,3
        ie=ie+1
c
c  .....evaluate edge bubble shape functions
        call shapHbe(Xi(3),Norder(ie),Nedge_orient(ie),
     .               nrdofHe,shapHb,dshapHb)
c
c  .....evaluate edge shape functions
        do l=1,nrdofHe
          k=k+1
          ShapH(k) = vshapt(i)*shapHb(l)
          DshapH(1:2,k) = dvshapt(1:2,i)*shapHb(l)
          DshapH(3,k)   = vshapt(i)*dshapHb(l)
        enddo
      enddo
c
c
c  ...mid-face nodes shape functions
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  ...loop through horizontal faces
      if=0
      do i=1,2
        if=if+1
c
c  .....compute the mid-face node kernel shape functions
        call shapHbt(Xi(1:2),Norder(9+if),Nface_orient(if),
     .               nrdofHt,shapHt,dshapHt)
c
c  .....evaluate the mid-face node shape functions
        do j=1,nrdofHt
          k=k+1
          ShapH(k) = shapHt(j)*vshap(if)
          DshapH(1:2,k) = dshapHt(1:2,j)*vshap(if)
          DshapH(3,k)   = shapHt(j)*dvshap(if)
        enddo
      enddo
c
c  ...loop through vertical faces
      do ii=1,3
        if=if+1
c
c  .....get the edge vertices specifying the local edge orientation
        iv1=TRIAN_EDGE_TO_VERT(1,ii) ; iv2=TRIAN_EDGE_TO_VERT(2,ii)
c
c  .....project Xi(1:2) onto the edge
        call proj_t2e(iv1,iv2,vshapt,dvshapt, te,dtedxi(1:2))
c
        call decode(Norder(9+if), nordh,nordv)
        select case(NFAXES(3,Nface_orient(if)))
c
c  .....the face axes have NOT been reversed
        case(0)
          call shapHke(te,nordh,NFAXES(1,Nface_orient(if)),
     .                 nrdofHk,shapHk,dshapHk)
          call shapHbe(Xi(3),nordv,NFAXES(2,Nface_orient(if)),
     .                 nrdofHb,shapHb,dshapHb)
          do j=1,nrdofHb
            do i=1,nrdofHk
              k=k+1
              ShapH(k) = vshapt(iv1)*vshapt(iv2)*shapHk(i)*shapHb(j)
              DshapH(1:2,k)
     .        = dvshapt(1:2,iv1)*vshapt(iv2)*shapHk(i)*shapHb(j)
     .        + vshapt(iv1)*dvshapt(1:2,iv2)*shapHk(i)*shapHb(j)
     .        + vshapt(iv1)*vshapt(iv2)*dshapHk(i)*dtedxi(1:2)*shapHb(j)
              DshapH(3,k) = vshapt(iv1)*vshapt(iv2)*shapHk(i)*dshapHb(j)
            enddo
          enddo
c
c  .....the face axes HAVE been reversed
        case(1)
          call shapHke(te,nordv,NFAXES(1,Nface_orient(if)),
     .                 nrdofHk,shapHk,dshapHk)
          call shapHbe(Xi(3),nordh,NFAXES(2,Nface_orient(if)),
     .                 nrdofHb,shapHb,dshapHb)
          do i=1,nrdofHk
            do j=1,nrdofHb
              k=k+1
              ShapH(k) = vshapt(iv1)*vshapt(iv2)*shapHk(i)*shapHb(j)
              DshapH(1:2,k)
     .        = dvshapt(1:2,iv1)*vshapt(iv2)*shapHk(i)*shapHb(j)
     .        + vshapt(iv1)*dvshapt(1:2,iv2)*shapHk(i)*shapHb(j)
     .        + vshapt(iv1)*vshapt(iv2)*dshapHk(i)*dtedxi(1:2)*shapHb(j)
              DshapH(3,k) = vshapt(iv1)*vshapt(iv2)*shapHk(i)*dshapHb(j)
            enddo
          enddo
        end select
      enddo
c
c
c  ...middle node shape functions
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      call decode(Norder(15), nordh,nordv)
      call shapHbt(Xi(1:2),nordh,0, nrdofHt,shapHt,dshapHt)
      call shapHbe(Xi(3),nordv,0, nrdofHb,shapHb,dshapHb)
      do j=1,nrdofHb
        do i=1,nrdofHt
          k=k+1
          ShapH(k) = shapHt(i)*shapHb(j)
          DshapH(1:2,k) = dshapHt(1:2,i)*shapHb(j)
          DshapH(3,k) = shapHt(i)*dshapHb(j)
        enddo
      enddo
c
c  ...save the number of dof
      Nrdof = k
c
c
      if (iprint.eq.1) then
        write(*,7001) Xi(1:3),Norder(1:15)
 7001   format('shapeHp: Xi = ',3e12.5,',  Norder = ',
     .          9i2,2x,2i2,3i3,2x,i3)
        write(*,7002)
 7002   format('         k,ShapH(k),DshapH(1:3,k) = ')
        do k=1,Nrdof
          write(*,7003) k, ShapH(k), DshapH(1:3,k)
 7003     format(i2,2x,e12.5,3x,3e12.5)
        enddo
        call pause
      endif
c
c
      end
c
c----------------------------------------------------------------------
c
c   routine name       - shapeHb
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 08
c
c   purpose            - routine evaluates shape functions for master
c                        brick
c
c   arguments :
c     in:
c         Xi           - master element coordinates
c         Norder       - order of approximation for the nodes
c         Nedge_orient - edge orientations
c         Nface_orient - face orientations
c     out:
c         Nrdof        - number of dof
c         Vshap        - values of shape functions
c         Dvshap       - values of derivatives of the shape functions
c                        wrt to master element coordinates
c
c----------------------------------------------------------------------
c
      subroutine shapeHb(Xi,Norder, Nedge_orient,Nface_orient,
     .                   Nrdof,ShapH,DshapH)
c
      use element_data , only : NFAXES,IJKV,IXIEDGE,IBLENDE,NBLENDE,
     .                                      IXIFACE,IBLENDF,NBLENDF
      use parameters   , only : MAXP,MAXbrickH
#include "syscom.blk"
c
      dimension Nedge_orient(12),Nface_orient(6),Norder(19)
      dimension Xi(3),ShapH(MAXbrickH),DshapH(3,MAXbrickH)
c
c
c  ...1D linear shape functions, first index indicates the function
c     number, second index the direction
      dimension vshap(2,3),dvshap(2,3)
c
c  ...1D bubble shape functions
      dimension shapH1(MAXP-1),dshapH1(MAXP-1),
     .          shapH2(MAXP-1),dshapH2(MAXP-1),
     .          shapH3(MAXP-1),dshapH3(MAXP-1)
c----------------------------------------------------------------------
c
c  ...calculate 1D linear shape functions
      do ixi=1,3
        vshap(1,ixi) = 1.d0 - Xi(ixi); dvshap(1,ixi) = -1.d0
        vshap(2,ixi) =        Xi(ixi); dvshap(2,ixi) =  1.d0
      enddo
c
      k=0
c
c----------------------------------------------------------------------
c  VERTEX SHAPE FUNCTIONS                                             |
c----------------------------------------------------------------------
      do iv=1,8
        k=k+1; i=IJKV(1,iv); j=IJKV(2,iv); kk=IJKV(3,iv)
        ShapH(k) = vshap(i,1)*vshap(j,2)*vshap(kk,3)
        DshapH(1,k) = dvshap(i,1)*vshap(j,2)*vshap(kk,3)
        DshapH(2,k) = vshap(i,1)*dvshap(j,2)*vshap(kk,3)
        DshapH(3,k) = vshap(i,1)*vshap(j,2)*dvshap(kk,3)
      enddo
c
c----------------------------------------------------------------------
c  EDGE SHAPE FUNCTIONS                                               |
c----------------------------------------------------------------------
      do ie=1,12
        ixi=IXIEDGE(ie)
        ibl1=IBLENDE(1,ie) ; ibl2=IBLENDE(2,ie)
        nv1 =NBLENDE(1,ie) ; nv2 =NBLENDE(2,ie)
        call shapHbe(Xi(ixi),Norder(ie),Nedge_orient(ie),
     .               nrdofHe,shapH1,dshapH1)
        do j=1,nrdofHe
          k=k+1
          ShapH(k) = shapH1(j)*vshap(nv1,ibl1)*vshap(nv2,ibl2)
          DshapH(ixi,k) = dshapH1(j)*vshap(nv1,ibl1)*vshap(nv2,ibl2)
          DshapH(ibl1,k) = shapH1(j)*dvshap(nv1,ibl1)*vshap(nv2,ibl2)
          DshapH(ibl2,k) = shapH1(j)*vshap(nv1,ibl1)*dvshap(nv2,ibl2)
        enddo
      enddo
c
c----------------------------------------------------------------------
c  FACE SHAPE FUNCTIONS                                               |
c----------------------------------------------------------------------
      do if=1,6
        ixi1 = IXIFACE(1,if); ixi2 = IXIFACE(2,if)
        ibl = IBLENDF(if); nv = NBLENDF(if)
        call decode(Norder(12+if), nordh,nordv)
        select case(NFAXES(3,Nface_orient(if)))
c
c  .....the face axes have NOT been reversed
        case(0)
          call shapHbe(Xi(ixi1),nordh,NFAXES(1,Nface_orient(if)),
     .                 nrdofH1,shapH1,dshapH1)
          call shapHbe(Xi(ixi2),nordv,NFAXES(2,Nface_orient(if)),
     .                 nrdofH2,shapH2,dshapH2)
          do j=1,nrdofH2
            do i=1,nrdofH1
              k=k+1
              ShapH(k) = shapH1(i)*shapH2(j)*vshap(nv,ibl)
              DshapH(ixi1,k) = dshapH1(i)*shapH2(j)*vshap(nv,ibl)
              DshapH(ixi2,k) = shapH1(i)*dshapH2(j)*vshap(nv,ibl)
              DshapH(ibl,k)  = shapH1(i)*shapH2(j)*dvshap(nv,ibl)
            enddo
          enddo
c
c  .....the face axes HAVE been reversed
        case(1)
          call shapHbe(Xi(ixi1),nordv,NFAXES(1,Nface_orient(if)),
     .                 nrdofH1,shapH1,dshapH1)
          call shapHbe(Xi(ixi2),nordh,NFAXES(2,Nface_orient(if)),
     .                 nrdofH2,shapH2,dshapH2)
c
c  .......the order of loops gets reversed...
          do i=1,nrdofH1
            do j=1,nrdofH2
              k=k+1
              ShapH(k) = shapH1(i)*shapH2(j)*vshap(nv,ibl)
              DshapH(ixi1,k) = dshapH1(i)*shapH2(j)*vshap(nv,ibl)
              DshapH(ixi2,k) = shapH1(i)*dshapH2(j)*vshap(nv,ibl)
              DshapH(ibl,k)  = shapH1(i)*shapH2(j)*dvshap(nv,ibl)
            enddo
          enddo
        end select
      enddo
c
c  ...calculate middle node shape functions
      call decode(Norder(19), nordh, nord3)
      call decode(nordh, nord1,nord2)
      call shapHbe(Xi(1),nord1,0, nrdofH1,shapH1,dshapH1)
      call shapHbe(Xi(2),nord2,0, nrdofH2,shapH2,dshapH2)
      call shapHbe(Xi(3),nord3,0, nrdofH3,shapH3,dshapH3)
      do kk=1,nrdofH3
        do j=1,nrdofH2
          do i=1,nrdofH1
            k=k+1
            ShapH(k) = shapH1(i)*shapH2(j)*shapH3(kk)
            DshapH(1,k) = dshapH1(i)*shapH2(j)*shapH3(kk)
            DshapH(2,k) = shapH1(i)*dshapH2(j)*shapH3(kk)
            DshapH(3,k) = shapH1(i)*shapH2(j)*dshapH3(kk)
          enddo
        enddo
      enddo
c
c  ...total number of dofs
      Nrdof=k
c
c
      end
c
c----------------------------------------------------------------------
c
c   routine name       - shapeHn
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 08
c
c   purpose            - evaluate H1 shape functions for a tetrahedral
c                        element
c
c   arguments :
c     in:
c         Nedge_orient - edge orientations
c         Nface_orient - face orientations
c         Xi           - master element coordinates
c         Norder       - order of approximation for the nodes
c     out:
c         Nrdof        - number of dof
c         ShapH        - values of shape functions
c         DshapH       - values of derivatives of the shape functions
c                        wrt to master element coordinates
c
c----------------------------------------------------------------------
c
      subroutine shapeHn(Xi,Norder, Nedge_orient,Nface_orient,
     .                   Nrdof,ShapH,DshapH)
c
      use parameters
      use element_data
#include "syscom.blk"
      common /cshapeHn/ iprint
c
      real*8 , dimension(3), intent(in) :: Xi
      dimension Nedge_orient(6), Nface_orient(4),Norder(11)
ccc      dimension Xi(3),ShapH(MAXbrickH),DshapH(3,MAXbrickH)
      dimension ShapH(MAXbrickH),DshapH(3,MAXbrickH)
c
c  ...derivatives of edge coordinate
      dimension dtedxi(3)
c
c  ...face coordinates
      dimension tf(2),dtfdxi(2,3)
c
c  ...edge kernel shape functions
      dimension shapHek(MAXP-1),dshapHek(MAXP-1)
c
c  ...mid-face node kernel shape functions
      dimension shapHtk(MAXmdltH),dshapHtk(2,MAXmdltH)
c
c  ...middle node kernel shape functions
      dimension shapHnk(MAXmdlnH),dshapHnk(3,MAXmdlnH)
c
c  ...blending function
      dimension dblend(1:3)
c
      iprint=0
c
c
      k=0
c
c  ...vertex linear shape functions (affine coordinates)
      k=k+1
      ShapH(k) = 1.d0 - (Xi(1)+Xi(2)+Xi(3))
      DshapH(1:3,k) = -1.d0
c
      do i=1,3
        k=k+1
        ShapH(k) = Xi(i)
        DshapH(1:3,k) = 0.d0; DshapH(i,k) = 1.d0
      enddo
c
c
c  ...mid-edge nodes shape functions
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c  ...loop through edges
      do ie=1,6
c
c  .....get the edge vertices specifying the local edge orientation
        iv1=TETRA_EDGE_TO_VERT(1,ie) ; iv2=TETRA_EDGE_TO_VERT(2,ie)
c
c  .....project Xi onto the edge
        call proj_n2e(iv1,iv2,ShapH,DshapH, te,dtedxi)
        if (iprint.eq.1) then
          write(*,7004) ie,te,dtedxi(1:3)
 7004     format('shapeHn: ie,te,dtedxi = ',i2,2x,f8.3,2x,3f8.3)
          call pause
        endif
c
c  .....evaluate edge kernel shape functions
        call shapHke(te,Norder(ie),Nedge_orient(ie),
     .               nrdofHe,shapHek,dshapHek)
c
c  .....blending function
        blend = ShapH(iv1)*ShapH(iv2)
        dblend(1:3) = DshapH(1:3,iv1)*ShapH(iv2)
     .              + ShapH(iv1)*DshapH(1:3,iv2)
c
c  .....evaluate edge shape functions
        do j=1,nrdofHe
          k=k+1
          ShapH(k) = blend*shapHek(j)
          DshapH(1:3,k) = dblend(1:3)*shapHek(j)
     .                  + blend*dshapHek(j)*dtedxi(1:3)
        enddo
      enddo
c
c  ...mid-face nodes shape functions
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     loop through faces
      do if=1,4
c
c  .....get the vertices for the face defining its local orientation
        iv1=TETRA_FACE_TO_VERT(1,if)
        iv2=TETRA_FACE_TO_VERT(2,if)
        iv3=TETRA_FACE_TO_VERT(3,if)
c
c  .....project Xi onto the face
        call proj_n2f(iv1,iv2,iv3,ShapH,DshapH, tf,dtfdxi)
c
c  .....compute the mid-face node kernel shape functions
        call shapHkt(tf,Norder(6+if),Nface_orient(if),
     .               nrdofHt,shapHtk,dshapHtk)
c
c  .....blending function
        blend = ShapH(iv1)*ShapH(iv2)*ShapH(iv3)
        dblend(1:3) = DshapH(1:3,iv1)*ShapH(iv2)*ShapH(iv3)
     .              + ShapH(iv1)*DshapH(1:3,iv2)*ShapH(iv3)
     .              + ShapH(iv1)*ShapH(iv2)*DshapH(1:3,iv3)
c
c  .....evaluate the mid-face node shape functions
        do j=1,nrdofHt
          k=k+1
          ShapH(k) = blend*shapHtk(j)
          do ivar=1,3
            DshapH(ivar,k) = dblend(ivar)*shapHtk(j)
     .    + blend*(dshapHtk(1,j)*dtfdxi(1,ivar)
     .            +dshapHtk(2,j)*dtfdxi(2,ivar))
          enddo
        enddo
      enddo
c
c  ...middle node shape functions
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  ...compute the middle node kernel shape functions
      call shapHkn(Xi,Norder(11), nrdofHn,shapHnk,dshapHnk)
      do j=1,nrdofHn
        k=k+1
        ShapH(k) = ShapH(1)*ShapH(2)*ShapH(3)*ShapH(4)*shapHnk(j)
        DshapH(1:3,k)
     .  = DshapH(1:3,1)*ShapH(2)*ShapH(3)*ShapH(4)*shapHnk(j)
     .  + ShapH(1)*DshapH(1:3,2)*ShapH(3)*ShapH(4)*shapHnk(j)
     .  + ShapH(1)*ShapH(2)*DshapH(1:3,3)*ShapH(4)*shapHnk(j)
     .  + ShapH(1)*ShapH(2)*ShapH(3)*DshapH(1:3,4)*shapHnk(j)
     .  + ShapH(1)*ShapH(2)*ShapH(3)*ShapH(4)*dshapHnk(1:3,j)
      enddo
c
c  ...save the number of dof
      Nrdof = k
c
c
      if (iprint.eq.1) then
        write(*,7001) Xi(1:3),Norder(1:11)
 7001   format('shapeHn: Xi = ',3e12.5,',  Norder = ',6i2,2x,4i2,2x,i2)
        write(*,7002)
 7002   format('         k,ShapH(k),DshapH(1:3,k) = ')
        do k=1,Nrdof
          write(*,7003) k, ShapH(k), DshapH(1:3,k)
 7003     format(i2,2x,e12.5,3x,3e12.5)
        enddo
        call pause
      endif
c
c
      end
c
c----------------------------------------------------------------------
c
c   routine name       - shapeHd
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 08
c
c   purpose            - a  3-D pyramid element routine
c                        evaluating shape functions and their derivati-
c                        ves
c
c   arguments :
c     in:
c         Xi           - master element coordinates
c         Norder       - order of approximation for the nodes
c         Nedge_orient - edge orientations
c         Nface_orient - face orientations
c     out:
c         Nrdof        - number of the element shape functions
c         ShapH        - values of shape functions
c         DshapH       - values of derivatives of the shape functions
c                        wrt to master element coordinates
c
c----------------------------------------------------------------------
c
      subroutine shapeHd(Xi,Norder, Nedge_orient,Nface_orient,
     .                   Nrdof,ShapH,DshapH)
c
      use element_data , only : NFAXES
      use control , only : GEOM_TOL
      use parameters
#include "syscom.blk"
      common /cshapeHd/ iprint
      common /cshapke/ iprint_shapke
c
      dimension Nedge_orient(8),Nface_orient(5),Norder(14)
      dimension Xi(3),ShapH(MAXpyramH),DshapH(3,MAXpyramH)
c
c  ...blending functions
      dimension dblend(3),blend_edge(4),dblend_edge(3,4)
c
c  ...projections
      dimension dtedxi(3),tf(2),dtfdxi(2,3)
c
c  ...1D bubble/kernel shape functions
      dimension shapH1(MAXP-1),dshapH1(MAXP-1),
     .          shapH2(MAXP-1),dshapH2(MAXP-1)
c
c  ...mid-face node kernel shape functions
      dimension shapHtk(MAXmdltH),dshapHtk(2,MAXmdltH)
c
      double precision z,z1
c----------------------------------------------------------------------
c
ccc      iprint=1
ccc      iprint_shapke = iprint
c
c  ...check that coordinates are inside reference pyramid
      if ((Xi(1)      .lt.-GEOM_TOL    ).or.
     .    (Xi(2)      .lt.-GEOM_TOL    ).or.
     .    (Xi(3)      .lt.-GEOM_TOL    ).or.
     .    (Xi(1)+Xi(3).gt.1.d0+GEOM_TOL).or.
     .    (Xi(2)+Xi(3).gt.1.d0+GEOM_TOL)     ) then
ccc        write(*,*)'shapeHd: point out of master pyramid!'
ccc        write(*,1000)Xi(1:3)
1000    format(' Xi = ',3(e12.5,2x))
ccc        call pause
      endif
c
      x = Xi(1); y = Xi(2); z = Xi(3)
c
c  ...pyramid shape functions cannot be evaluated at z=1
      if (abs(z-1.d0).lt.1.d-12) then
ccc        write(*,*) 'shapeHd: z = ',z
        z = 1.d0 - 1.d-12
      endif
c
c  ...auxiliary variables
      xz1 = 1.d0-x-z
      yz1 = 1.d0-y-z
      z1 = 1.d0-z
c
c  ...initialize dof counter
      k=0
c
c  ...first vertex shape function
      k=k+1
      ShapH(k) = xz1*yz1/z1
      DshapH(1,k) = -yz1/z1
      DshapH(2,k) = -xz1/z1
      DshapH(3,k) = (-yz1-xz1)/z1 + xz1*yz1/z1**2
c
c  ...second vertex shape function
      k=k+1
      ShapH(k) = x*yz1/z1
      DshapH(1,k) = yz1/z1
      DshapH(2,k) = -x/z1
      DshapH(3,k) = -x/z1 + x*yz1/z1**2
c
c  ...third vertex shape function
      k=k+1
      ShapH(k) = x*y/z1
      DshapH(1,k) = y/z1
      DshapH(2,k) = x/z1
      DshapH(3,k) = x*y/z1**2
c
c  ...fourth vertex shape function
      k=k+1
      ShapH(k) = xz1*y/z1
      DshapH(1,k) = -y/z1
      DshapH(2,k) = xz1/z1
      DshapH(3,k) = -y/z1 + xz1*y/z1**2
c
c  ...fifth vertex shape function
      k=k+1
      ShapH(k) = z
      DshapH(1,k) = 0.d0
      DshapH(2,k) = 0.d0
      DshapH(3,k) = 1.d0
c
c-----------------------------------------------------------------------
c
c  ...mid-edge node shape functions
      do ie=1,8
c
c  .....blending function
        select case(ie)
        case(1)
          blend = x*xz1*yz1/z1
          dblend(1) = (xz1-x)*yz1/z1
          dblend(2) = -x*xz1/z1
          dblend(3) = x*(-yz1-xz1)/z1 + x*xz1*yz1/z1**2
        case(2)
          blend = x*y*yz1/z1
          dblend(1) = y*yz1/z1
          dblend(2) = x*(yz1-y)/z1
          dblend(3) = -x*y/z1 + x*y*yz1/z1**2
        case(3)
          blend = x*xz1*y/z1
          dblend(1) = (xz1-x)*y/z1
          dblend(2) = x*xz1/z1
          dblend(3) = -x*y/z1 + x*xz1*y/z1**2
        case(4)
          blend = xz1*y*yz1/z1
          dblend(1) = -y*yz1/z1
          dblend(2) = xz1*(yz1-y)/z1
          dblend(3) = y*(-yz1-xz1)/z1 + xz1*y*yz1/z1**2
        case(5,6,7,8)
          iv = ie-4
          blend = ShapH(iv)*z
          dblend(1:2) = DshapH(1:2,iv)*z
          dblend(3) = DshapH(3,iv)*z + ShapH(iv)
        end select
        if (ie.le.4) then
          blend_edge(ie) = blend
          dblend_edge(1:3,ie) = dblend(1:3)
        endif
c
c  .....project x,y,z onto the edge
        call proj_d2e(Xi,ie,ShapH,DshapH, te,dtedxi)
c
c  .....evaluate the edge kernel functions taking into account
c       the orientation
        call shapHke(te,Norder(ie),Nedge_orient(ie),
     .               nrdofHe,shapH1,dshapH1)
c
c  .....evaluate the edge shape functions
        do j=1,nrdofHe
          k=k+1
          ShapH(k) = shapH1(j)*blend
          DshapH(1:3,k) = dshapH1(j)*dtedxi(1:3)*blend
     .                  + shapH1(j)*dblend(1:3)
        enddo
      enddo
c
c-----------------------------------------------------------------------
c
c  ...bottom mid-face node shape functions
      call decode(Norder(9), nordh,nordv)
      if ((nordh.ge.2).and.(nordv.ge.2)) then
c
c  .....blending function
        blend = ShapH(1)*ShapH(3)
        dblend(1:3) = DshapH(1:3,1)*ShapH(3) + ShapH(1)*DshapH(1:3,3)
c
c  .....project onto the bottom face
        xz = x/z1; yz = y/z1
c
        select case(NFAXES(3,Nface_orient(1)))
c
c  .....the face axes have NOT been reversed
        case(0)
          call shapHke(xz,nordh,NFAXES(1,Nface_orient(1)),
     .                 nrdofH1,shapH1,dshapH1)
          call shapHke(yz,nordv,NFAXES(2,Nface_orient(1)),
     .                 nrdofH2,shapH2,dshapH2)
          do j=1,nrdofH2
            do i=1,nrdofH1
c
c  ...........the (1-z) factor
              nexp = max(j-1,i-1)
              fact = z1**nexp
              if (nexp.eq.0) then
                dfact = 0.d0
              else
                dfact = -z1**(nexp-1)*nexp
              endif
              k=k+1
              ShapH(k) = shapH1(i)*shapH2(j)*fact
              DshapH(1,k) = dshapH1(i)/z1*shapH2(j)*fact
              DshapH(2,k) = shapH1(i)*dshapH2(j)/z1*fact
              DshapH(3,k) = (dshapH1(i)*x/z1**2*shapH2(j)
     .                     + shapH1(i)*dshapH2(j)*y/z1**2)*fact
     .                    + shapH1(i)*shapH2(j)*dfact
c
c  ...........blend
              DshapH(1:3,k) = DshapH(1:3,k)*blend + ShapH(k)*dblend(1:3)
              ShapH(k) = ShapH(k)*blend
            enddo
          enddo
c
c  .....the face axes HAVE been reversed
        case(1)
          call shapHke(xz,nordv,NFAXES(1,Nface_orient(1)),
     .                 nrdofH1,shapH1,dshapH1)
          call shapHke(yz,nordh,NFAXES(2,Nface_orient(1)),
     .                 nrdofH2,shapH2,dshapH2)
          do i=1,nrdofH1
            do j=1,nrdofH2
c
c  ...........the (1-z) factor
              nexp = max(j-1,i-1)
              fact = z1**nexp
              if (nexp.eq.0) then
                dfact = 0.d0
              else
                dfact = -z1**(nexp-1)*nexp
              endif
              k=k+1
              ShapH(k) = shapH1(i)*shapH2(j)*fact
              DshapH(1,k) = dshapH1(i)/z1*shapH2(j)*fact
              DshapH(2,k) = shapH1(i)*dshapH2(j)/z1*fact
              DshapH(3,k) = (dshapH1(i)*x/z1**2*shapH2(j)
     .                     + shapH1(i)*dshapH2(j)*y/z1**2)*fact
     .                    + shapH1(i)*shapH2(j)*dfact
c
c  ...........blend
              DshapH(1:3,k) = DshapH(1:3,k)*blend + ShapH(k)*dblend(1:3)
              ShapH(k) = ShapH(k)*blend
            enddo
          enddo
        end select
      endif
c
c-----------------------------------------------------------------------
c
c  ...lateral mid-face nodes shape functions
      do if=2,5
c
c  .....blending function
        blend = blend_edge(if-1)*z
        dblend(1:3) = dblend_edge(1:3,if-1)*z
        dblend(3) = dblend(3) + blend_edge(if-1)
c
c  .....project xyz onto the face
        call proj_d2f(Xi,if, tf,dtfdxi)
c
c  .....compute the mid-face node kernel shape functions, taking into
c       account the face orientation
        call shapHkt(tf,Norder(8+if),Nface_orient(if),
     .               nrdofHt,shapHtk,dshapHtk)
c
c  .....evaluate the mid-face node shape functions
        do j=1,nrdofHt
          k=k+1
          ShapH(k) = blend*shapHtk(j)
          do ivar=1,3
            DshapH(ivar,k) = dblend(ivar)*shapHtk(j)
     .    + blend*(dshapHtk(1,j)*dtfdxi(1,ivar)
     .            +dshapHtk(2,j)*dtfdxi(2,ivar))
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c
c  ...middle node shape functions
c
c  ...blending function
      blend = ShapH(1)*ShapH(3)*z
      dblend(1:3) = DshapH(1:3,1)*ShapH(3)*z + ShapH(1)*DshapH(1:3,3)*z
      dblend(3) = dblend(3) + ShapH(1)*ShapH(3)
c
c  ...project onto the bottom face
      xz = x/z1; yz = y/z1
      nord = Norder(14)
c
c  ...evaluate the kernel shape functions
      call shapHke(xz,nord,0, nrdofH1,shapH1,dshapH1)
      call shapHke(yz,nord,0, nrdofH2,shapH2,dshapH2)
c
c  ...evaluate the shape functions
      do kk=1,nord-1
c
c  .....evaluate the (1-z)**nexp factor
        nexp = kk-1
        fact = z1**nexp
        if (nexp.eq.0) then
          dfact = 0.d0
        else
          dfact = -z1**(nexp-1)*nexp
        endif
c
        do j=1,kk
        do i=1,kk
          k=k+1
          ShapH(k) = shapH1(i)*shapH2(j)*fact
          DshapH(1,k) = dshapH1(i)/z1*shapH2(j)*fact
          DshapH(2,k) = shapH1(i)*dshapH2(j)/z1*fact
          DshapH(3,k) = (dshapH1(i)*x/z1**2*shapH2(j)
     .                 + shapH1(i)*dshapH2(j)*y/z1**2)*fact
     .                + shapH1(i)*shapH2(j)*dfact
c
c  .......blend
          DshapH(1:3,k) = DshapH(1:3,k)*blend+shapH(k)*dblend(1:3)
          ShapH(k) = ShapH(k)*blend
        enddo
        enddo
      enddo
c
c
      Nrdof=k
c
      if (iprint.ge.1) then
        write(*,7001) Xi(1:3),Norder(1:14)
 7001   format('shapeHd: Xi = ',3e12.5,',  Norder = ',
     .          8i2,2x,i3,4i2,2x,i3)
        write(*,7002)
 7002   format('         k,ShapH(k),DshapH(1:3,k) = ')
        do k=1,Nrdof
          write(*,7003) k, ShapH(k), DshapH(1:3,k)
 7003     format(i3,2x,e12.5,3x,3e12.5)
        enddo
        call pause
      endif
c
c
      end

c----------------------------------------------------------------------
c
c   routine name       - shapeHt_0
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - evaluate triangle H^1 shape functions
c
c   arguments :
c     in:
c        Xi            - master element coordinates
c        Norder        - order of approximation for the nodes
c        Nedge_orient  - orientation for the edges
c     out:
c        Nrdof         - number of dof
c        ShapH         - values of shape functions
c        DshapH        - values of derivatives of the shape functions
c                        wrt to master element coordinates
c
c----------------------------------------------------------------------
c
      subroutine shapeHt_0(Xi,Norder,Nedge_orient, Nrdof,ShapH,DshapH)
c
      use parameters
      use element_data
      implicit none
c
      integer :: Norder(4),Nedge_orient(3),Nrdof
      real(8) :: Xi(2),ShapH(MAXtriaH),DshapH(2,MAXtriaH)
c
c  ...derivatives of edge coordinate
      real(8) :: dtdxi(2)
c
c  ...edge kernel shape functions
      real(8) :: shapHek(MAXP-1),dshapHek(MAXP-1)
c
      integer :: ie,iv,iv1,iv2,j,k,nrdofHb,nrdofHe
      real(8) :: t
c
c----------------------------------------------------------------------
c
c  ...calculate the affine coordinates (linear shape functions)
      ShapH(1) = 1.d0 -  (Xi(1)+Xi(2))
      DshapH(1,1) = -1.d0
      DshapH(2,1) = -1.d0
c
      ShapH(2) = Xi(1)
      DshapH(1,2) = 1.d0
      DshapH(2,2) = 0.d0
c
      ShapH(3) = Xi(2)
      DshapH(1,3) = 0.d0
      DshapH(2,3) = 1.d0
c
c  ...initiate shape function counter
      k=3
c
c  ...mid-edge nodes shape functions
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c  ...loop through edges
      do ie=1,3
c
c  .....get the edge vertices specifying the local edge orientation
        iv1=TRIAN_EDGE_TO_VERT(1,ie) ; iv2=TRIAN_EDGE_TO_VERT(2,ie)
c
c  .....project Xi onto the edge
        call proj_t2e(iv1,iv2,ShapH,DshapH, t,dtdxi)
c
c  .....evaluate edge kernel shape functions
        call shapHke(t,Norder(ie),Nedge_orient(ie),
     .               nrdofHe,shapHek,dshapHek)
c
c  .....evaluate edge shape functions
        do j=1,nrdofHe
          k=k+1
          ShapH(k) = ShapH(iv1)*ShapH(iv2)*shapHek(j)
          DshapH(1:2,k) = DshapH(1:2,iv1)*ShapH(iv2)*shapHek(j)
     .                  + ShapH(iv1)*DshapH(1:2,iv2)*shapHek(j)
     .                  + ShapH(iv1)*ShapH(iv2)*dshapHek(j)*dtdxi(1:2)
        enddo
      enddo
c
c  ...middle node shape functions
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  ...compute the middle node bubble shape functions
      call shapHbt(Xi,Norder(4),0, nrdofHb,ShapH(k+1),DshapH(1:2,k+1))
c
c  ...save the number of dof
      Nrdof = k+nrdofHb
c
c
      end subroutine shapeHt_0
c
c----------------------------------------------------------------------
c
c   routine name       - shapeHq_0
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - evaluate quad H^1 shape functions
c
c   arguments :
c     in:
c        Xi            - master element coordinates
c        Norder        - order of approximation for the nodes
c        Nedge_orient  - orientation for the edges
c     out:
c        Nrdof         - number of dof
c        ShapH         - values of shape functions
c        DshapH        - values of derivatives of the shape functions
c                        wrt to master element coordinates
c
c   required  routines -
c
c----------------------------------------------------------------------
c
      subroutine shapeHq_0(Xi,Norder,Nedge_orient, Nrdof,ShapH,DshapH)
c
      use parameters
      use element_data
      implicit none
c
      integer :: Norder(5),Nedge_orient(4),Nrdof
      real(8) :: Xi(2),ShapH(MAXquadH),DshapH(2,MAXquadH)
c
c  ...1D linear shape functions, first index indicates the function
c     number, second index the direction
      real(8) :: vshap(2,2),dvshap(2,2)
c
c  ...defines vertex shape functions; for vertex iv,
c     ijv(ixi,iv) = number of 1D linear shape function in Xi(ixi)
      integer :: ijv(2,4)
c
c  ...defines edge shape functions; for edge ie
c     ixiedg(ie) = edge coordinate
c     iblend(ie) = blending direction
c     nblend(ie) = number of blending linear shape function
      integer :: ixiedg(4),iblend(4),nblend(4)
c
c  ...1D bubble functions
      real(8) :: shapH1(MAXP-1),dshapH1(MAXP-1),
     .           shapH2(MAXP-1),dshapH2(MAXP-1)
c
      integer :: i,ibl,ie,ixi,iv,j,k,nv
      integer :: nordh,nordv,nrdofH1,nrdofH2,nrdofHe
c
      data ijv /1,1, 2,1, 2,2, 1,2/
      data ixiedg /1,2,1,2/
      data iblend /2, 1, 2, 1/
      data nblend /1, 2, 2, 1/
c
c----------------------------------------------------------------------
c
c  ...calculate 1D linear shape functions
      do ixi=1,2
        vshap(1,ixi) = 1.d0 - Xi(ixi); dvshap(1,ixi) = -1.d0
        vshap(2,ixi) =        Xi(ixi); dvshap(2,ixi) =  1.d0
      enddo
c
      k=0
c
c  ...calculate vertex shape functions
      do iv=1,4
        k=k+1; i=ijv(1,iv); j=ijv(2,iv)
        ShapH(k) = vshap(i,1)*vshap(j,2)
        DshapH(1,k) = dvshap(i,1)*vshap(j,2)
        DshapH(2,k) = vshap(i,1)*dvshap(j,2)
      enddo
c
c  ...calculate edge shape functions
      do ie=1,4
        ixi = ixiedg(ie); ibl = iblend(ie); nv = nblend(ie)
        call shapHbe(Xi(ixi),Norder(ie),Nedge_orient(ie),
     .               nrdofHe,shapH1,dshapH1)
        do j=1,nrdofHe
          k=k+1
          ShapH(k) = shapH1(j)*vshap(nv,ibl)
          DshapH(ixi,k) = dshapH1(j)*vshap(nv,ibl)
          DshapH(ibl,k) = shapH1(j)*dvshap(nv,ibl)
        enddo
      enddo
c
c  ...middle node shape functions
      call decode(Norder(5), nordh,nordv)
      call shapHbe(Xi(1),nordh,0, nrdofH1,shapH1,dshapH1)
      call shapHbe(Xi(2),nordv,0, nrdofH2,shapH2,dshapH2)
      do j=1,nrdofH2
        do i=1,nrdofH1
          k=k+1
          ShapH(k) = shapH1(i)*shapH2(j)
          DshapH(1,k) = dshapH1(i)*shapH2(j)
          DshapH(2,k) = shapH1(i)*dshapH2(j)
        enddo
      enddo
c
c  ...save the number of dof
      Nrdof = k
c
c
      end subroutine shapeHq_0

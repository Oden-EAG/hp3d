!----------------------------------------------------------------------
!
!   routine name       - setcnstr_trian_iso_hcurl
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine evaluates constraint coefficients for
!                        the triangular master element
!
!   arguments          - none
!
!----------------------------------------------------------------------
!
!    'big' parent nodes
!
!      *
!      * *
!      *   *
!      *     *
!      *       *
!      *         *
!    4 *           *  3
!      *             *
!      *       1       *
!      *                 *
!      *                   *
!      *                     *
!      *************************
!                  2
!
!      'small' constrained nodes
!      *
!      * *
!      *   *
!      *     *
!      *   3   *
!      *         *
!      ******7******
!      * *         * *
!      *   *   4   *   *
!      *     5     6     *
!      *   1   *   *   2   *
!      *         * *         *
!      *************************
!
!--------------------------------------------------------------------
!
   subroutine setcnstr_trian_iso_hcurl
!
      use parameters
      use constraints
      use element_data
!
      implicit none
!
!  ...use collocation points for lagrange shape functions h1
      integer, parameter :: npts  = 4*MAXtriaE
      integer, parameter :: lwork = MAXtriaE**2
!
      real(8) :: shapsma(2,MAXtriaE),curl(MAXtriaE),              &
                 shapbig(2,MAXtriaE),                             &
                 evalsma(npts,MAXtriaE),evalbig(npts,MAXtriaE),   &
                 work(lwork),r(MAXtriaE,MAXtriaE),xibig(2),xisma(2)
!
!  ...order and edge orientations
      integer :: norder(4),norient(3),nsize(2)
!
      real(8) :: dbds(2,2)
!
      integer :: i,ib,ic,ibcase,info,is,iscase
      integer :: ivoid,jvoid,kvoid
      integer :: j,k,nord,nc,nrb,nrdof,nrs,nspan
      integer :: ndof_i,ndof_j,ndof_mdlt,ndof_medg
!
#if HP3D_DEBUG
      real(8) :: val(2)
      integer :: ians,nel
      integer :: iprint
      iprint = 0
#endif
!
!  ...set up order and orientation
      nord = MAXP
      norder  = nord
      norient = 0
      nsize   = (/MAXP,MAXtriaE/)
!
      ivoid = 0
      call ndof_nod(MEDG,nord, ivoid, ndof_medg, jvoid, kvoid)
      call ndof_nod(MDLT,nord, ivoid, ndof_mdlt, jvoid, kvoid)
      nrdof = 3*ndof_medg + ndof_mdlt
!
!*********************************************************************
!
!  ...small triangle 1 -- small nodes 1 and 5:
!
!*********************************************************************
!
      ic =0; nspan = (2*nord);
      dbds    = 0.d0;   curl    = 0.d0;
      shapsma = 0.d0;   shapbig = 0.d0;
      evalsma = 0.d0;   evalbig = 0.d0;
!  ...derivative dx/dxi
      dbds(1,1) = 0.5d0;
      dbds(2,2) = 0.5d0
!
!  ...loop through collocation points:
      do i=1,nspan
!
        xibig(1) = (i-1)*0.5d0/(nspan-1)
        xisma(1) = (i-1)*1.0d0/(nspan-1)
!
        do j=1,nspan-(i-1)
!
          xibig(2) = 0.5d0/(nspan-1)*(j-1)
          xisma(2) = 1.0d0/(nspan-1)*(j-1)
!
!  .......find small shape functions at this point:
!          call shapeEt(                  &
!            xisma,norder,norient,ivoid,  &
!            shapsma(1:2,1:MAXtriaE),     &
!            curl(1:MAXtriaE))
          call shape2DETri(xisma,norder,norient,nsize,   &
                           ivoid,shapsma(1:2,:),curl)
!
!  .......find big shape functions at this point:
!          call shapeEt(                  &
!            xibig,norder,norient,ivoid,  &
!            shapbig(1:2,1:MAXtriaE),     &
!            curl(1:MAXtriaE))
          call shape2DETri(xibig,norder,norient,nsize,   &
                           ivoid,shapbig(1:2,:),curl)
!  .......for each component variable,
          if ((ic+2).gt.npts) then
            write(*,*) '1 setcnstr_trian_iso_hcurl: INCREASE WORKSPACE', &
              npts, (ic+2)
          endif
!
          do k=1,2
            evalbig(ic+k, 1:MAXtriaE) = shapbig(k,1:MAXtriaE)*dbds(k,1)  &
                                      + shapbig(k,1:MAXtriaE)*dbds(k,2)
            evalsma(ic+k, 1:MAXtriaE) = shapsma(k,1:MAXtriaE)
          enddo
          ic = ic + 2
        enddo
      enddo
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) 'ic, nrdof = ', ic, nrdof
        call pause
      endif
!
      if (iprint.eq.1) then
        write(*,*) 'setcnstr_trian_iso_hcurl: COLUMNWISE ....'
        do i=1,nc
          write(*,7003) evalsma(1:nc,i)
 7003     format(15(e12.5,2x))
        enddo
        call pause
      endif
#endif
!
      nc = ic;
      call dgels('N', nc, nrdof, nrdof,   &
                 evalsma, npts, evalbig, npts, work, lwork, info)
      if (info.ne.0) then
        write(*,*)'setcnstr_trian_iso_hcurl: HCURL DGELS INFO =',info
        call logic_error(FAILURE,__FILE__,__LINE__)
      endif
!!!      r = 2.d0*evalbig(1:nrdof,1:nrdof)
      r = evalbig(1:nrdof,1:nrdof)
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7002) r(1:nrdof,1:nrdof)
 7002   format('r=', 10e12.5)
      endif
#endif
!
!  ...degrees of freedom of small node 5:
!****************************************
      nrs = ndof_medg
      do i=1,ndof_medg
!
!  .....big node 2,3,4:
        do ib=2,4
          nrb = (ib-2)*ndof_medg
          do j=1,ndof_medg
            RRTE(ib,j,5,i) = r(nrs+i,nrb+j)
          enddo
        enddo
!
!  .....big node 1:
        nrb = 3*ndof_medg
        do j=1,ndof_mdlt
          RRTE(1,j,5,i) = r(nrs+i,nrb+j)
        enddo
      enddo
!
!  ...degrees of freedom of small node 1:
!*****************************************
      nrs = 3*ndof_medg
      do i=1,ndof_mdlt
!
!  .....big node 2,3,4:
        do ib=2,4
          nrb = (ib-2)*ndof_medg
          do j=1,ndof_medg
            RRTE(ib,j,1,i) = r(nrs+i,nrb+j)
          enddo
        enddo
!
!  .....big node 1:
        nrb = 3*ndof_medg
        do j=1,ndof_mdlt
          RRTE(1,j,1,i) = r(nrs+i,nrb+j)
        enddo
      enddo
!*********************************************************************
!
!  ...small triangle 2 -- small nodes 2 and 6:
!
!*********************************************************************
!
      ic =0; nspan = (2*nord);
      dbds    = 0.d0;   curl    = 0.d0;
      shapsma = 0.d0;   shapbig = 0.d0;
      evalsma = 0.d0;   evalbig = 0.d0;
!  ...derivative dx/dxi
      dbds(1,1) = 0.5d0;
      dbds(2,2) = 0.5d0
!
!  ...loop through collocation points:
      do i=1,nspan
!
        xibig(1) = (i-1)*0.5d0/(nspan-1) + 0.5d0
        xisma(1) = (i-1)*1.0d0/(nspan-1)
!
        do j=1,nspan-(i-1)
!
          xibig(2) = 0.5d0/(nspan-1)*(j-1)
          xisma(2) = 1.0d0/(nspan-1)*(j-1)
!
!  .......find small shape functions at this point:
!          call shapeEt(                  &
!            xisma,norder,norient,ivoid,  &
!            shapsma(1:2,1:MAXtriaE),     &
!            curl(1:MAXtriaE))
          call shape2DETri(xisma,norder,norient,nsize,   &
                           ivoid,shapsma(1:2,:),curl)
!
!  .......find big shape functions at this point:
!          call shapeEt(                  &
!            xibig,norder,norient,ivoid,  &
!            shapbig(1:2,1:MAXtriaE),     &
!            curl(1:MAXtriaE))
          call shape2DETri(xibig,norder,norient,nsize,   &
                           ivoid,shapbig(1:2,:),curl)
!
!  .......for each component variable,
          if ((ic+2).gt.npts) then
            write(*,*) '2 setcnstr_trian_iso_hcurl: INCREASE WORKSPACE',npts
          endif
!
          do k=1,2
            evalbig(ic+k, 1:MAXtriaE) = shapbig(k,1:MAXtriaE)*dbds(k,1)  &
                                      + shapbig(k,1:MAXtriaE)*dbds(k,2)
            evalsma(ic+k, 1:MAXtriaE) = shapsma(k,1:MAXtriaE)
          enddo
          ic = ic + 2
        enddo
      enddo
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) 'ic, nrdof = ', ic, nrdof
        call pause
      endif
!
      if (iprint.eq.1) then
        write(*,*) 'setcnstr_trian_iso_hcurl: COLUMNWISE ....'
        do i=1,nc
          write(*,7003) evalsma(1:nc,i)
        enddo
        call pause
      endif
#endif
!
      nc = ic;
      call dgels('N', nc, nrdof, nrdof,   &
                 evalsma, npts, evalbig, npts, work, lwork, info)
      if (info.ne.0) then
        write(*,*)'setcnstr_trian_iso_hcurl: HCURL DGELS INFO =',info
        call logic_error(FAILURE,__FILE__,__LINE__)
      endif
!!!      r = 2.d0*evalbig(1:nrdof,1:nrdof)
      r = evalbig(1:nrdof,1:nrdof)
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7002) r(1:nrdof,1:nrdof)
      endif
#endif
!
!  ...degrees of freedom of small node 6:
!****************************************
      nrs = 2*ndof_medg
      do i=1,ndof_medg
!
!  .....big node 2,3,4:
        do ib=2,4
          nrb = (ib-2)*ndof_medg
          do j=1,ndof_medg
            RRTE(ib,j,6,i) = r(nrs+i,nrb+j)
          enddo
        enddo
!
!  .....big node 1:
        nrb = 3*ndof_medg
        do j=1,ndof_mdlt
          RRTE(1,j,6,i) = r(nrs+i,nrb+j)
        enddo
      enddo
!
!  ...degrees of freedom of small node 2:
!*****************************************
      nrs = 3*ndof_medg
      do i=1,ndof_mdlt
!
!  .....big node 2,3,4:
        do ib=2,4
          nrb = (ib-2)*ndof_medg
          do j=1,ndof_medg
            RRTE(ib,j,2,i) = r(nrs+i,nrb+j)
          enddo
        enddo
!
!  .....big node 1:
        nrb = 3*ndof_medg
        do j=1,ndof_mdlt
          RRTE(1,j,2,i) = r(nrs+i,nrb+j)
        enddo
      enddo
!
!*********************************************************************
!
!  ...small triangle 3 -- small nodes 3 and 7:
!
!*********************************************************************
!
      ic =0; nspan = (2*nord);
      dbds    = 0.d0;   curl    = 0.d0;
      shapsma = 0.d0;   shapbig = 0.d0;
      evalsma = 0.d0;   evalbig = 0.d0;
!  ...derivative dx/dxi
      dbds(1,1) = 0.5d0;
      dbds(2,2) = 0.5d0
!
!  ...loop through collocation points:
      do i=1,nspan
!
        xibig(1) = (i-1)*0.5d0/(nspan-1)
        xisma(1) = (i-1)*1.0d0/(nspan-1)
!
        do j=1,nspan-(i-1)
!
          xibig(2) = 0.5d0/(nspan-1)*(j-1) + 0.5d0
          xisma(2) = 1.0d0/(nspan-1)*(j-1)
!
!  .......find small shape functions at this point:
!          call shapeEt(                  &
!            xisma,norder,norient,ivoid,  &
!            shapsma(1:2,1:MAXtriaE),     &
!            curl(1:MAXtriaE))
          call shape2DETri(xisma,norder,norient,nsize,   &
                           ivoid,shapsma(1:2,:),curl)
!
!  .......find big shape functions at this point:
!          call shapeEt(                  &
!            xibig,norder,norient,ivoid,  &
!            shapbig(1:2,1:MAXtriaE),     &
!            curl(1:MAXtriaE))
          call shape2DETri(xibig,norder,norient,nsize,   &
                           ivoid,shapbig(1:2,:),curl)
!
!  .......for each component variable,
          if ((ic+2).gt.npts) then
            write(*,*) '3 setcnstr_trian_iso_hcurl: INCREASE WORKSPACE',npts
          endif
!
          do k=1,2
            evalbig(ic+k, 1:MAXtriaE) = shapbig(k,1:MAXtriaE)*dbds(k,1)  &
                                      + shapbig(k,1:MAXtriaE)*dbds(k,2)
            evalsma(ic+k, 1:MAXtriaE) = shapsma(k,1:MAXtriaE)
          enddo
          ic = ic + 2
        enddo
      enddo
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) 'ic, nrdof = ', ic, nrdof
        call pause
      endif
!
      if (iprint.eq.1) then
        write(*,*) 'setcnstr_trian: COLUMNWISE ....'
        do i=1,nc
          write(*,7003) evalsma(1:nc,i)
        enddo
        call pause
      endif
#endif
!
      nc = ic;
      call dgels('N', nc, nrdof, nrdof,   &
                 evalsma, npts, evalbig, npts, work, lwork, info)
      if (info.ne.0) then
        write(*,*)'setcnstr_trian_iso_hcurl: HCURL DGELS INFO =',info
        call logic_error(FAILURE,__FILE__,__LINE__)
      endif
!!!      r = 2.d0*evalbig(1:nrdof,1:nrdof)
      r = evalbig(1:nrdof,1:nrdof)
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7002) r(1:nrdof,1:nrdof)
      endif
#endif
!
!  ...degrees of freedom of small node 7:
!****************************************
      nrs = 0
      do i=1,ndof_medg
!
!  .....big node 2,3,4:
        do ib=2,4
          nrb = (ib-2)*ndof_medg
          do j=1,ndof_medg
            RRTE(ib,j,7,i) = r(nrs+i,nrb+j)
          enddo
        enddo
!
!  .....big node 1:
        nrb = 3*ndof_medg
        do j=1,ndof_mdlt
          RRTE(1,j,7,i) = r(nrs+i,nrb+j)
        enddo
      enddo
!
!  ...degrees of freedom of small node 3:
!*****************************************
      nrs = 3*ndof_medg
      do i=1,ndof_mdlt
!
!  .....big node 2,3,4:
        do ib=2,4
          nrb = (ib-2)*ndof_medg
          do j=1,ndof_medg
            RRTE(ib,j,3,i) = r(nrs+i,nrb+j)
          enddo
        enddo
!
!  .....big node 1:
        nrb = 3*ndof_medg
        do j=1,ndof_mdlt
          RRTE(1,j,3,i) = r(nrs+i,nrb+j)
        enddo
      enddo
!
!*********************************************************************
!
!  ...small triangle 4 -- small nodes 4:
!
!*********************************************************************
!
      ic =0; nspan = (2*nord);
      dbds    = 0.d0;   curl    = 0.d0;
      shapsma = 0.d0;   shapbig = 0.d0;
      evalsma = 0.d0;   evalbig = 0.d0;
!  ...derivative dx/dxi
      dbds(1,1) = -0.5d0;
      dbds(2,2) = -0.5d0
!
!  ...loop through collocation points:
      do i=1,nspan
!
        xibig(1) = 0.5d0 - (i-1)*0.5d0/(nspan-1)
        xisma(1) =         (i-1)*1.0d0/(nspan-1)
!
        do j=1,nspan-(i-1)
!
          xibig(2) = 0.5d0 - 0.5d0/(nspan-1)*(j-1)
          xisma(2) =         1.0d0/(nspan-1)*(j-1)
!
!  .......find small shape functions at this point:
!          call shapeEt(                  &
!            xisma,norder,norient,ivoid,  &
!            shapsma(1:2,1:MAXtriaE),     &
!            curl(1:MAXtriaE))
          call shape2DETri(xisma,norder,norient,nsize,   &
                           ivoid,shapsma(1:2,:),curl)
!
!  .......find big shape functions at this point:
!          call shapeEt(                  &
!            xibig,norder,norient,ivoid,  &
!            shapbig(1:2,1:MAXtriaE),     &
!            curl(1:MAXtriaE))
          call shape2DETri(xibig,norder,norient,nsize,   &
                           ivoid,shapbig(1:2,:),curl)
!
!  .......for each component variable,
          if ((ic+2).gt.npts) then
            write(*,*) '4 setcnstr_trian_iso_hcurl: INCREASE WORKSPACE',npts
          endif
!
          do k=1,2
            evalbig(ic+k, 1:MAXtriaE) = shapbig(k,1:MAXtriaE)*dbds(k,1)  &
                                      + shapbig(k,1:MAXtriaE)*dbds(k,2)
            evalsma(ic+k, 1:MAXtriaE) = shapsma(k,1:MAXtriaE)
          enddo
          ic = ic + 2
        enddo
      enddo
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) 'ic, nrdof = ', ic, nrdof
        call pause
      endif
!
      if (iprint.eq.1) then
        write(*,*) 'setcnstr_trian_iso_hcurl: COLUMNWISE ....'
        do i=1,nc
          write(*,7003) evalsma(1:nc,i)
        enddo
        call pause
      endif
#endif
!
      nc = ic;
      call dgels('N', nc, nrdof, nrdof,   &
                 evalsma, npts, evalbig, npts, work, lwork, info)
      if (info.ne.0) then
        write(*,*)'setcnstr_trian_iso_hcurl: HCURL DGELS INFO =',info
        call logic_error(FAILURE,__FILE__,__LINE__)
      endif
!!!      r = -2.d0*evalbig(1:nrdof,1:nrdof)
      r = evalbig(1:nrdof,1:nrdof)
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7002) r(1:nrdof,1:nrdof)
      endif
#endif
!
!  ...degrees of freedom of small node 4:
!*****************************************
      nrs = 3*ndof_medg
      do i=1,ndof_mdlt
!
!  .....big node 2,3,4:
        do ib=2,4
          nrb = (ib-2)*ndof_medg
          do j=1,ndof_medg
            RRTE(ib,j,4,i) = r(nrs+i,nrb+j)
          enddo
        enddo
!
!  .....big node 1:
        nrb = 3*ndof_medg
        do j=1,ndof_mdlt
          RRTE(1,j,4,i) = r(nrs+i,nrb+j)
        enddo
      enddo
!
!-----------------------------------------------------------------------
!
!  ...clean up machine zeros...
      ibcase = 0; iscase = 0;
      do ib=1,4
        do is=1,7
!
          select case(ib)
          case(1);     ndof_i = ndof_mdlt; ibcase = 0
          case(2,3,4); ndof_i = ndof_medg; ibcase = 1
          end select
!
          select case(is)
          case(1,2,3,4); ndof_j = ndof_mdlt; iscase = 0
          case(5,6,7);   ndof_j = ndof_medg; iscase = 1
          end select
!
          do i=1,ndof_i
            do j=1,ndof_j
              if (abs(RRTE(ib,i,is,j)).lt.1.d-12) then
                RRTE(ib,i,is,j) = 0.d0
              endif
            enddo
          enddo
        enddo
      enddo
!
      return
!
!*********************************************************************
!
#if HP3D_DEBUG
!
!  ...begin testing
  777 continue
!
      write(*,*) 'setcnstr_trian_iso_hcurl: SET xibig '
      read(*,*) xibig(1:2)
!
!  ...shape functions of big element
!      call shapeEt(                   &
!        xibig,norder,norient,ivoid,   &
!        shapbig(1:2,1:MAXtriaE),      &
!        curl(1:MAXtriaE))
      call shape2DETri(xibig,norder,norient,nsize, &
                       ivoid,shapbig(1:2,:),curl)
!
      write(*,*) 'xibig = ',xibig
!
!  ...which small element is this:
!  ...coordinates in a small element:
!
      dbds = 0.d0
      dbds(1,1) = 0.5d0
      dbds(2,2) = 0.5d0
      if (xibig(1)+xibig(2).le.0.5d0) then
!
!  .....element 1 - nodes 1,4
        nel = 1
        xisma(1) = 2*xibig(1)
        xisma(2) = 2*xibig(2)
!
!
      elseif (xibig(1).ge.0.5d0) then
!
!  .....element 2 - nodes 2,5:
        nel = 2
        xisma(1) = (xibig(1)-0.5d0)*2.d0
        xisma(2) =  xibig(2)*2.d0
!
      elseif (xibig(2).gt.0.5d0) then
!
!  .....element 3 - nodes 3,6:
        nel = 3
        xisma(2) = (xibig(2)-0.5d0)*2.d0
        xisma(1) =  xibig(1)*2.d0
!
      else
!
!  .....element 4:
        nel = 4
        xisma(1) = (0.5d0-xibig(1))*2
        xisma(2) = (0.5d0-xibig(2))*2
        dbds(1:2,1:2) = -1.d0*dbds(1:2,1:2)
      endif

      write(*,*) 'setcnstr_trian_iso_hcurl: nel = ', nel
!
!  ...find small shape functions at this point:
!      call shapeEt(
!        xisma,norder,norient,ivoid,   &
!        shapsma(1:2,1:MAXtriaE),      &
!        curl(1:MAXtriaE))             &
      call shape2DETri(xisma,norder,norient,nsize, &
                       ivoid,shapsma(1:2,:),curl)
!
!*********************************************************************
!
!  ...verify if big shape functions are right combinations of small
!     ones:
!
!  ...loop through big shape functions - central node:
      write(*,*) 'testing = ', ndof_mdlt
      do i=1,ndof_mdlt
!
!  .....initiate value of the linear combination:
        val = 0.d0
!
        select case(nel)
        case(1)
!
!  .......small node 5:
          nrs = ndof_medg
          do j=1,ndof_medg
            val(1:2) = val(1:2)+2.d0*RRTE(1,i, 5,j)*shapsma(1:2,nrs+j)
          enddo
!
!  .......small node 1:
          nrs = 3*ndof_medg
          do j=1,ndof_mdlt
            val(1:2) = val(1:2)+2.d0*RRTE(1,i, 1,j)*shapsma(1:2,nrs+j)
          enddo
!
        case(2)
!
!  .......small node 6:
          nrs = 2*ndof_medg
          do j=1,ndof_medg
            val(1:2) = val(1:2)+2.d0*RRTE(1,i, 6,j)*shapsma(1:2,nrs+j)
          enddo
!
!  .......small node 2:
          nrs = 3*ndof_medg
          do j=1,ndof_mdlt
            val(1:2) = val(1:2)+2.d0*RRTE(1,i, 2,j)*shapsma(1:2,nrs+j)
          enddo
!
        case(3)
!
!  .......small node 7:
          nrs = 0
          do j=1,ndof_medg
            val(1:2) = val(1:2)+2.d0*RRTE(1,i, 7,j)*shapsma(1:2,nrs+j)
          enddo
!
!  .......small node 3:
          nrs = 3*ndof_medg
          do j=1,ndof_mdlt
            val(1:2) = val(1:2)+2.d0*RRTE(1,i, 3,j)*shapsma(1:2,nrs+j)
          enddo
!
        case(4)
!
!  ........small node 5:
           nrs = ndof_medg
           do j=1,ndof_medg
             val(1:2) = val(1:2)+2.d0*RRTE(1,i, 5,j)  &
               * shapsma(1:2,nrs+j)*(-1)**(j+1)
           enddo
!
!  ........small node 6:
           nrs = 2*ndof_medg
           do j=1,ndof_medg
             val(1:2) = val(1:2)+2.d0*RRTE(1,i, 6,j)  &
               * shapsma(1:2,nrs+j)*(-1)**(j+1)

           enddo
!
!  ........small node 7:
           nrs = 0
           do j=1,ndof_medg
             val(1:2) = val(1:2)+2.d0*RRTE(1,i, 7,j)  &
               * shapsma(1:2,nrs+j)*(-1)**(j+1)
           enddo
!
!  ........small node 4:
           nrs = 3*ndof_medg
           do j=1,ndof_mdlt
             val(1:2) = val(1:2)-2.d0*RRTE(1,i, 4,j)*shapsma(1:2,nrs+j)
           enddo
!
         end select
!
         nrb = 3*ndof_medg
!
         write(*,8002) i,shapbig(1:2,nrb+i), &
           abs(val(1:2)-shapbig(1:2,nrb+i))
 8002    format('setcnstr_trian_iso_hcurl: i,shapebig,difference = ',   &
           i3,2(2e12.5,2x))
!
!  ...end of loop through shape functions of big element
      enddo
!
      write(*,*) 'setcnstr_trian_iso_hcurl: CONTINUE ?(1/0)'
      read(*,*) ians
!
      if (ians.eq.1) go to 777
!
#endif
!
   end subroutine setcnstr_trian_iso_hcurl

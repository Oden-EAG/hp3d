c----------------------------------------------------------------------
c
c   routine name       - setcnstr_trian_iso
c
c----------------------------------------------------------------------
c
c   latest revision    - May 2020
c
c   purpose            - routine evaluates constraint coefficients for
c                        the triangular master element
c
c   usage              - call setcnstrt
c
c   arguments          - none
c
c----------------------------------------------------------------------
c
      subroutine setcnstr_trian_iso_hcurl
c
c    'big' parent nodes
c
c      *
c      * *
c      *   *
c      *     *
c      *       *
c      *         *
c    4 *           *  3
c      *             *
c      *       1       *
c      *                 *
c      *                   *
c      *                     *
c      *************************
c                  2
c
c      'small' constrained nodes
c      *
c      * *
c      *   *
c      *     *
c      *   3   *
c      *         *
c      ******7******
c      * *         * *
c      *   *   4   *   *
c      *     5     6     *
c      *   1   *   *   2   *
c      *         * *         *
c      *************************
c
c--------------------------------------------------------------------
c
      use parameters
      use constraints
      use element_data
#include "syscom.blk"
c  ...use collocation points for lagrange shape functions h1
      integer, parameter :: npts = 4*(MAXtriaE)
      integer, parameter :: lwork = MAXtriaE**2
      dimension shapsma(2,MAXtriaE),curl(MAXtriaE),
     .          shapbig(2,MAXtriaE),
     .          evalsma(npts,MAXtriaE),evalbig(npts,MAXtriaE),
     .          work(lwork),r(MAXtriaE,MAXtriaE),xibig(2),xisma(2)
c
c  ...order and edge orientations
      dimension norder(4),norient(3),nsize(2),val(2),dbds(2,2)
c
      iprint=0
c
c  ...set up order and orientation
      nord = MAXP
      norder  = nord
      norient = 0
      nsize   = (/MAXP,MAXtriaE/)
c
      ivoid = 0
      call ndof_nod(MEDG,nord, ivoid, ndof_medg, ivoid, ivoid)
      call ndof_nod(MDLT,nord, ivoid, ndof_mdlt, ivoid, ivoid)
      nrdof = 3*ndof_medg + ndof_mdlt
c
c*********************************************************************
c
c  ...small triangle 1 -- small nodes 1 and 5:
c
c*********************************************************************
c
      ic =0; nspan = (2*nord);
      dbds    = 0.d0;   curl    = 0.d0;
      shapsma = 0.d0;   shapbig = 0.d0;
      evalsma = 0.d0;   evalbig = 0.d0;
c  ...derivative dx/dxi
      dbds(1,1) = 0.5d0;
      dbds(2,2) = 0.5d0
c
c  ...loop through collocation points:
      do i=1,nspan
c
        xibig(1) = (i-1)*0.5d0/(nspan-1)
        xisma(1) = (i-1)*1.0d0/(nspan-1)
c
        do j=1,nspan-(i-1)
c
          xibig(2) = 0.5d0/(nspan-1)*(j-1)
          xisma(2) = 1.0d0/(nspan-1)*(j-1)
c
c  .......find small shape functions at this point:
c          call shapeEt(
c     .      xisma,norder,norient,ivoid,
c     .      shapsma(1:2,1:MAXtriaE),
c     .      curl(1:MAXtriaE))
          call shape2DETri(xisma,norder,norient,nsize,
     .                     ivoid,shapsma(1:2,:),curl)
c
c  .......find big shape functions at this point:
c          call shapeEt(
c     .      xibig,norder,norient,ivoid,
c     .      shapbig(1:2,1:MAXtriaE),
c     .      curl(1:MAXtriaE))
          call shape2DETri(xibig,norder,norient,nsize,
     .                     ivoid,shapbig(1:2,:),curl)
c  .......for each component variable,
          if ((ic+2).gt.npts) then
            write(*,*) '1 setcnstr_trian_iso_hcurl: INCREASE WORKSPACE',
     .        npts, (ic+2)
          endif
c
          do k=1,2
            evalbig(ic+k, 1:MAXtriaE) = shapbig(k,1:MAXtriaE)*dbds(k,1)
     .                                + shapbig(k,1:MAXtriaE)*dbds(k,2)
            evalsma(ic+k, 1:MAXtriaE) = shapsma(k,1:MAXtriaE)
          enddo
          ic = ic + 2
        enddo
      enddo
c
      if (iprint.eq.1) then
        write(*,*) 'ic, nrdof = ', ic, nrdof
        call pause
      endif
c
      if (iprint.eq.1) then
        write(*,*) 'setcnstr_trian: COLUMNWISE ....'
        do i=1,nc
          write(*,7003) evalsma(1:nc,i)
 7003     format(15(e12.5,2x))
        enddo
        call pause
      endif
c
      nc = ic;
      call dgels('N', nc, nrdof, nrdof,
     .           evalsma, npts, evalbig, npts, work, lwork, info)
      if (info.ne.0) then
        write(*,*)'setcnstr_trian_iso: HCURL DGELS INFO =',info
        call logic_error(FAILURE,
     .  __FILE__,__LINE__)
      endif
ccc      r = 2.d0*evalbig(1:nrdof,1:nrdof)
      r = evalbig(1:nrdof,1:nrdof)
      if (iprint.eq.1) then
        write(*,7002) r(1:nrdof,1:nrdof)
 7002   format('r=', 10e12.5)
      endif

c
c  ...degrees of freedom of small node 5:
c****************************************
      nrs = ndof_medg
      do i=1,ndof_medg
c
c  .....big node 2,3,4:
        do ib=2,4
          nrb = (ib-2)*ndof_medg
          do j=1,ndof_medg
            RRTE(ib,j,5,i) = r(nrs+i,nrb+j)
          enddo
        enddo
c
c  .....big node 1:
        nrb = 3*ndof_medg
        do j=1,ndof_mdlt
          RRTE(1,j,5,i) = r(nrs+i,nrb+j)
        enddo
      enddo
c
c  ...degrees of freedom of small node 1:
c*****************************************
      nrs = 3*ndof_medg
      do i=1,ndof_mdlt
c
c  .....big node 2,3,4:
        do ib=2,4
          nrb = (ib-2)*ndof_medg
          do j=1,ndof_medg
            RRTE(ib,j,1,i) = r(nrs+i,nrb+j)
          enddo
        enddo
c
c  .....big node 1:
        nrb = 3*ndof_medg
        do j=1,ndof_mdlt
          RRTE(1,j,1,i) = r(nrs+i,nrb+j)
        enddo
      enddo
c*********************************************************************
c
c  ...small triangle 2 -- small nodes 2 and 6:
c
c*********************************************************************
c
      ic =0; nspan = (2*nord);
      dbds    = 0.d0;   curl    = 0.d0;
      shapsma = 0.d0;   shapbig = 0.d0;
      evalsma = 0.d0;   evalbig = 0.d0;
c  ...derivative dx/dxi
      dbds(1,1) = 0.5d0;
      dbds(2,2) = 0.5d0
c
c  ...loop through collocation points:
      do i=1,nspan
c
        xibig(1) = (i-1)*0.5d0/(nspan-1) + 0.5d0
        xisma(1) = (i-1)*1.0d0/(nspan-1)
c
        do j=1,nspan-(i-1)
c
          xibig(2) = 0.5d0/(nspan-1)*(j-1)
          xisma(2) = 1.0d0/(nspan-1)*(j-1)
c
c  .......find small shape functions at this point:
c          call shapeEt(
c     .      xisma,norder,norient,ivoid,
c     .      shapsma(1:2,1:MAXtriaE),
c     .      curl(1:MAXtriaE))
          call shape2DETri(xisma,norder,norient,nsize,
     .                     ivoid,shapsma(1:2,:),curl)
c
c  .......find big shape functions at this point:
c          call shapeEt(
c     .      xibig,norder,norient,ivoid,
c     .      shapbig(1:2,1:MAXtriaE),
c     .      curl(1:MAXtriaE))
          call shape2DETri(xibig,norder,norient,nsize,
     .                     ivoid,shapbig(1:2,:),curl)
c
c  .......for each component variable,
          if ((ic+2).gt.npts) then
            write(*,*) '2 setcnstr_trian_iso_hcurl: INCREASE WORKSPACE',
     .        npts
          endif
c
          do k=1,2
            evalbig(ic+k, 1:MAXtriaE) = shapbig(k,1:MAXtriaE)*dbds(k,1)
     .                                + shapbig(k,1:MAXtriaE)*dbds(k,2)
            evalsma(ic+k, 1:MAXtriaE) = shapsma(k,1:MAXtriaE)
          enddo
          ic = ic + 2
        enddo
      enddo
c
      if (iprint.eq.1) then
        write(*,*) 'ic, nrdof = ', ic, nrdof
        call pause
      endif
c
      if (iprint.eq.1) then
        write(*,*) 'setcnstr_trian: COLUMNWISE ....'
        do i=1,nc
          write(*,7003) evalsma(1:nc,i)
        enddo
        call pause
      endif
c
      nc = ic;
      call dgels('N', nc, nrdof, nrdof,
     .           evalsma, npts, evalbig, npts, work, lwork, info)
      if (info.ne.0) then
        write(*,*)'setcnstr_trian_iso: HCURL DGELS INFO =',info
        call logic_error(FAILURE,
     .  __FILE__,__LINE__)
      endif
ccc      r = 2.d0*evalbig(1:nrdof,1:nrdof)
      r = evalbig(1:nrdof,1:nrdof)
      if (iprint.eq.1) then
        write(*,7002) r(1:nrdof,1:nrdof)
      endif

c
c  ...degrees of freedom of small node 6:
c****************************************
      nrs = 2*ndof_medg
      do i=1,ndof_medg
c
c  .....big node 2,3,4:
        do ib=2,4
          nrb = (ib-2)*ndof_medg
          do j=1,ndof_medg
            RRTE(ib,j,6,i) = r(nrs+i,nrb+j)
          enddo
        enddo
c
c  .....big node 1:
        nrb = 3*ndof_medg
        do j=1,ndof_mdlt
          RRTE(1,j,6,i) = r(nrs+i,nrb+j)
        enddo
      enddo
c
c  ...degrees of freedom of small node 2:
c*****************************************
      nrs = 3*ndof_medg
      do i=1,ndof_mdlt
c
c  .....big node 2,3,4:
        do ib=2,4
          nrb = (ib-2)*ndof_medg
          do j=1,ndof_medg
            RRTE(ib,j,2,i) = r(nrs+i,nrb+j)
          enddo
        enddo
c
c  .....big node 1:
        nrb = 3*ndof_medg
        do j=1,ndof_mdlt
          RRTE(1,j,2,i) = r(nrs+i,nrb+j)
        enddo
      enddo
c
c*********************************************************************
c
c  ...small triangle 3 -- small nodes 3 and 7:
c
c*********************************************************************
c
      ic =0; nspan = (2*nord);
      dbds    = 0.d0;   curl    = 0.d0;
      shapsma = 0.d0;   shapbig = 0.d0;
      evalsma = 0.d0;   evalbig = 0.d0;
c  ...derivative dx/dxi
      dbds(1,1) = 0.5d0;
      dbds(2,2) = 0.5d0
c
c  ...loop through collocation points:
      do i=1,nspan
c
        xibig(1) = (i-1)*0.5d0/(nspan-1)
        xisma(1) = (i-1)*1.0d0/(nspan-1)
c
        do j=1,nspan-(i-1)
c
          xibig(2) = 0.5d0/(nspan-1)*(j-1) + 0.5d0
          xisma(2) = 1.0d0/(nspan-1)*(j-1)
c
c  .......find small shape functions at this point:
c          call shapeEt(
c     .      xisma,norder,norient,ivoid,
c     .      shapsma(1:2,1:MAXtriaE),
c     .      curl(1:MAXtriaE))
          call shape2DETri(xisma,norder,norient,nsize,
     .                     ivoid,shapsma(1:2,:),curl)
c
c  .......find big shape functions at this point:
c          call shapeEt(
c     .      xibig,norder,norient,ivoid,
c     .      shapbig(1:2,1:MAXtriaE),
c     .      curl(1:MAXtriaE))
          call shape2DETri(xibig,norder,norient,nsize,
     .                     ivoid,shapbig(1:2,:),curl)
c
c  .......for each component variable,
          if ((ic+2).gt.npts) then
            write(*,*) '3 setcnstr_trian_iso_hcurl: INCREASE WORKSPACE',
     .        npts
          endif
c
          do k=1,2
            evalbig(ic+k, 1:MAXtriaE) = shapbig(k,1:MAXtriaE)*dbds(k,1)
     .                                + shapbig(k,1:MAXtriaE)*dbds(k,2)
            evalsma(ic+k, 1:MAXtriaE) = shapsma(k,1:MAXtriaE)
          enddo
          ic = ic + 2
        enddo
      enddo
c
      if (iprint.eq.1) then
        write(*,*) 'ic, nrdof = ', ic, nrdof
        call pause
      endif
c
      if (iprint.eq.1) then
        write(*,*) 'setcnstr_trian: COLUMNWISE ....'
        do i=1,nc
          write(*,7003) evalsma(1:nc,i)
        enddo
        call pause
      endif
c
      nc = ic;
      call dgels('N', nc, nrdof, nrdof,
     .           evalsma, npts, evalbig, npts, work, lwork, info)
      if (info.ne.0) then
        write(*,*)'setcnstr_trian_iso: HCURL DGELS INFO =',info
        call logic_error(FAILURE,
     .  __FILE__,__LINE__)
      endif
ccc      r = 2.d0*evalbig(1:nrdof,1:nrdof)
      r = evalbig(1:nrdof,1:nrdof)
      if (iprint.eq.1) then
        write(*,7002) r(1:nrdof,1:nrdof)
      endif

c
c  ...degrees of freedom of small node 7:
c****************************************
      nrs = 0
      do i=1,ndof_medg
c
c  .....big node 2,3,4:
        do ib=2,4
          nrb = (ib-2)*ndof_medg
          do j=1,ndof_medg
            RRTE(ib,j,7,i) = r(nrs+i,nrb+j)
          enddo
        enddo
c
c  .....big node 1:
        nrb = 3*ndof_medg
        do j=1,ndof_mdlt
          RRTE(1,j,7,i) = r(nrs+i,nrb+j)
        enddo
      enddo
c
c  ...degrees of freedom of small node 3:
c*****************************************
      nrs = 3*ndof_medg
      do i=1,ndof_mdlt
c
c  .....big node 2,3,4:
        do ib=2,4
          nrb = (ib-2)*ndof_medg
          do j=1,ndof_medg
            RRTE(ib,j,3,i) = r(nrs+i,nrb+j)
          enddo
        enddo
c
c  .....big node 1:
        nrb = 3*ndof_medg
        do j=1,ndof_mdlt
          RRTE(1,j,3,i) = r(nrs+i,nrb+j)
        enddo
      enddo
c
c*********************************************************************
c
c  ...small triangle 4 -- small nodes 4:
c
c*********************************************************************
c
      ic =0; nspan = (2*nord);
      dbds    = 0.d0;   curl    = 0.d0;
      shapsma = 0.d0;   shapbig = 0.d0;
      evalsma = 0.d0;   evalbig = 0.d0;
c  ...derivative dx/dxi
      dbds(1,1) = -0.5d0;
      dbds(2,2) = -0.5d0
c
c  ...loop through collocation points:
      do i=1,nspan
c
        xibig(1) = 0.5d0 - (i-1)*0.5d0/(nspan-1)
        xisma(1) =         (i-1)*1.0d0/(nspan-1)
c
        do j=1,nspan-(i-1)
c
          xibig(2) = 0.5d0 - 0.5d0/(nspan-1)*(j-1)
          xisma(2) =         1.0d0/(nspan-1)*(j-1)
c
c  .......find small shape functions at this point:
c          call shapeEt(
c     .      xisma,norder,norient,ivoid,
c     .      shapsma(1:2,1:MAXtriaE),
c     .      curl(1:MAXtriaE))
          call shape2DETri(xisma,norder,norient,nsize,
     .                     ivoid,shapsma(1:2,:),curl)
c
c  .......find big shape functions at this point:
c          call shapeEt(
c     .      xibig,norder,norient,ivoid,
c     .      shapbig(1:2,1:MAXtriaE),
c     .      curl(1:MAXtriaE))
          call shape2DETri(xibig,norder,norient,nsize,
     .                     ivoid,shapbig(1:2,:),curl)
c
c  .......for each component variable,
          if ((ic+2).gt.npts) then
            write(*,*) '4 setcnstr_trian_iso_hcurl: INCREASE WORKSPACE',
     .        npts
          endif
c
          do k=1,2
            evalbig(ic+k, 1:MAXtriaE) = shapbig(k,1:MAXtriaE)*dbds(k,1)
     .                                + shapbig(k,1:MAXtriaE)*dbds(k,2)
            evalsma(ic+k, 1:MAXtriaE) = shapsma(k,1:MAXtriaE)
          enddo
          ic = ic + 2
        enddo
      enddo
c
      if (iprint.eq.1) then
        write(*,*) 'ic, nrdof = ', ic, nrdof
        call pause
      endif
c
      if (iprint.eq.1) then
        write(*,*) 'setcnstr_trian: COLUMNWISE ....'
        do i=1,nc
          write(*,7003) evalsma(1:nc,i)
        enddo
        call pause
      endif
c
      nc = ic;
      call dgels('N', nc, nrdof, nrdof,
     .           evalsma, npts, evalbig, npts, work, lwork, info)
      if (info.ne.0) then
        write(*,*)'setcnstr_trian_iso: HCURL DGELS INFO =',info
        call logic_error(FAILURE,
     .  __FILE__,__LINE__)
      endif
ccc      r = -2.d0*evalbig(1:nrdof,1:nrdof)
      r = evalbig(1:nrdof,1:nrdof)
      if (iprint.eq.1) then
        write(*,7002) r(1:nrdof,1:nrdof)
      endif

c
c  ...degrees of freedom of small node 4:
c*****************************************
      nrs = 3*ndof_medg
      do i=1,ndof_mdlt
c
c  .....big node 2,3,4:
        do ib=2,4
          nrb = (ib-2)*ndof_medg
          do j=1,ndof_medg
            RRTE(ib,j,4,i) = r(nrs+i,nrb+j)
          enddo
        enddo
c
c  .....big node 1:
        nrb = 3*ndof_medg
        do j=1,ndof_mdlt
          RRTE(1,j,4,i) = r(nrs+i,nrb+j)
        enddo
      enddo
c
c-----------------------------------------------------------------------
c
c  ...clean up machine zeros...
      ibcase = 0; iscase = 0;
      do ib=1,4
        do is=1,7
c
          select case(ib)
          case(1);     ndof_i = ndof_mdlt; ibcase = 0
          case(2,3,4); ndof_i = ndof_medg; ibcase = 1
          end select
c
          select case(is)
          case(1,2,3,4); ndof_j = ndof_mdlt; iscase = 0
          case(5,6,7);   ndof_j = ndof_medg; iscase = 1
          end select
c
          do i=1,ndof_i
            do j=1,ndof_j
              if (abs(RRTE(ib,i,is,j)).lt.1.d-12) then
                RRTE(ib,i,is,j) = 0.d0
              endif
            enddo
          enddo
        enddo
      enddo
c
      return
c
c*********************************************************************
c
c  ...begin testing
  777 continue
c
      write(*,*) 'setcnstr_trian: SET xibig '
      read(*,*) xibig(1:2)
c
c  ...shape functions of big element
c      call shapeEt(
c     .  xibig,norder,norient,ivoid,
c     .  shapbig(1:2,1:MAXtriaE),
c     .  curl(1:MAXtriaE))
      call shape2DETri(xibig,norder,norient,nsize,
     .                 ivoid,shapbig(1:2,:),curl)
c
      write(*,*) 'xibig = ',xibig
c
c  ...which small element is this:
c  ...coordinates in a small element:
c
      dbds = 0.d0
      dbds(1,1) = 0.5d0
      dbds(2,2) = 0.5d0
      if (xibig(1)+xibig(2).le.0.5d0) then
c
c  .....element 1 - nodes 1,4
        nel = 1
        xisma(1) = 2*xibig(1)
        xisma(2) = 2*xibig(2)
c
c
      elseif (xibig(1).ge.0.5d0) then
c
c  .....element 2 - nodes 2,5:
        nel = 2
        xisma(1) = (xibig(1)-0.5d0)*2.d0
        xisma(2) =  xibig(2)*2.d0
c
      elseif (xibig(2).gt.0.5d0) then
c
c  .....element 3 - nodes 3,6:
        nel = 3
        xisma(2) = (xibig(2)-0.5d0)*2.d0
        xisma(1) =  xibig(1)*2.d0
c
      else
c
c  .....element 4:
        nel = 4
        xisma(1) = (0.5d0-xibig(1))*2
        xisma(2) = (0.5d0-xibig(2))*2
        dbds(1:2,1:2) = -1.d0*dbds(1:2,1:2)
      endif

      write(*,*) 'setcnstr_trian: nel = ', nel
c
c  ...find small shape functions at this point:
c      call shapeEt(
c     .  xisma,norder,norient,ivoid,
c     .  shapsma(1:2,1:MAXtriaE),
c     .  curl(1:MAXtriaE))
      call shape2DETri(xisma,norder,norient,nsize,
     .                 ivoid,shapsma(1:2,:),curl)
c
c*********************************************************************
c
c  ...verify if big shape functions are right combinations of small
c     ones:
c
c  ...loop through big shape functions - central node:
      write(*,*) 'testing = ', ndof_mdlt
      do i=1,ndof_mdlt
c
c  .....initiate value of the linear combination:
        val = 0.d0
c
        select case(nel)
        case(1)
c
c  .......small node 5:
          nrs = ndof_medg
          do j=1,ndof_medg
            val(1:2) = val(1:2)+2.d0*RRTE(1,i, 5,j)*shapsma(1:2,nrs+j)
          enddo
c
c  .......small node 1:
          nrs = 3*ndof_medg
          do j=1,ndof_mdlt
            val(1:2) = val(1:2)+2.d0*RRTE(1,i, 1,j)*shapsma(1:2,nrs+j)
          enddo
c
        case(2)
c
c  .......small node 6:
          nrs = 2*ndof_medg
          do j=1,ndof_medg
            val(1:2) = val(1:2)+2.d0*RRTE(1,i, 6,j)*shapsma(1:2,nrs+j)
          enddo
c
c  .......small node 2:
          nrs = 3*ndof_medg
          do j=1,ndof_mdlt
            val(1:2) = val(1:2)+2.d0*RRTE(1,i, 2,j)*shapsma(1:2,nrs+j)
          enddo
c
        case(3)
c
c  .......small node 7:
          nrs = 0
          do j=1,ndof_medg
            val(1:2) = val(1:2)+2.d0*RRTE(1,i, 7,j)*shapsma(1:2,nrs+j)
          enddo
c
c  .......small node 3:
          nrs = 3*ndof_medg
          do j=1,ndof_mdlt
            val(1:2) = val(1:2)+2.d0*RRTE(1,i, 3,j)*shapsma(1:2,nrs+j)
          enddo
c
        case(4)
c
c  ........small node 5:
           nrs = ndof_medg
           do j=1,ndof_medg
             val(1:2) = val(1:2)+2.d0*RRTE(1,i, 5,j)
     .         * shapsma(1:2,nrs+j)*(-1)**(j+1)
           enddo
c
c  ........small node 6:
           nrs = 2*ndof_medg
           do j=1,ndof_medg
             val(1:2) = val(1:2)+2.d0*RRTE(1,i, 6,j)
     .         * shapsma(1:2,nrs+j)*(-1)**(j+1)

           enddo
c
c  ........small node 7:
           nrs = 0
           do j=1,ndof_medg
             val(1:2) = val(1:2)+2.d0*RRTE(1,i, 7,j)
     .         * shapsma(1:2,nrs+j)*(-1)**(j+1)
           enddo
c
c  ........small node 4:
           nrs = 3*ndof_medg
           do j=1,ndof_mdlt
             val(1:2) = val(1:2)-2.d0*RRTE(1,i, 4,j)*shapsma(1:2,nrs+j)
           enddo
c
         end select
c
         nrb = 3*ndof_medg
c
         write(*,8002) i,shapbig(1:2,nrb+i),
     .     abs(val(1:2)-shapbig(1:2,nrb+i))
 8002    format('setcnstr_trian: i,shapebig,difference = ',
     .     i3,2(2e12.5,2x))
c
c   ...end of loop through shape functions of big element
       enddo
c
       write(*,*) 'setcnstr_trian: CONTINUE ?(1/0)'
       read(*,*) ians
c
       if (ians.eq.1) go to 777
c
       end subroutine setcnstr_trian_iso_hcurl

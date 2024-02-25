c----------------------------------------------------------------------
c
c   routine name       - setcnstr_trian_iso_h1
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - routine evaluates constraint coefficients for
c                        the triangular master element
c
c   arguments          - none
c
c----------------------------------------------------------------------
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
      subroutine setcnstr_trian_iso_h1
c
      use parameters
      use constraints
c
      implicit none
c
      real(8) :: shapsma(MAXtriaH,MAXtriaH),
     .           shapbig(MAXtriaH,MAXtriaH),
     .           void(2,MAXtriaH),
     .           r(MAXtriaH,MAXtriaH),xibig(2),xisma(2)
      integer :: ip(MAXtriaH)
c
c  ...order and edge orientations
      integer :: norder(4),norient(3),nsize(2)
c
      real(8) :: aux,val
      integer :: i,ipp,j,jbeg,jp,jpbeg,jpend,k,ians,icl,info
      integer :: n,nel,nrb,nord,nrdof,nrs
c
#if HP3D_DEBUG
      integer :: iprint
      iprint = 0
#endif
c
      nord = MAXP
c
c  ...set uniform order MAXP:
      norder  = nord
      norient = 0
      nsize   = (/MAXP,MAXtriaH/)
c
c  ...initiate the matrix of coefficients
cccc      RRTH = 0.d0
c
c*********************************************************************
c
c  ...small triangle 1 -- small nodes 1 and 5:
c
c*********************************************************************
c
      icl = 0
c
c  ...loop through collocation points:
      do i=1,nord+1
c
        xibig(1) = (i-1)*0.5d0/nord
        xisma(1) = (i-1)*1.0d0/nord
c
        do j=1,nord+1 -(i-1)
          icl = icl+1
c
          xibig(2) = 0.5d0/nord*(j-1)
          xisma(2) = 1.0d0/nord*(j-1)
c
c  .......find small shape functions at this point:
c          call shapeHt(xisma,norder,norient,
c     .                 nrdof,shapsma(1:MAXtriaH,icl),void)
          call shape2DHTri(xisma,norder,norient,nsize,
     .                     nrdof,shapsma(:,icl),void)
c
c  .......find big shape functions at this point:
c          call shapeHt(xibig,norder,norient,
c     .                 nrdof,shapbig(1:MAXtriaH,icl),void)
          call shape2DHTri(xibig,norder,norient,nsize,
     .                     nrdof,shapbig(:,icl),void)
        enddo
      enddo
c
c  ...transpose shapsma:
      do i=1,(nord+1)*(nord+2)/2
        do j=1,i
          aux = shapsma(i,j)
          shapsma(i,j) = shapsma(j,i)
          shapsma(j,i) = aux
c
          aux = shapbig(i,j)
          shapbig(i,j) = shapbig(j,i)
          shapbig(j,i) = aux
        enddo
      enddo
c
      n = (nord+1)*(nord+2)/2
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) 'setcnstr_trian_iso_h1: shapsma COLUMNWISE ....'
        do i=1,n
          write(*,7004) i
 7004     format('setcnstr_trian_iso_h1: i = ',i2)
          write(*,7003) shapsma(1:n,i)
 7003     format(10e12.5)
        enddo
        call pause
      endif
#endif
ccc      call decomp(n,na,shapsma,ip,iflag)
cccc
ccc      if (iflag.ne.0) then
ccc        write(*,*) 'setcnstr_trian: iflag = ',iflag
ccc        stop 1
ccc      endif
cccc
ccc      do i=1,n
ccc        call gauss2(n,na,shapsma,ip,shapbig(1:n,i), r(1:n,i))
ccc#if HP3D_DEBUG
ccc        if (iprint.eq.1) then
ccc          write(*,*) 'setcnstr_trian: r FOR i = ',i
ccc          write(*,7003) r(1:n,i)
ccc        endif
ccc#endif
ccc      enddo


      call dgetrf(n,n,shapsma,MAXtriaH,ip,info)
      if (info.ne.0) then
        write(*,*)'setcnstr_trian_iso_h1: H1 DGETRF RETURNED INFO =',
     .             info
        call logic_error(FAILURE,
     .  __FILE__,__LINE__)
      endif
c
c  ...right-hand side resolution
      call dlaswp(n,shapbig,n,1,MAXtriaH,ip,1)
      call dtrsm('L','L','N','U',n,n,1.d0,shapsma,MAXtriaH,shapbig,
     .            MAXtriaH)
      call dtrsm('L','U','N','N',n,n,1.d0,shapsma,MAXtriaH,shapbig,
     .            MAXtriaH)
      r = shapbig



c
c  ...degrees of freedom of small node 5:
c     **********************************

      nrs = 3+(nord-1)
      do i=1,nord-1
c
c  .....big node 2:
        nrb = 3
        do j=1,nord-1
          RRTH(2,j,5,i) = r(nrs+i,nrb+j)
        enddo
c
c  .....big node 3:
        nrb = 3+(nord-1)
        do j=1,nord-1
          RRTH(3,j,5,i) = r(nrs+i,nrb+j)
        enddo
c
c  .....big node 4:
        nrb = 3+2*(nord-1)
        do j=1,nord-1
          RRTH(4,j,5,i) = r(nrs+i,nrb+j)
        enddo
c
c  .....big node 1:
        nrb = 3*nord
        do j=1,(nord-1)*(nord-2)/2
          RRTH(1,j,5,i) = r(nrs+i,nrb+j)
        enddo
      enddo
c
c  ...degrees of freedom of small node 1:
c*****************************************
c
      nrs = 3*nord
      do i=1,(nord-1)*(nord-2)/2
c
c  .....big node 2:
        nrb = 3
        do j=1,nord-1
          RRTH(2,j,1,i) = r(nrs+i,nrb+j)
        enddo
c
c  .....big node 3:
        nrb = 3+(nord-1)
        do j=1,nord-1
          RRTH(3,j,1,i) = r(nrs+i,nrb+j)
        enddo
c
c  .....big node 4:
        nrb = 3+2*(nord-1)
        do j=1,nord-1
          RRTH(4,j,1,i) = r(nrs+i,nrb+j)
        enddo
c
c  .....big node 1:
        nrb = 3*nord
        do j=1,(nord-1)*(nord-2)/2
          RRTH(1,j,1,i) = r(nrs+i,nrb+j)
        enddo
c
       enddo
c
c*********************************************************************
c
c  ...small triangle 2 --- small nodes 2 and 6:
c
c*********************************************************************
c
c  ...loop through collocation points:
      icl = 0
      do i=1,nord+1
c
        xibig(1) = (i-1)*0.5d0/nord + 0.5d0
        xisma(1) = (i-1)*1.0d0/nord
c
        do j=1,nord+1 -(i-1)
          icl = icl+1
c
          xibig(2) = 0.5d0/nord*(j-1)
          xisma(2) = 1.0d0/nord*(j-1)
c
c  .......find small shape functions at this point:
c          call shapeHt(xisma,norder,norient,
c     .                 nrdof,shapsma(1:MAXtriaH,icl),void)
          call shape2DHTri(xisma,norder,norient,nsize,
     .                     nrdof,shapsma(:,icl),void)
c
c  .......find big   shape functions at this point:
c          call shapeHt(xibig,norder,norient,
c     .                 nrdof,shapbig(1:MAXtriaH,icl),void)
          call shape2DHTri(xibig,norder,norient,nsize,
     .                     nrdof,shapbig(:,icl),void)
        enddo
      enddo
c
c  ...transpose shapsma:
      do i=1,(nord+1)*(nord+2)/2
        do j=1,i
          aux = shapsma(i,j)
          shapsma(i,j) = shapsma(j,i)
          shapsma(j,i) = aux
c
          aux = shapbig(i,j)
          shapbig(i,j) = shapbig(j,i)
          shapbig(j,i) = aux
        enddo
      enddo
c
      n = (nord+1)*(nord+2)/2
ccc      call decomp(n,na,shapsma,ip,Iflag)
ccc      if (iflag.eq.1)then
ccc        write(*,*) 'setcnstr_trian: WARNING 2 !! iflag = ',iflag
ccc        call pause
ccc      endif
ccc      do i=1,n
ccc        call gauss2(n,na,shapsma,ip,shapbig(1:n,i), r(1:n,i))
ccc      enddo
c
c  ...decompose the matrix
      call dgetrf(n,n,shapsma,MAXtriaH,ip,info)
      if (info.ne.0) then
        write(*,*)'setcnstr_trian_iso_h1: H1 DGETRF RETURNED INFO =',
     .             info
        call logic_error(FAILURE,
     . __FILE__,__LINE__)
      endif
c
c  ...right-hand side resolution
      call dlaswp(n,shapbig,n,1,MAXtriaH,ip,1)
      call dtrsm('L','L','N','U',n,n,1.d0,shapsma,MAXtriaH,shapbig,
     .            MAXtriaH)
      call dtrsm('L','U','N','N',n,n,1.d0,shapsma,MAXtriaH,shapbig,
     .            MAXtriaH)
      r = shapbig
c
c
c  ...degrees of freedom of small node 6:
c     **********************************
c
      nrs = 3+2*(nord-1)
      do i=1,nord-1
c
c  .....big node 2:
        nrb = 3
        do j=1,nord-1
          RRTH(2,j,6,i) = r(nrs+i,nrb+j)
        enddo
c
c  .....big node 3:
        nrb = 3+(nord-1)
        do j=1,nord-1
          RRTH(3,j,6,i) = r(nrs+i,nrb+j)
        enddo
c
c  .....big node 4:
        nrb = 3+2*(nord-1)
        do j=1,nord-1
          RRTH(4,j,6,i) = r(nrs+i,nrb+j)
        enddo
c
c  .....big node 1:
        nrb = 3*nord
        do j=1,(nord-1)*(nord-2)/2
          RRTH(1,j,6,i) = r(nrs+i,nrb+j)
        enddo
c
      enddo
c
c  ...degrees of freedom of small node 2:
c*****************************************
c
      nrs = 3*nord
      do i=1,(nord-1)*(nord-2)/2
c
c  .....big node 2:
        nrb = 3
        do j=1,nord-1
          RRTH(2,j,2,i) = r(nrs+i,nrb+j)
        enddo
c
c  .....big node 3:
        nrb = 3+(nord-1)
        do j=1,nord-1
          RRTH(3,j,2,i) = r(nrs+i,nrb+j)
        enddo
c
c  .....big node 4:
        nrb = 3+2*(nord-1)
        do j=1,nord-1
          RRTH(4,j,2,i) = r(nrs+i,nrb+j)
        enddo
c
c  .....big node 1:
        nrb = 3*nord
        do j=1,(nord-1)*(nord-2)/2
          RRTH(1,j,2,i) = r(nrs+i,nrb+j)
        enddo
c
      enddo
c
c*********************************************************************
c
c  ...small triangle 3 --- small nodes 3 and 7:
c
c*********************************************************************
c
c  ...loop through collocation points:
      icl = 0
      do i=1,nord+1
c
        xibig(1) = (i-1)*0.5d0/nord
        xisma(1) = (i-1)*1.0d0/nord
c
        do j=1,nord+1 -(i-1)
          icl = icl + 1
c
          xibig(2) = 0.5d0/nord*(j-1) + 0.5d0
          xisma(2) = 1.0d0/nord*(j-1)
c
c  .......find small shape functions at this point:
c          call shapeHt(xisma,norder,norient,
c     .                 nrdof,shapsma(1:MAXtriaH,icl),void)
          call shape2DHTri(xisma,norder,norient,nsize,
     .                     nrdof,shapsma(:,icl),void)
c
c  .......find big   shape functions at this point:
c          call shapeHt(xibig,norder,norient,
c     .                 nrdof,shapbig(1:MAXtriaH,icl),void)
          call shape2DHTri(xibig,norder,norient,nsize,
     .                     nrdof,shapbig(:,icl),void)
        enddo
c
      enddo
c
c  ...transpose shapsma:
      do i=1,(nord+1)*(nord+2)/2
        do j=1,i
          aux = shapsma(i,j)
          shapsma(i,j) = shapsma(j,i)
          shapsma(j,i) = aux
c
          aux = shapbig(i,j)
          shapbig(i,j) = shapbig(j,i)
          shapbig(j,i) = aux
        enddo
      enddo
c
      n = (nord+1)*(nord+2)/2
ccc      call decomp(n,na,shapsma,ip,iflag)
ccc      if (iflag.eq.1)then
ccc        write(*,*) 'setcnstr_trian: WARNING 3 !! iflag = ',iflag
ccc        call pause
ccc      endif
cccc
ccc      do i=1,n
ccc        call gauss2(n,na,shapsma,ip,shapbig(1:n,i),r(1:n,i))
ccc      enddo
c
c  ...decompose the matrix
      call dgetrf(n,n,shapsma,MAXtriaH,ip,info)
      if (info.ne.0) then
        write(*,*)'setcnstr_trian_iso_h1: H1 DGETRF RETURNED INFO =',
     .             info
        call logic_error(FAILURE,
     .  __FILE__,__LINE__)
      endif
c
c  ...right-hand side resolution
      call dlaswp(n,shapbig,n,1,MAXtriaH,ip,1)
      call dtrsm('L','L','N','U',n,n,1.d0,shapsma,MAXtriaH,shapbig,
     .            MAXtriaH)
      call dtrsm('L','U','N','N',n,n,1.d0,shapsma,MAXtriaH,shapbig,
     .            MAXtriaH)
      r = shapbig
c
c
c  ...degrees of freedom of small node 7:
c     **********************************
c
      nrs = 3
      do i=1,nord-1
c
c  .....big node 2:
        nrb = 3
        do j=1,nord-1
          RRTH(2,j,7,i) = r(nrs+i,nrb+j)
        enddo
c
c  .....big node 3:
        nrb = 3+(nord-1)
        do j=1,nord-1
          RRTH(3,j,7,i) = r(nrs+i,nrb+j)
        enddo
c
c  .....big node 4:
        nrb = 3+2*(nord-1)
        do j=1,nord-1
          RRTH(4,j,7,i) = r(nrs+i,nrb+j)
        enddo
c
c  .....big node 1:
        nrb = 3*nord
        do j=1,(nord-1)*(nord-2)/2
          RRTH(1,j,7,i) = r(nrs+i,nrb+j)
        enddo
c
      enddo
c
c  ...degrees of freedom of small node 3:
c*****************************************
c
      nrs = 3*nord
      do i=1,(nord-1)*(nord-2)/2
c
c  .....big node 2:
        nrb = 3
        do j=1,nord-1
          RRTH(2,j,3,i) = r(nrs+i,nrb+j)
        enddo
c
c  .....big node 3:
        nrb = 3+(nord-1)
        do j=1,nord-1
          RRTH(3,j,3,i) = r(nrs+i,nrb+j)
        enddo
c
c  .....big node 4:
        nrb = 3+2*(nord-1)
        do j=1,nord-1
          RRTH(4,j,3,i) = r(nrs+i,nrb+j)
        enddo
c
c  .....big node 1:
        nrb = 3*nord
        do j=1,(nord-1)*(nord-2)/2
          RRTH(1,j,3,i) = r(nrs+i,nrb+j)
        enddo
c
      enddo
c
c*********************************************************************
c
c  ...small triangle 4 --- small node 4:
c
c*********************************************************************

c  ...loop through collocation points:
      icl = 0
      do i=1,nord+1
c
        xibig(1) = 0.5d0 - (i-1)*0.5d0/nord
        xisma(1) =         (i-1)*1.0d0/nord
c
        do j=1,nord+1 -(i-1)
          icl = icl + 1
c
          xibig(2) = 0.5d0 - 0.5d0/nord*(j-1)
          xisma(2) =         1.0d0/nord*(j-1)
c
c  .......find small shape functions at this point:
c          call shapeHt(xisma,norder,norient,
c     .                 nrdof,shapsma(1:MAXtriaH,icl),void)
          call shape2DHTri(xisma,norder,norient,nsize,
     .                     nrdof,shapsma(:,icl),void)
c
c  .......find big   shape functions at this point:
c          call shapeHt(xibig,norder,norient,
c     .                 nrdof,shapbig(1:MAXtriaH,icl),void)
          call shape2DHTri(xibig,norder,norient,nsize,
     .                     nrdof,shapbig(:,icl),void)
        enddo
c
      enddo
c
c  ...transpose shapsma:
      do i=1,(nord+1)*(nord+2)/2
        do j=1,i
          aux = shapsma(i,j)
          shapsma(i,j) = shapsma(j,i)
          shapsma(j,i) = aux
c
          aux = shapbig(i,j)
          shapbig(i,j) = shapbig(j,i)
          shapbig(j,i) = aux
        enddo
      enddo
c
      n = (nord+1)*(nord+2)/2
ccc      call decomp(n,na,shapsma,ip,iflag)
ccc      if (iflag.eq.1)then
ccc        write(*,*) 'setcnstr_trian: WARNING 4 !! iflag = ',iflag
ccc        call pause
ccc      endif
cccc
ccc      do i=1,n
ccc        call gauss2(n,na,shapsma,ip,shapbig(1:n,i),r(1:n,i))
ccc      enddo
c
c  ...decompose the matrix
      call dgetrf(n,n,shapsma,MAXtriaH,ip,info)
      if (info.ne.0) then
        write(*,*)'setcnstr_trian_iso_h1: H1 DGETRF RETURNED INFO =',
     .             info
        call logic_error(FAILURE,
     .  __FILE__,__LINE__)
      endif
c
c  ...right-hand side resolution
      call dlaswp(n,shapbig,n,1,MAXtriaH,ip,1)
      call dtrsm('L','L','N','U',n,n,1.d0,shapsma,MAXtriaH,shapbig,
     .            MAXtriaH)
      call dtrsm('L','U','N','N',n,n,1.d0,shapsma,MAXtriaH,shapbig,
     .            MAXtriaH)
      r = shapbig
c
c
c  ...degrees of freedom of small node 4:
c*****************************************
c
      nrs = 3*nord
      do i=1,(nord-1)*(nord-2)/2
c
c  .....big node 2:
        nrb = 3
        do j=1,nord-1
          RRTH(2,j,4,i) = r(nrs+i,nrb+j)
        enddo
c
c  .....big node 3:
        nrb = 3+(nord-1)
        do j=1,nord-1
          RRTH(3,j,4,i) = r(nrs+i,nrb+j)
        enddo
c
c  .....big node 4:
        nrb = 3+2*(nord-1)
        do j=1,nord-1
          RRTH(4,j,4,i) = r(nrs+i,nrb+j)
        enddo
c
c  .....big node 1:
        nrb = 3*nord
        do j=1,(nord-1)*(nord-2)/2
          RRTH(1,j,4,i) = r(nrs+i,nrb+j)
        enddo
c
      enddo
c
c-----------------------------------------------------------------------
c
c  ...clean up machine zeros...
c
c  ...loop through possible orders of approximations
      do n=2,nord
c
c  .....loop through parent nodes
        do ipp=1,4
c
c  .......beginning and ending dof
          select case(ipp)
          case(1)
            jpbeg = (n-2)*(n-3)/2 + 1; jpend = (n-1)*(n-2)/2
          case(2,3,4)
            jpbeg = n-1; jpend = n-1
          end select
c
c  .......loop through the parent dof
          do jp=jpbeg,jpend
c
c  .........loop through constrained nodes
            do i=1,7
c
c  ...........beginning of dof to be zero out
              select case(i)
              case(1,2,3,4)
                jbeg = (n-1)*(n-2)/2+1
              case(5,6,7)
                jbeg = n
              end select
c
c  ...........loop through the constrained node dof
              do j=jbeg,(nord-2)*(nord-1)/2
                if (abs(RRTH(ipp,jp,i,j)).gt.1.d-12) then
                  write(*,7001) n,ipp,jp,i,j,RRTH(ipp,jp,i,j)
 7001             format('setcnstr_trian_iso_h1: n,ipp,jp,i,j,RRTH = ',
     .                    5i3,2x,e12.5)
                  stop 1
                endif
                RRTH(ipp,jp,i,j) = 0.d0
              enddo
            enddo
          enddo
        enddo
      enddo
c
c
      return
c
c*********************************************************************
c
#if HP3D_DEBUG
c
c  ...begin testing
  777 continue
c
      write(*,*) 'setcnstr_trian_iso_h1: SET xibig '
      read(*,*) xibig(1:2)
c
c  ...shape functions of big element
c      call shapeHt(xibig,norder,norient,
c     .             nrdof,shapbig(1:MAXtriaH,1),void)
      call shape2DHTri(xibig,norder,norient,nsize,
     .                 nrdof,shapbig(:,1),void)
      write(*,*) 'xibig = ',xibig
      write(*,*) 'shapbig = '
      do k=1,nrdof
        write(*,*) k,shapbig(k,1)
      enddo
c
c  ...which small element is this:
c  ...coordinates in a small element:
c
      if (xibig(1)+xibig(2).le.0.5d0) then
c
c  .....element 1 - nodes 1,4
        nel = 1
        xisma(1) = 2*xibig(1)
        xisma(2) = 2*xibig(2)
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
      endif
c
c  ...find small shape functions at this point:
c      call shapeHt(xisma,norder,norient,
c     .             nrdof,shapsma(1:MAXtriaH,1),void)
      call shape2DHTri(xisma,norder,norient,nsize,
     .                 nrdof,shapsma(:,1),void)
c
c*********************************************************************
c
c  ...verify if big shape functions are right combinations of small
c     ones:
c
c  ...loop through big shape functions - central node:
      nrb = 3*nord
      do i=1,(nord-1)*(nord-2)/2
c
c  .....initiate value of the linear combination:
        val = 0.d0
c
        select case(nel)
        case(1)
c
c  .......small node 5:
          nrs = 3 + (nord-1)
          do j=1,nord-1
            val = val + RRTH(1,i, 5,j)*shapsma(nrs+j,1)
          enddo
c
c  .......small node 1:
          nrs = 3*nord
          do j=1,(nord-1)*(nord-2)/2
            val = val + RRTH(1,i, 1,j)*shapsma(nrs+j,1)
          enddo
c
        case(2)
c
c  .......small node 6:
          nrs = 3 + (nord-1)*2
          do j=1,nord-1
            val = val + RRTH(1,i, 6,j)*shapsma(nrs+j,1)
          enddo
c
c  .......small node 2:
          nrs = 3*nord
          do j=1,(nord-1)*(nord-2)/2
            val = val + RRTH(1,i, 2,j)*shapsma(nrs+j,1)
          enddo
c
        case(3)
c
c  .......small node 7:
          nrs = 3
          do j=1,nord-1
            val = val + RRTH(1,i, 7,j)*shapsma(nrs+j,1)
          enddo
c
c  .......small node 3:
          nrs = 3*nord
          do j=1,(nord-1)*(nord-2)/2
            val = val + RRTH(1,i, 3,j)*shapsma(nrs+j,1)
          enddo
c
        case(4)
c
c  ........small node 5:
           nrs = 3 + (nord-1)
           do j=1,nord-1
             val = val + RRTH(1,i, 5,j)*shapsma(nrs+j,1)*(-1)**(j+1)
           enddo
c
c  ........small node 6:
           nrs = 3 + 2*(nord-1)
           do j=1,nord-1
             val = val + RRTH(1,i, 6,j)*shapsma(nrs+j,1)*(-1)**(j+1)
           enddo
c
c  ........small node 7:
           nrs = 3
           do j=1,nord-1
             val = val + RRTH(1,i, 7,j)*shapsma(nrs+j,1)*(-1)**(j+1)
           enddo
c
c  ........small node 4:
           nrs = 3*nord
           do j=1,(nord-1)*(nord-2)/2
             val = val + RRTH(1,i, 4,j)*shapsma(nrs+j,1)
           enddo
c
         end select
c
c
         write(*,7002) i,shapbig(nrb+i,1),abs(val-shapbig(nrb+i,1))
 7002    format('setcnstr_trian_iso_h1: i,shapebig,difference = '
     . ,i3,2e12.5)

c
c  ...end of loop through shape functions of big element
      enddo
c
      write(*,*) 'setcnstr_trian_iso_h1: CONTINUE ?(1/0)'
      read(*,*) ians
      if (ians.eq.1) go to 777
c
#endif
c
      end subroutine setcnstr_trian_iso_h1

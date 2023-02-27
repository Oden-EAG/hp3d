c----------------------------------------------------------------------
c
c   routine name       - setcnstr_trian_aniso
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
      subroutine setcnstr_trian_aniso_h1(Iref)
c
c     'big' parent nodes
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
c      Iref = 2; 'small' constrained nodes
c      *
c      * *
c      *   *
c      *     *
c      *       *
c      *         *
c      *           *
c      * *      1    *
c      *   *           *
c      *     ?           *
c      *   ?   *           *
c      *         *           *
c      *************************
c
c      Iref = 3; 'small' constrained nodes
c      *
c      * *
c      *   *
c      *     *
c      *       *
c      *         *
c      *         * *
c      *    1    *   *
c      *         *     *
c      *         ?   ?   *
c      *         *         *
c      *         *           *
c      *************************
c
c      Iref = 4; 'small' constrained nodes
c      *
c      * *
c      *   *
c      *     *
c      *   ?   *
c      *         *
c      *           *
c      * * * ? * * * *
c      *               *
c      *        1        *
c      *                   *
c      *                     *
c      *************************
c
c--------------------------------------------------------------------
c
      use parameters
      use constraints
#include "syscom.blk"
c
      dimension shapsma(MAXquadH,MAXquadH),
     .          shapbig(MAXquadH,MAXquadH),
     .          void(2,MAXquadH),ip(MAXquadH),
     .          r(MAXquadH,MAXquadH),xibig(2),xisma(2)
c
c  ...order and edge orientations
      dimension norder(5),norient(4)
c
      iprint=0
c
      if (MAXP.ne.MAXP) then
        write(*,7005) MAXP,MAXP
 7005   format('setcnstrt: INCONSISTENT MAXP,MAXP = ',2i3)
        stop 1
      endif
      nord = MAXP
c
c  ...set uniform order MAXP:
      norder(1:4) = nord
      norder(5)   = nord*10 + nord
c
c  ...set orientations
      norient = 0
c
c*********************************************************************
c
      icl = 0
c
c  ...loop through collocation points:
      do j = 1,nord+1
        do i = 1,nord+1
          xisma(1) = (i-1)*1.0d0/nord
          xisma(2) = (j-1)*1.0d0/nord
          call map_quad(Iref,xisma(1:2), xibig(1:2))
c
          icl = icl + 1
c
c  .......find small shape functions at this point:
c          call shapeHq(xisma,norder(1:5),norient(1:4),
c     .                 nrdof,shapsma(1:MAXquadH,icl),void)
          call shape2DH(QUAD,xisma,norder,norient,
     .                  nrdof,shapsma(:,icl),void)
c
c  .......find big shape functions at this point:
c          call shapeHt(xibig,norder(1:4),norient(1:3),
c     .                 nrdof,shapbig(1:MAXtriaH,icl),void)
          call shape2DH(TRIA,xibig,norder,norient,
     .                  nrdof,shapbig(:,icl),void)
        enddo
      enddo
c
c  ...transpose shapsma:
      do i=1,(nord+1)**2
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
      n = (nord+1)**2
      na = MAXquadH + 1
      if (iprint.eq.1) then
        write(*,*) 'setcnstr_trian: shapsma COLUMNWISE ....'
        do i=1,n
          write(*,7004) i
 7004     format('setcnstr_trian: i = ',i2)
          write(*,7003) shapsma(1:n,i)
 7003     format(10e12.5)
        enddo
        call pause
      endif
ccc      call decomp(n,na,shapsma,ip,iflag)
cccc
ccc      if (iflag.ne.0) then
ccc        write(*,*) 'setcnstr_trian: iflag = ',iflag
ccc        stop 1
ccc      endif
cccc
ccc      do i = 1,n
ccc        call gauss2(n,na,shapsma,ip,shapbig(1:n,i), r(1:n,i))
ccc        if (iprint.eq.1) then
ccc          write(*,*) 'setcnstr_trian: r FOR i = ',i
ccc          write(*,7003) r(1:n,i)
ccc        endif
ccc      enddo

c  ...decompose the matrix
      call dgetrf(n,n,shapsma,MAXquadH,ip,info)
      if (info.ne.0) then
        write(*,*)'setcnstr_trian_aniso: H1 DGETRF RETURNED INFO =',info
        call logic_error(FAILURE,
     .  __FILE__,__LINE__)
      endif
c
c  ...right-hand side resolution
      call dlaswp(n,shapbig,n,1,MAXquadH,ip,1)
      call dtrsm('L','L','N','U',n,n,1.d0,shapsma,MAXquadH,shapbig,
     .            MAXquadH)
      call dtrsm('L','U','N','N',n,n,1.d0,shapsma,MAXquadH,shapbig,
     .            MAXquadH)
      r = shapbig


c
c  ...degrees of freedom of small node 1:
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      nrs = 4*nord
c
      do i = 1,(nord-1)**2
c
c  .....big node 2,3,4
        do ie=2,4
          nrb = 3 + (nord-1)*(ie-2)
          do j = 1,nord-1
            RRQH(Iref,ie,j,i) = r(nrs+i,nrb+j)
          enddo
        enddo
c
c  .....big node 1 (bubbles):
        nrb = 3*nord
        do j = 1,(nord-1)*(nord-2)/2
          RRQH(Iref,1,j,i) = r(nrs+i,nrb+j)
        enddo
      enddo
c
c  ...skip the cleaning
      return
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
c  .........loop through constrained dof
            do i2=1,nord-1
c
c  ...........beginning of dof to be zero out
              if (i2.le.n-1) then
                i1beg = n
              else
                i1beg = 1
              endif
              do i1=i1beg,nord-1
c
c  .............constrained dof is
                j = (i2-1)*(nord-1)+i1
ccc              if (abs(get_rrqh(Iref, n, ipp,jp,j)).gt.1.d-12) then
                if (abs(RRQH(Iref,ipp,jp,j)).gt.1.d-12) then
                  write(*,7001) n,ipp,jp,i1,i2,RRQH(Iref,ipp,jp,j)
 7001             format('setcnstr_trian_aniso: n,ipp,jp,i1,i2,RRQH = ',
     .                    5i3,2x,e12.5)
                  stop 1
                endif
                RRQH(Iref,ipp,jp,j) = 0.d0
              enddo
            enddo
          enddo
        enddo
      enddo
c
ccc      call test_trian_aniso(Iref)
c
      end




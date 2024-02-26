!----------------------------------------------------------------------
!
!   routine name       - setcnstr_trian_aniso_h1
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine evaluates constraint coefficients for
!                        the triangular master element
!
!----------------------------------------------------------------------
!
!     'big' parent nodes
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
!      Iref = 2; 'small' constrained nodes
!      *
!      * *
!      *   *
!      *     *
!      *       *
!      *         *
!      *           *
!      * *      1    *
!      *   *           *
!      *     ?           *
!      *   ?   *           *
!      *         *           *
!      *************************
!
!      Iref = 3; 'small' constrained nodes
!      *
!      * *
!      *   *
!      *     *
!      *       *
!      *         *
!      *         * *
!      *    1    *   *
!      *         *     *
!      *         ?   ?   *
!      *         *         *
!      *         *           *
!      *************************
!
!      Iref = 4; 'small' constrained nodes
!      *
!      * *
!      *   *
!      *     *
!      *   ?   *
!      *         *
!      *           *
!      * * * ? * * * *
!      *               *
!      *        1        *
!      *                   *
!      *                     *
!      *************************
!
!--------------------------------------------------------------------
!
      subroutine setcnstr_trian_aniso_h1(Iref)
!
      use parameters
      use constraints
!
      implicit none
!
      integer, intent(in) :: Iref
!
      real(8) :: shapsma(MAXquadH,MAXquadH), &
                 shapbig(MAXquadH,MAXquadH), &
                 void(2,MAXquadH),           &
                 r(MAXquadH,MAXquadH),xibig(2),xisma(2)
      integer :: ip(MAXquadH)
!
!  ...order and edge orientations
      integer :: norder(5),norient(4)
!
      real(8) :: aux
      integer :: i,ie,i1,i2,i1beg,ipp,jp,jpbeg,jpend,j,icl
      integer :: n,na,nord,nrb,nrdof,nrs,info
!
#if HP3D_DEBUG
      integer :: iprint
      iprint = 0
#endif
!
      nord = MAXP
!
!  ...set uniform order MAXP:
      norder(1:4) = nord
      norder(5) = nord*10 + nord
!
!  ...set orientations
      norient = 0
!
!*********************************************************************
!
      icl = 0
!
!  ...loop through collocation points:
      do j = 1,nord+1
        do i = 1,nord+1
          xisma(1) = (i-1)*1.0d0/nord
          xisma(2) = (j-1)*1.0d0/nord
          call map_quad(Iref,xisma(1:2), xibig(1:2))
!
          icl = icl + 1
!
!  .......find small shape functions at this point:
!          call shapeHq(xisma,norder(1:5),norient(1:4),  &
!                        nrdof,shapsma(1:MAXquadH,icl),void)
          call shape2DH(QUAD,xisma,norder,norient, &
                         nrdof,shapsma(:,icl),void)
!
!  .......find big shape functions at this point:
!          call shapeHt(xibig,norder(1:4),norient(1:3),  &
!                        nrdof,shapbig(1:MAXtriaH,icl),void)
          call shape2DH(TRIA,xibig,norder,norient, &
                         nrdof,shapbig(:,icl),void)
        enddo
      enddo
!
!  ...transpose shapsma:
      do i=1,(nord+1)**2
        do j=1,i
          aux = shapsma(i,j)
          shapsma(i,j) = shapsma(j,i)
          shapsma(j,i) = aux
!
          aux = shapbig(i,j)
          shapbig(i,j) = shapbig(j,i)
          shapbig(j,i) = aux
        enddo
      enddo
!
      n = (nord+1)**2
      na = MAXquadH + 1
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) 'setcnstr_trian_aniso_h1: shapsma COLUMNWISE ....'
        do i=1,n
          write(*,7004) i
 7004     format('setcnstr_trian_aniso_h1: i = ',i2)
          write(*,7003) shapsma(1:n,i)
 7003     format(10e12.5)
        enddo
        call pause
      endif
#endif
!!!      call decomp(n,na,shapsma,ip,iflag)
!!!!
!!!      if (iflag.ne.0) then
!!!        write(*,*) 'setcnstr_trian: iflag = ',iflag
!!!        stop 1
!!!      endif
!!!!
!!!      do i = 1,n
!!!        call gauss2(n,na,shapsma,ip,shapbig(1:n,i), r(1:n,i))
!!!#if HP3D_DEBUG
!!!        if (iprint.eq.1) then
!!!          write(*,*) 'setcnstr_trian: r FOR i = ',i
!!!          write(*,7003) r(1:n,i)
!!!        endif
!!!#endif
!!!      enddo

!  ...decompose the matrix
      call dgetrf(n,n,shapsma,MAXquadH,ip,info)
      if (info.ne.0) then
        write(*,*)'setcnstr_trian_aniso_h1: H1 DGETRF RETURNED INFO =',info
        call logic_error(FAILURE,__FILE__,__LINE__)
      endif
!
!  ...right-hand side resolution
      call dlaswp(n,shapbig,n,1,MAXquadH,ip,1)
      call dtrsm('L','L','N','U',n,n,1.d0,shapsma,MAXquadH,shapbig,MAXquadH)
      call dtrsm('L','U','N','N',n,n,1.d0,shapsma,MAXquadH,shapbig,MAXquadH)
      r = shapbig


!
!  ...degrees of freedom of small node 1:
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      nrs = 4*nord
!
      do i = 1,(nord-1)**2
!
!  .....big node 2,3,4
        do ie=2,4
          nrb = 3 + (nord-1)*(ie-2)
          do j = 1,nord-1
            RRQH(Iref,ie,j,i) = r(nrs+i,nrb+j)
          enddo
        enddo
!
!  .....big node 1 (bubbles):
        nrb = 3*nord
        do j = 1,(nord-1)*(nord-2)/2
          RRQH(Iref,1,j,i) = r(nrs+i,nrb+j)
        enddo
      enddo
!
!  ...skip the cleaning
!!!      return
!
!-----------------------------------------------------------------------
!
!  ...clean up machine zeros...
!
!  ...loop through possible orders of approximations
      do n=2,nord
!
!  .....loop through parent nodes
        do ipp=1,4
!
!  .......beginning and ending dof
          select case(ipp)
          case(1)
            jpbeg = (n-2)*(n-3)/2 + 1; jpend = (n-1)*(n-2)/2
          case(2,3,4)
            jpbeg = n-1; jpend = n-1
          end select
!
!  .......loop through the parent dof
          do jp=jpbeg,jpend
!
!  .........loop through constrained dof
            do i2=1,nord-1
!
!  ...........beginning of dof to be zero out
              if (i2.le.n-1) then
                i1beg = n
              else
                i1beg = 1
              endif
              do i1=i1beg,nord-1
!
!  .............constrained dof is
                j = (i2-1)*(nord-1)+i1
!!!              if (abs(get_rrqh(Iref, n, ipp,jp,j)).gt.1.d-12) then
                if (abs(RRQH(Iref,ipp,jp,j)).gt.1.d-12) then
                  write(*,7001) n,ipp,jp,i1,i2,RRQH(Iref,ipp,jp,j)
 7001             format('setcnstr_trian_aniso_h1: ', &
                         /'n,ipp,jp,i1,i2,RRQH = ',5i3,2x,e12.5)
                  stop 1
                endif
                RRQH(Iref,ipp,jp,j) = 0.d0
              enddo
            enddo
          enddo
        enddo
      enddo
!
!!!      call test_trian_aniso_h1(Iref)
!
      end subroutine setcnstr_trian_aniso_h1

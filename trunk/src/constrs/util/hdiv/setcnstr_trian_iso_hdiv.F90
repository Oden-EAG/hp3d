!----------------------------------------------------------------------
!
!   routine name       - setcnstr_trian_iso_hdiv
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine evaluates constraint coefficients for
!                        the triangular master element and H(div) (L2)
!                        constrained approximation
!
!   arguments          - none
!
!----------------------------------------------------------------------
!
!    'big' parent node (just one)
!
!      *
!      * *
!      *   *
!      *     *
!      *       *
!      *         *
!      *           *
!      *             *
!      *       1       *
!      *                 *
!      *                   *
!      *                     *
!      *************************
!
!
!      'small' constrained nodes
!      *
!      * *
!      *   *
!      *     *
!      *   3   *
!      *         *
!      *************
!      * *         * *
!      *   *   4   *   *
!      *     *     *     *
!      *   1   *   *   2   *
!      *         * *         *
!      *************************
!
!--------------------------------------------------------------------
!
   subroutine setcnstr_trian_iso_hdiv
!
      use parameters
      use constraints
      use element_data
!
      implicit none
!
!  ...big(parent) and small element shape functions, collocation points
      real(8) :: shapsma(MAXtriaQ,MAXtriaQ), &
                 shapbig(MAXtriaQ,MAXtriaQ), &
                 xibig(2),xisma(2)
      integer :: ip(MAXtriaQ)
!
      integer :: norder(4),nsize(2)
!
      real(8) :: aux
      integer :: i,j,icl,info
      integer :: n,nord,nel,nrdofQ
!
      real(8), parameter :: eps = 1.0d-13
!
#if HP3D_DEBUG
      real(8) :: diff,val
      integer :: i1,i2,ians,mp,np
#endif
!
!  ...set up order
      nord   = MAXP
      norder = (/1,1,1,nord/)
      nsize  = (/MAXP,MAXtriaQ/)
!
!  ...loop through small triangles
      do nel=1,4
!
!  .....loop through collocation points
        icl=0
        do i=1,nord
          select case(nel)
          case(1)
            xibig(1) = (i-1)*0.5d0/(nord-1)
            xisma(1) = (i-1)*1.0d0/(nord-1)
          case(2)
            xibig(1) = (i-1)*0.5d0/(nord-1) + 0.5d0
            xisma(1) = (i-1)*1.0d0/(nord-1)
          case(3)
            xibig(1) = (i-1)*0.5d0/(nord-1)
            xisma(1) = (i-1)*1.0d0/(nord-1)
          case(4)
            xibig(1) = 0.5d0 - (i-1)*0.5d0/(nord-1)
            xisma(1) =         (i-1)*1.0d0/(nord-1)
          end select
!
          do j=1,nord -(i-1)
            icl = icl+1
            select case(nel)
            case(1)
              xibig(2) = 0.5d0/(nord-1)*(j-1)
              xisma(2) = 1.0d0/(nord-1)*(j-1)
            case(2)
              xibig(2) = 0.5d0/(nord-1)*(j-1)
              xisma(2) = 1.0d0/(nord-1)*(j-1)
            case(3)
              xibig(2) = 0.5d0/(nord-1)*(j-1) + 0.5d0
              xisma(2) = 1.0d0/(nord-1)*(j-1)
            case(4)
              xibig(2) = 0.5d0 - 0.5d0/(nord-1)*(j-1)
              xisma(2) =         1.0d0/(nord-1)*(j-1)
            end select
!
!  .........compute small shape functions at this point:
!            call shapeQt(xisma,(/1,1,1,nord/), &
!                              nrdofQ,shapsma(1:MAXtriaQ,icl))
            call shape2DQTri(xisma,norder,nsize, &
                             nrdofQ,shapsma(:,icl))
!
!  .........Piola transform (divide by jacobian)
            shapsma(1:MAXtriaQ,icl) = shapsma(1:MAXtriaQ,icl)/.25d0
!
!  .........compute big shape functions at this point:
!            call shapeQt(xibig,(/1,1,1,nord/), &
!                              nrdofQ,shapbig(1:MAXtriaQ,icl))
            call shape2DQTri(xibig,norder,nsize, &
                             nrdofQ,shapbig(:,icl))
          enddo
        enddo
        if (icl.ne.nrdofQ) then
          write(*,*) 'setcnstr_trian_iso_hdiv: icl,nrdofQ = ',icl,nrdofQ
          stop 1
        endif
!
!  .....transpose both matrices (due to the way, we interface with linear solver)
        do i=1,nrdofQ
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
!  .....solve for the constrained approximation coefficients
        n=nrdofQ
        call dgetrf(n,n,shapsma,MAXtriaQ,ip,info)
        if (info.ne.0) then
          write(*,*)'setcnstr_trian_iso_hdiv: DGETRF INFO =',info
          call logic_error(FAILURE, __FILE__,__LINE__)
        endif
!
!  .....right-hand side resolution
        call dlaswp(n,shapbig,n,1,MAXtriaQ,ip,1)
        call dtrsm('L','L','N','U',n,n,1.d0,shapsma,MAXtriaQ,shapbig,MAXtriaQ)
        call dtrsm('L','U','N','N',n,n,1.d0,shapsma,MAXtriaQ,shapbig,MAXtriaQ)
!
!  .....save coefficients cleaning machine zeros
        do i=1,nrdofQ
          do j=1,nrdofQ
            if (abs(shapbig(i,j)).gt.1.d-12) then
              RRTQ(j,nel,i) = shapbig(i,j)
            else
              RRTQ(j,nel,i) = 0.d0
            endif
          enddo
        enddo
!
!  ...end of loop through small triangles
      enddo
!
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
      write(*,*) 'setcnstr_trian_iso_hdiv: SET xibig '
      read(*,*) xibig(1:2)
      write(*,*) 'xibig = ',xibig
!
!  ...shape functions of big element
!     call shapeQt(xibig,(/1,1,1,nord/), nrdofQ,shapbig(1:MAXtriaQ,1))
      call shape2DQTri(xibig,norder,nsize, nrdofQ,shapbig(:,1))
!
!  ...determine which small element is this:
!  ...coordinates in a small element:
      if (xibig(1)+xibig(2).le.0.5d0) then
!
!  .....element 1
        nel = 1
        xisma(1) = 2*xibig(1)
        xisma(2) = 2*xibig(2)
!
!
      elseif (xibig(1).ge.0.5d0) then
!
!  .....element 2
        nel = 2
        xisma(1) = (xibig(1)-0.5d0)*2.d0
        xisma(2) =  xibig(2)*2.d0
!
      elseif (xibig(2).gt.0.5d0) then
!
!  .....element 3
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
      endif

      write(*,*) 'setcnstr_trian_iso_hdiv: nel = ', nel
!
!  ...find small shape functions at this point:
!     call shapeQt(xisma,(/1,1,1,nord/), nrdofQ,shapsma(1:MAXtriaQ,1))
      call shape2DQTri(xisma,norder,nsize, nrdofQ,shapsma(:,1))
!
!*********************************************************************
!
!  ...verify if big shape functions are right combinations of small
!     ones:
!
!  ...loop through big shape functions (central node)
!
!  ...loop through polynomial orders
      i=0
      do np=0,nord-1
      do i1=0,np
        i=i+1
        write(*,7010) i
 7010   format('BIG ELEMENT SHAPE FUNTION ',i3)
        write(*,7020) nel
 7020   format('nel = ',i2,' COEFFICIENTS = ')
        write(*,7030) RRTQ(i,nel,1:nrdofQ)
 7030   format(10e12.5)
!
!  .....initiate value of the linear combination:
        val = 0.d0
        j=0
        do mp=0,nord-1
        do i2=0,mp
          j=j+1
          if ((mp.gt.np).and.(dabs(RRTQ(i,nel,j)) > eps)) then
            write(*,7040) np,i1,mp,i2,RRTQ(i,nel,j)
 7040       format('np,i1,mp,i2,RRTQ(i,nel,j) = ',4i4,2x,e12.5)
            stop 1
          endif
          val  = val + RRTQ(i,nel,j)*shapsma(j,1)/.25d0
        enddo
        enddo
        diff = abs(val-shapbig(i,1))
        if (diff.gt.1.d-12) then
          write(*,*) 'i,shapbig(i,1),val,diff = ', &
                      i,shapbig(i,1),val,diff
          call pause
        endif
      enddo
      enddo
!
      write(*,*) 'setcnstr_trian_iso_hdiv: CONTINUE ?(1/0)'
      read(*,*) ians
!
      if (ians.eq.1) go to 777
!
#endif
!
   end subroutine setcnstr_trian_iso_hdiv

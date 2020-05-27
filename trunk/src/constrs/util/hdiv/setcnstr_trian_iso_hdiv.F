c----------------------------------------------------------------------
c
c   routine name       - setcnstr_trian_iso_div
c
c----------------------------------------------------------------------
c
c   latest revision    - May 2020
c
c   purpose            - routine evaluates constraint coefficients for
c                        the triangular master element and H(div) (L2)
c                        constrained approximation
c
c   arguments          - none
c
c----------------------------------------------------------------------
c
      subroutine setcnstr_trian_iso_hdiv
c
c    'big' parent node (just one)
c
c      *
c      * *
c      *   *
c      *     *
c      *       *
c      *         *
c      *           *
c      *             *
c      *       1       *
c      *                 *
c      *                   *
c      *                     *
c      *************************
c
c
c      'small' constrained nodes
c      *
c      * *
c      *   *
c      *     *
c      *   3   *
c      *         *
c      *************
c      * *         * *
c      *   *   4   *   *
c      *     *     *     *
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
c
c  ...big(parent) and small element shape functions, collocation points
      dimension shapsma(MAXtriaQ,MAXtriaQ),
     .          shapbig(MAXtriaQ,MAXtriaQ),
     .          ip(MAXtriaQ),xibig(2),xisma(2)
c
      dimension norder(4),nsize(2)
c
      iprint=0
c
c  ...set up order
      nord   = MAXP
      norder = (/1,1,1,nord/)
      nsize  = (/MAXP,MAXtriaQ/)
c
c  ...loop through small triangles
      do nel=1,4
c
c  .....loop through collocation points
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
c
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
c
c  .........compute small shape functions at this point:
c            call shapeQt(xisma,(/1,1,1,nord/),
c     .                        nrdofQ,shapsma(1:MAXtriaQ,icl))
            call shape2DQTri(xisma,norder,nsize,
     .                       nrdofQ,shapsma(:,icl))
c
c  .........Piola transform (divide by jacobian)
            shapsma(1:MAXtriaQ,icl) = shapsma(1:MAXtriaQ,icl)/.25d0
c
c  .........compute big shape functions at this point:
c            call shapeQt(xibig,(/1,1,1,nord/),
c     .                        nrdofQ,shapbig(1:MAXtriaQ,icl))
            call shape2DQTri(xibig,norder,nsize,
     .                       nrdofQ,shapbig(:,icl))
          enddo
        enddo
        if (icl.ne.nrdofQ) then
          write(*,*) ' setcnstr_trian_iso_div: icl,nrdofQ = ',icl,nrdofQ
          stop 1
        endif
c
c  .....transpose both matrices (due to the way, we interface with linear solver)
        do i=1,nrdofQ
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
c  .....solve for the constrained approximation coefficients
        n=nrdofQ
        call dgetrf(n,n,shapsma,MAXtriaQ,ip,info)
        if (info.ne.0) then
          write(*,*)'setcnstr_trian_iso_hdiv: DGETRF INFO =',info
          call logic_error(FAILURE,
     .    __FILE__,__LINE__)
        endif
c
c  .....right-hand side resolution
        call dlaswp(n,shapbig,n,1,MAXtriaQ,ip,1)
        call dtrsm('L','L','N','U',n,n,1.d0,shapsma,MAXtriaQ,shapbig,
     .              MAXtriaQ)
        call dtrsm('L','U','N','N',n,n,1.d0,shapsma,MAXtriaQ,shapbig,
     .              MAXtriaQ)
c
c  .....save coefficients cleaning machine zeros
        do i=1,nrdofQ
          do j=1,nrdofQ
            if (abs(shapbig(i,j)).gt.1.d-12) then
              RRTQ(j,nel,i) = shapbig(i,j)
            else
              RRTQ(j,nel,i) = 0.d0
            endif
          enddo
        enddo
c
c  ...end of loop through small triangles
      enddo
c
c
      return
c
c*********************************************************************
c
c  ...begin testing
  777 continue
c
      write(*,*) 'setcnstr_trian_iso_hdiv: SET xibig '
      read(*,*) xibig(1:2)
      write(*,*) 'xibig = ',xibig
c
c  ...shape functions of big element
c     call shapeQt(xibig,(/1,1,1,nord/), nrdofQ,shapbig(1:MAXtriaQ,1))
      call shape2DQTri(xibig,norder,nsize, nrdofQ,shapbig(:,1))
c
c  ...determine which small element is this:
c  ...coordinates in a small element:
      if (xibig(1)+xibig(2).le.0.5d0) then
c
c  .....element 1
        nel = 1
        xisma(1) = 2*xibig(1)
        xisma(2) = 2*xibig(2)
c
c
      elseif (xibig(1).ge.0.5d0) then
c
c  .....element 2
        nel = 2
        xisma(1) = (xibig(1)-0.5d0)*2.d0
        xisma(2) =  xibig(2)*2.d0
c
      elseif (xibig(2).gt.0.5d0) then
c
c  .....element 3
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

      write(*,*) 'setcnstr_trian: nel = ', nel
c
c  ...find small shape functions at this point:
c     call shapeQt(xisma,(/1,1,1,nord/), nrdofQ,shapsma(1:MAXtriaQ,1))
      call shape2DQTri(xisma,norder,nsize, nrdofQ,shapsma(:,1))
c
c*********************************************************************
c
c  ...verify if big shape functions are right combinations of small
c     ones:
c
c  ...loop through big shape functions (central node)
c
c  ...loop through polynomial orders
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
c
c  .....initiate value of the linear combination:
        val = 0.d0
        j=0
        do mp=0,nord-1
        do i2=0,mp
          j=j+1
          if ((mp.gt.np).and.(RRTQ(i,nel,j).ne.0.d0)) then
            write(7040) np,i1,mp,i2,RRTQ(i,nel,j)
 7040       format('np,i1,mp,i2,RRTQ(i,nel,j) = ',4i4,2x,e12.5)
            stop 1
          endif
          val  = val + RRTQ(i,nel,j)*shapsma(j,1)/.25d0
        enddo
        enddo
        diff = abs(val-shapbig(i,1))
        if (diff.gt.1.d-12) then
          write(*,*) 'i,shapbig(i,1),val,diff = ',
     .                i,shapbig(i,1),val,diff
          call pause
        endif
      enddo
      enddo
c
       write(*,*) 'setcnstr_trian: CONTINUE ?(1/0)'
       read(*,*) ians
c
       if (ians.eq.1) go to 777
c
       end subroutine setcnstr_trian_iso_hdiv

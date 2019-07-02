c----------------------------------------------------------------------
c
c   routine name       - setcnstr_edge
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 09
c
c   purpose            - routine sets coefficients for the constrained
c                        approximation on the master edge
c
c   arguments          - none
c
c   required  routines -
c
c----------------------------------------------------------------------
c
c     parent nodes:
c
c     |----------------->------------------|
c    big 2           big 1                big 3
c
c
c     constrained nodes:
c
c     |--------->-------|--------->--------|
c             sm.1     sm.3     sm.2
c
c
c---------------------------------------------------------------------
c
      subroutine setcnstr_edge
c
      use parameters
      use constraints
#include "syscom.blk"
c
c  ...values of small elements shape functions, big element
c     shape functions, derivatives (not used)
      dimension shapsma(MAXP+2,MAXP+1),shapbig(MAXP+1,MAXP+1),
     .          void(MAXP+1),
c
c  ...pivoting array, work array
     .          ip(MAXP+1),r(MAXP+1,MAXP+1)
c
      iprint=0
c
c  ...clear:
      RRRH = 0.d0
c
      nord = MAXP
c
c---------------------------------------------------------------------
c
c  ...small element 1 - small nodes 1 and 3:
c
c---------------------------------------------------------------------
c
      icl = 0
c
c  ...loop through collocation points:
      do i=1,nord+1
        icl = icl+1
c
c  .....local coordinates for big and small elements
        xibig = (i-1)*0.5d0/nord
        xisma = (i-1)*1.0d0/nord
c
c  .....find small shape functions at this point:
        call shape1H(xisma,nord,0, shapsma(1:MAXP+1,icl),void)
c
c  .....find big shape functions at this point:
        call shape1H(xibig,nord,0, shapbig(1:MAXP+1,icl),void)
      enddo
c
c  ...transpose shapsma and shapbig
      do i=1,nord+1
        do j=1,i
c
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
c  ...solve the collocation system to determine the constrained
c     approximation coefficients
c
c  ...number of unknowns
      n = nord+1
c
c  ...invert the collocation matrix
      if (iprint.eq.1) then
        write(*,7010)
 7010   format('setcnstr_edge: COLLOCATION MATRIX = ')
        do i=1,n
          write(*,7011) i,shapsma(i,1:n)
 7011     format('i = ',i3,10e12.5)
        enddo
        call pause
      endif
      call decomp(n,MAXP+2,shapsma,ip,iflag)
      if (iflag.ne.0) then
        write(*,7002) iflag
 7002   format('setcnstr_edge: TROUBLE WITH INVERTING THE COLLOCATION',
     .         ' MATRIX, iflag = ',i3)
        call pause
      endif
c
c  ...right-hand side resolution
      do i=1,n
        call gauss2(n,MAXP+2,shapsma,ip,shapbig(1:n,i),r(1:n,i))
      enddo
c
c  ...degrees of freedom of small node 3:
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c  ...big node 2:
      RRRH(2,1,3,1) = r(2,1)
c
c  ...big node 3:
      RRRH(3,1,3,1) = r(2,2)
c
c  ...big node 1:
      nrb = 2
      do j=1,nord-1
        RRRH(1,j,3,1) = r(2,nrb+j)
      enddo
c
c  ...degrees of freedom of small node 1:
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      nrs = 2
      do i=1,nord-1
c
c  .....big node 1:
        nrb = 2
        do j=1,nord-1
          RRRH(1,j,1,i) = r(nrs+i,nrb+j)
        enddo
      enddo
c
c---------------------------------------------------------------------
c
c  ...small element 2 --- small nodes 2 and 3:
c
c---------------------------------------------------------------------
c
      icl = 0
c
c  ...loop through collocation points:
      do i=1,nord+1
        icl = icl+1
c
c  .....big and small element coordinates of the collocation point
        xibig = (i-1)*0.5d0/nord + 0.5d0
        xisma = (i-1)*1.0d0/nord
c
c  .....find small shape functions at this point:
        call shape1H(xisma,nord,0, shapsma(1:MAXP+1,icl),void)
c
c  .....find big shape functions at this point:
        call shape1H(xibig,nord,0, shapbig(1:MAXP+1,icl),void)
      enddo
c
c  ...transpose shapsma and shapbig
      do i=1,nord+1
        do j=1,i
c
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
c  ...solve the collocation system
      n = nord+1
      call decomp(n,MAXP+2,shapsma,ip,Iflag)
      if (iflag.ne.0) then
        write(*,7002) iflag
        call pause
      endif
c
      do i=1,n
        call gauss2(n,MAXP+2,shapsma,ip,shapbig(1:n,i),r(1:n,i))
      enddo
c
c  ...degrees of freedom of small node 2:
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      nrs = 2
      do i=1,nord-1
c
c  .....big node 1:
        nrb = 2
        do j=1,nord-1
          RRRH(1,j,2,i) = r(nrs+i,nrb+j)
        enddo
      enddo
c
c-----------------------------------------------------------------------
c
c  ...clean up machine zeros...
c
c  ...loop through possible orders of approximations
      do n=2,nord
        jpbeg = n-1; jpend = n-1
c
c  .....loop through the parent dof
        do jp=jpbeg,jpend
c
c  .......loop through constrained middle nodes
          do i=1,2
c
c  .........beginning of dof to be zero out
            jbeg = n
c
c  .........loop through the constrained node dof
            do j=jbeg,nord-1
              if (abs(RRRH(1,jp,i,j)).gt.1.d-12) then
                write(*,7001) n,jp,i,j,RRRH(ipp,jp,i,j)
 7001           format('setcnstr_edge: n,jp,i,j,RRTH = ',
     .                  4i3,2x,e12.5)
                stop 1
              endif
              RRRH(1,jp,i,j) = 0.d0
            enddo
          enddo
        enddo
      enddo
c
      return
c
c---------------------------------------------------------------------
c
c  ...this is a test which should be run when installing the code
c     on a new machine
c
   10 write(*,*) 'setcnstr_edge: SET BIG ELEMENT COORDINATE'
      read(*,*) xibig
c
c  ...shape functions of the big element
      call shape1H(xibig,nord,0, shapbig(1:MAXP+1,1),void)
c
c  ...determine which small element is this, and find the small
c     element coordinate of the point
      if (xibig.le.0.5d0) then
c
c  .....element 1 - verifying small nodes 1,3
        nel = 1
        xisma = 2*xibig
      else
c
c  .....element 2 - verifying small nodes 2,3
        nel = 2
        xisma = (xibig-0.5)*2
      endif
c
c  ...find small shape functions at this point:
      call shape1H(xisma,nord,0, shapsma(1:MAXP+1,1),void)
c
c  ...verify if the values of the big shape functions match
c     combinations of the small shape functions.......................
c
c  ...loop through the big shape functions - central node:
      nrb = 2 ! skip vertex shape functions
      do i=1,nord-1
c
c  .....initiate value of the linear combination:
        val = 0
        select case(nel)
        case(1)
c
c  .......small node 3:
          nrs = 1
          do j=1,1
            val = val + RRRH(1,i, 3,j)*shapsma(nrs+j,1)
          enddo
c
c  .......small node 1:
          nrs = 2
          do j=1,nord-1
            val = val + RRRH(1,i, 1,j)*shapsma(nrs+j,1)
          enddo
c
        case(2)
c
c  .......small node 3:
          nrs = 0
          do j=1,1
            val = val + RRRH(1,i, 3,j)*shapsma(nrs+j,1)
          enddo
c
c  .......small node 2:
          nrs = 2
          do j=1,nord-1
            val = val + RRRH(1,i, 2,j)*shapsma(nrs+j,1)
          enddo
        end select
c
        write(*,7004) i,shapbig(nrb+i,1),abs(val-shapbig(nrb+i,1))
 7004   format('setcnstr_edge: MIDDLE NODE SHAPE FUNCTION i = ',i2,
     .         ' VALUE = ',e12.5, ' DIFFERENCE = ',e12.5)
c
c  ...end of loop through middle node big shape functions
      enddo
c
c  ...left big vertex node
      select case(nel)
      case(1)
        val = RRRH(2,1, 3,1)*shapsma(2,1)
     .      + 1.d0         *shapsma(1,1)
      case(2)
        val = RRRH(2,1, 3,1)*shapsma(1,1)
      end select
      write(*,7005) shapbig(1,1),abs(val-shapbig(1,1))
 7005 format('setcnstr_edge: LEFT VERTEX NODE SHAPE FUNCTION '
     .       ' VALUE = ',e12.5, ' DIFFERENCE = ',e12.5)
c
c  ...right big vertex node
      select case(nel)
      case(1)
        val = RRRH(3,1, 3,1)*shapsma(2,1)
      case(2)
        val = RRRH(3,1, 3,1)*shapsma(1,1)
     .      + 1.0d0        *shapsma(2,1)
      end select
      write(*,7006) shapbig(2,1),abs(val-shapbig(2,1))
 7006 format('setcnstr_edge: RIGHT VERTEX NODE SHAPE FUNCTION '
     .       ' VALUE = ',e12.5, ' DIFFERENCE = ',e12.5)
c
      write(*,*) 'setcnstr: CONTINUE ?(1/0)'
      read(*,*) ians
      if (ians.eq.1) go to 10
c
c
      end subroutine











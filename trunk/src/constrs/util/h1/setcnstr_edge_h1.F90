!----------------------------------------------------------------------
!
!   routine name       - setcnstr_edge_h1
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine sets coefficients for the constrained
!                        approximation on the master edge
!
!   arguments          - none
!
!----------------------------------------------------------------------
!
!     parent nodes:
!
!     |----------------->------------------|
!    big 2           big 1                big 3
!
!
!     constrained nodes:
!
!     |--------->-------|--------->--------|
!             sm.1     sm.3     sm.2
!
!
!---------------------------------------------------------------------
!
      subroutine setcnstr_edge_h1
!
      use parameters
      use constraints
!
      implicit none
!
!  ...values of small elements shape functions, big element
!     shape functions, derivatives (not used)
      real(8) :: shapsma(MAXP+1,MAXP+1),shapbig(MAXP+1,MAXP+1),   &
                 void(MAXP+1)
!
!  ...pivoting array, work array
      integer :: ip(MAXP+1)
      real(8) :: r(MAXP+1,MAXP+1)
!
      real(8) :: aux,xibig,xisma
      integer :: i,j,icl,info,jbeg,jp,jpbeg,jpend
      integer :: n,nord,ndofh,nrb,nrs
!
#if HP3D_DEBUG
      real(8) :: val
      integer :: ians,nel
      integer :: iprint
      iprint = 0
#endif
!
!  ...clear:
      RRRH = 0.d0
!
      nord = MAXP
!
!---------------------------------------------------------------------
!
!  ...small element 1 - small nodes 1 and 3:
!
!---------------------------------------------------------------------
!
      icl = 0
!
!  ...loop through collocation points:
      do i=1,nord+1
        icl = icl+1
!
!  .....local coordinates for big and small elements
        xibig = (i-1)*0.5d0/nord
        xisma = (i-1)*1.0d0/nord
!
!  .....find small shape functions at this point:
        call shape1DH(xisma,nord, ndofh,shapsma(1:MAXP+1,icl),void)
!
!  .....find big shape functions at this point:
        call shape1DH(xibig,nord, ndofh,shapbig(1:MAXP+1,icl),void)
      enddo
!
!  ...transpose shapsma and shapbig
      do i=1,nord+1
        do j=1,i
!
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
!  ...solve the collocation system to determine the constrained
!     approximation coefficients
!
!  ...number of unknowns
      n = nord+1
!
!  ...invert the collocation matrix
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7010)
 7010   format('setcnstr_edge_h1: COLLOCATION MATRIX = ')
        do i=1,n
          write(*,7011) i,shapsma(i,1:n)
 7011     format('i = ',i3,10e12.5)
        enddo
        call pause
      endif
#endif
      call dgesv(n,MAXP+1,shapsma,MAXP+1,ip,shapbig,MAXP+1,info)
      if (info.ne.0) then
        write(*,*)'setcnstr_edge_h1: H1 DGETRS RETURNED INFO =',info
        call logic_error(FAILURE,__FILE__,__LINE__)
      endif
!
      r = shapbig
!
!  ...degrees of freedom of small node 3:
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  ...big node 2:
      RRRH(2,1,3,1) = r(2,1)
!
!  ...big node 3:
      RRRH(3,1,3,1) = r(2,2)
!
!  ...big node 1:
      nrb = 2
      do j=1,nord-1
        RRRH(1,j,3,1) = r(2,nrb+j)
      enddo
!
!  ...degrees of freedom of small node 1:
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      nrs = 2
      do i=1,nord-1
!
!  .....big node 1:
        nrb = 2
        do j=1,nord-1
          RRRH(1,j,1,i) = r(nrs+i,nrb+j)
        enddo
      enddo
!
!---------------------------------------------------------------------
!
!  ...small element 2 --- small nodes 2 and 3:
!
!---------------------------------------------------------------------
!
      icl = 0
!
!  ...loop through collocation points:
      do i=1,nord+1
        icl = icl+1
!
!  .....big and small element coordinates of the collocation point
        xibig = (i-1)*0.5d0/nord + 0.5d0
        xisma = (i-1)*1.0d0/nord
!
!  .....find small shape functions at this point:
        call shape1DH(xisma,nord, ndofh,shapsma(1:MAXP+1,icl),void)
!
!  .....find big shape functions at this point:
        call shape1DH(xibig,nord, ndofh,shapbig(1:MAXP+1,icl),void)
      enddo
!
!  ...transpose shapsma and shapbig
      do i=1,nord+1
        do j=1,i
!
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
!  ...solve the collocation system
      n = nord+1
      call dgesv(n,MAXP+1,shapsma,MAXP+1,ip,shapbig,MAXP+1,info)
      if (info.ne.0) then
        write(*,*)'setcnstr_edge_h1: H1 DGETRS RETURNED INFO =',info
        call logic_error(FAILURE,__FILE__,__LINE__)
      endif
!
      r = shapbig
!
!  ...degrees of freedom of small node 2:
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      nrs = 2
      do i=1,nord-1
!
!  .....big node 1:
        nrb = 2
        do j=1,nord-1
          RRRH(1,j,2,i) = r(nrs+i,nrb+j)
        enddo
      enddo
!
!-----------------------------------------------------------------------
!
!  ...clean up machine zeros...
!
!  ...loop through possible orders of approximations
      do n=2,nord
        jpbeg = n-1; jpend = n-1
!
!  .....loop through the parent dof
        do jp=jpbeg,jpend
!
!  .......loop through constrained middle nodes
          do i=1,2
!
!  .........beginning of dof to be zero out
            jbeg = n
!
!  .........loop through the constrained node dof
            do j=jbeg,nord-1
              if (abs(RRRH(1,jp,i,j)).gt.1.d-12) then
                write(*,7001) n,jp,i,j,RRRH(1,jp,i,j)
 7001           format('setcnstr_edge_h1: n,jp,i,j,RRTH = ', &
                        4i3,2x,e12.5)
                stop 1
              endif
              RRRH(1,jp,i,j) = 0.d0
            enddo
          enddo
        enddo
      enddo
!
      return
!
!---------------------------------------------------------------------
!
#if HP3D_DEBUG
!
!  ...this is a test which should be run when installing the code
!     on a new machine
!
   10 write(*,*) 'setcnstr_edge_h1: SET BIG ELEMENT COORDINATE'
      read(*,*) xibig
!
!  ...shape functions of the big element
      call shape1DH(xibig,nord, ndofh,shapbig(1:MAXP+1,1),void)
!
!  ...determine which small element is this, and find the small
!     element coordinate of the point
      if (xibig.le.0.5d0) then
!
!  .....element 1 - verifying small nodes 1,3
        nel = 1
        xisma = 2*xibig
      else
!
!  .....element 2 - verifying small nodes 2,3
        nel = 2
        xisma = (xibig-0.5)*2
      endif
!
!  ...find small shape functions at this point:
      call shape1DH(xisma,nord, ndofh,shapsma(1:MAXP+1,1),void)
!
!  ...verify if the values of the big shape functions match
!     combinations of the small shape functions.......................
!
!  ...loop through the big shape functions - central node:
      nrb = 2 ! skip vertex shape functions
      do i=1,nord-1
!
!  .....initiate value of the linear combination:
        val = 0
        select case(nel)
        case(1)
!
!  .......small node 3:
          nrs = 1
          do j=1,1
            val = val + RRRH(1,i, 3,j)*shapsma(nrs+j,1)
          enddo
!
!  .......small node 1:
          nrs = 2
          do j=1,nord-1
            val = val + RRRH(1,i, 1,j)*shapsma(nrs+j,1)
          enddo
!
        case(2)
!
!  .......small node 3:
          nrs = 0
          do j=1,1
            val = val + RRRH(1,i, 3,j)*shapsma(nrs+j,1)
          enddo
!
!  .......small node 2:
          nrs = 2
          do j=1,nord-1
            val = val + RRRH(1,i, 2,j)*shapsma(nrs+j,1)
          enddo
        end select
!
        write(*,7004) i,shapbig(nrb+i,1),abs(val-shapbig(nrb+i,1))
 7004   format('setcnstr_edge_h1: MIDDLE NODE SHAPE FUNCTION i = ',i2, &
               ' VALUE = ',e12.5, ' DIFFERENCE = ',e12.5)
!
!  ...end of loop through middle node big shape functions
      enddo
!
!  ...left big vertex node
      select case(nel)
      case(1)
        val = RRRH(2,1, 3,1)*shapsma(2,1) &
            + 1.d0          *shapsma(1,1)
      case(2)
        val = RRRH(2,1, 3,1)*shapsma(1,1)
      end select
      write(*,7005) shapbig(1,1),abs(val-shapbig(1,1))
 7005 format('setcnstr_edge_h1: LEFT VERTEX NODE SHAPE FUNCTION ',
              ' VALUE = ',e12.5, ' DIFFERENCE = ',e12.5)
!
!  ...right big vertex node
      select case(nel)
      case(1)
        val = RRRH(3,1, 3,1)*shapsma(2,1)
      case(2)
        val = RRRH(3,1, 3,1)*shapsma(1,1) &
            + 1.0d0         *shapsma(2,1)
      end select
      write(*,7006) shapbig(2,1),abs(val-shapbig(2,1))
 7006 format('setcnstr_edge_h1: RIGHT VERTEX NODE SHAPE FUNCTION ', &
             ' VALUE = ',e12.5, ' DIFFERENCE = ',e12.5)
!
      write(*,*) 'setcnstr_edge_h1: CONTINUE ?(1/0)'
      read(*,*) ians
      if (ians.eq.1) go to 10
!
#endif
!
      end subroutine setcnstr_edge_h1

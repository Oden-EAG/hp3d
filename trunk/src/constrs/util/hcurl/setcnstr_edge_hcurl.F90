!----------------------------------------------------------------------
!
!   routine name       - setcnstr_edge_hcurl
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
   subroutine setcnstr_edge_hcurl
!
      use parameters
      use constraints
!
      implicit none
!
!  ...values of small elements shape functions, big element
!     shape functions, derivatives (not used)
      real(8) :: shapsma(MAXP,MAXP),shapbig(MAXP,MAXP)
!
!  ...pivoting array, work array
      integer :: ip(MAXP)
      real(8) :: r(MAXP,MAXP)
!
      real(8) :: aux,xibig,xisma
      integer :: i,j,icl,is,info
      integer :: nord,nvoid
!
#if HP3D_DEBUG
      real(8) :: val
      integer :: ians,nel
      integer :: iprint
      iprint = 0
#endif
!
!  ...clear:
      RRRE = 0.d0
!
      nord = MAXP
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) 'Compute constraints wrt p = ', nord
      end if
#endif
!
!---------------------------------------------------------------------
!
!  ...small nodes 1:
!
!---------------------------------------------------------------------
!
      icl = 0
!
!  ...loop through collocation points:
      do i=1,nord
        icl = icl+1
!
!  .....local coordinates for big and small elements
        xibig = (i-1)*0.5d0/(nord-1)
        xisma = (i-1)*1.0d0/(nord-1)
!
!  .....find small shape functions at this point:
        call TraceEshapeE(xisma,nord,0, nvoid,shapsma(1:MAXP,icl))
!
!  .....find big shape functions at this point:
        call TraceEshapeE(xibig,nord,0, nvoid,shapbig(1:MAXP,icl))
      enddo
!
!  ...transpose shapsma and shapbig
      do i=1,nord
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
!  ...invert the collocation matrix
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7010)
 7010   format('setcnstr_edge_hcurl: COLLOCATION MATRIX = ')
        do i=1,nord
          write(*,7011) i,shapsma(i,1:nord)
 7011     format('i = ',i3,10e12.5)
        enddo
        call pause
      endif
#endif
      call dgesv(nord,nord,shapsma,MAXP,ip,shapbig,MAXP,info)
      if (info.ne.0) then
        write(*,*)'setcnstr_edge_hcurl: HCURL DGESV RETURNED INFO =',info
        call logic_error(FAILURE,__FILE__,__LINE__)
      endif
!
      r = shapbig
!
!  ...degrees of freedom of small node 1:
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do i=1,nord
!
!  .....big node 1:
        do j=1,nord
          RRRE(1,j,1,i) = r(i,j)
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
      do i=1,nord
        icl = icl+1
!
!  .....big and small element coordinates of the collocation point
        xibig = (i-1)*0.5d0/(nord-1) + 0.5d0
        xisma = (i-1)*1.0d0/(nord-1)
!
!  .....find small shape functions at this point:
        call TraceEshapeE(xisma,nord,0, nvoid,shapsma(1:MAXP,icl))
!
!  .....find big shape functions at this point:
        call TraceEshapeE(xibig,nord,0, nvoid,shapbig(1:MAXP,icl))
      enddo
!
!  ...transpose shapsma and shapbig
      do i=1,nord
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
      call dgesv(nord,nord,shapsma,MAXP,ip,shapbig,MAXP,info)
      if (info.ne.0) then
        write(*,*)'setcnstr_edge_hcurl: HCURL DGESV RETURNED INFO =',info
        call logic_error(FAILURE,__FILE__,__LINE__)
      endif
!
      r = shapbig
!
!  ...degrees of freedom of small node 2:
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do i=1,nord
!
!  .....big node 1:
        do j=1,nord
          RRRE(1,j,2,i) = r(i,j)
        enddo
      enddo
!
!-----------------------------------------------------------------------
!
!  ...clean up machine zeros...
      do is=1,2
        do j=1,nord
          do i=1,nord
!  .........multiply jacobian
            RRRE(1,i,is,j) = 0.5d0*RRRE(1,i,is,j)
            if (abs(RRRE(1,i,is,j)).lt.1.d-12) then
              RRRE(1,i,is,j) = 0.d0
            endif
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
   10 write(*,*) 'setcnstr_edge_hcurl: SET BIG ELEMENT COORDINATE'
      read(*,*) xibig
!
!  ...shape functions of the big element
      call TraceEshapeE(xibig,nord,0, nvoid,shapbig(1:MAXP,1))
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
      call TraceEshapeE(xisma,nord,0, nvoid,shapsma(1:MAXP,1))
!
!  ...verify if the values of the big shape functions match
!     combinations of the small shape functions.......................
!
!  ...loop through the big shape functions - central node:
      do i=1,nord
!
!  .....initiate value of the linear combination:
        val = 0
        do j=1,nord
!  .......tricky part... need to think more about this piola transforms
          val = val + 2.d0*RRRE(1,i, nel,j)*shapsma(j,1)
        enddo
!
        write(*,7004) i,shapbig(i,1),abs(val-shapbig(i,1))
 7004   format('setcnstr_edge_hcurl: MIDDLE NODE SHAPE FUNCTION i = ', &
                i2,' VALUE = ',e12.5, ' DIFFERENCE = ',e12.5)
!
!  ...end of loop through middle node big shape functions
      enddo
!
      write(*,*) 'setcnstr_edge_hcurl: CONTINUE ?(1/0)'
      read(*,*) ians
      if (ians.eq.1) goto 10
!
#endif
!
   end subroutine setcnstr_edge_hcurl

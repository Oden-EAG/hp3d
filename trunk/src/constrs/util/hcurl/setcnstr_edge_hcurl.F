c----------------------------------------------------------------------
c
c   routine name       - setcnstr_edge_hcurl
c
c----------------------------------------------------------------------
c
c   latest revision    - May 12
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
      subroutine setcnstr_edge_hcurl
c
      use parameters
      use constraints
#include "syscom.blk"
c
c  ...values of small elements shape functions, big element
c     shape functions, derivatives (not used) 
      dimension shapsma(MAXP,MAXP),shapbig(MAXP,MAXP), 
     .          void(MAXP), 
c
c  ...pivoting array, work array
     .          ip(MAXP),r(MAXP,MAXP)
c
      iprint=0
c
c  ...clear:
      RRRE = 0.d0
c
      nord = MAXP
      if (iprint.eq.1) then
        write(*,*) 'Compute constraints wrt p = ', nord
      end if
c
c---------------------------------------------------------------------
c
c  ...small nodes 1:
c
c---------------------------------------------------------------------
c
      icl = 0 
c
c  ...loop through collocation points:
      do i=1,nord
        icl = icl+1
c
c  .....local coordinates for big and small elements
        xibig = (i-1)*0.5d0/(nord-1)
        xisma = (i-1)*1.0d0/(nord-1)
c
c  .....find small shape functions at this point:
        call shapEbe(xisma,nord,0, void, shapsma(1:MAXP,icl))
c
c  .....find big shape functions at this point:
        call shapEbe(xibig,nord,0, void, shapbig(1:MAXP,icl))
      enddo
c
c  ...transpose shapsma and shapbig
      do i=1,nord
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
c  ...invert the collocation matrix
      if (iprint.eq.1) then
        write(*,7010)
 7010   format('setcnstr_edge: COLLOCATION MATRIX = ')
        do i=1,nord
          write(*,7011) i,shapsma(i,1:nord)
 7011     format('i = ',i3,10e12.5)
        enddo
        call pause
      endif

      call dgesv(nord,nord,shapsma,MAXP,ip,shapbig,MAXP,info)
      if (info.ne.0) then
        write(*,*)'setcnstr_edge: HCURL DGESV RETURNED INFO =',info
        call logic_error(FAILURE,
     .    __FILE__,__LINE__)
      endif
c
      r = shapbig
c
c  ...degrees of freedom of small node 1:
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do i=1,nord
c
c  .....big node 1:
        do j=1,nord
          RRRE(1,j,1,i) = r(i,j)            
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
      do i=1,nord
        icl = icl+1
c
c  .....big and small element coordinates of the collocation point
        xibig = (i-1)*0.5d0/(nord-1) + 0.5d0
        xisma = (i-1)*1.0d0/(nord-1)
c
c  .....find small shape functions at this point:
        call shapEbe(xisma,nord,0, void, shapsma(1:MAXP,icl))
c
c  .....find big shape functions at this point:
        call shapEbe(xibig,nord,0, void, shapbig(1:MAXP,icl))
      enddo
c
c  ...transpose shapsma and shapbig
      do i=1,nord
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
      call dgesv(nord,nord,shapsma,MAXP,ip,shapbig,MAXP,info)
      if (info.ne.0) then
        write(*,*)'setcnstr_edge: HCURL DGESV RETURNED INFO =',info
        call logic_error(FAILURE,
     .    __FILE__,__LINE__)
      endif
c
      r = shapbig
c
c  ...degrees of freedom of small node 2:
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do i=1,nord
c
c  .....big node 1:
        do j=1,nord
          RRRE(1,j,2,i) = r(i,j)            
        enddo
      enddo
c
c-----------------------------------------------------------------------
c
c  ...clean up machine zeros...
      do is=1,2
        do j=1,nord
          do i=1,nord
c  .........multiply jacobian            
            RRRE(1,i,is,j) = 0.5d0*RRRE(1,i,is,j)
            if (abs(RRRE(1,i,is,j)).lt.1.d-12) then
              RRRE(1,i,is,j) = 0.d0
            endif
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
      call shapEbe(xibig,nord,0, void, shapbig(1:MAXP,1))
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
      call shapEbe(xisma,nord,0, void, shapsma(1:MAXP,1))
c
c  ...verify if the values of the big shape functions match
c     combinations of the small shape functions.......................
c
c  ...loop through the big shape functions - central node:
      do i=1,nord
c
c  .....initiate value of the linear combination:
        val = 0
        do j=1,nord
c  .......tricky part... need to think more about this piola transforms
          val = val + 2.d0*RRRE(1,i, nel,j)*shapsma(j,1)
        enddo
c
        write(*,7004) i,shapbig(i,1),abs(val-shapbig(i,1))
 7004   format('setcnstr_edge: MIDDLE NODE SHAPE FUNCTION i = ',i2,
     .         ' VALUE = ',e12.5, ' DIFFERENCE = ',e12.5)
c
c  ...end of loop through middle node big shape functions
      enddo
c
      write(*,*) 'setcnstr: CONTINUE ?(1/0)'
      read(*,*) ians
      if (ians.eq.1) go to 10
c
c
      end subroutine

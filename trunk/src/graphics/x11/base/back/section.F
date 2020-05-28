c----------------------------------------------------------------------
c
c   routine name       - section
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Nov 06
c
c   purpose            - routine displays the solution along
c                        a line
c
c   required  routines -
c
c----------------------------------------------------------------------
c
      subroutine section
c
      use graphmod
#include "syscom.blk"
c
c  ...origin of systme of coordinates
      dimension origin(2)
      character*5 comment5
      character*9 comment9
      character*10 comment10
      character*20 comment20
c
c  ...point and direction
      dimension xp(3),vect(3)
c
c  ...work array
      parameter (na=1000) !array dimension
      parameter (nrval=2) !number of quantities to be drawn
      double precision ,allocatable :: array(:,:)
c
c----------------------------------------------------------------------
c
      iprint=0
c
      write(*,*) 'section: SET xp,vect'
      xp(1:3) = .0000000001
      vect = (/1.d0,0.d0,0.d0/)
      allocate(array(4,na),  STAT=ic )
      call compute_section(xp,vect,nrval,na, nrsegm,array)
c
c  ...select window
      write(*,*) 'PLEASE SELECT KIND OF WINDOW FOR GRAPHICS'
      write(*,*) '1-B&W, 2-Gray16, 3-Gray256, 4-Color16, 5-Color256'
      read(*,*) iwin
      call initwin(iwin)
c
      nblack = NPCOL(2)
c
c----------------------------------------------------------------------
c
c  ...step 1: Determine the maximum and minimum values of
c             the solution
c
      y_max = -1.d30
      y_min =  1.d30
c
c  ...loop through the points
      do l=1,nrsegm
        y_max = max(y_max,array(3,l))
        y_min = min(y_min,array(3,l))
        y_max = max(y_max,array(4,l))
        y_min = min(y_min,array(4,l))
      enddo
      x_min = array(1,1)
      x_max = array(2,nrsegm)
   10 continue
      write(*,7045) x_min,x_max,y_min,y_max
 7045 format('section: x_min,x_max,y_min,y_max = ',2(2e12.5,2x))
      dx = x_max - x_min
      dy = max(y_max - y_min,1.d-5)
      dx = dx*1.1d0
      dy = dy*1.1d0
c
c  ...open the window
      call selwin(iwin)
c
c----------------------------------------------------------------------
c
c  ...step 2: set up coordinates of the origin and draw horizontal and
c             vertical axes
c
      origin(1) = 2.d0*RMARGIN
      origin(2) = 2.d0*RMARGIN
      xlen = XLENGTH - 2.d0*RMARGIN
      ylen = YLENGTH - 2.d0*RMARGIN
c
c  ...zero elevation
      y1 = origin(2)+ (0.d0 - y_min)/dy*ylen
      if ((y1.ge.origin(2)).and.(y1.le.origin(2)+ylen)) then
        call drawline(origin(1),y1,
     .                origin(1)+xlen,y1,nblack)
      else
        call drawline(origin(1),origin(2),
     .                origin(1)+xlen,origin(2),nblack)
      endif
      call drawline(origin(1),origin(2),
     .              origin(1),origin(2)+ylen,nblack)
c
c
c
c----------------------------------------------------------------------
c
c  ...step 3: draw the graph
c
c  ...loop through the components
      nrdraw = 1
      do i=1,Nrdraw
c
c  .....loop through segments
        do n=1,nrsegm
          x1 = origin(1) + (array(1,n) - x_min)/dx*xlen
          y1 = origin(2) + (array(3,n) - y_min)/dy*ylen
          x2 = origin(1) + (array(2,n) - x_min)/dx*xlen
          y2 = origin(2) + (array(4,n) - y_min)/dy*ylen
          call draw_line(x1,y1,x2,y2,nblack,i-1)
c
c  .....end of loop through segments
        enddo
      enddo
c
      write(*,*) 'section: NEW RANGE,COMMENTS,EXIT ? (1/2/0)'
      read(*,*) ians
c
      select case(ians)
      case(0)
c
c  .....close the window
        call closwind(1)
        return
      case(1)
        write(*,*) 'section: SET y_min,y_max '
        read(*,*) y_min, y_max
        go to 10
      case(2)
c
        write(*,*) 'section: SET x-axis description'
        read(*,*) comment10
        x = origin(1) + xlen
ccc        y = origin(2) + 3*RMARGIN/8.d0
        y = origin(2) + (0.d0 - y_min)/dy*ylen + RMARGIN/8.d0
        call symbol(x,y,10.d0,comment10,0.d0,10,1)
c
        write(*,*) 'section: SET y-axis description'
        read(*,*) comment10
        write(*,*) comment10
        x = origin(1) + 3*RMARGIN/8.d0
        y = origin(2) + ylen + RMARGIN/4.d0
        call symbol(x,y,10.d0,comment10,0.d0,10,1)
c
c  .....marks on the horizontal axis
        write(*,*) 'section: SET nr_marks'
        read(*,*) nr_marks
        do i=0,nr_marks
          x1 = origin(1) + float(i)/float(nr_marks)*xlen
          y1 = origin(2) + (0.d0 - y_min)/dy*ylen
          x2 = x1
          y2 = y1 - RMARGIN/2.d0
          call drawline(x1,y1,x2,y2,nblack)
          x = x1 + RMARGIN/8.d0
          y = y2
c
c  .......find out the value to write out
          xx = x_min + (x_max-x_min)/nr_marks*i
          write(*,*) 'xx = ',xx
          write(comment9,1002) xx
 1001     format(i5)
          call symbol(x,y,10.d0,comment9,0.d0,9,1)
        enddo
c
c  .....marks on the vertical axis
        nr_marks=10
        do i=0,nr_marks
          x1 = origin(1) - RMARGIN*.75
          y1 = origin(2) + float(i)/float(nr_marks)*ylen
          x2 = origin(1)
          y2 = y1
          call drawline(x1,y1,x2,y2,nblack)
c
c  .......find out the value to write out
          val = y_min + (y_max - y_min)/nr_marks*i
          x = origin(1) - 2*RMARGIN
          y = y2 + RMARGIN/4.d0
          write(comment9,1002) val
 1002     format(f9.3)
          call symbol(x,y,10.d0,comment9,0.d0,9,1)
        enddo
c
c  .....write out additional comments
        write(*,*) 'rates: ADDITIONAL COMMENTS ?(1/0)'
        read(*,*) ians
        if (ians.eq.1) then
          x1 =  origin(1)+ .7d0*xlen
          x2 = x1 + RMARGIN
          do nplot=1,Nrdraw
            y1 = YLENGTH - 3.*RMARGIN - nplot*RMARGIN
            y2 = y1
            call draw_line(x1,y1,x2,y2,nblack,nplot-1)
            write(*,7003) nplot
 7003       format('rates: nplot = ',i2,
     .             ' WRITE THE COMMENT (max 20 characters)')
            read(*,7004) comment20
 7004       format(a20)
            x = x2 + RMARGIN/8.
            y = y1 - RMARGIN/8.
            call symbol(x,y,20.d0,comment20,0.d0,20,1)
          enddo
        endif
        x = origin(1)+ .7d0*xlen
        y = YLENGTH - 2.*RMARGIN
c
c  .....close the window
        call closwind(1)
c
        return
      end select
c
 999  deallocate(array)
c
c
      end

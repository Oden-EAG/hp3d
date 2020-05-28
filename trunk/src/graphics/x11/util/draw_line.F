c----------------------------------------------------------------------
c
c   routine name       - draw_line
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Sep 01
c
c   purpose            - routine draws various types of lines
c
c   arguments:
c     in:
c             X1,Y2,X2,Y2 - coordinates of endpoints of a line
c             Ic       - color to draw the line with
c             Itype    - type of a line
c
c   required  routines -
c
c----------------------------------------------------------------------
c
      subroutine draw_line(X1,Y1,X2,Y2,Ic,Itype)
c
c
      use graphmod
#include "syscom.blk"
c
c  ...compute the line increment
      dl = RMARGIN/16.d0
c
      if (Itype.eq.0) then
c
c  .....solid line case
        call drawline(X1,Y1,X2,Y2,Ic)
        return
      endif
c
c  ...compute the increments in the line and space
      select case(Itype)
c
        case(1)
        ds = dl
c
        case(2)
        ds = 3.d0*dl
c
        case(3)
        ds = 7.d0*dl
      end select
c
      xl1 = X1
      yl1 = Y1
      do iter=1,1000
        dsl = sqrt((X2-xl1)**2+(Y2-yl1)**2)
        if (ds.le.dsl) then
          xl2 = xl1 + (X2-xl1)/dsl*ds
          yl2 = yl1 + (y2-yl1)/dsl*ds
          call drawline(xl1,yl1,xl2,yl2,Ic)
          xl1 = xl2
          yl1 = yl2
        else
          call drawline(xl1,yl1,X2,Y2,Ic)
          exit
        endif
c
        dsl = sqrt((X2-xl1)**2+(Y2-yl1)**2)
        if (ds.le.dsl) then
          xl1 = xl1 + (X2-xl1)/dsl*ds
          yl1 = yl1 + (y2-yl1)/dsl*ds
        else
          exit
        endif
      enddo
c
c
      return
      end





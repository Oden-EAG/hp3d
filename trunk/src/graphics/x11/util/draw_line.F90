#if HP3D_USE_X11

!----------------------------------------------------------------------
!
!   routine name       - draw_line
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine draws various types of lines
!
!   arguments:
!     in:
!             (X1,Y1)  - coordinates of endpoints of a line
!             (X2,Y2)
!             Ic       - color to draw the line with
!             Itype    - type of a line
!
!----------------------------------------------------------------------
!
      subroutine draw_line(X1,Y1,X2,Y2,Ic,Itype)
!
      use graphmod
!
      implicit none
!
      real(8) :: X1,Y1,X2,Y2
      integer :: Ic,Itype
!
      real(8) :: dl,ds,dsl,xl1,xl2,yl1,yl2
      integer :: iter
!
!  ...compute the line increment
      dl = RMARGIN/16.d0
!
      if (Itype.eq.0) then
!
!  .....solid line case
        call drawline(X1,Y1,X2,Y2,Ic)
        return
      endif
!
!  ...compute the increments in the line and space
      select case(Itype)
!
        case(1)
        ds = dl
!
        case(2)
        ds = 3.d0*dl
!
        case(3)
        ds = 7.d0*dl
      end select
!
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
!
        dsl = sqrt((X2-xl1)**2+(Y2-yl1)**2)
        if (ds.le.dsl) then
          xl1 = xl1 + (X2-xl1)/dsl*ds
          yl1 = yl1 + (y2-yl1)/dsl*ds
        else
          exit
        endif
      enddo
!
      end subroutine draw_line

#endif

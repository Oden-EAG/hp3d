c----------------------------------------------------------------------
c
c   routine name       - dpborder
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Mar 02
c
c   purpose            - routine displays the border of the window
c                        and the color table
c
c   usage              - call dpborder(Ibf)
c
c   arguments 
c     in: 
c                  Ibf - a flag whether to display color table
c
c   required  routines - none
c
c----------------------------------------------------------------------
c
      subroutine dpborder(Ibf)
c
      use graphmod
#include "syscom.blk"
c
c  ...screen coordinates of up to four points, components of
c     a 3D unit vector, components of the same vector after
c     transformation to screen coordinates
      dimension xgr(2,4),coord(3),cgr(3)
c
      iprint=0
      if (iprint.eq.1) then
        write(*,*) 'dpborder: Ibf = ',Ibf
      endif
c
      nblack=NPCOL(2)
c
c  ...fill in borders with white color...
      ncol=NPCOL(1)
      c0 = 0.d0
      c1x = rwindl
      c1y = rwindh
      c2 = 0.9d0*rmargin
      c3 = c2
      c4 = xlength+1.1d0*rmargin
      c5 = ylength+1.1d0*rmargin
      xgr(1,1)=c0
      xgr(2,1)=c0
      xgr(1,2)=c1x
      xgr(2,2)=c0
      xgr(1,3)=c1x
      xgr(2,3)=c3
      xgr(1,4)=c0
      xgr(2,4)=c3
      call fillpoly(4,xgr,ncol,ncol)
      xgr(1,1)=c0
      xgr(2,1)=c0
      xgr(1,2)=c2
      xgr(2,2)=c0
      xgr(1,3)=c2
      xgr(2,3)=c1y
      xgr(1,4)=c0
      xgr(2,4)=c1y
      call fillpoly(4,xgr,ncol,ncol)
      xgr(1,1)=c1x
      xgr(2,1)=c1y
      xgr(1,2)=c1x
      xgr(2,2)=c0
      xgr(1,3)=c4
      xgr(2,3)=c0
      xgr(1,4)=c4
      xgr(2,4)=c1y
      call fillpoly(4,xgr,ncol,ncol)
      xgr(1,1)=c1x
      xgr(2,1)=c1y
      xgr(1,2)=c0
      xgr(2,2)=c1y
      xgr(1,3)=c0
      xgr(2,3)=c5
      xgr(1,4)=c1x
      xgr(2,4)=c5
      call fillpoly(4,xgr,ncol,ncol)
c
c  ...specify the box dimensions for the system of coordinates icon
      xgr(1,1)=rwindl-rmargin
      xgr(2,1)=rmargin
      XY_BOX(1,1,0) = xgr(1,1) - rmargin/2.d0
      XY_BOX(2,1,0) = xgr(2,1) - rmargin/2.d0
      XY_BOX(1,2,0) = xgr(1,1) + rmargin/2.d0
      XY_BOX(2,2,0) = xgr(2,1) - rmargin/2.d0
      XY_BOX(1,3,0) = xgr(1,1) + rmargin/2.d0
      XY_BOX(2,3,0) = xgr(2,1) + rmargin/2.d0
      XY_BOX(1,4,0) = xgr(1,1) - rmargin/2.d0
      XY_BOX(2,4,0) = xgr(2,1) + rmargin/2.d0
c
c  ...draw the icon for the system of coordinates...
      nacol=NPCOL(2)
      nscol=NPCOL(2)
      coord(1)=1.d0
      coord(2)=0.d0
      coord(3)=0.d0
      call trobs(coord,cgr)
      xgr(1,2)=xgr(1,1)+cgr(1)*0.9d0*rmargin
      xgr(2,2)=xgr(2,1)+cgr(2)*0.9d0*rmargin
      call drawline(xgr(1,1),xgr(2,1),xgr(1,2),xgr(2,2),nacol)
      height=1.
      angle=0.
      call symbol(xgr(1,2),xgr(2,2),height,'x',angle,1,nscol)
      coord(1)=0.d0
      coord(2)=1.d0
      coord(3)=0.d0
      call trobs(coord,cgr)
      xgr(1,2)=xgr(1,1)+cgr(1)*0.9d0*rmargin
      xgr(2,2)=xgr(2,1)+cgr(2)*0.9d0*rmargin
      call drawline(xgr(1,1),xgr(2,1),xgr(1,2),xgr(2,2),nacol)
      height=1.
      angle=0.
      call symbol(xgr(1,2),xgr(2,2),height,'y',angle,1,nscol)
      coord(1)=0.d0
      coord(2)=0.d0
      coord(3)=1.d0
      call trobs(coord,cgr)
      xgr(1,2)=xgr(1,1)+cgr(1)*0.9d0*rmargin
      xgr(2,2)=xgr(2,1)+cgr(2)*0.9d0*rmargin
      call drawline(xgr(1,1),xgr(2,1),xgr(1,2),xgr(2,2),nacol)
      height=1.
      angle=0.
      call symbol(xgr(1,2),xgr(2,2),height,'z',angle,1,nscol)
c
c  ...draw the color table
      if (Ibf.eq.0) then
c
c  .....loop through the color boxes of p-approximation
        vbox=(rwindh-4.d0*rmargin)/8.d0
        do i=1,8
          xgr(1,1) = c1x-1.5d0*rmargin
          xgr(2,1) = 3.d0*rmargin+(i-1)*vbox
          xgr(1,2) = c1x-0.5d0*rmargin
          xgr(2,2) = xgr(2,1)
          xgr(1,3) = xgr(1,2)
          xgr(2,3) = xgr(2,1)+vbox
          xgr(1,4) = xgr(1,1) 
          xgr(2,4) = xgr(2,3)
          call fillpoly(4,xgr,NPCOL(i+2),NPCOL(i+2))
          do j=1,2
            XY_BOX(j,1:4,i) = xgr(j,1:4)
          enddo
        enddo
c
c  .....loop one more time and redraw the lines
        vbox=(rwindh-4.d0*rmargin)/8.d0
        do i=1,9
          xgr(1,1) = c1x-1.5d0*rmargin
          xgr(2,1) = 3.d0*rmargin+(i-1)*vbox
          xgr(1,2) = c1x-0.5d0*rmargin
          xgr(2,2) = xgr(2,1)
          call drawline(xgr(1,1),xgr(2,1),xgr(1,2),xgr(2,2),nblack)
        enddo
        xgr(2,1) = 3.d0*rmargin
        xgr(2,2) = 3.d0*rmargin+8*vbox
        call drawline(xgr(1,1),xgr(2,1),xgr(1,1),xgr(2,2),nblack)
        call drawline(xgr(1,2),xgr(2,1),xgr(1,2),xgr(2,2),nblack)
c
      elseif (ibf.gt.0) then
c
c  .....loop through the color boxes of isolines
        vbox=(rwindh-4.d0*rmargin)/float(NR_COLORS-10)*1.001d0
        do i=1,NR_COLORS-10
          xgr(1,1) = c1x-1.5d0*rmargin
          xgr(2,1) = 3.d0*rmargin+(i-1)*vbox
          xgr(1,2) = c1x-0.5d0*rmargin
          xgr(2,2) = xgr(2,1)
          xgr(1,3) = xgr(1,2)
          xgr(2,3) = xgr(2,1)+vbox
          xgr(1,4) = xgr(1,1) 
          xgr(2,4) = xgr(2,3)
          call fillpoly(4,xgr,NPCOL(i+10),NPCOL(i+10))
        enddo
      endif
c
      return
      end

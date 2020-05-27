c-----------------------------------------------------------------------c
c   routine name       - dpallcu
c
c-----------------------------------------------------------------------c
c   computer           - machine independent
c
c   latest revision    - Dec, 01
c
c   purpose            - routine draws all blocks' edges
c
c   usage              - call dpallcu
c
c   arguments          -none
c
c   required  routines -
c
c-----------------------------------------------------------------------c
      subroutine dpallcu
c
      use graphmod
      use GMP
c
#include "syscom.blk"
c
      dimension dxdxi(3)
      dimension xcoord(3,2),xobs(3,2)
c
      data small,one /1.d-15,1.d0/
c
      ncol=npcol(2)
c
c  ...loop through curves
      do 150 nc=1,NRCURVE
c
c  ...get first point and transform its coordinates
        xi = 0.d0
        call curve(nc,xi, xcoord(1,2),dxdxi)
        call trobs(xcoord(1,2),xobs(1,2))
c
c  ...loop through subsegments
        do 40 i=1,NRSUB
          do 10 j=1,3
            xobs(j,1) = xobs(j,2)
   10     continue
          xi = xi + DX
          call curve(nc,xi, xcoord(1,2),dxdxi)
          call trobs(xcoord(1,2),xobs(1,2))
c
c  .......rescale
          do 30 k=1,2
          do 30 ivar=1,2
            XY(ivar,k) = (xobs(ivar,k)-XCIM(ivar))/DIMIM*SIZE
     .                   + XCWIN(ivar)
   30     continue
c
c  .......draw the segment of line
          call drawline(xy(1,1),xy(2,1),xy(1,2),xy(2,2),ncol)
c
c  ...end of the loop through segments
   40   continue
  150 continue
c
      return
      end




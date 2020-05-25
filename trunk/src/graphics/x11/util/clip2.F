      subroutine clip2(Trc,Imod)
c
c----------------------------------------------------------------------
c
c   routine name       - clip2
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - jan, 1992
c
c   purpose            - routine performs 2D clipping of visible
c                        triangles
c
c   usage              - call clip2(Trc,Imod)
c
c   arguments 
c            in :  Trc - triangle's vertices coordinates
c           out : Imod - clipping indicator
c
c   required  routines - none
c
c----------------------------------------------------------------------
c
      use graphmod
#include "syscom.blk"
c
      dimension Trc(3,3)
c
      Imod = 0
c
c  ...get triangle's vertices coordinates, rescale them and
c     if any vertex far from the window ignore the triangle
      do k=1,3
         do ivar=1,2
            xy(ivar,k) = (Trc(ivar,k)-XCIM(ivar))/DIMIM
            if(xy(ivar,k).lt.-1.5d0) return
            if(xy(ivar,k).gt.1.5d0) return
         enddo
      enddo
   10 continue
c
c  .....if visible 
      Imod = 1 
c         
      return
      end



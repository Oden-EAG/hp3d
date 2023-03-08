c----------------------------------------------------------------------
c
c   routine name       - clip2
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - routine performs 2D clipping of visible
c                        triangles
c
c   arguments
c            in :  Trc - triangle's vertices coordinates
c           out : Imod - clipping indicator
c
c----------------------------------------------------------------------
      subroutine clip2(Trc, Imod)
c
      use graphmod, only: XY,DIMIM,XCIM
c
      implicit none
c
      real(8), intent(in)  :: Trc(3,3)
      integer, intent(out) :: Imod
c
      integer :: k,ivar
c
      Imod = 0
c
c  ...get triangle's vertices coordinates, rescale them and
c     if any vertex far from the window ignore the triangle
      do k=1,3
         do ivar=1,2
            XY(ivar,k) = (Trc(ivar,k)-XCIM(ivar))/DIMIM
            if(xy(ivar,k).lt.-1.5d0) return
            if(xy(ivar,k).gt. 1.5d0) return
         enddo
      enddo
c
c  .....if visible
      Imod = 1
c
      end subroutine clip2

#if HP3D_USE_X11

!----------------------------------------------------------------------
!
!   routine name       - clip2
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine performs 2D clipping of visible
!                        triangles
!
!   arguments
!            in :  Trc - triangle's vertices coordinates
!           out : Imod - clipping indicator
!
!----------------------------------------------------------------------
   subroutine clip2(Trc, Imod)
!
      use graphmod, only: XY,DIMIM,XCIM
!
      implicit none
!
      real(8), intent(in)  :: Trc(3,3)
      integer, intent(out) :: Imod
!
      integer :: k,ivar
!
      Imod = 0
!
!  ...get triangle's vertices coordinates, rescale them and
!     if any vertex far from the window ignore the triangle
      do k=1,3
         do ivar=1,2
            XY(ivar,k) = (Trc(ivar,k)-XCIM(ivar))/DIMIM
            if(xy(ivar,k).lt.-1.5d0) return
            if(xy(ivar,k).gt. 1.5d0) return
         enddo
      enddo
!
!  .....if visible
      Imod = 1
!
   end subroutine clip2

#endif

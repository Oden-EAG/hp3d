#if HP3D_USE_X11

!-----------------------------------------------------------------------
!
!   routine name       - trobs
!
!-----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine transforms cartesian coordinates
!                        of a point into observer's coordinates; the
!                        first two correspond to the projection plane,
!                        the third to the direction of projection;
!                        the third coordinate decides about the order
!                        in which objects are displayed
!
!   arguments :
!     in:
!              Xcart   - 3D cartesian coordinates of a point
!     out:
!              Xobs    - 2D observer's coordinates of a point
!
!-----------------------------------------------------------------------
!
   subroutine trobs(Xcart, Xobs)
!
      use graphmod, only: RMTR
!
      implicit none
!
      real(8), intent(in)  :: Xcart(3)
      real(8), intent(out) :: Xobs(3)
!
      real(8) :: s
      integer :: i,j
!
      do i=1,3
        s = 0.d0
        do j=1,3
          s = s + RMTR(i,j)*Xcart(j)
        enddo
        Xobs(i) = s
      enddo
!
   end subroutine trobs

#endif

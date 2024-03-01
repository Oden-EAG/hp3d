#if HP3D_USE_X11

!-----------------------------------------------------------------------!
!   routine name       - dpallcu
!
!-----------------------------------------------------------------------!
!
!   latest revision    - Mar 2023
!
!   purpose            - routine draws all blocks' edges
!
!-----------------------------------------------------------------------!
      subroutine dpallcu
!
      use graphmod
      use GMP
!
      implicit none
!
      real(8) :: xi,dxdxi(3)
      real(8) :: xcoord(3,2),xobs(3,2)
!
      integer :: i,ivar,j,k,ncol,nc
!
      real(8) :: small,one
      data small,one /1.d-15,1.d0/
!
      ncol=npcol(2)
!
!  ...loop through curves
      do nc=1,NRCURVE
!
!  .....get first point and transform its coordinates
        xi = 0.d0
        call curve(nc,xi, xcoord(1,2),dxdxi)
        call trobs(xcoord(1,2),xobs(1,2))
!
!  .....loop through subsegments
        do i=1,NRSUB
          do j=1,3
            xobs(j,1) = xobs(j,2)
          enddo
          xi = xi + DX
          call curve(nc,xi, xcoord(1,2),dxdxi)
          call trobs(xcoord(1,2),xobs(1,2))
!
!  .......rescale
          do k=1,2; do ivar=1,2
            XY(ivar,k)=(xobs(ivar,k)-XCIM(ivar))/DIMIM*SIZE+XCWIN(ivar)
          enddo; enddo
!
!  .......draw the segment of line
          call drawline(xy(1,1),xy(2,1),xy(1,2),xy(2,2),ncol)
!
!  .....end of the loop through segments
        enddo
!  ...end of the loop through curves
      enddo
!
!
      end subroutine dpallcu

#endif

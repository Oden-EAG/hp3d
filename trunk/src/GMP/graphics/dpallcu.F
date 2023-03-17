c-----------------------------------------------------------------------c
c   routine name       - dpallcu
c
c-----------------------------------------------------------------------c
c
c   latest revision    - Mar 2023
c
c   purpose            - routine draws all blocks' edges
c
c-----------------------------------------------------------------------c
      subroutine dpallcu
c
      use graphmod
      use GMP
c
      implicit none
c
      real(8) :: xi,dxdxi(3)
      real(8) :: xcoord(3,2),xobs(3,2)
c
      integer :: i,ivar,j,k,ncol,nc
c
      real(8) :: small,one
      data small,one /1.d-15,1.d0/
c
      ncol=npcol(2)
c
c  ...loop through curves
      do nc=1,NRCURVE
c
c  .....get first point and transform its coordinates
        xi = 0.d0
        call curve(nc,xi, xcoord(1,2),dxdxi)
        call trobs(xcoord(1,2),xobs(1,2))
c
c  .....loop through subsegments
        do i=1,NRSUB
          do j=1,3
            xobs(j,1) = xobs(j,2)
          enddo
          xi = xi + DX
          call curve(nc,xi, xcoord(1,2),dxdxi)
          call trobs(xcoord(1,2),xobs(1,2))
c
c  .......rescale
          do k=1,2; do ivar=1,2
            XY(ivar,k)=(xobs(ivar,k)-XCIM(ivar))/DIMIM*SIZE+XCWIN(ivar)
          enddo; enddo
c
c  .......draw the segment of line
          call drawline(xy(1,1),xy(2,1),xy(1,2),xy(2,2),ncol)
c
c  .....end of the loop through segments
        enddo
c  ...end of the loop through curves
      enddo
c
c
      end subroutine dpallcu

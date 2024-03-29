#if HP3D_USE_X11

!----------------------------------------------------------------------
!
!   latest revision     - Mar 2023
!
!   purpose             - routine displays a manifold in 2D space
!
!   argument in : Iwind - type of the window
!
!   required  routines  - curve, openwind, clrwind, drawline, closwind
!
!----------------------------------------------------------------------
!
   subroutine object2(Iwind)
!
      use graphmod
      use GMP
      implicit none
!
      integer :: Iwind
!
      real(8) :: x1(3),dxdx(3),x2(3)
      real(8) :: x(3),xmin(2),xmax(2)
!
      integer :: i,ichoice,ivar,j,nbgcol,ncol,nc,nrsub1
      real(8) :: dimobj,q,q1,rncp1,rncp2,rnsc,xc1,xc2,xi
!
      real(8) :: bigp,bign,small
      data bigp,bign,small /1.d30,-1.d30,1.d-14/
!
!  ...set colour
!     ...for curves  ...
      ncol = npcol(8)
!     ...and for background
      nbgcol = npcol(1)
!
!  ...set up number of line segments for drawing a curve
      nrsub1=2
      NRSUB = 4
!      write(*,*)'NRSUB=',NRSUB
      DX=1.d0/NRSUB
!
!  ...determine bounds for the picture
      do ivar=1,2
        xmin(ivar) = bigp
        xmax(ivar) = bign
      enddo
!
!  ...loop through curves
      do nc = 1,NRCURVE
!
!  .....loop through the points on the curve
        xi = -DX
        do i=1,NRSUB+1
          xi = xi + DX
          call curve(nc,xi, x,dxdx)
          do j=1,2
            xmin(j) = min(xmin(j),x(j))
            xmax(j) = max(xmax(j),x(j))
          enddo
        enddo
      enddo
      do j=1,2
        DIMOB(j) = (xmax(j) - xmin(j))/2.d0
        XCENTR(j) = (xmax(j) + xmin(j))/2.d0
      enddo
        XCIM(1) = XCENTR(1)
        XCIM(2) = XCENTR(2)
        xc1 = XCENTR(1)
        xc2 = XCENTR(2)
!
!  .....decide whether to scale by x or y axes
        q = DIMOB(1)/DIMOB(2)
        q1 = xlength/ylength
        if (q.gt.q1) then
          DIMIM = DIMOB(1)
          SIZE = xlength/2.d0
        else
          DIMIM = DIMOB(2)
          SIZE = ylength/2.d0
        endif
!
        dimobj = DIMIM
        XCWIN(1) = rmargin + xlength/2.d0
        XCWIN(2) = rmargin + ylength/2.d0
!
!  ...in an infinite loop
  100 continue
!
!  ...display the object
      call selwin(Iwind)
!
!  ...loop through curves
      do nc = 1,NRCURVE
!
!  .....loop through the points on the curve
        xi = 0.d0
        call curve(nc,xi, x2,dxdx)
        do i=1,NRSUB
          do j=1,2
            x1(j) = x2(j)
          enddo
          xi = xi + DX
          call curve(nc,xi, x2,dxdx)
!
!  .......scale the points
          do  j=1,2
            XY(j,1) = (x1(j)-XCIM(j))/DIMIM*SIZE + XCWIN(j)
            XY(j,2) = (x2(j)-XCIM(j))/DIMIM*SIZE + XCWIN(j)
          enddo
!          write(*,*)"XY(1,1),XY(2,1)",XY(1,1),XY(2,1),XY(1,2),XY(2,2)
          call drawline(XY(1,1),XY(2,1),XY(1,2),XY(2,2),ncol)
        enddo
!
!  ...end of the loop through curves
      enddo
!
!  ...close window
      write(*,*) 'PLEASE CLICK THE MOUSE INSIDE THE GRAPHICS'
      write(*,*) 'WINDOW TO CONTINUE ...'
      call closwind(iwindnum)
!
      write(*,*) 'SELECT OPTION :'
      write(*,*) '0 - EXIT'
      write(*,*) '1 - CHANGE THE DEGREE OF COARSENESS FOR OVALS'
      write(*,*) '2 - CHANGE THE CENTRAL POINT OF THE IMAGE'
      write(*,*) '3 - RESCALE THE IMAGE'
      write(*,*) '4 - DISPLAY THE WHOLE OBJECT'
      read(*,*) ichoice
!
      if (ichoice.eq.0) then
        return
      endif
!
      if (ichoice.eq.1) then
        write(*,*)'SET NEW DEGREE OF COARSENESS FOR OVALS.'
        write(*,*)'MUST BE BETWEEN 1 AND 4'
        write(*,*)'CURRENT VALUE IS ', nrsub1
        read(*,*) nrsub1
        if (nrsub1.gt.4) then
          write(*,*)'HAVE CHANGED TO 4'
          nrsub1 = 4
        endif
        if (nrsub1.lt.1) then
          write(*,*)'HAVE CHANGED TO 1'
          nrsub1 = 1
        endif
        NRSUB=2**nrsub1
        DX=1.d0/NRSUB
      endif
      if (ichoice.eq.2) then
        write(*,*) 'SET NEW COORDINATES OF THE CENTRAL POINT.'
        write(*,*) 'MUST BE BETWEEN -1 AND 1 (WRT ACTUAL WINDOW).'
        read(*,*) rncp1,rncp2
        XCIM(1)=XCIM(1)+rncp1*DIMIM
        XCIM(2)=XCIM(2)+rncp2*DIMIM
      endif
      if (ichoice.eq.3) then
        write(*,*) 'SET MAGNIFICATION FACTOR'
        read(*,*) rnsc
        rnsc=dabs(rnsc)
        DIMIM=DIMIM/rnsc
      endif
      if (ichoice.eq.4) then
        DIMIM = dimobj
        XCIM(1) = xc1
        XCIM(2) = xc2
      endif
!
      goto 100
!
   end subroutine object2

#endif

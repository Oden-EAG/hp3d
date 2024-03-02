#if HP3D_USE_X11

!----------------------------------------------------------------------
!
!   routine name       - dpblno
!
!----------------------------------------------------------------------
!
!   latest revision    - Mar 2023
!
!   purpose            - routine displays the nicknames of all blocks
!
!   usage              - call dpblno
!
!----------------------------------------------------------------------
!
   subroutine dpblno
!
      use graphmod
      use GMP
!
      implicit none
!
      real(8) :: xcoord(3,8),dxdxi(3,2)
      real(8) :: xav(3),xavob(3),xi(2)
      character(3) :: text
!
      real(8) :: angle,height,xtext,ytext
      integer :: i,iv,ivar,inv,j,k,nc,ncol,nick,np,ntb,nt1
!
!  ...set colour
      ncol=npcol(2)
!
!  ...loop through triangular blocks
      do ntb = 1, NRPRISM
        nick=ntb*10+1
!
!  ...check if not invisible
        do inv=1,NRINVBL
          if(nick.eq.IGINV(inv)) goto 90
        enddo
!
!  .....bilinear blocks
        if (PRISMS(ntb)%Type.eq.'LiLiPrism') then

!
!  .....get the bases, next edges and then the vertices
          iv=0
          do i=1,2
            nt1 = abs(PRISMS(ntb)%FigNo(i)/10)
            do k=1,3
              nc = TRIANGLES(nt1)%EdgeNo(k)
              if (nc.gt.0) then
                np = CURVES(nc)% EndPoNo(1)
              else
                np = CURVES(-nc)% EndPoNo(2)
              endif
              iv=iv+1
              call pointr(np, xcoord(1,iv))
            enddo
          enddo
!
!  .....compute the average
          do i = 1,3
            xav(i) = 0.d0
            do j = 1,6
              xav(i) = xav(i) + xcoord(i,j)
            enddo
            xav(i) = xav(i)/6.d0
          enddo
        endif
!
!  .....spherical blocks
        if (PRISMS(ntb)%Type.eq. 'TriaShell') then
!
!  .....get the bases, next its central points
          xi(1) = 0.33333333333333d0
          xi(2) = 0.33333333333333d0
          do i=1,2
            nt1 = abs(PRISMS(ntb)%FigNo(i)/10)
            call trian(nt1,xi,xcoord(1,i),dxdxi)
          enddo
          do i = 1,3
            xav(i) = (xcoord(i,2)+xcoord(i,1))/2.d0
          enddo
        endif
!
!  .....transform to observer's system
        call trobs(xav,xavob)
!
!  .....rescale
        do ivar=1,2
          XY(ivar,1) = (xavob(ivar)-XCIM(ivar))/DIMIM*SIZE &
                       + XCWIN(ivar)
        enddo
!
!  .....and finally write the block nickname
        write(text,'(i3)') ntb*10+1
        xtext = xy(1,1)
        ytext = xy(2,1)
        height=1.
        angle=0.
        call symbol(xtext,ytext,height,text,angle,3,ncol)
   90 enddo
!
!  ...loop through rectangular blocks

      do ntb = 1, NRHEXAS
        nick=ntb*10+2
!
!  ...check if not invisible
        do inv=1,NRINVBL
          if(nick.eq.IGINV(inv)) goto 190
        enddo
!
!  .....get the bases, next edges and then the vertices
        iv=0
        do i=1,2
          nt1 = abs(HEXAS(ntb)%FigNo(i)/10)
          do k=1,4
            nc = RECTANGLES(nt1)%EdgeNo(k)
            if (nc.gt.0) then
              np = CURVES(nc)% EndPoNo(1)
            else
              np = CURVES(-nc)% EndPoNo(2)
            endif
            iv=iv+1
            call pointr(np, xcoord(1,iv))
          enddo
        enddo
!
!  .....compute the average
        do i = 1,3
          xav(i) = 0.d0
          do j = 1,8
            xav(i) = xav(i) + xcoord(i,j)
          enddo
          xav(i) = xav(i)/8.d0
        enddo
!
!  .....transform
        call trobs(xav,xavob)
!
!  .....rescale
        do ivar=1,2
          XY(ivar,1) = (xavob(ivar)-XCIM(ivar))/DIMIM*SIZE &
                       + XCWIN(ivar)
        enddo
!
!  .....write the block nickname
        write(text,'(i3)') ntb*10+2
        xtext = xy(1,1)
        ytext = xy(2,1)
        height=1.
        angle=0.
        call symbol(xtext,ytext,height,text,angle,3,ncol)
     190 enddo
!
!
   end subroutine dpblno

#endif

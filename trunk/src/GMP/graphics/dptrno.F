c----------------------------------------------------------------------
c
c   routine name       - dptrno
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - routine displays the nicknames of
c                        all 2D blocksc
c
c----------------------------------------------------------------------
c
      subroutine dptrno
c
      use graphmod
      use GMP
c
      implicit none
c
      real(8) :: xcoord(3,8), dxdxi(3,2)
      real(8) :: xav(3),xavob(3), xi(2)
      character(3) :: text
c
      real(8) :: angle,height,xtext,ytext
      integer :: i,inv,iv,ivar,j,k,nc,ncol,nick,np,nt1,ntb,ntv
c
c  ...set colour
      ncol=npcol(2)
c
c  ...loop through triangular blocks
      do 90 ntb = 1, NRPRISM
        nick=ntb*10+1
c
c  ...check if not invisible
        do 10 inv=1,NRINVBL
          if(nick.eq.IGINV(inv)) go to 90
   10     continue
c
c  .....bilinear blocks
        if (PRISMS(ntb)%Type.eq.'LiLiPrism') then
c
c  .....get the bases, next edges and then the vertices
          iv=0
          do 30 i=1,2
            nt1 = abs(PRISMS(ntb)%FigNo(i)/10)
            do 20 k=1,3
              nc = TRIANGLES(nt1)%EdgeNo(k)
              if (nc.gt.0) then
                np = CURVES(nc)% EndPoNo(1)
              else
                np = CURVES(-nc)% EndPoNo(2)
              endif
              iv=iv+1
              call pointr(np, xcoord(1,iv))
   20       continue
   30     continue
c
c  .....compute the average
          do 80 i = 1,3
            xav(i) = 0.d0
            do 70 j = 1,6
              xav(i) = xav(i) + xcoord(i,j)
  70        continue
            xav(i) = xav(i)/6.d0
  80      continue
        endif
c
c  .....spherical blocks
        if (PRISMS(ntb)%Type.eq. 'TriaShell') then
c
c  .....get the bases, next its central points
          xi(1) = 0.33333333333333d0
          xi(2) = 0.33333333333333d0
          do 230 i=1,2
            nt1 = abs(PRISMS(ntb)%FigNo(i)/10)
            call trian(nt1,xi,xcoord(1,i),dxdxi)
  230     continue
          do 280 i = 1,3
            xav(i) = (xcoord(i,2)+xcoord(i,1))/2.d0
  280     continue
        endif
c
c  .....transform to observer's system
        call trobs(xav,xavob)
c
c  .....rescale
        do 95 ivar=1,2
          XY(ivar,1) = (xavob(ivar)-XCIM(ivar))/DIMIM*SIZE
     .                 + XCWIN(ivar)
   95   continue
c
c  .....and finally write the block nickname
        write(text,'(i3)') ntb*10+1
        xtext = xy(1,1)
        ytext = xy(2,1)
        height=1.
        angle=0.
        call symbol(xtext,ytext,height,text,angle,3,ncol)
   90 continue
c
c  ...loop through rectangular blocks

      do 190 ntb = 1, NRHEXAS
        nick=ntb*10+2
c
c  ...check if not invisible
        do 110 inv=1,NRINVBL
          if(nick.eq.IGINV(inv)) go to 190
  110     continue
c
c  .....get the bases, next edges and then the vertices
        iv=0
        do 130 i=1,2
          nt1 = abs(HEXAS(ntb)%FigNo(i)/10)
          do 120 k=1,4
            nc = RECTANGLES(nt1)%EdgeNo(k)
            if (nc.gt.0) then
              np = CURVES(nc)% EndPoNo(1)
            else
              np = CURVES(-nc)% EndPoNo(2)
            endif
            iv=iv+1
            call pointr(np, xcoord(1,iv))
  120     continue
  130   continue
c
c  .....compute the average
        do 180 i = 1,3
          xav(i) = 0.d0
          do 170 j = 1,8
            xav(i) = xav(i) + xcoord(i,j)
  170      continue
          xav(i) = xav(i)/8.d0
  180    continue
c
c  .....transform
        call trobs(xav,xavob)
c
c  .....rescale
        do 195 ivar=1,2
          XY(ivar,1) = (xavob(ivar)-XCIM(ivar))/DIMIM*SIZE
     .                 + XCWIN(ivar)
  195     continue
c
c  .....write the block nickname
        write(text,'(i3)') ntb*10+2
        xtext = xy(1,1)
        ytext = xy(2,1)
        height=1.
        angle=0.
        call symbol(xtext,ytext,height,text,angle,3,ncol)
  190 continue
c
      end subroutine dptrno

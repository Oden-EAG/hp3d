c-----------------------------------------------------------------------
c
c   routine name       - dpvisid
c
c-----------------------------------------------------------------------
c
c   latest revision    - Jan 08
c
c   purpose            - routine displays triangles from the list RGTR
c                        in order stored in the file IGTR
c                        (DP-display,VI-visible,SID-sides)
c
c   arguments :
c     in:
c                Nf    = 10*(number of levels for contour plot)+
c                        flag whether to display elements numbers 
c
c-----------------------------------------------------------------------
c
      subroutine dpvisid(Nf)
c
      use graphmod
c
#include "syscom.blk"
c
      dimension xyn(2,3)
      character*5 text
c
      i3(j) = j-(j-1)/3*3
c   
      nick=0; nickno=0
c
c  ...decode flag
      call decode(Nf, numlev,nof)
c
c  ...decide whether to display edges
      if (numlev.gt.0) then
c        write(*,*) 'DISPLAY EDGES ?  (1/0)'
c        read(*,*) iedflag
         iedflag = 1
      endif
c
c  ...loop through all visible triangles
      do ilo=1,NRVISTR
c
c  .....get the triangle number
        itr=IGTR(ilo)
c
c  .....store triangle symbol
        itrtemp=IGTRCU(itr)
c
c  .....decode nickname to be displayed
        nick=IGTRCU(itr)/1000
        IGTRCU(itr)=IGTRCU(itr)-1000*nick
c
c  .....when displaying p-approximation colors
        if (numlev.eq.0) then
c
c  .......decode node number
          if (Nof.eq.2) nickno=IGTRNO(itr)
c
c  .......set color
          ncol=IGTRCU(itr)/10
          if ((ncol.lt.1).or.(ncol.gt.8)) then
            write(*,*) 'lsvisid: ncol = ',ncol
            call pause
          endif
          if (ncol.eq.0) then
            ncol=npcol(2)
          else
            ncol=npcol(2+ncol)
          endif
c
c  .......rescale
          do k=1,3
            do ivar=1,2
              XY(ivar,k) = RGTR((itr-1)*9+(k-1)*3+ivar)
              XY(ivar,k) = (XY(ivar,k)-XCIM(ivar))/DIMIM*SIZE 
     .                     + XCWIN(ivar)
              xyn(ivar,k) = XY(ivar,k)
            enddo
          enddo
c
c  .......fill the subtriangle
          call fillpoly(3,XY,ncol,ncol)
c
c  .......clear the buffer
          call flushx
c
c  .....when displaying solution values
        else
c
c  .......find entries in matrix of small triangles
          nfstr=IGTRNO(itr)/100
          nstr=IGTRNO(itr)-100*nfstr
c
c  .......in a loop over small triangles
          do i=1,nstr
c              
c  .........set color
            ncol=IGSTR(nfstr+i)
            if (ncol.ge.0) then
c               
c  ...........interpolate the color scale for the right color index
              aux = float(ncol)/float(numlev+1)
cldem 04.03.02
ccc              ncol=NPCOL(11+int(aux*(NR_COLORS-10)))
              ncol=NPCOL(11+int(aux*(NR_COLORS-11)))
            endif
c                    
c  .........rescale
            do k=1,3
              do ivar=1,2
                isub = (nfstr+i-1)*9+(k-1)*3+ivar
                XY(ivar,k) = RGTR(isub)
                XY(ivar,k) = (XY(ivar,k)-XCIM(ivar))/DIMIM*SIZE 
     .                       + XCWIN(ivar)
                if (i.eq.1) xyn(ivar,k) = XY(ivar,k)
              enddo
            enddo
c                      
            if (IDISPLAY_TYPE.eq.1) then
              ncol_line=npcol(2)
c               
c  ...........decode the edge number to draw
              nsid = iabs(ncol)
              nsid1 = i3(nsid+1)
ccc              write(*,*) 'dpvisid: ncol = ',ncol
              if (ncol.ne.0) 
     .        call drawline(XY(1,nsid),XY(2,nsid),
     .                      XY(1,nsid1),XY(2,nsid1),ncol_line)
            else
c                
c  ...........fill the subtriangle
              call fillpoly(3,XY,ncol,ncol)
            endif
c                    
c  .........clear the buffer
            call flushx
c              
c  .......end of loop through subtriangles
          enddo
        endif
c
c  .....set color
        ncol=npcol(2)
c
c  .....write block number when required
        if (Nof.ne.0) then
          shiftx=0.017d0*xlength
          shifty=0.004d0*ylength
          if (((nick.ne.0).and.(Nof.eq.1)).or.
     .       ((nickno.ne.0).and.(Nof.eq.2))) then
            if (Nof.eq.2) nick=nickno
            ncar=5
            write(text,'(i5)') nick
C             write(*,*) 'nick = ',nick
            xtext = (xyn(1,1)+xyn(1,2)+xyn(1,3))/3.d0-shiftx
            ytext = (xyn(2,1)+xyn(2,2)+xyn(2,3))/3.d0-shifty
            height=0.5
ccc            height=1.
            angle=0.
            call symbol(xtext,ytext,height,text,angle,ncar,ncol)
c
c  .........clear the buffer
            call flushx
c
          endif
        endif 
c
c  .....check whether subsequent edges are to be displayed
c       and draw them when needed
        if ((numlev.eq.0).or.(iedflag.eq.1)) then
          iflag = IGTRCU(itr)-10*(IGTRCU(itr)/10)
          iflag1 = iflag/4
          iflag2 = (iflag-4*iflag1)/2
          iflag3 = iflag-4*iflag1-2*iflag2
ccc          write(*,*) 'dpvisid: iflag1,iflag2,iflag3,ncol = ',
ccc     .                         iflag1,iflag2,iflag3,ncol
          if (iflag1.eq.1) 
     .       call drawline(xyn(1,1),xyn(2,1),xyn(1,2),xyn(2,2),ncol)
          if (iflag2.eq.1) 
     .       call drawline(xyn(1,3),xyn(2,3),xyn(1,2),xyn(2,2),ncol)
          if (iflag3.eq.1)
     .       call drawline(xyn(1,1),xyn(2,1),xyn(1,3),xyn(2,3),ncol)
c
c  .......clear the buffer
          call flushx
c
        endif
c
c  .....restore element symbol
        IGTRCU(itr)=itrtemp
c
c  ...end of loop through traingles
      enddo
c
c
      return
      end

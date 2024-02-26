#if HP3D_USE_X11

!-----------------------------------------------------------------------
!
!   routine name       - dpvisid
!
!-----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine displays triangles from the list RGTR
!                        in order stored in the file IGTR
!                        (DP-display,VI-visible,SID-sides)
!
!   arguments :
!     in:
!                Nf    = 10*(number of levels for contour plot)+
!                        flag whether to display elements numbers
!
!-----------------------------------------------------------------------
      subroutine dpvisid(Nf)
!
      use graphmod
!
      implicit none
!
      integer, intent(in) :: Nf
!
      real(8) :: xyn(2,3)
!
      real(8) :: angle,aux,height,shiftx,shifty,xtext,ytext
      integer :: nick,nickno,numlev,nof,iedflag,ilo,itr
      integer :: i,j,iflag,iflag1,iflag2,iflag3,isub,itrtemp
      integer :: ivar,k,ncar,ncol,ncol_line,nfstr,nsid,nsid1,nstr
!
      character(5) :: text
!
      nick=0; nickno=0
!
!  ...decode flag
      call decode(Nf, numlev,nof)
!
!  ...decide whether to display edges
      if (numlev.gt.0) then
!        write(*,*) 'DISPLAY EDGES ?  (1/0)'
!        read(*,*) iedflag
         iedflag = 1
      endif
!
!  ...loop through all visible triangles
      do ilo=1,NRVISTR
!
!  .....get the triangle number
        itr=IGTR(ilo)
!
!  .....store triangle symbol
        itrtemp=IGTRCU(itr)
!
!  .....decode nickname to be displayed
        nick=IGTRCU(itr)/1000
        IGTRCU(itr)=IGTRCU(itr)-1000*nick
!
!  .....when displaying p-approximation colors
        if (numlev.eq.0) then
!
!  .......decode node number
          if (Nof.eq.2) nickno=IGTRNO(itr)
!
!  .......set color
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
!
!  .......rescale
          do k=1,3
            do ivar=1,2
              XY(ivar,k) = RGTR((itr-1)*9+(k-1)*3+ivar)
              XY(ivar,k) = (XY(ivar,k)-XCIM(ivar))/DIMIM*SIZE &
                           + XCWIN(ivar)
              xyn(ivar,k) = XY(ivar,k)
            enddo
          enddo
!
!  .......fill the subtriangle
          call fillpoly(3,XY,ncol,ncol)
!
!  .......clear the buffer
          call flushx
!
!  .....when displaying solution values
        else
!
!  .......find entries in matrix of small triangles
          nfstr=IGTRNO(itr)/100
          nstr=IGTRNO(itr)-100*nfstr
!
!  .......in a loop over small triangles
          do i=1,nstr
!
!  .........set color
            ncol=IGSTR(nfstr+i)
            if (ncol.ge.0) then
!
!  ...........interpolate the color scale for the right color index
              aux = float(ncol)/float(numlev+1)
! ldem 04.03.02
!!!              ncol=NPCOL(11+int(aux*(NR_COLORS-10)))
              ncol=NPCOL(11+int(aux*(NR_COLORS-11)))
            endif
!
!  .........rescale
            do k=1,3
              do ivar=1,2
                isub = (nfstr+i-1)*9+(k-1)*3+ivar
                XY(ivar,k) = RGTR(isub)
                XY(ivar,k) = (XY(ivar,k)-XCIM(ivar))/DIMIM*SIZE &
                             + XCWIN(ivar)
                if (i.eq.1) xyn(ivar,k) = XY(ivar,k)
              enddo
            enddo
!
            if (IDISPLAY_TYPE.eq.1) then
              ncol_line=npcol(2)
!
!  ...........decode the edge number to draw
              nsid = abs(ncol)
              j = nsid+1
              nsid1 = j-(j-1)/3*3
!!!              write(*,*) 'dpvisid: ncol = ',ncol
              if (ncol.ne.0) &
              call drawline(XY(1,nsid),XY(2,nsid), &
                            XY(1,nsid1),XY(2,nsid1),ncol_line)
            else
!
!  ...........fill the subtriangle
              call fillpoly(3,XY,ncol,ncol)
            endif
!
!  .........clear the buffer
            call flushx
!
!  .......end of loop through subtriangles
          enddo
        endif
!
!  .....set color
        ncol=npcol(2)
!
!  .....write block number when required
        if (Nof.ne.0) then
          shiftx=0.017d0*xlength
          shifty=0.004d0*ylength
          if (((nick.ne.0).and.(Nof.eq.1)).or. &
             ((nickno.ne.0).and.(Nof.eq.2))) then
            if (Nof.eq.2) nick=nickno
            ncar=5
            write(text,'(i5)') nick
!             write(*,*) 'nick = ',nick
            xtext = (xyn(1,1)+xyn(1,2)+xyn(1,3))/3.d0-shiftx
            ytext = (xyn(2,1)+xyn(2,2)+xyn(2,3))/3.d0-shifty
            height=0.5
!!!            height=1.
            angle=0.
            call symbol(xtext,ytext,height,text,angle,ncar,ncol)
!
!  .........clear the buffer
            call flushx
!
          endif
        endif
!
!  .....check whether subsequent edges are to be displayed
!       and draw them when needed
        if ((numlev.eq.0).or.(iedflag.eq.1)) then
          iflag = IGTRCU(itr)-10*(IGTRCU(itr)/10)
          iflag1 = iflag/4
          iflag2 = (iflag-4*iflag1)/2
          iflag3 = iflag-4*iflag1-2*iflag2
!!!          write(*,*) 'dpvisid: iflag1,iflag2,iflag3,ncol = ', &
!!!                               iflag1,iflag2,iflag3,ncol
          if (iflag1.eq.1) &
             call drawline(xyn(1,1),xyn(2,1),xyn(1,2),xyn(2,2),ncol)
          if (iflag2.eq.1) &
             call drawline(xyn(1,3),xyn(2,3),xyn(1,2),xyn(2,2),ncol)
          if (iflag3.eq.1) &
             call drawline(xyn(1,1),xyn(2,1),xyn(1,3),xyn(2,3),ncol)
!
!  .......clear the buffer
          call flushx
!
        endif
!
!  .....restore element symbol
        IGTRCU(itr)=itrtemp
!
!  ...end of loop through traingles
      enddo
!
!
      end subroutine dpvisid

#endif

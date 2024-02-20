#if HP3D_USE_X11

C******************************************************
C
C THE IMPORTANT ROUTINES FOR THE POSTSCRIPT FILES
C
C wvopps - to open a ps file
C wvpscl(clip(4)) - to set clipping region
C wvpsdl(X1,Y1,X2,Y2,IC) - to draw a line
C wvpsfp(N,Xc,Icf,Icl) - to fill a polygon
C wvpspx(X,Y,IC) - to draw a pixel
C wvpste(X1,Y1,H,MESS,A,N,ICOL) - to write text
C wvpsout(x,y,icode) - to write ps parameters into buffer
C flushps - to rewrite contents os the buffer into output ps file
C wvclps - to close the ps file
C
C OTHER PARAMETERS ARE DISCUSSED IN THE ROUTINES
C
C********************************************************

        subroutine wvopps
C---------------------------------------------------------
C      AUTHOR(S)        - Satish Chavva,Tadeusz Liszka and others
C
C      version  -   1998
C
C    purpose - to write header of ps file
C
C========================================================
cwr07.08.99
cwr03.14.02
c
      implicit real (a-h,o-y)
#include "wvglob.blk"
#include "wvpost.blk"
#include "wvscrn.blk"
c
      character(10) name
c
  110 write(*,*)  'Give name for Postscript file (0 - no postscript)'
   10 read(*,1011) name

 1011 format(a10)
      if (name(1:1).le.' ') goto 10
      if (name(1:1).eq.'0') then
          iflag=0
      else
          ifile=79
          open(ifile,file=name,status='new',err=100)
          write(*,*)' PostScript send to ',name
          write(ifile,1000) name
c        Check if color Postscript required

          write(*,*) 'B/W or  Color Postscript (0 - B/W  1-Color)'
          read(*,*) icolval

 1000     format('%!  TL PostScript file, name ',a10/
     &        '%%BoundingBox: 36 36 518 756',/
     &        '0 setlinewidth',/
     &        '/n {newpath moveto} def',/
     &        '/c {lineto closepath} def',/
     &        '/m {moveto} def',/
     &        '/d {lineto} def',/
     &        '/f {setgray fill} def',/
     &        '/s {0 setgray stroke} def',/
     &        '/rf {setrgbcolor fill} def',/
     &        '/rs {setrgbcolor stroke} def',/
     &        '/gs {gsave} def',/
     &        '/gr {grestore} def',/
     &        '/fo { /Times-Roman findfont } def ',
     &        '/sc { scalefont setfont } def'/
     &        '50 750 translate -90 rotate fo 12 sc newpath')
          ilbuf=0
      endif
      ilepsp=1
c      write(*,*) 'NO POSTSCRIPT GENERATED, iflag=',iflag
      return
  100 continue
      write(*,*) ' File ',name,' exists !!'
      write(*,*) ' I will not overwrite it !!'
      goto 110
      end


c-----------------------------------------------------------------
      subroutine wvpscl(Clip)
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    -
c
c   purpose            - to set clipping region for the ps page
c
c
c----------------------------------------------------------------------
c
      dimension clip(4)
c
#include "wvglob.blk"
#include "wvpost.blk"
#include "wvscrn.blk"
c
      iprint=0
      if (iprint.eq.1) then
        write(*,*) 'WVPSCL: Clip = ',Clip
      endif
c
      if (IFILE.eq.0) return
c        call flushps
      call wvpsout(dum,dum,9)
      call wvpsout(Clip(1),Clip(2),1)
      call wvpsout(Clip(1),Clip(4),2)
      call wvpsout(Clip(3),Clip(4),2)
      call wvpsout(Clip(3),Clip(2),4)
      call wvpsout(dum,dum,7)
      if (iprint.eq.1) then
        write(*,*) 'WVPSCL: call pause'
        call pause
      endif
c
      return
      end

C========================================================
        SUBROUTINE WVPSDL(X1,Y1,X2,Y2,IC)
C========================================================
C
C  author - Satish Chavva and many otheres
C
C  version - 1998
C
C  purpose - to draw a line
C
C  arguments - X1,Y1,X2,Y2-coordinates of the two endpoints
C              IC - line color
C---------------------------------------------------------
c
#include "wvglob.blk"
#include "wvpost.blk"
#include "wvscrn.blk"
c
        iprint=0
c
        if (ifile.eq.0) return
        call wvpsout(x1,y1,1)
        call wvpsout(x2,y2,4)
        call wvpsout(dum,dum,6)
        RETURN
        END
      subroutine wvpsfp(N,Xc,Icf,Icl)
c
c----------------------------------------------------------------------
c
c   routine name       - wvpsfp
c
c----------------------------------------------------------------------
c   author  - Satish Chavva and many others
c
c   computer           - machine independent
c
c   latest revision    -  98
c
c   purpose            - routine fills a polygon with a color
c                        in postcript
c                        It is assumed that the first and last
c                        vertex point coincide with each other !!
c
c   usage              - call wvpsfp(N,Xc,Icf,Icl)
c
c   arguments :
c     in:
c               N      - number of vertices
c               Coor   - coordinates of the vertices
c               Icf    - field color number
c               Icl    - line color number
c
c   required  routines - wvpsfp,wvwarn
c
c----------------------------------------------------------------------
c
#include "wvglob.blk"
#include "wvpost.blk"
#include "wvscrn.blk"
c
      real Xc(2,N)
c
      if (IFILE.eq.0) return
      if (Icf.ge.0) then
c
c  .....begin a newpath with the first point coordinates
        call wvpsout(Xc(1,1),Xc(2,1),1)
c
c  .....loop through the consequtive vertices of the polygon
        do 10 i=2,N-1
c
c  .......line to the vertex point
          call wvpsout(Xc(1,i),Xc(2,i),2)
   10   continue
c
c  .....line to the last vertex point and close the path
        call wvpsout(Xc(1,n),Xc(2,n),4)
c
c  .....fill up the polygon with a grey pattern
        call wvpsout(dum,float(icf),5)
      endif
c
c  ...draw the boundary line if a different color
      if (Icl.ne.Icf) then
        call wvpsout(Xc(1,1),Xc(2,1),1)
        do 20 i=2,n-1
          call wvpsout(Xc(1,i),Xc(2,i),2)
   20   continue
        call wvpsout(Xc(1,n),Xc(2,n),4)
        call wvpsout(dum,dum,6)
      endif
c
c
      return
      end

C========================================================
c     draw pixel in postcript
C========================================================


        SUBROUTINE WVPSPX(X,Y,IC)
c
#include "wvglob.blk"
#include "wvpost.blk"
#include "wvscrn.blk"
c
        iprint=0

        call wvpsout(x,y,1)
        call wvpsout(x,y,4)
        call wvpsout(dum,dum,6)
        RETURN
        END


C========================================================
c     draw text in postscript
C========================================================

        SUBROUTINE WVPSTE(X1,Y1,H,MESS,A,N,ICOL)
        integer H,A
c        CHARACTER MESS*80
        CHARACTER MESS*(*)
c
#include "wvglob.blk"
#include "wvpost.blk"
#include "wvscrn.blk"
c
        if (ifile.eq.0) return
        call wvpsout(x1,y1,3)
        k=12
        if (h.ge.16) k=15
        call wvpsout(float(k),dum,8)
        call flushps
        write(ifile,1000) ' (',mess(1:n),') show'
 1000   format(a)
        RETURN
        END


C-------------------------------------------------------------------

      subroutine wvpsout(x,y,Icode)
c
c=====================================================================
c     author - Satish Chavva and many others
c
c     version - 1998
c
c=================================== post script dump to file ========
c       send real x,y (in inches) to postcript file
c       icode values 1:n 2:d 3:m 4:c 5:f 6:s 7:clip 8:font, 9:initclip
c       explanation of codes as in flushps:
c          1 - n   x,y,newpath,moveto
c          2 - d   x,y,lineto
c          3 - m   x,y,moveto
c          4 - c   x,y,lineto,closepath
c          5 - f   y,setgray,fill
c          6 - s   0,setgray,stroke
c          7 - clip
c          8 - font  == kx,font
c          9 - initclip
c       text printing is inside wvpstx routine
c       initialization - in wvopps
c
c=====================================================================
c
#include "wvglob.blk"
#include "wvpost.blk"
#include "wvscrn.blk"
c
      data ymax /0.d0/
      save ymax
c
      iprint=2

c      icolval=1
c
      if (ifile.eq.0) return
      ilbuf=ilbuf+1
c
c     size is in inches, output in points (1/72")
      scale=72.0

      if (icode.eq.5) then
c         set grayscale colors
         if (icolval.eq.0) then
          color = 1.0 - y / 9.0
          if (color .lt.0.0) color = 1.0 - (y-16)/50.0
          if ((iprint.eq.2).and.(y.gt.ymax)) then
            ymax=y
            write(*,*) 'WVPSOUT: y,color = ',y,color
            call pause
          endif
          bufor(1,ilbuf)=color
         elseif (icolval.eq.1) then
          bufor(1,ilbuf)=y
         endif
      elseif (icode.eq.8) then
          bufor(1,ilbuf)=x
      elseif (icode.le.4) then
          bufor(1,ilbuf)=x*scale
          bufor(2,ilbuf)=y*scale
c      make sure they are in correct limits
          do 20 i=1,2
              bufor(i,ilbuf)=max(-999.0,min(9990.0,bufor(i,ilbuf)))
   20     continue
      endif
c
      if (iprint.eq.1) then
        write(*,*) 'WVPSOUT: ilbuf,BUFOR(*,ilbuf) = '
        write(*,*) ilbuf,(BUFOR(ii,ilbuf),ii=1,2)
      endif
c
      ibufc(ilbuf)=icode
      if (ilbuf.eq.5) call flushps
      return
      end subroutine
c
c-----------------------------------------------------------------------
c
      subroutine flushps
c
c-----------------------------------------------------------------------
c
c   routine name       - flushps
c
c-----------------------------------------------------------------------
c   author  - Satish Chavva and many others
c
c   computer           - machine independent
c
c   latest revision    - 1998
c
c   purpose            - routine sends poscript buffer as one line to
c                        output
c
c   usage              - call flushps
c
c   arguments          - none
c
c   required  routines -
c
c-----------------------------------------------------------------------
c
      use graphmod
c
#include "wvglob.blk"
#include "wvpost.blk"
#include "wvscrn.blk"
c
      character line*100,text*20,code(4)*1,lineout*100
c
      data code/'n','d','m','c'/
c
      iprint=0

       vmaxval=float(ired(1))

      if (IFILE.eq.0) then
        ILBUF = 0
        return
      endif
c
      if (ILBUF.ne.0) then
        line=' '
        linelen=0
        do 100 i=1, ILBUF

          if (IBUFC(i).ne.0) then

            if (IBUFC(i).lt.5) then

              if (iprint.eq.1) then
         write(*,*) 'FLUSHPS: BUFOR(1,i),BUFOR(2,i),code(IBUFC(i)) = '
                write(*,1010) BUFOR(1,i),BUFOR(2,i),code(IBUFC(i))
              endif
              write(ifile,1010)
     .        BUFOR(1,i),BUFOR(2,i),code(IBUFC(i))

 1010         format(2f8.1,1x,1a)
              lentext=18
            endif
            if (IBUFC(i).eq.5) then
             if (icolval.eq.1) then
              i7=int(BUFOR(1,i))+1
              ir1=ired(i7)
              ig1=igreen(i7)
              ib1=iblue(i7)
c    here we make use of the color values and scale them to values
c    0 and 1
              r1=float(ir1)/vmaxval
              g1=float(ig1)/vmaxval
              b1=float(ib1)/vmaxval
             endif
            endif
            if (IBUFC(i).eq.5) then
             if (icolval.eq.0) then
              write(ifile,1020) BUFOR(1,i)

 1020         format(f5.2,' f')
              lentext=7
             elseif (icolval.eq.1) then
              write(ifile,1024) r1,g1,b1
 1024         format(3f6.3,1x,' rf')
              lentext=19
             endif
            elseif (IBUFC(i).eq.6) then

              write(ifile,*)' s'
              lentext=2

            elseif (IBUFC(i).eq.7) then
              write(ifile,*)' clip '
              lentext=5
            elseif (IBUFC(i).eq.8) then
              k=BUFOR(1,i)
              write(ifile,1030) k
c              write(*,1030) k
c              write(*,*)"3"
c              read(*,*)

 1030         format(' fo',i3,' sc')
              lentext=9
            elseif (IBUFC(i).eq.9) then
              write(ifile,1040)
c              write(*,1040)
c              read(*,*)

 1040         format(' initclip')
              lentext=9
            else
              lentext=0
            endif
c
            if (lentext.gt.0) then
              line(linelen+1:linelen+lentext)=text(1:lentext)
              linelen=linelen+lentext
            endif
           endif
  100     continue
c
c.......here one may reduce size of output by:
c.......skipping repeated blanks
c.......skipping .0 in floating point numbers
        ifrom=1
        ito=2
        lineout=' '
c
  200   continue
        if (line(ifrom:ifrom+1).eq.'  ') then
          ifrom = ifrom+1
cc        elseif (line(ifrom:ifrom+2).eq.'.0 ') then
cc          ifrom = ifrom +2
cc        elseif (line(ifrom:ifrom+2).eq.'.00') then
cc          ifrom = ifrom +3
        elseif (line(ifrom:ifrom+1).eq.'. ') then
          ifrom = ifrom +1
        else
          lineout(ito:ito)=line(ifrom:ifrom)
          ito = ito +1
          ifrom = ifrom+1
        endif
c
        if (ifrom.lt.linelen) goto 200
        lineout(ito:ito)=line(linelen:linelen)
c
c.......initial character in line is always blank
        if(ito.gt.100)then
c           write(*,*)"Hi Cam"
c           read(*,*)
        endif
c        write(IFILE,1000) lineout(3:ito)
c        write(*,1000) lineout(3:ito)
c              write(*,*)"4"
c        read(*,*)

 1000   format(a)
        ILBUF=0
      endif
c
      return
      end
c========================================================
        subroutine wvclps
c
#include "wvglob.blk"
#include "wvpost.blk"
#include "wvscrn.blk"
c
        ilepsp=0
        do i=1,10
          igpsp(i)=0
        enddo
        if (ifile.eq.0) return
        call flushps
        write(ifile,1000)
 1000   format(' stroke showpage'/)
        close(ifile)
        write(*,*)' PostScript file closed'
        return
        end subroutine

#endif

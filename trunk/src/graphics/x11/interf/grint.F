#if HP3D_USE_X11

c  GRAPHICAL INTERFACE FOR X-WINDOWS AND POSTCRIPT FILES
c
c  THE AVALIBLE BASIC ROUTINES ARE AS FOLLOWS
c
c  CLOSWIND(Iwin)- to end pictures and set the (ps off)
c  OPENWIND(Iwtype,Iwsize,Iwin)- initialzing the graphics
c  SELWIND(Iwin) - to begin pictures( gets the dimensions of the object
c      to be drawn and returns the dimensions of the image space; also
c      rescales and clears the window)
c  DRAWLINE(X1,Y1,X2,Y2,Ilcol) - to draw a line
c  DRAWPOLY(N,Coor,Ilcol) - to draw a polygon
c  FILLPOLY(Nin,Coor,Ifcol,Ilcol) - to draw a filled polygon
c  SYMBOL(X1,Y1,Height,Mess,Angle,N,Icol)- to write text
c
c  AUXILLARY ROUTINES
c  SETPOST(Ionoff)- to set 1-on/0-off ps echo
c  CLRWIND
c
c  PARAMETERS
c    Iwtype - window type
c    Ixsize,Iysize - x,y window dimensions(in pixels)
c    Iback,Ilcol,Ifcol - color numbers (Ifcol to fill)
c    Mess - charcter array of length N
c    Height,Angle - parameters of text(for ps only!)
c    X1,Y1 (or array Coor) -x,y coordinates of point(s)
c
c    Further Information in the Routines
c
c  Colors are as Follows
c        1-white 2-black
c        3-10 - sparse spectrum of colors
c        11-100 - dense spectrum of colors
c
c
c
c
c**************************************************************

      subroutine closwind(Iwin)
c
c-----------------------------------------------------------------------c
c   routine name       - closwind
c
c-----------------------------------------------------------------------c
c   author - Satish C. and many others
c
c   computer           - machine independent
c
c   latest revision    - 1998
c
c   purpose            - routine closes a window
c
c   usage              - call closwind(Iwin)
c
c   arguments :
c     in:
c          Iwin        - number of the window to close
c
c   required  routines - selwind,setpost,wverr,wvcldf,wvclps,wvtext,
c                        wvwarn
c
c----------------------------------------------------------------------
c
#include "wvglob.blk"
c
      if(IOPEN.ne.1) call wverr(2)
c
      if (Iwin.lt.0.or.Iwin.gt.9) then
        call wvwarn(1)
      else
        ic=iwin+1
        if (LWTYPE(ic).lt.0) return
c
c.......check for postscript
        call setpost(0)
        if (LWTYPE(ic).eq.2) then
          ILEPSP=ILEPSP-1
          if (ILEPSP.eq.0) call wvclps
        endif
c
c.......check for graph files
        if (LWTYPE(ic).eq.1) then
          ILEPF=ILEPF-1
c          if (ILEPF.eq.0) call wvcldf
        endif
        if (LWTYPE(ic).eq.0) then
          ILESCR=ILESCR-1
          if (ILESCR.eq.0) call wvtext
        endif
        LWTYPE(ic)=-1
        if (Iwin+1.eq.ic) then
cccc        if (iwind+1.eq.ic) then
c
c.........switch to any active window
          do 10 i=1,10
            if (LWTYPE(i).ge.0) then
              k=i-1
              call selwind(k)
              return
            endif
   10     continue
        endif
      endif
c
      return
      end
c----------------------------------------------------------------------
c
c   routine name       - openwind
c
c----------------------------------------------------------------------
c   author - Satish C. and many others
c
c   computer           - machine independent
c
c   latest revision    - 1998
c
c   purpose            - routine opens a graphic window
c
c   usage              - call openwind(Iwtype,Iwsize,Iwin)
c
c   arguments :
c   in:
c       Iwtype = 0     - screen
c                1     - disk file (not used)
c                2     - poscript file
c       Iwsize(4)      - Xmin,Ymin,Xmax,Ymax
c       Iwin           - number of the window = 0,1,...,9
c
c   out:
c       Iwsize         - adjusted size of the window
c
c   required routines  - closewind,newscale,scaleon,selwind,setpost,
c                        wverr,wvgrph,wvopdf,wvopps
c
c---------------------------------------------------------------------
c
      subroutine openwind(Iwtype,Iwsize,Iwin)
c
      integer(4) Iwtype,Iwsize(4),Iwin
c
#include "wvglob.blk"
#include "wvscrn.blk"
c
c     dumpfile maxx and maxy
c
      real(8) range(4)
c
      data range /0d0,0d0,1d0,1d0/
c
      iprint=0
      if (iprint.eq.1) then
        write(*,*) 'openwind: IOPEN = ',IOPEN
      endif
c
      if (IOPEN.ne.1.or.Iwtype.le.-1) then
c
c.......initialize window system
        do 10 i=1,10
          IGPSP(i)=0
          LWTYPE(i)=-1
   10   continue
        IOPEN=1
        ILEPSP=0
        ILESCR=0
        ILEDF=0
        MAXCOL=0
        if (Iwtype.lt.0) return
      endif
c
      if (Iwin.lt.0.or.Iwin.gt.9) call wverr(1)
      kwin=Iwin+1
c
      if (iprint.eq.1) then
        write(*,*) 'openwind: kwin,LWTYPE(kwin) = ',kwin,LWTYPE(kwin)
      endif
c
      if (LWTYPE(kwin).ge.0) then
c
c.......clear if window in use
        call closwind(Iwin)
      endif
c
      do i=1,4
        LWSIZE(kwin,i)=max0(0,Iwsize(i))
      enddo
c
      if (Iwtype.eq.2) then
        ILEPSP=ILEPSP+1
        if (Iwsize(3).eq.0) then
          Iwsize(3)=10
          Iwsize(4)=8
        endif
        LWSIZE(kwin,3)=min(Iwsize(3),10)
        LWSIZE(kwin,4)=min(Iwsize(4),8)
        if (ILEPSP.eq.1) call wvopps
      endif
c
ccc        write(*,*)'done with calling the file'
      if (Iwtype.eq.1) then
        ILEDF=ILEDF+1
c        ilepf=ilepf+1
c        if (ILEDF.eq.1) call wvopdf
c        if (ilepf .eq.1) call wvopdf
        if (Iwsize(3).eq.0) then
          Iwsize(3)=maxdfx
          Iwsize(4)=maxdfy
        endif
        LWSIZE(kwin,3)=min(Iwsize(3),MAXDFX)
        LWSIZE(kwin,4)=min(Iwsize(4),MAXDFY)
      endif
c
      if (Iwtype.eq.0) then
        ILESCR=ILESCR+1
c
        if (iprint.eq.1) then
          write(*,*) 'openwind: ILESCR = ',ILESCR
        endif
c
        if (ILESCR.eq.1) call wvgrph
ccc        write(*,*) 'openwind: AFTER CALL TO wvgrph'
        if (Iwsize(3).eq.0) then
          Iwsize(3)=MAXSCX
          Iwsize(4)=MAXSCY
        endif
ccc        write(*,*)'in openwind  MAXSCX,MAXSCY = ', MAXSCX,MAXSCY
        LWSIZE(kwin,3)=min(Iwsize(3),MAXSCX)
        LWSIZE(kwin,4)=min(Iwsize(4),MAXSCY)
      endif
c
c.....adjust Iwsize on output
      do 45 i=1,4
        Iwsize(i)=LWSIZE(kwin,i)
   45 continue
      LWTYPE(kwin)=Iwtype
      if (iprint.eq.1) then
        write(*,*)'in openwind after label 45 kwin,LWTYPE(kwin) = ',
     .             kwin,LWTYPE(kwin)
      endif
c
      call selwind(Iwin)
c
c.....set default scale
      call newscale(range)
      call scaleon(0)
      if (Iwtype.eq.0) then
        open(99,file='LASER',status='old',err=123)
        close(99)
        call setpost(1)
      endif
c
  123 continue
c
      return
      end
      subroutine selwind(Iwin)
c
c-----------------------------------------------------------------------
c
c   routine name       - selwind
c
c-----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - dec 90
c
c   purpose            - routine selects a graphic window to use
c
c   usage              - call selwind(Iwin)
c
c   arguments :
c     in:
c             Iwin     - number of the window
c
c   required  routines - wvwarn,wverr,wvpscl,wvsccl
c
c-----------------------------------------------------------------------
c
#include "wvglob.blk"
#include "wvscrn.blk"
c
      dimension iclip(4),clip(4)
c
      iprint=0
c
      if (IOPEN.ne.1) call wverr(2)
      if (Iwin.lt.0.or.Iwin.gt.9) then
        call wvwarn(1)
      else
        ICURRWIN = Iwin + 1
        do 10 i=1,4
          iclip(i) = LWSIZE(ICURRWIN,i)
   10   continue
        if (LWTYPE(ICURRWIN).eq.0) call wvsccl(iclip)
c        if (LWTYPE(ICURRWIN).eq.1) call wvdfcl(iclip)
        if (LWTYPE(ICURRWIN).eq.2) then
          do 20 i=1,4
            clip(i)=LWSIZE(ICURRWIN,i)
   20     continue
          if (iprint.eq.1) then
            write(*,*) 'SELWIND: clip = ',clip
          endif
          call wvpscl(clip)
        endif
        if (IGPSP(ICURRWIN).eq.1) then
          do 30 i=1,4
            clip(i)=iclip(i)*10.0/MAXSCX
   30     continue
          if (iprint.eq.1) then
            write(*,*) 'SELWIND: clip = ',clip
          endif
          call wvpscl(clip)
        endif
      endif
c
      return
      end
      subroutine  drawline(X1,Y1,X2,Y2,Ilcol)
c
c---------------------------------------------------------------------
c
c   routine name       - drawline
c
c---------------------------------------------------------------------
c   author     - Satish C. and many others
c
c   computer           - machine independent
c
c   latest revision    - 1998
c
c   purpose            - routine draws a line
c
c   usage              - call drawline(X1,Y1,X2,Y2,Ilcol)
c
c   arguments :
c     in:
c          X1,Y1,X2,Y2 - coordinates of two points
c          Ilcol       - color of the line
c     out:
c
c   required  routines -
c
c---------------------------------------------------------------------
c
      real(8) X1,X2,Y1,Y2
      integer sx1,sx2,sy1,sy2
c
#include "wvglob.blk"
#include "wvscrn.blk"
c
      il = Ilcol
c
      if (LWTYPE(ICURRWIN).le.1) then
          if (LSCALE(ICURRWIN).eq.0) then
              sx1=X1
              sx2=X2
              sy2=Y2
              sy1=Y1
          else
              sx1=xgscal(ICURRWIN,1)*x1+xgscal(ICURRWIN,3)
              sx2=xgscal(ICURRWIN,1)*x2+xgscal(ICURRWIN,3)
              sy1=xgscal(ICURRWIN,2)*y1+xgscal(ICURRWIN,4)
              sy2=xgscal(ICURRWIN,2)*y2+xgscal(ICURRWIN,4)
          endif
          if (LWTYPE(ICURRWIN).eq.0) call wvscdl(sx1,sy1,sx2,sy2,il)
c          if (LWTYPE(ICURRWIN).eq.1) call wvdfdl(sx1,sy1,sx2,sy2,il)
      endif
      if (LWTYPE(ICURRWIN).eq.2) then
          if (LSCALE(ICURRWIN).eq.0) then
              xx1=X1
              xx2=X2
              xy2=Y2
              xy1=Y1
              call wvpsdl(xx1,xy1,xx2,xy2,il)
          else
              xx1=xgscal(ICURRWIN,1)*x1+xgscal(ICURRWIN,3)
              xx2=xgscal(ICURRWIN,1)*x2+xgscal(ICURRWIN,3)
              xy1=xgscal(ICURRWIN,2)*y1+xgscal(ICURRWIN,4)
              xy2=xgscal(ICURRWIN,2)*y2+xgscal(ICURRWIN,4)
              call wvpsdl(xx1,xy1,xx2,xy2,il)
          endif
      endif
      if (IGPSP(ICURRWIN).eq.1) then
        if (LSCALE(ICURRWIN).eq.0) then
          xx1=x1*10.0/MAXSCX
          xx2=x2*10.0/MAXSCX
          xy2=y2*10.0/MAXSCX
          xy1=y1*10.0/MAXSCX
          call wvpsdl(xx1,xy1,xx2,xy2,il)
        else
          xx1=(XGSCAL(ICURRWIN,1)*X1+xgscal(ICURRWIN,3))*10.0/MAXSCX
          xx2=(XGSCAL(ICURRWIN,1)*X2+xgscal(ICURRWIN,3))*10.0/MAXSCX
          xy1=(XGSCAL(ICURRWIN,2)*Y1+xgscal(ICURRWIN,4))*10.0/MAXSCX
          xy2=(XGSCAL(ICURRWIN,2)*Y2+xgscal(ICURRWIN,4))*10.0/MAXSCX
          call wvpsdl(xx1,xy1,xx2,xy2,il)
        endif
      endif
c
      return
      end

      subroutine drawpoly(N,Coor,Ilcol)
c
c-----------------------------------------------------------------------c
c   routine name       - drawpoly
c
c-----------------------------------------------------------------------c
c   author  - Satish C. and many others
c
c   computer           - machine independent
c
c   latest revision    - 1998
c
c   purpose            - routine draws a polygon
c
c   usage              - call drawpoly(N,Coor,Ilcol)
c
c   arguments :
c     in:
c               N      - number of the vertices
c               Coor   - coordinates of the vertices
c               Ilcol  - line color number
c
c   required  routines - fillpoly
c
c-----------------------------------------------------------------------c
      real(8) Coor(2,N)
c
      call fillpoly(N,Coor,-1,Ilcol)
c
      return
      end

      subroutine fillpoly(Nin,Coor,Ifcol,Ilcol)
c
c-----------------------------------------------------------------------
c
c   routine name       - fillpoly
c
c-----------------------------------------------------------------------
c   author   - Satish C. and many others
c
c   computer           - machine independent
c
c   latest revision    - 1998
c
c   purpose            - routine fills a polygon with a color
c
c   usage              - call fillpoly(Nin,Coor,Ifcol,Ilcol)
c
c   arguments :
c     in:
c               Nin    - number of the vertices
c               Coor   - coordinates of the vertices
c               Ifcol  - field color number
c               Ilcol  - line color number
c
c   required  routines - wvpsfp,wvwarn
c
c-----------------------------------------------------------------------
c
      real(8) Coor(2,Nin)
c
      integer scoor(2,65)
      real xcoor(2,65)
c
#include "wvscrn.blk"
#include "wvglob.blk"
c
      if (Nin.gt.65) call wvwarn(2)
      n=min(65,Nin)
      if=Ifcol
      il=Ilcol
c
      if (LWTYPE(ICURRWIN).le.1) then
        if (LSCALE(ICURRWIN).eq.0) then
          do 17 i=1,n
            scoor(1,i)=coor(1,i)
            scoor(2,i)=coor(2,i)
   17     continue
        else
          do 20 i=1,n
            scoor(1,i)=XGSCAL(ICURRWIN,1)*Coor(1,i)+
     .                 XGSCAL(ICURRWIN,3)
            scoor(2,i)=XGSCAL(ICURRWIN,2)*Coor(2,i)+
     .                 XGSCAL(ICURRWIN,4)
   20     continue
        endif
        if (LWTYPE(ICURRWIN).eq.0) call wvscfp(n,scoor,il,if)
c        if (LWTYPE(ICURRWIN).eq.1) call wvdffp(n,scoor,il,if)
      endif
c
      if (LWTYPE(ICURRWIN).eq.2.or.igpsp(ICURRWIN).eq.1) then
        xmult = 1.0
        if (IGPSP(ICURRWIN).eq.1) xmult = 10.0/maxscx
        if (LSCALE(ICURRWIN).eq.0) then
          do 30 i=1,n
            xcoor(1,i)=Coor(1,i)*xmult
            xcoor(2,i)=Coor(2,i)*xmult
   30     continue
        else
          do 40 i=1,n
            xcoor(1,i)=(XGSCAL(ICURRWIN,1)*Coor(1,i)+
     .                 XGSCAL(ICURRWIN,3))*xmult
            xcoor(2,i)=(XGSCAL(ICURRWIN,2)*Coor(2,i)+
     .                 XGSCAL(ICURRWIN,4))*xmult
   40     continue
        endif
c..ldem        call wvpsfp(n,xcoor,il,if)
        call wvpsfp(n,xcoor,if,il)
      endif
c
      return
      end

      subroutine symbol(X1,Y1,Height,Mess,Angle,N,Icol)
c
c-----------------------------------------------------------------------c
c   routine name       - symbol
c
c-----------------------------------------------------------------------c
c   computer           - machine independent
c
c   latest revision    - dec 90
c
c   purpose            - routine draws a text or char marker
c
c   usage              - call symbol(X1,Y1,Height,Mess,Angle,N,Icol)
c
c   arguments :
c     in:
c            X1,Y1     - position to draw
c            Height    - height of the symbol
c            Angle     - angle of drawing
c            Mess      - character(s) to draw
c            N         - number of characters in Mess
c            Icol      - color to use
c     out:
c
c   required  routines -
c
c-----------------------------------------------------------------------c
      real(8) X1,Y1,Height,Angle
      integer h,a,sx1,sy1
      character Mess*(*)
c
#include "wvglob.blk"
#include "wvscrn.blk"
c
c      write(*,*) 'in symbol'
c      write(*,*) 'Mess',Mess
      nn = N
      ic = Icol
      h = Height*16
      a = Angle
c
      if (LWTYPE(ICURRWIN).le.1) then
        if (LSCALE(ICURRWIN).eq.0) then
          sx1=x1
          sy1=y1
        else
          sx1=XGSCAL(ICURRWIN,1)*x1+XGSCAL(ICURRWIN,3)
          sy1=XGSCAL(ICURRWIN,2)*y1+XGSCAL(ICURRWIN,4)
        endif
        if (LWTYPE(ICURRWIN).eq.0) call wvscte(sx1,sy1,h,Mess,a,nn,ic)
c        if (LWTYPE(ICURRWIN).eq.1) call wvdfte(sx1,sy1,h,Mess,a,nn,ic)
      endif
c      write(*,*) 'symbol1'
c
      if (LWTYPE(ICURRWIN).eq.2) then
        if (LSCALE(ICURRWIN).eq.0) then
           xx1 = x1
           xy1 = y1
        else
           xx1=XGSCAL(ICURRWIN,1)*x1+XGSCAL(ICURRWIN,3)
           xy1=XGSCAL(ICURRWIN,2)*y1+XGSCAL(ICURRWIN,4)
        endif
        call wvpste(xx1,xy1,h,Mess,a,nn,ic)
      endif
c      write(*,*) 'symbol2'
c
      if (IGPSP(ICURRWIN).eq.1) THEN
        if (LSCALE(ICURRWIN).eq.0) then
          xx1 = x1*10.0/MAXSCX
          xy1 = y1*10.0/MAXSCX
        else
          xx1=(XGSCAL(ICURRWIN,1)*x1+XGSCAL(ICURRWIN,3))*10.0/MAXSCX
c          xx2=(XGSCAL(ICURRWIN,1)*x2+XGSCAL(ICURRWIN,3))*10.0/MAXSCX
          xx2=(XGSCAL(ICURRWIN,1)*x1+XGSCAL(ICURRWIN,3))*10.0/MAXSCX
        endif
c        write(*,*) 'before wvpste'
        call wvpste(xx1,xy1,h,Mess,a,nn,ic)
      endif
c      write(*,*) 'symbol3'
c
      return
      end

      subroutine setpost(ionoff)
c
c-----------------------------------------------------------------------c
c   routine name       - setpost
c
c-----------------------------------------------------------------------c
c   computer           - machine independent
c
c   latest revision    - dec 90
c
c   purpose            - routine turns on/off postscript echo for
c                        a window (not for the postscripts windows
c                        only)
c
c   usage              - call setpost(Ionoff)
c
c   arguments :
c     in:
c            Ionoff    = 0  no poscript echo (close file)
c                      = 1  poscript echo (new file)
c
c   required  routines - wverr
c
c-----------------------------------------------------------------------c
#include "wvglob.blk"
c
      if(IOPEN.ne.1) call wverr(2)
      if (LWTYPE(ICURRWIN).eq.2) return
      if (Ionoff.eq.0) then
        if (IGPSP(ICURRWIN).eq.0) return
        IGPSP(ICURRWIN)=0
        ILEPSP=ILEPSP-1
        if (ILEPSP.eq.0) call wvclps
      else
        if (IGPSP(ICURRWIN).eq.1) return
        IGPSP(ICURRWIN)=1
        ILEPSP=ILEPSP+1
        if (ILEPSP.eq.1) call wvopps
        call selwind(ICURRWIN-1)
      endif
c
      return
      end
      subroutine clrwind(Iback)
c
c-----------------------------------------------------------------------c
c   routine name       - clrwind
c
c-----------------------------------------------------------------------c
c   computer           - machine independent
c
c   latest revision    - dec 90
c
c   purpose            - clears a window
c
c   usage              - call clrwind(Iback)
c
c   arguments :
c     in:
c            Iback     - color number for the background
c
c   required  routines -
c
c-----------------------------------------------------------------------c
#include "wvglob.blk"
c
      real(8) xy(2,4)
c
      is = LSCALE(ICURRWIN)
      xy(1,1)=LWSIZE(ICURRWIN,1)
      xy(2,1)=LWSIZE(ICURRWIN,2)
      xy(1,2)=LWSIZE(ICURRWIN,1)
      xy(2,2)=LWSIZE(ICURRWIN,4)-1
      xy(1,3)=LWSIZE(ICURRWIN,3)-1
      xy(2,3)=LWSIZE(ICURRWIN,4)-1
      xy(1,4)=LWSIZE(ICURRWIN,3)-1
      xy(2,4)=LWSIZE(ICURRWIN,2)
      LSCALE(ICURRWIN)=0
      call fillpoly(4,xy,iback,iback)
      LSCALE(ICURRWIN)=is
c
      return
      end


      subroutine scaleon(Ionoff)
c
c-----------------------------------------------------------------------c
c   routine name       - scaleon
c
c-----------------------------------------------------------------------c
c   computer           - machine independent
c
c   latest revision    - dec 90
c
c   purpose            - routine turns scaling on/off
c
c   usage              - call scaleon(Ionoff)
c
c   arguments :
c     in:
c          Ionoff      = 0  no scaling
c                      = 1  scale all data
c
c   required  routines - wverr
c
c-----------------------------------------------------------------------c
#include "wvglob.blk"
c
      if (IOPEN.ne.1) call wverr(2)
      LSCALE(ICURRWIN) = Ionoff
c
      return
      end
      subroutine newscale(Range)
c
c-----------------------------------------------------------------------c
c   routine name       - newscale
c
c-----------------------------------------------------------------------c
c   computer           - machine independent
c
c   latest revision    - dec 90
c
c   purpose            - routine sets offset and scale factors
c
c   usage              - call newscale(Range)
c
c   arguments :
c     in:
c            Range     - xmin,ymin,xmax,ymax
c
c   required  routines -
c
c-----------------------------------------------------------------------c
      real(8) Range(4)
c
#include "wvglob.blk"
c
      if (IOPEN.ne.1) call wverr(2)
      do 10 i=1,2
        XGSCAL(ICURRWIN,i)  =(LWSIZE(ICURRWIN,i+2)-LWSIZE(ICURRWIN,i))/
     .                       (Range(2+i)-Range(i))
        XGSCAL(ICURRWIN,i+2)=LWSIZE(ICURRWIN,i)-
     .                       Range(i)*XGSCAL(ICURRWIN,i)
10    continue
c
      call scaleon(1)
c
      return
      end

      subroutine initcm
c
c-----------------------------------------------------------------------c
c   routine name       - initcm
c
c-----------------------------------------------------------------------c
c   computer           - machine independent
c
c   latest revision    - dec 90
c
c   purpose            - routine initiates colormaps for xwindow driver
c
c   usage              - call initcm
c
c   arguments          - none
c
c   required  routines -
c
c-----------------------------------------------------------------------c
c-----------------------------------------------------------------------c     defined colors are 0,1 - black,white
c      2,..,15 - colors for mesh plot
c     16,..,66 - 50 colors in the rampe sequence:
c                yellow, green, cyan, blue, magenta, red, yellow, white
c       exact definition of rampe is given in IRAM array below
c---------------------------------------------------------------------- c
c.....set of full rampe of colors
c
      dimension iram(4,8)
      data iram /  0,200,255,  0, 130,  0,255,  0,
     .           220,  0,255,255, 320,  0, 60,250,
     .           520,255,  0,255, 630,255,  0,  0,
     .           870,255,255,0,  1000,255,255,255 /
c
      call xgsetcm(0,  0,  0,  0)
      call xgsetcm(1,255,255,255)
c
c##c     for b/w machine we have only 2 colors
c##      if (NCOLOR.eq.0) return
c##
c##      if (NGREY.eq.1) then
c##          ncol = 14
c##          if (Iset.gt.10 ) ncol = Iset
c##          ndcol = 255 / (ncol+1)
c##          do 30 i = 1,ncol
c##              irgb = 255-ndcol*i
c##              call xgsetcm(nshft+i,irgb,irgb,irgb)
c##   30     continue
c##          return
c##      endif
c##c
c##      if (ISET.EQ.1) THEN
c
c         color map for mesh plot
c
      call xgsetcm(2, 255,255,0)
      call xgsetcm(3, 255,162,0)
      call xgsetcm(4, 255,0,  0)
      call xgsetcm(5, 200,55, 230)
      call xgsetcm(6, 0,  0,  255)
      call xgsetcm(7, 0,  155,255)
      call xgsetcm(8, 0,  255,0)
      call xgsetcm(9, 140,140,140)
      call xgsetcm(10,255,255,255)
      call xgsetcm(11,140,140,140)
      call xgsetcm(12,255,255,255)
      call xgsetcm(13,140,140,140)
      call xgsetcm(14,0,  0,  0)
      call xgsetcm(15,200,200,200)
c
c##      ENDIF
c##
c##      IF (ISET.EQ.2) THEN
c##c
c##c         color map for contour plot
c##c
c##         call xgsetcm(ncolshft+2, 255,255,0)
c##         call xgsetcm(ncolshft+3, 255,235,45)
c##         call xgsetcm(ncolshft+4, 255,195,75)
c##         call xgsetcm(ncolshft+5, 255,155,90)
c##         call xgsetcm(ncolshft+6, 255,105,75)
c##         call xgsetcm(ncolshft+7, 255,70, 75)
c##         call xgsetcm(ncolshft+8, 255,0,  25)
c##  	      call xgsetcm(ncolshft+9, 210,0,  85)
c##         call xgsetcm(ncolshft+10,200,55, 150)
c##         call xgsetcm(ncolshft+11,155,55, 220)
c##         call xgsetcm(ncolshft+12,55, 25, 230)
c##         call xgsetcm(ncolshft+13,0,  0,  255)
c##         call xgsetcm(ncolshft+14,0,  100,255)
c##         call xgsetcm(ncolshft+15,200,200,200)
c##
c##      ENDIF
c##
c##      if (Iset.gt.13) then
c##c
c##c         soft ramp of colors
c##c
      nrcol = 0
      nrrow = 1
      iset = 51

      do 10 i=0,iset
        if (nrcol.ge.iram(1,nrrow+1)) nrrow = nrrow + 1
        alfa = float(nrcol-iram(1,nrrow))/
     .        (iram(1,nrrow+1)-iram(1,nrrow))
        ir = alfa*iram(2,nrrow+1)+(1-alfa)*iram(2,nrrow)
        ig = alfa*iram(3,nrrow+1)+(1-alfa)*iram(3,nrrow)
        ib = alfa*iram(4,nrrow+1)+(1-alfa)*iram(4,nrrow)
        call xgsetcm(16+i,ir,ig,ib)
        nrcol = nrcol + 1000/iset
   10 continue
c
c##      endif
c
      return
      end


c ********************************************************************
c
      subroutine flushx

c********************************************************************

c
c
      return
      end



      subroutine wvsccl(Iclip)
c
c-----------------------------------------------------------------------
c
c   routine name       - wvsccl
c
c-----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    -
c
c   purpose            -
c
c   usage              - call wvsccl(Iclip)
c
c   arguments :
c     in:
c              Iclip   -
c     out:
c
c   required  routines -
c
c-----------------------------------------------------------------------
c
#include "rastcom.blk"
c
      dimension Iclip(4)
c
      iprint=0
c
      if (iprint.eq.1) then
        write(*,*)'IN WVSCCL IFILE = ',IFILE
      endif
c
      if (IFILE.eq.1000) then
c
c  .....set clip for tektronics
        return
      endif
c
      if (IFILE.eq.1001) then
c
c  .....set clip for X-windows
        call xgsetclip(iclip(1),iclip(2),iclip(3),iclip(4))
        return
      endif
c
      return
      end


C========================================================
        SUBROUTINE WVSCDL(X1,Y1,X2,Y2,IC)
C========================================================
        common/rastcom/ifile,idump,lastcol,lastx,lasty
        integer X1,X2,Y1,Y2

        if (ifile.eq.1001) then
c           X-window
            call xgsetcol(ic)
            call xgdrawl(x1,y1,x2,y2)
            return
        endif
        RETURN
        END
C========================================================
        SUBROUTINE WVSCFP(N,XC,ICF,ICL)
C========================================================
        common/rastcom/ifile,idump,lastcol,lastx,lasty
        integer XC(*),ix(100),iy(100)
        if (ifile.eq.1001) then
c           fillpoly for X
            call xgsetcol(icf)
            do 5 i=1,n
                ix(i)=xc(2*i-1)
                iy(i)=xc(2*i)
    5       continue
            call xgdrfpoly(n,ix,iy)
            if (icf.ne.icl) then
                call xgsetcol(icl)
                call xgdrpoly(n,ix,iy)
            endif
            return
        endif
        RETURN
        END
C========================================================
        SUBROUTINE WVSCPX(X,Y,IC)
C=======================================================
        integer X,Y
        common/rastcom/ifile,idump,lastcol,lastx,lasty
        if (ifile.eq.1001) then
            call xgsetcol(ic)
            call xgdrawl(x,y,x,y)
            return
        endif
        RETURN
        END
C========================================================
        SUBROUTINE WVSCTE(X1,Y1,H,MESS,A,N,ICOL)
C========================================================
        integer X1,Y1,H,A
        CHARACTER MESS*(*),text*100
        common/rastcom/ifile,idump,lastcol,lastx,lasty
        if (ifile.eq.1001) then
c           text for X
            call xgsetcol(Icol)
            call xgsettl(x1,y1)
            text(1:n)=mess(1:n)
            text(n+1:n+1)=char(0)
            call xgdrtxt(text,n)
            return
        endif
        RETURN
        END
      subroutine wvwarn(N)
c
c----------------------------------------------------------------------
c
c   routine name       - wvwarn
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - feb, 92
c
c   purpose            - routine issues warning messages
c
c   usage              - call wvwarn(N)
c
c   arguments :
c     in:
c               N      - error flag
c
c   required  routines - wverr
c
c----------------------------------------------------------------------
c
      common /wvwar1/ IERRNO
c
      IERRNO=IERRNO+1
c
c     if warning messages are desired add appropriate
c     write statements here
c
      if (IERRNO.lt.100) return
c
      write(*,*) 'WVWARN: More warnings than 100.'
      if (n.eq.0) write(*,*) 'Last message : Clip error'
      if (n.eq.1) write(*,*) 'Last message : Window number'
      if (n.eq.2) write(*,*) 'Last message : FILLPOLY size'
      call wverr(0)
c
      end

      subroutine wverr(N)
c
c-----------------------------------------------------------------------c
c   routine name       - wverr
c
c-----------------------------------------------------------------------c
c   computer           - machine independent
c
c   latest revision    - dec 90
c
c   purpose            - error message routine
c
c   usage              - call wverr(N)
c
c   arguments :
c     in:
c                  N   - error flag
c     out:
c
c   required  routines - closwind
c
c-----------------------------------------------------------------------c
      do 10 i=0,9
        call closwind(i)
   10 continue
c
      write(*,*) 'GRAPHICS ERROR : '
      if (n.eq.0) write(*,*) ' END OF PROGRAM '
      if (n.eq.1) write(*,*) ' WRONG WINDOW NUMBER'
      if (n.eq.2) write(*,*) ' INITIALIZE GRAPHIC'
c
      stop 1
      end
c---------------------------------------------------------------------
c
c   routine name       - wvgrph
c
c---------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - dec 90
c
c   purpose            - sets screen to graphic mode and sets screen
c                        sizes and max color in /wvscrn/
c
c   usage              - call wvgrph
c
c   arguments          - none
c
c   required  routines - inittek,xgclear,xginfo,xgnxtev,xgopw,xgsetcol
c
c---------------------------------------------------------------------
c
      subroutine wvgrph
      use graphmod, only: IWINDL, IWINDH
c
#include "rastcom.blk"
#include "wvscrn.blk"
c
      integer xglength
      character(len=1024) :: res_string
c
      iprint=0
      if (iprint.eq.1) then
        write(*,*) 'wvgrph: IFILE = ',IFILE
      endif
c
      if (IFILE.eq.1001) then
c
c.......another entry into X
  110   continue
        if (xglength().gt.0) then
          call xgnxtev(i1,i2,i3,i4)
          go to 110
        endif
        if (xglength().ne.0) write(*,*) 'Events !!'
        call xginfo(MAXSCX,MAXSCY,MAXCOL,idum1,idum2)
        call xgsetcol(0)
        call xgclear
        return
      endif
c
c---------------------------------------------------------------------
c
c  test for file T4014 in current dir
c  if it exists then open Tektronics graph
c
c      open (99,file='T4014',status='old',err=91)
c
c  file exists -
c      close (99)
c
c  open Tektronics - B-W draw only - no fill
c      ifile=1000
c      MAXSCX=4000
c      MAXSCY=3000
c      MAXCOL=1
c      call inittek
c      return
c
c---------------------------------------------------------------------
c
c   91 continue
c
c  open X-window connection - hardwired
ccc      call xgopw(3,'GRAPHICS'//char(0),'-geometry'//char(0),
ccc     .    '700x500'//char(0),char(0),char(0),char(0))
ccc      call xgopw(3,'GRAPHICS'//char(0),'-geometry'//char(0),
ccc     .    '1000x600'//char(0),char(0),char(0),char(0))
ccc      call xgopw(3,'GRAPHICS'//char(0),'-geometry'//char(0),
ccc     .    '667x400'//char(0),char(0),char(0),char(0))
      write(res_string, fmt="(i0,'x',i0)") IWINDL, IWINDH ! jz 1/21/2013
      call xgopw(3,'GRAPHICS'//char(0),'-geometry'//char(0),
     .    trim(res_string)//char(0),char(0),char(0),char(0))
c     .    '1334x800'//char(0),char(0),char(0),char(0))
  100   call xgnxtev(i1,i2,i3,i4)
      if (xglength().gt.0) goto 100
      if (xglength().gt.0) write(*,*) ' EVENTS !'
      call xginfo(MAXSCX,MAXSCY,MAXCOL,idum1,idum2)
c
c  initialize color map
      call initcm
      IFILE=1001
c
c---------------------------------------------------------------------
c
      return
      end
C========================================================
        subroutine wvtext
C======================================================
c               performs QUIT
        integer xglength
        common/rastcom/ifile,idump,lastcol,lastx,lasty
        common/wvscrn/MAXSCx,MAXSCy,maxcol,MAXDFx,MAXDFy
c          if it was opened....
        if (maxcol.ne.0) then
               if (ifile.eq.1001) then
c               flush and wait for any event
                call xginfo(MAXSCX,MAXSCY,MAXCOL,id1,id2)
                write(*,*)'Ev queue :',xglength()
   90           if (xglength().gt.0) then
                    call xgnxtev(i1,i2,i3,i4)
                    goto 90
                endif
                if (xglength().ne.0) write(*,*)'EV!'
  100           call xgnxtev(i1,i2,i3,i4)
                if (xglength().gt.0) goto 100
                if (xglength().ne.0) write(*,*)'Ev!'
c               do nothing - no close for xwindow
            endif
            write(*,*) 'Graphics closed'
        endif
        maxcol=0
        return
        end

#endif

#if HP3D_USE_X11

!  GRAPHICAL INTERFACE FOR X-WINDOWS AND POSTCRIPT FILES
!
!  THE AVALIBLE BASIC ROUTINES ARE AS FOLLOWS
!
!  CLOSWIND(Iwin)- to end pictures and set the (ps off)
!  OPENWIND(Iwtype,Iwsize,Iwin)- initialzing the graphics
!  SELWIND(Iwin) - to begin pictures( gets the dimensions of the object
!      to be drawn and returns the dimensions of the image space; also
!      rescales and clears the window)
!  DRAWLINE(X1,Y1,X2,Y2,Ilcol) - to draw a line
!  DRAWPOLY(N,Coor,Ilcol) - to draw a polygon
!  FILLPOLY(Nin,Coor,Ifcol,Ilcol) - to draw a filled polygon
!  SYMBOL(X1,Y1,Height,Mess,Angle,N,Icol)- to write text
!
!  AUXILLARY ROUTINES
!  SETPOST(Ionoff)- to set 1-on/0-off ps echo
!  CLRWIND
!
!  PARAMETERS
!    Iwtype - window type
!    Ixsize,Iysize - x,y window dimensions(in pixels)
!    Iback,Ilcol,Ifcol - color numbers (Ifcol to fill)
!    Mess - charcter array of length N
!    Height,Angle - parameters of text(for ps only!)
!    X1,Y1 (or array Coor) -x,y coordinates of point(s)
!
!    Further Information in the Routines
!
!  Colors are as Follows
!        1-white 2-black
!        3-10 - sparse spectrum of colors
!        11-100 - dense spectrum of colors
!
!
!
!
!**************************************************************

   subroutine closwind(Iwin)
      use wvglob
      implicit none
!
!-----------------------------------------------------------------------!
!   routine name       - closwind
!
!-----------------------------------------------------------------------!
!   author - Satish C. and many others
!
!   computer           - machine independent
!
!   latest revision    - 1998
!
!   purpose            - routine closes a window
!
!   usage              - call closwind(Iwin)
!
!   arguments :
!     in:
!          Iwin        - number of the window to close
!
!   required  routines - selwind,setpost,wverr,wvcldf,wvclps,wvtext,
!                        wvwarn
!
!----------------------------------------------------------------------
!
!
      integer :: Iwin
      integer :: ic,i,k,ilepf
!
      if(IOPEN.ne.1) call wverr(2)
!
      if (Iwin.lt.0.or.Iwin.gt.9) then
        call wvwarn(1)
      else
        ic=iwin+1
        if (LWTYPE(ic).lt.0) return
!
!.......check for postscript
        call setpost(0)
        if (LWTYPE(ic).eq.2) then
          ILEPSP=ILEPSP-1
          if (ILEPSP.eq.0) call wvclps
        endif
!
!.......check for graph files
        if (LWTYPE(ic).eq.1) then
          ilepf=ilepf-1
!          if (ilepf.eq.0) call wvcldf
        endif
        if (LWTYPE(ic).eq.0) then
          ILESCR=ILESCR-1
          if (ILESCR.eq.0) call wvtext
        endif
        LWTYPE(ic)=-1
        if (Iwin+1.eq.ic) then
!!!!        if (iwind+1.eq.ic) then
!
!.........switch to any active window
          do 10 i=1,10
            if (LWTYPE(i).ge.0) then
              k=i-1
              call selwind(k)
              return
            endif
   10     continue
        endif
      endif
!
      return
   end subroutine closwind
!
!----------------------------------------------------------------------
!
!   routine name       - openwind
!
!----------------------------------------------------------------------
!   author - Satish C. and many others
!
!   computer           - machine independent
!
!   latest revision    - 1998
!
!   purpose            - routine opens a graphic window
!
!   usage              - call openwind(Iwtype,Iwsize,Iwin)
!
!   arguments :
!   in:
!       Iwtype = 0     - screen
!                1     - disk file (not used)
!                2     - poscript file
!       Iwsize(4)      - Xmin,Ymin,Xmax,Ymax
!       Iwin           - number of the window = 0,1,...,9
!
!   out:
!       Iwsize         - adjusted size of the window
!
!   required routines  - closewind,newscale,scaleon,selwind,setpost,
!                        wverr,wvgrph,wvopdf,wvopps
!
!---------------------------------------------------------------------
!
   subroutine openwind(Iwtype,Iwsize,Iwin)
      use wvglob
      use wvscrn
      implicit none
!
      integer(4) :: Iwtype,Iwsize(4),Iwin
!
!     dumpfile maxx and maxy
!
      integer :: i,kwin
      integer :: iprint
!
      real(8) :: range(4)
      range = (/0d0,0d0,1d0,1d0/)
!
      iprint=0
      if (iprint.eq.1) then
        write(*,*) 'openwind: IOPEN = ',IOPEN
      endif
!
      if (IOPEN.ne.1.or.Iwtype.le.-1) then
!
!.......initialize window system
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
!
      if (Iwin.lt.0.or.Iwin.gt.9) call wverr(1)
      kwin=Iwin+1
!
      if (iprint.eq.1) then
        write(*,*) 'openwind: kwin,LWTYPE(kwin) = ',kwin,LWTYPE(kwin)
      endif
!
      if (LWTYPE(kwin).ge.0) then
!
!.......clear if window in use
        call closwind(Iwin)
      endif
!
      do i=1,4
        LWSIZE(kwin,i)=max0(0,Iwsize(i))
      enddo
!
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
!
!!!        write(*,*)'done with calling the file'
      if (Iwtype.eq.1) then
        ILEDF=ILEDF+1
!        ilepf=ilepf+1
!        if (ILEDF.eq.1) call wvopdf
!        if (ilepf .eq.1) call wvopdf
        if (Iwsize(3).eq.0) then
          Iwsize(3)=MAXDFX
          Iwsize(4)=MAXDFY
        endif
        LWSIZE(kwin,3)=min(Iwsize(3),MAXDFX)
        LWSIZE(kwin,4)=min(Iwsize(4),MAXDFY)
      endif
!
      if (Iwtype.eq.0) then
        ILESCR=ILESCR+1
!
        if (iprint.eq.1) then
          write(*,*) 'openwind: ILESCR = ',ILESCR
        endif
!
        if (ILESCR.eq.1) call wvgrph
!!!        write(*,*) 'openwind: AFTER CALL TO wvgrph'
        if (Iwsize(3).eq.0) then
          Iwsize(3)=MAXSCX
          Iwsize(4)=MAXSCY
        endif
!!!        write(*,*)'in openwind  MAXSCX,MAXSCY = ', MAXSCX,MAXSCY
        LWSIZE(kwin,3)=min(Iwsize(3),MAXSCX)
        LWSIZE(kwin,4)=min(Iwsize(4),MAXSCY)
      endif
!
!.....adjust Iwsize on output
      do 45 i=1,4
        Iwsize(i)=LWSIZE(kwin,i)
   45 continue
      LWTYPE(kwin)=Iwtype
      if (iprint.eq.1) then
        write(*,*)'in openwind after label 45 kwin,LWTYPE(kwin) = ', &
                   kwin,LWTYPE(kwin)
      endif
!
      call selwind(Iwin)
!
!.....set default scale
      call newscale(range)
      call scaleon(0)
      if (Iwtype.eq.0) then
        open(99,file='LASER',status='old',err=123)
        close(99)
        call setpost(1)
      endif
!
  123 continue
!
      return
   end subroutine openwind
!
!
   subroutine selwind(Iwin)
      use wvglob
      use wvscrn
      implicit none
!
!-----------------------------------------------------------------------
!
!   routine name       - selwind
!
!-----------------------------------------------------------------------
!
!   computer           - machine independent
!
!   latest revision    - dec 90
!
!   purpose            - routine selects a graphic window to use
!
!   usage              - call selwind(Iwin)
!
!   arguments :
!     in:
!             Iwin     - number of the window
!
!   required  routines - wvwarn,wverr,wvpscl,wvsccl
!
!-----------------------------------------------------------------------
!
      integer :: Iwin
!
      integer :: iclip(4),clip(4)
!
      integer :: i
      integer :: iprint
      iprint=0
!
      if (IOPEN.ne.1) call wverr(2)
      if (Iwin.lt.0.or.Iwin.gt.9) then
        call wvwarn(1)
      else
        ICURRWIN = Iwin + 1
        do 10 i=1,4
          iclip(i) = LWSIZE(ICURRWIN,i)
   10   continue
        if (LWTYPE(ICURRWIN).eq.0) call wvsccl(iclip)
!        if (LWTYPE(ICURRWIN).eq.1) call wvdfcl(iclip)
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
            clip(i)=int(iclip(i)*10.d0/MAXSCX)
   30     continue
          if (iprint.eq.1) then
            write(*,*) 'SELWIND: clip = ',clip
          endif
          call wvpscl(clip)
        endif
      endif
!
      return
   end subroutine selwind
!
!
   subroutine drawline(X1,Y1,X2,Y2,Ilcol)
      use wvglob
      use wvscrn
      implicit none
!
!---------------------------------------------------------------------
!
!   routine name       - drawline
!
!---------------------------------------------------------------------
!   author     - Satish C. and many others
!
!   computer           - machine independent
!
!   latest revision    - 1998
!
!   purpose            - routine draws a line
!
!   usage              - call drawline(X1,Y1,X2,Y2,Ilcol)
!
!   arguments :
!     in:
!          X1,Y1,X2,Y2 - coordinates of two points
!          Ilcol       - color of the line
!     out:
!
!   required  routines -
!
!---------------------------------------------------------------------
!
      real(8) :: X1,X2,Y1,Y2
      integer :: Ilcol
!
      real(8) :: xx1,xx2,xy1,xy2
      integer :: sx1,sx2,sy1,sy2
      integer :: il
!
      il = Ilcol
!
      if (LWTYPE(ICURRWIN).le.1) then
          if (LSCALE(ICURRWIN).eq.0) then
              sx1=int(X1)
              sx2=int(X2)
              sy2=int(Y2)
              sy1=int(Y1)
          else
              sx1=int(XGSCAL(ICURRWIN,1)*x1+XGSCAL(ICURRWIN,3))
              sx2=int(XGSCAL(ICURRWIN,1)*x2+XGSCAL(ICURRWIN,3))
              sy1=int(XGSCAL(ICURRWIN,2)*y1+XGSCAL(ICURRWIN,4))
              sy2=int(XGSCAL(ICURRWIN,2)*y2+XGSCAL(ICURRWIN,4))
          endif
          if (LWTYPE(ICURRWIN).eq.0) call wvscdl(sx1,sy1,sx2,sy2,il)
!          if (LWTYPE(ICURRWIN).eq.1) call wvdfdl(sx1,sy1,sx2,sy2,il)
      endif
      if (LWTYPE(ICURRWIN).eq.2) then
          if (LSCALE(ICURRWIN).eq.0) then
              xx1=X1
              xx2=X2
              xy2=Y2
              xy1=Y1
              call wvpsdl(xx1,xy1,xx2,xy2,il)
          else
              xx1=XGSCAL(ICURRWIN,1)*x1+XGSCAL(ICURRWIN,3)
              xx2=XGSCAL(ICURRWIN,1)*x2+XGSCAL(ICURRWIN,3)
              xy1=XGSCAL(ICURRWIN,2)*y1+XGSCAL(ICURRWIN,4)
              xy2=XGSCAL(ICURRWIN,2)*y2+XGSCAL(ICURRWIN,4)
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
          xx1=(XGSCAL(ICURRWIN,1)*X1+XGSCAL(ICURRWIN,3))*10.0/MAXSCX
          xx2=(XGSCAL(ICURRWIN,1)*X2+XGSCAL(ICURRWIN,3))*10.0/MAXSCX
          xy1=(XGSCAL(ICURRWIN,2)*Y1+XGSCAL(ICURRWIN,4))*10.0/MAXSCX
          xy2=(XGSCAL(ICURRWIN,2)*Y2+XGSCAL(ICURRWIN,4))*10.0/MAXSCX
          call wvpsdl(xx1,xy1,xx2,xy2,il)
        endif
      endif
!
      return
   end subroutine drawline
!
!
   subroutine drawpoly(N,Coor,Ilcol)
      implicit none
!
!-----------------------------------------------------------------------!
!   routine name       - drawpoly
!
!-----------------------------------------------------------------------!
!   author  - Satish C. and many others
!
!   computer           - machine independent
!
!   latest revision    - 1998
!
!   purpose            - routine draws a polygon
!
!   usage              - call drawpoly(N,Coor,Ilcol)
!
!   arguments :
!     in:
!               N      - number of the vertices
!               Coor   - coordinates of the vertices
!               Ilcol  - line color number
!
!   required  routines - fillpoly
!
!-----------------------------------------------------------------------!
      integer :: N
      real(8) :: Coor(2,N)
      integer :: Ilcol
!
      call fillpoly(N,Coor,-1,Ilcol)
!
      return
   end subroutine drawpoly
!
!
   subroutine fillpoly(Nin,Coor,Ifcol,Ilcol)
      use wvscrn
      use wvglob
      implicit none
!
!-----------------------------------------------------------------------
!
!   routine name       - fillpoly
!
!-----------------------------------------------------------------------
!   author   - Satish !. and many others
!
!   computer           - machine independent
!
!   latest revision    - 1998
!
!   purpose            - routine fills a polygon with a color
!
!   usage              - call fillpoly(Nin,Coor,Ifcol,Ilcol)
!
!   arguments :
!     in:
!               Nin    - number of the vertices
!               Coor   - coordinates of the vertices
!               Ifcol  - field color number
!               Ilcol  - line color number
!
!   required  routines - wvpsfp,wvwarn
!
!-----------------------------------------------------------------------
!
      integer :: Nin
      real(8) :: Coor(2,Nin)
      integer :: Ifcol,Ilcol
!
      integer :: scoor(2,65)
      real(8) :: xcoor(2,65)
      integer :: i,if,il,n
      real(8) :: xmult
!
      if (Nin.gt.65) call wvwarn(2)
      n=min(65,Nin)
      if=Ifcol
      il=Ilcol
!
      if (LWTYPE(ICURRWIN).le.1) then
        if (LSCALE(ICURRWIN).eq.0) then
          do 17 i=1,n
            scoor(1,i)=int(Coor(1,i))
            scoor(2,i)=int(Coor(2,i))
   17     continue
        else
          do 20 i=1,n
            scoor(1,i)=int(XGSCAL(ICURRWIN,1)*Coor(1,i)+ &
                           XGSCAL(ICURRWIN,3))
            scoor(2,i)=int(XGSCAL(ICURRWIN,2)*Coor(2,i)+ &
                           XGSCAL(ICURRWIN,4))
   20     continue
        endif
        if (LWTYPE(ICURRWIN).eq.0) call wvscfp(n,scoor,il,if)
!        if (LWTYPE(ICURRWIN).eq.1) call wvdffp(n,scoor,il,if)
      endif
!
      if (LWTYPE(ICURRWIN).eq.2.or.IGPSP(ICURRWIN).eq.1) then
        xmult = 1.d0
        if (IGPSP(ICURRWIN).eq.1) xmult = 10.0/maxscx
        if (LSCALE(ICURRWIN).eq.0) then
          do 30 i=1,n
            xcoor(1,i)=Coor(1,i)*xmult
            xcoor(2,i)=Coor(2,i)*xmult
   30     continue
        else
          do 40 i=1,n
            xcoor(1,i)=(XGSCAL(ICURRWIN,1)*Coor(1,i)+ &
                        XGSCAL(ICURRWIN,3))*xmult
            xcoor(2,i)=(XGSCAL(ICURRWIN,2)*Coor(2,i)+ &
                        XGSCAL(ICURRWIN,4))*xmult
   40     continue
        endif
!..ldem        call wvpsfp(n,xcoor,il,if)
        call wvpsfp(n,xcoor,if,il)
      endif
!
      return
   end subroutine fillpoly
!
!
   subroutine symbol(X1,Y1,Height,Mess,Angle,N,Icol)
      use wvglob
      use wvscrn
      implicit none
!
!-----------------------------------------------------------------------!
!   routine name       - symbol
!
!-----------------------------------------------------------------------!
!   computer           - machine independent
!
!   latest revision    - dec 90
!
!   purpose            - routine draws a text or char marker
!
!   usage              - call symbol(X1,Y1,Height,Mess,Angle,N,Icol)
!
!   arguments :
!     in:
!            X1,Y1     - position to draw
!            Height    - height of the symbol
!            Angle     - angle of drawing
!            Mess      - character(s) to draw
!            N         - number of characters in Mess
!            Icol      - color to use
!     out:
!
!   required  routines -
!
!-----------------------------------------------------------------------!
      real(8)   :: X1,Y1,Height,Angle
      character :: Mess*(*)
      integer   :: N,Icol
!
      real(8) :: xx1,xy1,xx2
      integer :: h,a,sx1,sy1
      integer :: nn,ic
!
!      write(*,*) 'in symbol'
!      write(*,*) 'Mess',Mess
      nn = N
      ic = Icol
      h = int(Height*16)
      a = int(Angle)
!
      if (LWTYPE(ICURRWIN).le.1) then
        if (LSCALE(ICURRWIN).eq.0) then
          sx1=int(X1)
          sy1=int(Y1)
        else
          sx1=int(XGSCAL(ICURRWIN,1)*X1+XGSCAL(ICURRWIN,3))
          sy1=int(XGSCAL(ICURRWIN,2)*Y1+XGSCAL(ICURRWIN,4))
        endif
        if (LWTYPE(ICURRWIN).eq.0) call wvscte(sx1,sy1,h,Mess,a,nn,ic)
!        if (LWTYPE(ICURRWIN).eq.1) call wvdfte(sx1,sy1,h,Mess,a,nn,ic)
      endif
!      write(*,*) 'symbol1'
!
      if (LWTYPE(ICURRWIN).eq.2) then
        if (LSCALE(ICURRWIN).eq.0) then
           xx1 = X1
           xy1 = Y1
        else
           xx1=XGSCAL(ICURRWIN,1)*X1+XGSCAL(ICURRWIN,3)
           xy1=XGSCAL(ICURRWIN,2)*Y1+XGSCAL(ICURRWIN,4)
        endif
        call wvpste(xx1,xy1,h,Mess,a,nn,ic)
      endif
!      write(*,*) 'symbol2'
!
      if (IGPSP(ICURRWIN).eq.1) then
        if (LSCALE(ICURRWIN).eq.0) then
          xx1 = X1*10.0/MAXSCX
          xy1 = Y1*10.0/MAXSCX
        else
          xx1=(XGSCAL(ICURRWIN,1)*x1+XGSCAL(ICURRWIN,3))*10.0/MAXSCX
!          xx2=(XGSCAL(ICURRWIN,1)*x2+XGSCAL(ICURRWIN,3))*10.0/MAXSCX
          xx2=(XGSCAL(ICURRWIN,1)*x1+XGSCAL(ICURRWIN,3))*10.0/MAXSCX
        endif
!        write(*,*) 'before wvpste'
        call wvpste(xx1,xy1,h,Mess,a,nn,ic)
      endif
!      write(*,*) 'symbol3'
!
      return
   end subroutine symbol
!
!
   subroutine setpost(Ionoff)
      use wvglob
      implicit none
!
!-----------------------------------------------------------------------!
!   routine name       - setpost
!
!-----------------------------------------------------------------------!
!   computer           - machine independent
!
!   latest revision    - dec 90
!
!   purpose            - routine turns on/off postscript echo for
!                        a window (not for the postscripts windows
!                        only)
!
!   usage              - call setpost(Ionoff)
!
!   arguments :
!     in:
!            Ionoff    = 0  no poscript echo (close file)
!                      = 1  poscript echo (new file)
!
!   required  routines - wverr
!
!-----------------------------------------------------------------------!
!
      integer :: Ionoff
!
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
!
      return
   end subroutine setpost
!
!
   subroutine clrwind(Iback)
      use wvglob
      implicit none
!
!-----------------------------------------------------------------------!
!   routine name       - clrwind
!
!-----------------------------------------------------------------------!
!   computer           - machine independent
!
!   latest revision    - dec 90
!
!   purpose            - clears a window
!
!   usage              - call clrwind(Iback)
!
!   arguments :
!     in:
!            Iback     - color number for the background
!
!   required  routines -
!
!-----------------------------------------------------------------------!
!
      integer :: Iback
!
      integer :: is
      real(8) :: xy(2,4)
!
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
      call fillpoly(4,xy,Iback,Iback)
      LSCALE(ICURRWIN)=is
!
      return
   end subroutine clrwind
!
!
   subroutine scaleon(Ionoff)
      use wvglob
      implicit none
!
!-----------------------------------------------------------------------!
!   routine name       - scaleon
!
!-----------------------------------------------------------------------!
!   computer           - machine independent
!
!   latest revision    - dec 90
!
!   purpose            - routine turns scaling on/off
!
!   usage              - call scaleon(Ionoff)
!
!   arguments :
!     in:
!          Ionoff      = 0  no scaling
!                      = 1  scale all data
!
!   required  routines - wverr
!
!-----------------------------------------------------------------------!
!
      integer :: Ionoff
!
      if (IOPEN.ne.1) call wverr(2)
      LSCALE(ICURRWIN) = Ionoff
!
      return
   end subroutine scaleon
!
!
   subroutine newscale(Range)
      use wvglob
      implicit none
!
!-----------------------------------------------------------------------!
!   routine name       - newscale
!
!-----------------------------------------------------------------------!
!   computer           - machine independent
!
!   latest revision    - dec 90
!
!   purpose            - routine sets offset and scale factors
!
!   usage              - call newscale(Range)
!
!   arguments :
!     in:
!            Range     - xmin,ymin,xmax,ymax
!
!   required  routines -
!
!-----------------------------------------------------------------------!
      real(8) :: Range(4)
      integer :: i
!
      if (IOPEN.ne.1) call wverr(2)
      do 10 i=1,2
        XGSCAL(ICURRWIN,i)  =(LWSIZE(ICURRWIN,i+2)-LWSIZE(ICURRWIN,i))/ &
                             (Range(2+i)-Range(i))
        XGSCAL(ICURRWIN,i+2)=LWSIZE(ICURRWIN,i)-   &
                             Range(i)*XGSCAL(ICURRWIN,i)
10    continue
!
      call scaleon(1)
!
      return
   end subroutine newscale
!
!
   subroutine initcm
      implicit none
!
!-----------------------------------------------------------------------!
!   routine name       - initcm
!
!-----------------------------------------------------------------------!
!   computer           - machine independent
!
!   latest revision    - dec 90
!
!   purpose            - routine initiates colormaps for xwindow driver
!
!   usage              - call initcm
!
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!     defined colors are 0,1 - black,white
!      2,..,15 - colors for mesh plot
!     16,..,66 - 50 colors in the rampe sequence:
!                yellow, green, cyan, blue, magenta, red, yellow, white
!       exact definition of rampe is given in IRAM array below
!---------------------------------------------------------------------- !
!.....set of full rampe of colors
!
      integer :: nrcol,nrrow,i,iset,ir,ig,ib
      real(8) :: alfa
!
      integer, parameter :: iram(4,8) = reshape(         &
              (/    0,200,255,  0 ,  130,  0,255,  0 ,   &
                  220,  0,255,255 ,  320,  0, 60,250 ,   &
                  520,255,  0,255 ,  630,255,  0,  0 ,   &
                  870,255,255,  0 , 1000,255,255,255  /) &
               ,(/4,8/))
!
      call xgsetcm(0,  0,  0,  0)
      call xgsetcm(1,255,255,255)
!
!##!     for b/w machine we have only 2 colors
!##      if (NCOLOR.eq.0) return
!##
!##      if (NGREY.eq.1) then
!##          ncol = 14
!##          if (Iset.gt.10 ) ncol = Iset
!##          ndcol = 255 / (ncol+1)
!##          do 30 i = 1,ncol
!##              irgb = 255-ndcol*i
!##              call xgsetcm(nshft+i,irgb,irgb,irgb)
!##   30     continue
!##          return
!##      endif
!##!
!##      if (ISET.eq.1) then
!
!         color map for mesh plot
!
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
!
!##      endif
!##
!##      if (ISET.eq.2) then
!##!
!##!         color map for contour plot
!##!
!##         call xgsetcm(ncolshft+2, 255,255,0)
!##         call xgsetcm(ncolshft+3, 255,235,45)
!##         call xgsetcm(ncolshft+4, 255,195,75)
!##         call xgsetcm(ncolshft+5, 255,155,90)
!##         call xgsetcm(ncolshft+6, 255,105,75)
!##         call xgsetcm(ncolshft+7, 255,70, 75)
!##         call xgsetcm(ncolshft+8, 255,0,  25)
!##  	      call xgsetcm(ncolshft+9, 210,0,  85)
!##         call xgsetcm(ncolshft+10,200,55, 150)
!##         call xgsetcm(ncolshft+11,155,55, 220)
!##         call xgsetcm(ncolshft+12,55, 25, 230)
!##         call xgsetcm(ncolshft+13,0,  0,  255)
!##         call xgsetcm(ncolshft+14,0,  100,255)
!##         call xgsetcm(ncolshft+15,200,200,200)
!##
!##      endif
!##
!##      if (Iset.gt.13) then
!##!
!##!         soft ramp of colors
!##!
      nrcol = 0
      nrrow = 1
      iset = 51

      do 10 i=0,iset
        if (nrcol.ge.iram(1,nrrow+1)) nrrow = nrrow + 1
        alfa = float(nrcol-iram(1,nrrow))/   &
              (iram(1,nrrow+1)-iram(1,nrrow))
        ir = int(alfa*iram(2,nrrow+1)+(1-alfa)*iram(2,nrrow))
        ig = int(alfa*iram(3,nrrow+1)+(1-alfa)*iram(3,nrrow))
        ib = int(alfa*iram(4,nrrow+1)+(1-alfa)*iram(4,nrrow))
        call xgsetcm(16+i,ir,ig,ib)
        nrcol = nrcol + 1000/iset
   10 continue
!
!##      endif
!
      return
   end subroutine initcm
!
!
! ********************************************************************
   subroutine flushx
        implicit none
        return
   end subroutine flushx
!
!
! ********************************************************************
   subroutine wvsccl(Iclip)
      use rastcom
      implicit none
!
!-----------------------------------------------------------------------
!
!   routine name       - wvsccl
!
!-----------------------------------------------------------------------
!
!   computer           - machine independent
!
!   latest revision    -
!
!   purpose            -
!
!   usage              - call wvsccl(Iclip)
!
!   arguments :
!     in:
!              Iclip   -
!     out:
!
!   required  routines -
!
!-----------------------------------------------------------------------
!
      integer :: Iclip(4)
!
      integer :: iprint
      iprint=0
!
      if (iprint.eq.1) then
        write(*,*)'IN WVSCCL IFILE = ',IFILE
      endif
!
      if (IFILE.eq.1000) then
!
!  .....set clip for tektronics
        return
      endif
!
      if (IFILE.eq.1001) then
!
!  .....set clip for X-windows
        call xgsetclip(iclip(1),iclip(2),iclip(3),iclip(4))
        return
      endif
!
      return
   end subroutine wvsccl
!
!
!========================================================
   subroutine WVSCDL(X1,Y1,X2,Y2,IC)
      use rastcom
      implicit none
!========================================================
        integer :: X1,X2,Y1,Y2,IC

        if (IFILE.eq.1001) then
!           X-window
            call xgsetcol(ic)
            call xgdrawl(x1,y1,x2,y2)
            return
        endif
        return
   end subroutine WVSCDL
!
!
!========================================================
   subroutine WVSCFP(N,XC,ICF,ICL)
      use rastcom
      implicit none
!========================================================
        integer :: N,XC(*),ICF,ICL
        integer i,ix(100),iy(100)
!
        if (IFILE.eq.1001) then
!           fillpoly for X
            call xgsetcol(ICF)
            do 5 i=1,N
                ix(i)=XC(2*i-1)
                iy(i)=XC(2*i)
    5       continue
            call xgdrfpoly(N,ix,iy)
            if (ICF.ne.ICL) then
                call xgsetcol(ICL)
                call xgdrpoly(N,ix,iy)
            endif
            return
        endif
        return
   end subroutine WVSCFP
!
!
!========================================================
   subroutine WVSCPX(X,Y,IC)
      use rastcom
      implicit none
!=======================================================
        integer X,Y,IC
        if (IFILE.eq.1001) then
            call xgsetcol(IC)
            call xgdrawl(X,Y,X,Y)
            return
        endif
        return
   end subroutine WVSCPX
!
!
!========================================================
   subroutine WVSCTE(X1,Y1,H,MESS,A,N,ICOL)
      use rastcom
      implicit none
!========================================================
        integer   :: X1,Y1,H
        character :: MESS*(*)
        integer   :: A,N,ICOL
!
        character(len=100) :: text
!
        if (IFILE.eq.1001) then
!           text for X
            call xgsetcol(ICOL)
            call xgsettl(X1,Y1)
            text(1:N)=mess(1:N)
            text(N+1:N+1)=char(0)
            call xgdrtxt(text,N)
            return
        endif
        return
   end subroutine WVSCTE
!
!
   subroutine wvwarn(N)
        implicit none
!
!----------------------------------------------------------------------
!
!   routine name       - wvwarn
!
!----------------------------------------------------------------------
!
!   computer           - machine independent
!
!   latest revision    - feb, 92
!
!   purpose            - routine issues warning messages
!
!   usage              - call wvwarn(N)
!
!   arguments :
!     in:
!               N      - error flag
!
!   required  routines - wverr
!
!----------------------------------------------------------------------
!
      integer :: N
!
      integer, save :: IERRNO = 0
!
      IERRNO=IERRNO+1
!
!     if warning messages are desired add appropriate
!     write statements here
!
      if (IERRNO.lt.100) return
!
      write(*,*) 'WVWARN: More warnings than 100.'
      if (N.eq.0) write(*,*) 'Last message : Clip error'
      if (N.eq.1) write(*,*) 'Last message : Window number'
      if (N.eq.2) write(*,*) 'Last message : FILLPOLY size'
      call wverr(0)
!
   end subroutine wvwarn


   subroutine wverr(N)
        implicit none
!
!-----------------------------------------------------------------------!
!   routine name       - wverr
!
!-----------------------------------------------------------------------!
!   computer           - machine independent
!
!   latest revision    - dec 90
!
!   purpose            - error message routine
!
!   usage              - call wverr(N)
!
!   arguments :
!     in:
!                  N   - error flag
!     out:
!
!   required  routines - closwind
!
!-----------------------------------------------------------------------!
        integer :: N
        integer :: i
!
        do 10 i=0,9
          call closwind(i)
     10 continue
!
        write(*,*) 'GRAPHICS ERROR : '
        if (N.eq.0) write(*,*) ' END OF PROGRAM '
        if (N.eq.1) write(*,*) ' WRONG WINDOW NUMBER'
        if (N.eq.2) write(*,*) ' INITIALIZE GRAPHIC'
!
        stop 1
   end subroutine wverr
!---------------------------------------------------------------------
!
!   routine name       - wvgrph
!
!---------------------------------------------------------------------
!
!   computer           - machine independent
!
!   latest revision    - dec 90
!
!   purpose            - sets screen to graphic mode and sets screen
!                        sizes and max color in /wvscrn/
!
!   usage              - call wvgrph
!
!   arguments          - none
!
!   required  routines - inittek,xgclear,xginfo,xgnxtev,xgopw,xgsetcol
!
!---------------------------------------------------------------------
!
   subroutine wvgrph
      use graphmod, only: IWINDL, IWINDH
      use rastcom
      use wvscrn
      implicit none
!
        integer             :: xglength
        character(len=1024) :: res_string
!
        integer :: i1,i2,i3,i4,idum1,idum2
!
        integer :: iprint
        iprint=0
        if (iprint.eq.1) then
          write(*,*) 'wvgrph: IFILE = ',IFILE
        endif
!
        if (IFILE.eq.1001) then
!
!.......another entry into X
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
!
!---------------------------------------------------------------------
!
!  test for file T4014 in current dir
!  if it exists then open Tektronics graph
!
!      open (99,file='T4014',status='old',err=91)
!
!  file exists -
!      close (99)
!
!  open Tektronics - B-W draw only - no fill
!      IFILE=1000
!      MAXSCX=4000
!      MAXSCY=3000
!      MAXCOL=1
!      call inittek
!      return
!
!---------------------------------------------------------------------
!
!   91 continue
!
!  open X-window connection - hardwired
!!!      call xgopw(3,'GRAPHICS'//char(0),'-geometry'//char(0),   &
!!!          '700x500'//char(0),char(0),char(0),char(0))
!!!      call xgopw(3,'GRAPHICS'//char(0),'-geometry'//char(0),   &
!!!          '1000x600'//char(0),char(0),char(0),char(0))
!!!      call xgopw(3,'GRAPHICS'//char(0),'-geometry'//char(0),   &
!!!          '667x400'//char(0),char(0),char(0),char(0))
        write(res_string, fmt="(i0,'x',i0)") IWINDL, IWINDH ! jz 1/21/2013
        call xgopw(3,'GRAPHICS'//char(0),'-geometry'//char(0),   &
            trim(res_string)//char(0),char(0),char(0),char(0))
!     .    '1334x800'//char(0),char(0),char(0),char(0))
  100   call xgnxtev(i1,i2,i3,i4)
        if (xglength().gt.0) goto 100
        if (xglength().gt.0) write(*,*) ' EVENTS !'
        call xginfo(MAXSCX,MAXSCY,MAXCOL,idum1,idum2)
!
!  initialize color map
        call initcm
        IFILE=1001
!
!---------------------------------------------------------------------
!
        return
   end subroutine wvgrph
!
!========================================================
   subroutine wvtext
      use rastcom
      use wvscrn
      implicit none
!======================================================
!
!               performs QUIT
        integer :: xglength
        integer :: i1,i2,i3,i4,id1,id2
!
!          if it was opened....
        if (MAXCOL.ne.0) then
               if (IFILE.eq.1001) then
!               flush and wait for any event
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
!               do nothing - no close for xwindow
            endif
            write(*,*) 'Graphics closed'
        endif
        MAXCOL=0
        return
   end subroutine wvtext

#endif

#if HP3D_USE_X11

!******************************************************
!
! THE IMPORTANT ROUTINES FOR THE POSTSCRIPT FILES
!
! wvopps - to open a ps file
! wvpscl(clip(4)) - to set clipping region
! wvpsdl(X1,Y1,X2,Y2,IC) - to draw a line
! wvpsfp(N,Xc,Icf,Icl) - to fill a polygon
! wvpspx(X,Y,IC) - to draw a pixel
! wvpste(X1,Y1,H,MESS,A,N,ICOL) - to write text
! wvpsout(x,y,icode) - to write ps parameters into buffer
! flushps - to rewrite contents os the buffer into output ps file
! wvclps - to close the ps file
!
! OTHER PARAMETERS ARE DISCUSSED IN THE ROUTINES
!
!********************************************************

   subroutine wvopps
      use wvglob
      use wvpost
      use wvscrn
      implicit none
!---------------------------------------------------------
!      AUTHOR(S)        - Satish Chavva,Tadeusz Liszka and others
!
!      version  -   1998
!
!    purpose - to write header of ps file
!
!========================================================
! wr07.08.99
! wr03.14.02
!
      character(10) :: name
      integer       :: iflag
!
  110 write(*,*)  'Give name for Postscript file (0 - no postscript)'
   10 read(*,1011) name

 1011 format(a10)
      if (name(1:1).le.' ') goto 10
      if (name(1:1).eq.'0') then
          iflag=0
      else
          IFILE=79
          open(IFILE,file=name,status='new',err=100)
          write(*,*)' PostScript send to ',name
          write(IFILE,1000) name
!        Check if color Postscript required

          write(*,*) 'B/W or  Color Postscript (0 - B/W  1-Color)'
          read(*,*) ICOLVAL

 1000     format('%!  TL PostScript file, name ',a10/    &
              '%%BoundingBox: 36 36 518 756',/           &
              '0 setlinewidth',/                         &
              '/n {newpath moveto} def',/                &
              '/c {lineto closepath} def',/              &
              '/m {moveto} def',/                        &
              '/d {lineto} def',/                        &
              '/f {setgray fill} def',/                  &
              '/s {0 setgray stroke} def',/              &
              '/rf {setrgbcolor fill} def',/             &
              '/rs {setrgbcolor stroke} def',/           &
              '/gs {gsave} def',/                        &
              '/gr {grestore} def',/                     &
              '/fo { /Times-Roman findfont } def ',      &
              '/sc { scalefont setfont } def'/           &
              '50 750 translate -90 rotate fo 12 sc newpath')
          ILBUF=0
      endif
      ILEPSP=1
!      write(*,*) 'NO POSTSCRIPT GENERATED, iflag=',iflag
      return
  100 continue
      write(*,*) ' File ',name,' exists !!'
      write(*,*) ' I will not overwrite it !!'
      goto 110
   end subroutine wvopps


!-----------------------------------------------------------------
   subroutine wvpscl(Clip)
      use wvglob
      use wvpost
      use wvscrn
      implicit none
!
!----------------------------------------------------------------------
!
!   computer           - machine independent
!
!   latest revision    -
!
!   purpose            - to set clipping region for the ps page
!
!
!----------------------------------------------------------------------
!
      real(8) :: Clip(4)
!
      real(8) :: dum1,dum2
!
      integer :: iprint
      iprint=0
      if (iprint.eq.1) then
        write(*,*) 'WVPSCL: Clip = ',Clip
      endif
!
      if (IFILE.eq.0) return
!        call flushps
      call wvpsout(dum1,dum2,9)
      call wvpsout(Clip(1),Clip(2),1)
      call wvpsout(Clip(1),Clip(4),2)
      call wvpsout(Clip(3),Clip(4),2)
      call wvpsout(Clip(3),Clip(2),4)
      call wvpsout(dum1,dum2,7)
!
      if (iprint.eq.1) then
        write(*,*) 'WVPSCL: call pause'
        call pause
      endif
!
      return
   end subroutine wvpscl

!========================================================
   subroutine WVPSDL(X1,Y1,X2,Y2,IC)
      use wvglob
      use wvpost
      use wvscrn
      implicit none
!========================================================
!
!  author - Satish Chavva and many others
!
!  version - 1998
!
!  purpose - to draw a line
!
!  arguments - X1,Y1,X2,Y2-coordinates of the two endpoints
!              IC - line color
!---------------------------------------------------------
!
        real(8) :: X1,Y1,X2,Y2
        integer :: IC
!
        real(8) :: dum1,dum2
        integer :: iprint
        iprint=0
!
        if (IFILE.eq.0) return
        call wvpsout(x1,y1,1)
        call wvpsout(x2,y2,4)
        call wvpsout(dum1,dum2,6)
        return
   end subroutine WVPSDL


   subroutine wvpsfp(N,Xc,Icf,Icl)
      use wvglob
      use wvpost
      use wvscrn
      implicit none
!
!----------------------------------------------------------------------
!
!   routine name       - wvpsfp
!
!----------------------------------------------------------------------
!   author  - Satish Chavva and many others
!
!   computer           - machine independent
!
!   latest revision    -  98
!
!   purpose            - routine fills a polygon with a color
!                        in postcript
!                        It is assumed that the first and last
!                        vertex point coincide with each other !!
!
!   usage              - call wvpsfp(N,Xc,Icf,Icl)
!
!   arguments :
!     in:
!               N      - number of vertices
!               Coor   - coordinates of the vertices
!               Icf    - field color number
!               Icl    - line color number
!
!   required routines  - wvpsfp,wvwarn
!
!----------------------------------------------------------------------
!
      integer :: N
      real(8) :: Xc(2,N)
      integer :: Icf,Icl
!
      integer :: i
      real(8) :: dum1,dum2
!
      if (IFILE.eq.0) return
      if (Icf.ge.0) then
!
!  .....begin a newpath with the first point coordinates
        call wvpsout(Xc(1,1),Xc(2,1),1)
!
!  .....loop through the consequtive vertices of the polygon
        do i=2,N-1
!
!  .......line to the vertex point
          call wvpsout(Xc(1,i),Xc(2,i),2)
        enddo
!
!  .....line to the last vertex point and close the path
        call wvpsout(Xc(1,n),Xc(2,n),4)
!
!  .....fill up the polygon with a grey pattern
        call wvpsout(dum1,real(Icf,8),5)
      endif
!
!  ...draw the boundary line if a different color
      if (Icl.ne.Icf) then
        call wvpsout(Xc(1,1),Xc(2,1),1)
        do i=2,n-1
          call wvpsout(Xc(1,i),Xc(2,i),2)
        enddo
        call wvpsout(Xc(1,n),Xc(2,n),4)
        call wvpsout(dum1,dum2,6)
      endif
!
!
      return
   end subroutine wvpsfp

!========================================================
!     draw pixel in postcript
!========================================================
   subroutine WVPSPX(X,Y,IC)
      use wvglob
      use wvpost
      use wvscrn
      implicit none
!
        real(8) :: X,Y
        integer :: IC
!
        real(8) :: dum1,dum2
        integer :: iprint
        iprint=0

        call wvpsout(X,Y,1)
        call wvpsout(X,Y,4)
        call wvpsout(dum1,dum2,6)
        return
   end subroutine WVPSPX


!========================================================
!     draw text in postscript
!========================================================

   subroutine WVPSTE(X1,Y1,H,MESS,A,N,ICOL)
      use wvglob
      use wvpost
      use wvscrn
      implicit none
!
        real(8)   :: X1,Y1
        integer   :: H,A
        character :: MESS*(*)
        integer   :: N,ICOL
!
        integer   :: k
        real(8)   :: dum
!
        if (IFILE.eq.0) return
        call wvpsout(X1,Y1,3)
        k=12
        if (H.ge.16) k=15
        call wvpsout(real(k,8),dum,8)
        call flushps
        write(IFILE,1000) ' (',mess(1:N),') show'
 1000   format(A)
        return
   end subroutine WVPSTE


!-------------------------------------------------------------------

   subroutine wvpsout(X,Y,Icode)
      use wvglob
      use wvpost
      use wvscrn
      implicit none
!
!=====================================================================
!     author - Satish Chavva and many others
!
!     version - 1998
!
!=================================== post script dump to file ========
!       send real x,y (in inches) to postcript file
!       icode values 1:n 2:d 3:m 4:c 5:f 6:s 7:clip 8:font, 9:initclip
!       explanation of codes as in flushps:
!          1 - n   x,y,newpath,moveto
!          2 - d   x,y,lineto
!          3 - m   x,y,moveto
!          4 - c   x,y,lineto,closepath
!          5 - f   y,setgray,fill
!          6 - s   0,setgray,stroke
!          7 - clip
!          8 - font  == kx,font
!          9 - initclip
!       text printing is inside wvpstx routine
!       initialization - in wvopps
!
!=====================================================================
!
      real(8) :: X,Y
      integer :: Icode
!
      real(8) :: color,scale
      integer :: i,ii,iprint
!
      real(8), save :: ymax = 0.d0
!
      iprint=2

!      ICOLVAL=1
!
      if (IFILE.eq.0) return
      ILBUF=ILBUF+1
!
!     size is in inches, output in points (1/72")
      scale=72.d0

      if (Icode.eq.5) then
!         set grayscale colors
         if (ICOLVAL.eq.0) then
          color = 1.d0 - Y / 9.d0
          if (color .lt.0.0) color = 1.0 - (Y-16.d0)/50.d0
          if ((iprint.eq.2).and.(Y.gt.ymax)) then
            ymax=Y
            write(*,*) 'WVPSOUT: Y,color = ',Y,color
            call pause
          endif
          BUFOR(1,ILBUF)=color
         elseif (ICOLVAL.eq.1) then
          BUFOR(1,ILBUF)=Y
         endif
      elseif (Icode.eq.8) then
          BUFOR(1,ILBUF)=x
      elseif (Icode.le.4) then
          BUFOR(1,ILBUF)=X*scale
          BUFOR(2,ILBUF)=Y*scale
!      make sure they are in correct limits
          do i=1,2
              BUFOR(i,ILBUF)=max(-999.d0,min(9990.d0,BUFOR(i,ILBUF)))
          enddo
      endif
!
      if (iprint.eq.1) then
        write(*,*) 'WVPSOUT: ILBUF,BUFOR(*,ILBUF) = '
        write(*,*) ILBUF,(BUFOR(ii,ILBUF),ii=1,2)
      endif
!
      IBUFC(ILBUF)=Icode
      if (ILBUF.eq.5) call flushps
      return
   end subroutine wvpsout
!
!-----------------------------------------------------------------------
!
   subroutine flushps
!
!-----------------------------------------------------------------------
!
!   routine name       - flushps
!
!-----------------------------------------------------------------------
!   author  - Satish Chavva and many others
!
!   computer           - machine independent
!
!   latest revision    - 1998
!
!   purpose            - routine sends poscript buffer as one line to
!                        output
!
!   usage              - call flushps
!
!-----------------------------------------------------------------------
!
      use graphmod
      use wvglob
      use wvpost
      use wvscrn
      implicit none
!
      character(len=100) :: line,lineout
      character(len=20)  :: text
      integer   :: i,i7,ir1,ig1,ib1,ifrom,ito
      real(8)   :: r1,g1,b1
      integer   :: iprint,k,linelen,lentext
      real(8)   :: vmaxval
!
      character(len=1) :: code(4)
      code = (/'n','d','m','c'/)
!
      iprint=0

      vmaxval=float(IRED(1))

      if (IFILE.eq.0) then
        ILBUF = 0
        return
      endif
!
      if (ILBUF.ne.0) then
        line=' '
        linelen=0
        do i=1, ILBUF

          if (IBUFC(i).ne.0) then

            if (IBUFC(i).lt.5) then

              if (iprint.eq.1) then
                write(*,*) 'FLUSHPS: BUFOR(1,i),BUFOR(2,i),code(IBUFC(i)) = '
                       write(*,1010) BUFOR(1,i),BUFOR(2,i),code(IBUFC(i))
              endif
              write(IFILE,1010) &
              BUFOR(1,i),BUFOR(2,i),code(IBUFC(i))

 1010         format(2f8.1,1x,1a)
              lentext=18
            endif
            if (IBUFC(i).eq.5) then
             if (ICOLVAL.eq.1) then
              i7=int(BUFOR(1,i))+1
              ir1=IRED(i7)
              ig1=igreen(i7)
              ib1=iblue(i7)
!    here we make use of the color values and scale them to values
!    0 and 1
              r1=real(ir1,8)/vmaxval
              g1=real(ig1,8)/vmaxval
              b1=real(ib1,8)/vmaxval
             endif
            endif
            if (IBUFC(i).eq.5) then
             if (ICOLVAL.eq.0) then
              write(IFILE,1020) BUFOR(1,i)

 1020         format(f5.2,' f')
              lentext=7
             elseif (ICOLVAL.eq.1) then
              write(IFILE,1024) r1,g1,b1
 1024         format(3f6.3,1x,' rf')
              lentext=19
             endif
            elseif (IBUFC(i).eq.6) then

              write(IFILE,*)' s'
              lentext=2

            elseif (IBUFC(i).eq.7) then
              write(IFILE,*)' clip '
              lentext=5
            elseif (IBUFC(i).eq.8) then
              k=int(BUFOR(1,i))
              write(IFILE,1030) k
!              write(*,1030) k
!              write(*,*)"3"
!              read(*,*)

 1030         format(' fo',i3,' sc')
              lentext=9
            elseif (IBUFC(i).eq.9) then
              write(IFILE,1040)
!              write(*,1040)
!              read(*,*)

 1040         format(' initclip')
              lentext=9
            else
              lentext=0
            endif
!
            if (lentext.gt.0) then
              line(linelen+1:linelen+lentext)=text(1:lentext)
              linelen=linelen+lentext
            endif
           endif
        enddo
!
!.......here one may reduce size of output by:
!.......skipping repeated blanks
!.......skipping .0 in floating point numbers
        ifrom=1
        ito=2
        lineout=' '
!
  200   continue
        if (line(ifrom:ifrom+1).eq.'  ') then
          ifrom = ifrom+1
!!        elseif (line(ifrom:ifrom+2).eq.'.0 ') then
!!          ifrom = ifrom+2
!!        elseif (line(ifrom:ifrom+2).eq.'.00') then
!!          ifrom = ifrom+3
        elseif (line(ifrom:ifrom+1).eq.'. ') then
          ifrom = ifrom+1
        else
          lineout(ito:ito)=line(ifrom:ifrom)
          ito = ito+1
          ifrom = ifrom+1
        endif
!
        if (ifrom.lt.linelen) goto 200
        lineout(ito:ito)=line(linelen:linelen)
!
!.......initial character in line is always blank
        if(ito.gt.100) then
!           write(*,*)"Hi Cam"
!           read(*,*)
        endif
!        write(IFILE,1000) lineout(3:ito)
!        write(*,1000) lineout(3:ito)
!           write(*,*)"4"
!        read(*,*)

 1000   format(a)
        ILBUF=0
      endif
!
      return
   end subroutine flushps


!========================================================
   subroutine wvclps
      use wvglob
      use wvpost
      use wvscrn
      implicit none
!
        integer :: i
!
        ILEPSP=0
        do i=1,10
          IGPSP(i)=0
        enddo
        if (IFILE.eq.0) return
        call flushps
        write(IFILE,1000)
 1000   format(' stroke showpage'/)
        close(IFILE)
        write(*,*)' PostScript file closed'
        return
   end subroutine wvclps

#endif

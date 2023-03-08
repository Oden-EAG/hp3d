c----------------------------------------------------------------------
c
c   routine name       - object3
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Sep 07
c
c   purpose            - routine displays a 3D manifold
c
c
c   arguments :
c     in:
c               Iwind  - window number
c
c   required  routines -
c
c----------------------------------------------------------------------
c
      subroutine object3(Iwind)
c
      use GMP
      use graphmod

#include "syscom.blk"
c
c
      dimension xcl(3),igv(20),im(2),xxy(2)
c
      data bigp,bign,small /1.d30,-1.d30,1.d-14/
      data pi /3.14159265358979312d0/
c
c  ...set initial clipping plane for sections
      CLPL(1) = 0.d0
      CLPL(2) = 0.d0
      CLPL(3) = 0.d0
      CLPL(4) = 0.d0
c
c  ...initialize array of invisible blocks
      NRINVBL = 0
ccc      IGINV(1) = 0
c
c  ...initiate the number of curvilinear blocks
      NRCURVBL=0
ccc      eta(1:2) = 0.5d0
ccc      do nr=1,NRRECTA
ccc        if (RECTANGLES(nr)%BlockNo(2).eq.0) then
ccc          call recta_linear(nr,eta, x,dxdeta)
ccc          if (x(3).gt.-29.d0) then
ccc            NRCURVBL = NRCURVBL + 1
ccc            NLINBLOCKS(NRCURVBL) = abs(RECTANGLES(nr)%BlockNo(1))
ccc          endif
ccc        endif
ccc      enddo
c
c  ...set initial scaling constants
      DIMIM = bigp
      XCIM(1) = 0.d0
      XCIM(2) = 0.d0
c
c  ...set up number of line segments for drawing a curve
      nrsub1 = 2
      NRSUB = 2**nrsub1
      DX=1.d0/NRSUB
c
c
c  ...determine data for projection...

      theta = pi/5.d0+pi
      psi   = pi/6.d0

C       theta = pi/4.d0
C       psi = pi/4.d0
      RN(1) = cos(psi)*cos(theta)
      RN(2) = cos(psi)*sin(theta)
      RN(3) = sin(psi)
c
c  ...set all domains to be displayed
      if (NRDOMAIN.gt.MAXNRDOMAIN) then
        write(*,*) 'object3: INCREASE MAXNRDOMAIN'
        stop 1
      endif
      NDOMAIN(1:NRDOMAIN)=1
c
c  ...select the window
      call selwin(Iwind)
      write(*,*) '...PLEASE WAIT, PREPARING IMAGE'
c
c  ...create global image
      call cartobs
c
c  ...create the list of triangles to display
      call lsvistr
      XCIM(1) = XCENTR(1)
      XCIM(2) = XCENTR(2)
c
c  ...decide whether to scale by x or y axes
      q = DIMOB(1)/DIMOB(2)
      q1 = xlength/ylength
      if (q.gt.q1) then
        DIMIM = DIMOB(1)
        SIZE = xlength/2.d0
      else
        DIMIM = DIMOB(2)
        SIZE = ylength/2.d0
      endif
      XCWIN(1) = rmargin + xlength/2.d0
      XCWIN(2) = rmargin + ylength/2.d0
c
c  ...display image
      call dpvistr(0)
c
c  ...display menu in infinite loop
   10 continue
c
c  ...close window
      call dpborder(-1)
      write(*,*) 'PLEASE CLICK THE MOUSE INSIDE THE GRAPHICS'
      write(*,*) 'WINDOW TO CONTINUE ...'
      call closwind(iwindnum)
c
      write(*,*) 'SELECT OPTION :'
      write(*,*) '0 - EXIT'
      write(*,*) '1 - CHANGE THE POINT OF VIEW (WITH RESCALING)'
      write(*,*) '2 - CHANGE THE POINT OF VIEW (WITHOUT RESCALING)'
      write(*,*) '3 - CHANGE THE CENTRAL POINT OF THE IMAGE AND',
     .           ' RESCALE THE IMAGE'
      write(*,*) '41 - MAKE ONLY FEW BLOCKS INVISIBLE (MOUSE)'
      write(*,*) '42 - MAKE SOME BLOCKS INVISIBLE (MOUSE)'
      write(*,*) '5  - DRAW A CROSS-SECTION'
      write(*,*) '6  - DISPLAY THE WHOLE OBJECT'
      write(*,*) '7  - DISPLAY EDGES AND NUMBERS OF ALL BLOCKS'
      write(*,*) '8  - MAKE SOME BLOCKS INVISIBLE (MOUSE CONTROLLED)'
      write(*,*) '9 - MAKE ONLY FEW BLOCKS VISIBLE'
      write(*,*) '10 - MAKE ALL BLOCKS VISIBLE'
      write(*,*) '11 - CHANGE THE DEGREE OF COARSENESS FOR OVALS'
      write(*,*) '12 - UPGRADE THE BLOCK TO CURVILINEAR'
      write(*,*) '13 - SELECT SUBDOMAINS TO BE DISPLAYED'
      read(*,*) ichoice
c
      select case(ichoice)
c
      case(0)
        return
c
      case(8)
        write(*,*) 'DISPLAY BLOCK NUMBERS ? (0-NO/1-YES)'
        read(*,*) iblno
        if(iblno.eq.1) then
          call selwin(Iwind)
          call dpvistr(0)
          if(MANDIM.eq.3) call dpblno
          if(MANDIM.eq.2) call dptrno
          write(*,*) 'PLEASE CLICK THE MOUSE INSIDE THE GRAPHICS'
          write(*,*) 'WINDOW TO CONTINUE ...'
          call closwind(iwindnum)
        endif
        write(*,*) 'GIVE THE NUMBER OF BLOCKS TO MAKE INVISIBLE'
        read(*,*) inickbl
        write(*,*) 'GIVE NICKNAMES OF BLOCKS TO MAKE INVISIBLE'
        do  ibl=1,inickbl
          read(*,*) nickbl
          NRINVBL = NRINVBL + 1
          IGINV(NRINVBL) = nickbl
        enddo
c
      case(41)
        write(*,*) 'DISPLAY BLOCK NUMBERS ? (0-NO/1-YES)'
        read(*,*) iblno
        if (iblno.eq.1) then
          call selwin(iwind)
          call dpvistr(0)
          if(MANDIM.eq.3) call dpblno
          if(MANDIM.eq.2) call dptrno
          write(*,*) 'PLEASE CLICK THE MOUSE INSIDE THE GRAPHICS'
          write(*,*) 'WINDOW TO CONTINUE ...'
          call closwind(1)
        endif
        icount = 0
        do
          write(*,*) 'CLICK AT THE BLOCK TO BE DISPLAYED'
          call xmousepos(im(1),im(2))
          xxy(1) = (im(1) - XCWIN(1))/SIZE*DIMIM + XCIM(1)
          xxy(2) = (im(2) - XCWIN(2))/SIZE*DIMIM + XCIM(2)
          call findno(xxy, nel)
          nel = abs(nel)
          if (nel.gt.0) then
            icount = icount +1
            igv(icount) = nel
          else
            write(*,*)'OBJECT3: END OF SELECTION'
            go to 410
          endif
        enddo
 410    continue
        NRINVBL = 1
        do nel = 1,NRPRISM
          nick=10*nel+1
          ivis=0
          do ibl=1,icount
            if (nick.eq.igv(ibl)) ivis=1
          enddo
          if (ivis.eq.0) then
            NRINVBL=NRINVBL+1
            IGINV(NRINVBL) = nick
          endif
        enddo
        do nel = 1,NRHEXAS
          nick=10*nel+2
          ivis=0
          do ibl=1,icount
            if (nick.eq.igv(ibl)) ivis=1
          enddo
          if (ivis.eq.0) then
            NRINVBL=NRINVBL+1
            IGINV(NRINVBL) = nick
          endif
        enddo
c
      case(42)
        write(*,*) 'DISPLAY BLOCK NUMBERS ? (0-NO/1-YES)'
        read(*,*) iblno
        if (iblno.eq.1) then
          call selwin(Iwind)
          call dpvistr(0)
          if (MANDIM.eq.3) call dpblno
          if (MANDIM.eq.2) call dptrno
          write(*,*) 'PLEASE CLICK THE MOUSE INSIDE THE GRAPHICS'
          write(*,*) 'WINDOW TO CONTINUE ...'
          call closwind(iwindnum)
        endif
        icount = 0
        igv = 0
        do
          write(*,*) 'CLICK AT THE BLOCK TO REMOVE'
          call xmousepos(im(1),im(2))
          xxy(1) = (im(1) - XCWIN(1))/SIZE*DIMIM + XCIM(1)
          xxy(2) = (im(2) - XCWIN(2))/SIZE*DIMIM + XCIM(2)
          call findno(xxy, nel)
          nel = abs(nel)
          if (nel.gt.0) then
            icount = icount +1
            igv(icount) = nel
          else
            go to 420
          endif
        enddo
 420    continue
        do i=1, icount
          NRINVBL = NRINVBL + 1
          IGINV(NRINVBL) = igv(i)
        enddo
c
      case(9)
        write(*,*) 'DISPLAY BLOCK NUMBERS ? (0-NO/1-YES)'
        read(*,*) iblno
        if(iblno.eq.1) then
          call selwin(iwind)
          call dpvistr(0)
          if(MANDIM.eq.3) call dpblno
          if(MANDIM.eq.2) call dptrno
          write(*,*) 'PLEASE CLICK THE MOUSE INSIDE THE GRAPHICS'
          write(*,*) 'WINDOW TO CONTINUE ...'
          call closwind(1)
        endif
        write(*,*) 'GIVE THE NUMBER OF BLOCKS TO STAY VISIBLE'
        read(*,*) inickbl
        write(*,*) 'GIVE NICKNAMES OF BLOCKS TO STAY VISIBLE'
        do ibl=1,inickbl
          read(*,*) nickbl
          igv(ibl) = nickbl
        enddo
        NRINVBL = 1
        do npri=1,NRPRISM
          nick=10*npri+1
          call locate(nick,igv,inickbl, loc)
          if (loc.eq.0) then
            NRINVBL=NRINVBL+1
            IGINV(NRINVBL) = nick
          endif
        enddo
        do nh = 1,NRHEXAS
          nick=10*nh+2
          call locate(nick,igv,inickbl, loc)
          if (loc.eq.0) then
            NRINVBL=NRINVBL+1
            IGINV(NRINVBL) = nick
          endif
        enddo
        do ntet = 1,NRTETRA
          nick=10*ntet+3
          call locate(nick,igv,inickbl, loc)
          if (loc.eq.0) then
            NRINVBL=NRINVBL+1
            IGINV(NRINVBL) = nick
          endif
        enddo
        do npyr = 1,NRPYRAM
          nick=10*npyr+4
          call locate(nick,igv,inickbl, loc)
          if (loc.eq.0) then
            NRINVBL=NRINVBL+1
            IGINV(NRINVBL) = nick
          endif
        enddo
c
      case(1,2)
        write(*,6024) psi/pi*180,theta/pi*180,RN(1:3)
 6024   format('CURRENT psi,theta,RN = ',2f6.1, 2x, 3f8.3)
        write(*,*) 'SET NEW psi,theta IN DEGREES'
        read(*,*)  psi,theta
        psi = psi/180*pi
        theta = theta/180*pi
        RN(1) = cos(psi)*cos(theta)
        RN(2) = cos(psi)*sin(theta)
        RN(3) = sin(psi)
        write(*,6025) psi/pi*180,theta/pi*180,RN(1:3)
 6025   format('NEW psi,theta,RN = ',2f6.1, 2x, 3f8.3)
        call cartobs
        if (ichoice.eq.1) DIMIM = bigp
c
      case(11)
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
c
      case(3)
        write(*,*) 'CLICK AT THE NEW CENTRAL POINT'
        call xmousepos(im(1),im(2))
        XCIM(1) = (im(1) - XCWIN(1))/SIZE*DIMIM + XCIM(1)
        XCIM(2) = (im(2) - XCWIN(2))/SIZE*DIMIM + XCIM(2)
        write(*,*) 'XCIM = ',XCIM
        write(*,*) 'SET MAGNIFICATION FACTOR'
        read(*,*) rnsc
        rnsc=dabs(rnsc)
        DIMIM=DIMIM/rnsc
c
      case(5)
        write(*,*) 'THE EXTREME VALUES OF COORDINATES '
        write(*,*) 'ARE (IN ORIGINAL SYSTEM OXYZ) :'
        write(*,*) 'X - ',XEX(5),XEX(6)
        write(*,*) 'Y - ',XEX(1),XEX(2)
        write(*,*) 'Z - ',XEX(3),XEX(4)
        write(*,*) 'SET NORMAL TO THE SECTION PLANE'
        read(*,*) CLPL(1),CLPL(2),CLPL(3)
        write(*,*) 'SET COORDINATES OF A POINT ',
     .             'BELONGING TO THE SECTION PLANE'
        read(*,*) xcl(1),xcl(2),xcl(3)
        CLPL(4)=-CLPL(1)*xcl(1)-CLPL(2)*xcl(2)-CLPL(3)*xcl(3)
c
      case(6)
        DIMIM = bigp
        CLPL(1) = 0.d0
        CLPL(2) = 0.d0
        CLPL(3) = 0.d0
        CLPL(4) = 0.d0
c
      case(10)
        NRINVBL = 1
        IGINV(1) = 0
c
      case(12)
        do
          write(*,*) 'CLICK AT THE BLOCK TO BE DISPLAYED'
          call xmousepos(im(1),im(2))
          xxy(1) = (im(1) - XCWIN(1))/SIZE*DIMIM + XCIM(1)
          xxy(2) = (im(2) - XCWIN(2))/SIZE*DIMIM + XCIM(2)
          call findno(xxy, nh)
          nh = abs(nh)
          if (nh.gt.0) then
            NRCURVBL = NRCURVBL + 1
            NLINBLOCKS(NRCURVBL) = nh*10+2
          else
            write(*,*)'OBJECT3: END OF SELECTION'
            exit
          endif
        enddo
c
      case(13)
        write(*,*) 'SELECT WHICH DOMAINS TO DISPLAY'
        write(*,*) 'NRDOMAIN = ',NRDOMAIN
        read(*,*) NDOMAIN(1:NRDOMAIN)
        write(*,7038) NDOMAIN(1:NRDOMAIN)
 7038   format('object3: NDOMAIN = ',20i2)
c
      end select
c
c-----------------------------------------------------------------------
c
      write(*,*) '...PLEASE WAIT, PREPARING IMAGE'
      if (ichoice.ne.7) call lsvistr
c
c  ...rescale when necessary
      if ((ichoice.eq.1).or.(ichoice.eq.6)) then
        XCIM(1) = XCENTR(1)
        XCIM(2) = XCENTR(2)
c
c  .....decide whether to scale by x or y axes
        q = DIMOB(1)/DIMOB(2)
        q1 = xlength/ylength
        if (q.gt.q1) then
          DIMIM = DIMOB(1)
          SIZE = xlength/2.d0
        else
          DIMIM = DIMOB(2)
          SIZE = ylength/2.d0
        endif
c
      endif
c
c  ...open the window again
      call selwin(Iwind)
c
      if (ichoice.ne.7) then
        write(*,*) 'DISPLAY BLOCK NUMBERS ? (0-NO/1-YES)'
        read(*,*) iblno
        call dpvistr(iblno)
      else
        call dpvistr(0)
        if (MANDIM.eq.3) call dpblno
        if (MANDIM.eq.2) call dptrno
        call dpallcu
      endif
c
      goto 10
c
c
      end







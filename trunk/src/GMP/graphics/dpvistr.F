c----------------------------------------------------------------------
c
c   routine name       - dpvistr
c
c----------------------------------------------------------------------
c
c   latest revision    - Jan 08
c
c   purpose            - routine displays triangles from the RGTR
c                        list in order stored in IGTR list
c
c   argument
c     in :
c                  Nof - flag whether to display block numbers or not
c
c   required  routines - none
c
c----------------------------------------------------------------------
c
      subroutine dpvistr(Nof)
c
      use graphmod
c
#include "syscom.blk"
c
c
c  ...flags for drawing boundary segments
      dimension nplotb(3)
      character*5 text
c
      iprint=0
c
c  ...define colors
      nblack  = npcol(2)
      ngreen  = npcol(5)
      nyellow = npcol(6)
      norange = npcol(7)
      nred    = npcol(8)
c
c  ...loop through all visible triangles
      do i=1,NRVISTR
c
c  .....get the triangle number
        itr=IGTR(i)
c
c  .....rescale
        do k=1,3
          do ivar=1,2
            XY(ivar,k) = RGTR((itr-1)*9+(k-1)*3+ivar)
            XY(ivar,k) = (XY(ivar,k)-XCIM(ivar))/DIMIM*SIZE
     .                   + XCWIN(ivar)
          enddo
        enddo
c
        iblno = IGTRCU(itr)/10000
        IGTRCU(itr) = IGTRCU(itr)-10000*iblno
c
c  .....set color
        icolor = IGTRCU(itr)/1000
        select case(icolor)
        case(1)
          ncol = ngreen
        case(2)
          ncol = nyellow
        case(3)
          ncol = norange
        case(4)
          ncol = nred
        case default
          write(*,*) 'itr,IGTRCU(itr) = ',itr,IGTRCU(itr)
          write(*,*) 'dpvistr: UNKNOWN COLOR '
cccc          stop 1
        end select
c
c  .....fill the subtriangle
        call fillpoly(3,XY,ncol,ncol)
c
c  .....clear the buffer
        call flushx
c
c  .....check whether subsequent edges are to be displayed
c       and draw them when needed
        call decode(IGTRCU(itr), nvoid,nickb)
        call decodg(nickb,2,3, nplotb)
c
        if (nplotb(1).eq.1)
     .    call drawline(XY(1,1),XY(2,1),XY(1,2),XY(2,2),nblack)
        if (nplotb(2).eq.1)
     .    call drawline(XY(1,2),XY(2,2),XY(1,3),XY(2,3),nblack)
        if (nplotb(3).eq.1)
     .    call drawline(XY(1,3),XY(2,3),XY(1,1),XY(2,1),nblack)
c
c  .....write block number when required
        shiftx=0.017d0*xlength
        shifty=0.004d0*ylength
        if ((Nof.eq.1).and.(iblno.ne.0)) then
          ncar=5
          write(text,'(i5)') iblno
          xtext = (xy(1,1)+xy(1,2)+xy(1,3))/3.d0 - shiftx
          ytext = (xy(2,1)+xy(2,2)+xy(2,3))/3.d0 - shifty
          height=0.5
          angle=0.
          call symbol(xtext,ytext,height,text,angle,ncar,nblack)
        endif
c
c  ...end of loop through small triangles
      enddo
c
      return
      end
